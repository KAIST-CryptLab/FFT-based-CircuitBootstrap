use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use rand::Rng;
use tfhe::core_crypto::prelude::*;
use auto_base_conv::{
    aes_instances::*, automorphism::gen_all_auto_keys, byte_array_to_mat, generate_scheme_switching_key, generate_vec_keyed_lut_accumulator, get_he_state_error, he_add_round_key, he_mix_columns_precomp, he_shift_rows, he_sub_bytes_8_to_24_by_patched_wwlp_cbs, he_sub_bytes_by_patched_wwlp_cbs, keygen_pbs_with_glwe_ds, keyswitch_lwe_ciphertext_by_glwe_keyswitch, known_rotate_keyed_lut, Aes128Ref, BLOCKSIZE_IN_BIT, BLOCKSIZE_IN_BYTE, BYTESIZE, NUM_ROUNDS
};

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(1000);
    targets =
        criterion_benchmark_aes,
);
criterion_main!(benches);

fn criterion_benchmark_aes(c: &mut Criterion) {
    let mut group = c.benchmark_group("aes evaluation by patched WWL+ circuit bootstrapping");

    let param_list = [
        (*AES_TIGHT, "tight"),
    ];

    for (param, id) in param_list.iter() {
        let lwe_dimension = param.lwe_dimension();
        let lwe_modular_std_dev = param.lwe_modular_std_dev();
        let glwe_dimension = param.glwe_dimension();
        let polynomial_size = param.polynomial_size();
        let glwe_modular_std_dev = param.glwe_modular_std_dev();
        let pbs_base_log = param.pbs_base_log();
        let pbs_level = param.pbs_level();
        let glwe_ds_base_log = param.glwe_ds_base_log();
        let glwe_ds_level = param.glwe_ds_level();
        let common_polynomial_size = param.common_polynomial_size();
        let fft_type_ds = param.fft_type_ds();
        let auto_base_log = param.auto_base_log();
        let auto_level = param.auto_level();
        let fft_type_auto = param.fft_type_auto();
        let ss_base_log = param.ss_base_log();
        let ss_level = param.ss_level();
        let cbs_base_log = param.cbs_base_log();
        let cbs_level = param.cbs_level();
        let log_lut_count = param.log_lut_count();
        let ciphertext_modulus = param.ciphertext_modulus();

        // Set random generators and buffers
        let mut boxed_seeder = new_seeder();
        let seeder = boxed_seeder.as_mut();

        let mut secret_generator = SecretRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed());
        let mut encryption_generator = EncryptionRandomGenerator::<ActivatedRandomGenerator>::new(seeder.seed(), seeder);

        // Generate keys
        let (
            lwe_sk,
            glwe_sk,
            lwe_sk_after_ks,
            fourier_bsk,
            fourier_ksk,
        ) = keygen_pbs_with_glwe_ds(
            lwe_dimension,
            glwe_dimension,
            polynomial_size,
            lwe_modular_std_dev,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            glwe_ds_base_log,
            glwe_ds_level,
            common_polynomial_size,
            fft_type_ds,
            ciphertext_modulus,
            &mut secret_generator,
            &mut encryption_generator,
        );
        let fourier_bsk = fourier_bsk.as_view();

        let ss_key = generate_scheme_switching_key(
            &glwe_sk,
            ss_base_log,
            ss_level,
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );
        let ss_key = ss_key.as_view();

        let auto_keys = gen_all_auto_keys(
            auto_base_log,
            auto_level,
            fft_type_auto,
            &glwe_sk,
            glwe_modular_std_dev,
            &mut encryption_generator,
        );

        // ======== Plain ========
        let mut rng = rand::thread_rng();
        let mut key = [0u8; BLOCKSIZE_IN_BYTE];
        for i in 0..BLOCKSIZE_IN_BYTE {
            key[i] = rng.gen_range(0..=u8::MAX);
        }

        let aes = Aes128Ref::new(&key);
        let round_keys = aes.get_round_keys();

        let mut message = [0u8; BLOCKSIZE_IN_BYTE];
        for i in 0..16 {
            message[i] = rng.gen_range(0..=255);
        }
        let correct_output = byte_array_to_mat(aes.encrypt_block(message));

        // ======== HE ========
        let mut he_round_keys = Vec::<LweCiphertextListOwned<u64>>::with_capacity(NUM_ROUNDS + 1);
        for r in 0..=NUM_ROUNDS {
            let mut lwe_list_rk = LweCiphertextList::new(
                0u64,
                fourier_bsk.output_lwe_dimension().to_lwe_size(),
                LweCiphertextCount(BLOCKSIZE_IN_BIT),
                ciphertext_modulus,
            );

            let rk = PlaintextList::from_container((0..BLOCKSIZE_IN_BIT).map(|i| {
                let byte_idx = i / BYTESIZE;
                let bit_idx = i % BYTESIZE;
                let round_key_byte = round_keys[r][byte_idx];
                let round_key_bit = (round_key_byte & (1 << bit_idx)) >> bit_idx;
                (round_key_bit as u64) << 63
            }).collect::<Vec<u64>>());
            encrypt_lwe_ciphertext_list(
                &lwe_sk,
                &mut lwe_list_rk,
                &rk,
                glwe_modular_std_dev,
                &mut encryption_generator,
            );

            he_round_keys.push(lwe_list_rk);
        }

        let mut he_state = LweCiphertextList::new(
            0u64,
            fourier_bsk.output_lwe_dimension().to_lwe_size(),
            LweCiphertextCount(BLOCKSIZE_IN_BIT),
            ciphertext_modulus,
        );
        let mut he_state_mult_by_2 = LweCiphertextList::new(
            0u64,
            fourier_bsk.output_lwe_dimension().to_lwe_size(),
            LweCiphertextCount(BLOCKSIZE_IN_BIT),
            ciphertext_modulus,
        );
        let mut he_state_mult_by_3 = LweCiphertextList::new(
            0u64,
            fourier_bsk.output_lwe_dimension().to_lwe_size(),
            LweCiphertextCount(BLOCKSIZE_IN_BIT),
            ciphertext_modulus,
        );
        let mut he_state_ks = LweCiphertextList::new(
            0u64,
            lwe_sk_after_ks.lwe_dimension().to_lwe_size(),
            LweCiphertextCount(BLOCKSIZE_IN_BIT),
            ciphertext_modulus,
        );

        let vec_keyed_sbox_round_1 = generate_vec_keyed_lut_accumulator(
            aes.get_keyed_sbox(0),
            u64::BITS as usize - 1,
            &glwe_sk,
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );
        let vec_keyed_sbox_round_1_mult_by_2 = generate_vec_keyed_lut_accumulator(
            aes.get_keyed_sbox_mult_by_2(0),
            u64::BITS as usize - 1,
            &glwe_sk,
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );
        let vec_keyed_sbox_round_1_mult_by_3 = generate_vec_keyed_lut_accumulator(
            aes.get_keyed_sbox_mult_by_3(0),
            u64::BITS as usize - 1,
            &glwe_sk,
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );

        // Bench
        he_state.as_mut().fill(0u64);
        for (bit_idx, mut he_bit) in he_state.iter_mut().enumerate() {
            let byte_idx = bit_idx / 8;
            let pt = (message[byte_idx] & (1 << bit_idx)) >> bit_idx;
            *he_bit.get_mut_body().data += (pt as u64) << 63;
        }

        { // r = 1
            // Keyed LUT
            group.bench_function(
                BenchmarkId::new(
                    format!("Round 1 Keyed-LUT"),
                    id,
                ),
                |b| b.iter(|| {
                    known_rotate_keyed_lut(
                        black_box(message),
                        black_box(&vec_keyed_sbox_round_1),
                        black_box(&mut he_state),
                    );
                    known_rotate_keyed_lut(
                        black_box(message),
                        black_box(&vec_keyed_sbox_round_1_mult_by_2),
                        black_box(&mut he_state_mult_by_2),
                    );
                    known_rotate_keyed_lut(
                        black_box(message),
                        black_box(&vec_keyed_sbox_round_1_mult_by_3),
                        black_box(&mut he_state_mult_by_3),
                    );
                }),
            );

            // ShiftRows, MixColumns, AddRoundKey
            group.bench_function(
                BenchmarkId::new(
                    "Ruond 1 ShiftRows, MixColumns and AddRoundKey",
                    id,
                ),
                |b| b.iter(|| {
                    // ShiftRows
                    he_shift_rows(black_box(&mut he_state));
                    he_shift_rows(black_box(&mut he_state_mult_by_2));
                    he_shift_rows(black_box(&mut he_state_mult_by_3));

                    // MixColumns
                    he_mix_columns_precomp(
                        black_box(&mut he_state),
                        black_box(&he_state_mult_by_2),
                        black_box(&he_state_mult_by_3),
                    );

                    // AddRoundKey
                    he_add_round_key(
                        black_box(&mut he_state),
                        black_box(&he_round_keys[1]),
                    );
                }),
            );
        }

        for r in 2..NUM_ROUNDS {
            // LWE KS
            group.bench_function(
                BenchmarkId::new(
                    format!("Round {r} LWE Keyswitching"),
                    id,
                ),
                |b| b.iter(|| {
                    for (lwe, mut lwe_ks) in he_state.iter().zip(he_state_ks.iter_mut()) {
                        keyswitch_lwe_ciphertext_by_glwe_keyswitch(
                            black_box(&lwe),
                            black_box(&mut lwe_ks),
                            black_box(&fourier_ksk),
                        );
                    }
                })
            );

            // SubBytes
            group.bench_function(
                BenchmarkId::new(
                    format!("Round {r} SubBytes"),
                    id,
                ),
                |b| b.iter(|| {
                    he_sub_bytes_8_to_24_by_patched_wwlp_cbs(
                        black_box(&he_state_ks),
                        black_box(&mut he_state),
                        black_box(&mut he_state_mult_by_2),
                        black_box(&mut he_state_mult_by_3),
                        black_box(fourier_bsk),
                        black_box(&auto_keys),
                        black_box(ss_key),
                        black_box(cbs_base_log),
                        black_box(cbs_level),
                        black_box(log_lut_count),
                    );
                })
            );

            // Linear
            group.bench_function(
                BenchmarkId::new(
                    format!("Round {r} ShiftRows, MixColumns, AddRoundKey"),
                    id,
                ),
                |b| b.iter(|| {
                    // ShiftRows
                    he_shift_rows(black_box(&mut he_state));
                    he_shift_rows(black_box(&mut he_state_mult_by_2));
                    he_shift_rows(black_box(&mut he_state_mult_by_3));

                    // MixColumns
                    he_mix_columns_precomp(
                        black_box(&mut he_state),
                        black_box(&he_state_mult_by_2),
                        black_box(&he_state_mult_by_3),
                    );

                    // AddRoundKey
                    he_add_round_key(
                        black_box(&mut he_state),
                        black_box(&he_round_keys[r]),
                    );
                })
            );
        }

        // LWE KS
        group.bench_function(
            BenchmarkId::new(
                format!("Final Round LWE Keyswitching"),
                id,
            ),
            |b| b.iter(|| {
                for (lwe, mut lwe_ks) in he_state.iter().zip(he_state_ks.iter_mut()) {
                    keyswitch_lwe_ciphertext_by_glwe_keyswitch(
                        black_box(&lwe),
                        black_box(&mut lwe_ks),
                        black_box(&fourier_ksk),
                    );
                }
            }));

        // SubBytes
        group.bench_function(
            BenchmarkId::new(
                format!("Final Round SubBytes"),
                id,
            ),
            |b| b.iter(|| {
                he_sub_bytes_by_patched_wwlp_cbs(
                    black_box(&he_state_ks),
                    black_box(&mut he_state),
                    black_box(fourier_bsk),
                    black_box(&auto_keys),
                    black_box(ss_key),
                    black_box(cbs_base_log),
                    black_box(cbs_level),
                    black_box(log_lut_count),
                );
            })
        );

        group.bench_function(
            BenchmarkId::new(
                format!("Final Round ShiftRows, AddRoundKey"),
                id,
            ),
            |b| b.iter(|| {
                // ShiftRows
                he_shift_rows(&mut he_state);

                // AddRoundKey
                he_add_round_key(&mut he_state, &he_round_keys[NUM_ROUNDS]);
            })
        );

        let (_, max_err) = get_he_state_error(&he_state, correct_output, &lwe_sk);

        println!(
            "n: {}, N: {}, k: {}, l_pbs: {}, B_pbs: 2^{}, l_cbs: {}, B_cbs: 2^{}
B_ds: 2^{}, l_ds: {},
l_auto: {}, B_auto: 2^{}, l_ss: {}, B_ss: 2^{}, log_lut_count: {},
max err: {:.2} bits",
            lwe_dimension.0, polynomial_size.0, glwe_dimension.0, pbs_level.0, pbs_base_log.0, cbs_level.0, cbs_base_log.0,
            glwe_ds_base_log.0, glwe_ds_level.0,
            auto_level.0, auto_base_log.0, ss_level.0, ss_base_log.0, log_lut_count.0,
            (max_err as f64).log2(),
        );
        println!();
    }
}