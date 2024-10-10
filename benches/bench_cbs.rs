use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use tfhe::core_crypto::prelude::*;
use tfhe::core_crypto::fft_impl::fft64::c64;
use auto_base_conv::{blind_rotate_for_msb, convert_to_ggsw_after_blind_rotate, gen_all_auto_keys, generate_scheme_switching_key, get_max_err_ggsw_bit, glwe_ciphertext_clone_from, glwe_ciphertext_monic_monomial_div, keygen_pbs, wopbs_instance::*};


criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(1000);
    targets =
        criterion_benchmark_improved_wopbs,
);
criterion_main!(benches);


#[allow(unused)]
fn criterion_benchmark_improved_wopbs(c: &mut Criterion) {
    let mut group = c.benchmark_group("wopbs");

    let param_list = [
        (*BITWISE_CBS_CMUX1, 1, "CMUX1"),
        (*BITWISE_CBS_CMUX2, 1, "CMUX2"),
        (*BITWISE_CBS_CMUX3, 1, "CMUX3"),
    ];

    for (param, extract_size, id) in param_list.iter() {
        let lwe_dimension = param.lwe_dimension();
        let lwe_modular_std_dev = param.lwe_modular_std_dev();
        let glwe_dimension = param.glwe_dimension();
        let polynomial_size = param.polynomial_size();
        let glwe_modular_std_dev = param.glwe_modular_std_dev();
        let pbs_base_log = param.pbs_base_log();
        let pbs_level = param.pbs_level();
        let ks_base_log = param.ks_base_log();
        let ks_level = param.ks_level();
        let auto_base_log = param.auto_base_log();
        let auto_level = param.auto_level();
        let auto_fft_type = param.fft_type_auto();
        let ss_base_log = param.ss_base_log();
        let ss_level = param.ss_level();
        let cbs_base_log = param.cbs_base_log();
        let cbs_level = param.cbs_level();
        let log_lut_count = param.log_lut_count();
        let ciphertext_modulus = param.ciphertext_modulus();
        let message_size = param.message_size();

        let extract_size = *extract_size;
        let glwe_size = glwe_dimension.to_glwe_size();

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
            bsk,
            ksk,
        ) = keygen_pbs(
            lwe_dimension,
            glwe_dimension,
            polynomial_size,
            lwe_modular_std_dev,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            ks_base_log,
            ks_level,
            &mut secret_generator,
            &mut encryption_generator,
        );
        let bsk = bsk.as_view();

        let ksk = allocate_and_generate_new_lwe_keyswitch_key(
            &lwe_sk,
            &lwe_sk_after_ks,
            ks_base_log,
            ks_level,
            lwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );

        let auto_keys = gen_all_auto_keys(
            auto_base_log,
            auto_level,
            auto_fft_type,
            &glwe_sk,
            glwe_modular_std_dev,
            &mut encryption_generator,
        );

        let ss_key = generate_scheme_switching_key(
            &glwe_sk,
            ss_base_log,
            ss_level,
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );
        let ss_key = ss_key.as_view();

        // Set input LWE ciphertext
        let msg = (1 << message_size) - 1;
        let lwe = allocate_and_encrypt_new_lwe_ciphertext(
            &lwe_sk,
            Plaintext(msg << (u64::BITS as usize - message_size)),
            glwe_modular_std_dev,
            ciphertext_modulus,
            &mut encryption_generator,
        );

        // Bench
        let fft = Fft::new(polynomial_size);
        let fft = fft.as_view();

        let mut ggsw_list_out = GgswCiphertextList::new(0u64, glwe_size, polynomial_size, cbs_base_log, cbs_level, GgswCiphertextCount(message_size), ciphertext_modulus);

        let mut buf = LweCiphertext::new(u64::ZERO, lwe.lwe_size(), ciphertext_modulus);
        buf.as_mut().clone_from_slice(lwe.as_ref());

        for (idx, mut ggsw_chunk) in ggsw_list_out.chunks_exact_mut(extract_size).enumerate() {
            let mut acc_glev = GlweCiphertextList::new(u64::ZERO, glwe_size, polynomial_size, GlweCiphertextCount(cbs_level.0), ciphertext_modulus);

            group.bench_function(
                BenchmarkId::new(
                    format!("Idx {idx} Extr. + Refr."),
                    id,
                ),
                |b| b.iter(|| {
                    let mut lwe_extract = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);
                    lwe_ciphertext_cleartext_mul(&mut lwe_extract, &buf, Cleartext(1u64 << (message_size - extract_size * (idx + 1))));

                    let mut lwe_extract_ks = LweCiphertext::new(u64::ZERO, ksk.output_lwe_size(), ciphertext_modulus);
                    keyswitch_lwe_ciphertext(&ksk, &lwe_extract, &mut lwe_extract_ks);

                    blind_rotate_for_msb(
                        black_box(&lwe_extract_ks),
                        black_box(&mut acc_glev),
                        black_box(bsk),
                        black_box(log_lut_count),
                        black_box(cbs_base_log),
                        black_box(cbs_level),
                        black_box(extract_size),
                        black_box(ciphertext_modulus),
                    );
                }),
            );

            let mut fourier_ggsw_chunk_out = FourierGgswCiphertextList::new(
                vec![c64::default();
                    extract_size * polynomial_size.to_fourier_polynomial_size().0
                        * glwe_size.0
                        * glwe_size.0
                        * cbs_level.0
                ],
                extract_size,
                glwe_size,
                polynomial_size,
                cbs_base_log,
                cbs_level,
            );
            let mut buf_next = LweCiphertext::new(u64::ZERO, buf.lwe_size(), ciphertext_modulus);

            group.bench_function(
                BenchmarkId::new(
                    format!("Idx {idx} Conv."),
                    id,
                ),
                |b| b.iter(|| {

                    let log_scale = u64::BITS as usize - message_size + idx * extract_size;
                    let acc_plaintext = PlaintextList::from_container((0..polynomial_size.0).map(|i| {
                        if i < (1 << extract_size) {
                            if (i >> (extract_size - 1)) == 0 {
                                (i << log_scale) as u64
                            } else {
                                (((1 << (extract_size - 1)) + ((1 << extract_size) - 1 - i)) << log_scale) as u64
                            }
                        } else {
                            u64::ZERO
                        }
                    }).collect::<Vec<u64>>());
                    let acc_id = allocate_and_trivially_encrypt_new_glwe_ciphertext(
                        glwe_size,
                        &acc_plaintext,
                        ciphertext_modulus,
                    );

                    let mut ct0 = GlweCiphertext::new(u64::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
                    let mut ct1 = GlweCiphertext::new(u64::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
                    glwe_ciphertext_clone_from(&mut ct0, &acc_id);

                    for (i, (mut ggsw, mut fourier_ggsw)) in ggsw_chunk.iter_mut()
                    .zip(fourier_ggsw_chunk_out.as_mut_view().into_ggsw_iter())
                    .enumerate()
                    {
                        convert_to_ggsw_after_blind_rotate(
                            &acc_glev,
                            &mut ggsw,
                            extract_size - i - 1,
                            &auto_keys,
                            ss_key,
                            ciphertext_modulus,
                        );

                        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);

                        glwe_ciphertext_monic_monomial_div(&mut ct1, &ct0, MonomialDegree(1 << i));
                        cmux_assign(&mut ct0, &mut ct1, &fourier_ggsw);
                    }

                    let mut lwe_extract = LweCiphertext::new(u64::ZERO, lwe.lwe_size(), ciphertext_modulus);
                    extract_lwe_sample_from_glwe_ciphertext(&ct0, &mut lwe_extract, MonomialDegree(0));

                    lwe_ciphertext_sub(&mut buf_next, &buf, &lwe_extract);
                }),
            );
            buf.as_mut().clone_from_slice(buf_next.as_ref());
        }

        println!(
"n: {}, N: {}, k: {}, l_pbs: {}, B_pbs: 2^{}, l_auto: {}, B_auto: 2^{}, l_ss: {}, B_ss: 2^{}, l_cbs: {}, B_cbs: 2^{}",
            lwe_dimension.0, polynomial_size.0, glwe_dimension.0, pbs_level.0, pbs_base_log.0, auto_level.0, auto_base_log.0, ss_level.0, ss_base_log.0, cbs_level.0, cbs_base_log.0,
        );
        for (bit_idx, ggsw) in ggsw_list_out.iter().enumerate() {
            let extract_bit = (msg & (1 << bit_idx)) >> bit_idx;
            let correct_val = if bit_idx % extract_size == extract_size - 1 {
                extract_bit
            } else {
                let mask_idx = (bit_idx / extract_size) * extract_size + (extract_size - 1);
                let mask_bit = (msg & (1 << mask_idx)) >> mask_idx;
                extract_bit ^ mask_bit
            };
            let err = get_max_err_ggsw_bit(&glwe_sk, ggsw, correct_val);
            println!("[{bit_idx}] {:.3} bits", (err as f64).log2());
        }

    }
}
