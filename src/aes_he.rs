use aligned_vec::ABox;
use tfhe::{
    shortint::prelude::*,
    core_crypto::{
        prelude::{*, CiphertextModulus},
        fft_impl::fft64::{
            c64,
            crypto::{
                bootstrap::FourierLweBootstrapKeyView,
                ggsw::{
                    FourierGgswCiphertextListView,
                    FourierGgswCiphertextListMutView,
                },
            },
        },
    },
};
use std::collections::HashMap;
use crate::{aes_ref::*, convert_lwe_to_glwe_const, ggsw_conv::*, lwe_preprocessing, trace_assign, utils::*, AutomorphKey};

#[inline]
pub fn he_add_round_key<Scalar, StateCont, RkCont>(
    he_state: &mut LweCiphertextList<StateCont>,
    he_round_key: &LweCiphertextList<RkCont>,
) where
    Scalar: UnsignedInteger,
    StateCont: ContainerMut<Element=Scalar>,
    RkCont: Container<Element=Scalar>,
{
    lwe_ciphertext_list_add_assign(he_state, he_round_key.as_view());
}

pub fn he_sub_bytes_by_patched_wwlp_cbs<Scalar, InputCont, OutputCont>(
    he_state_input: &LweCiphertextList<InputCont>,
    he_state_output: &mut LweCiphertextList<OutputCont>,
    fourier_bsk: FourierLweBootstrapKeyView,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: FourierGgswCiphertextListView,
    ggsw_base_log: DecompositionBaseLog,
    ggsw_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
) where
    Scalar: UnsignedTorus + CastInto<usize> + CastFrom<usize>,
    InputCont: Container<Element = Scalar>,
    OutputCont: ContainerMut<Element = Scalar>,
{
    for (input_byte, mut output_byte) in he_state_input.chunks_exact(BYTESIZE)
        .zip(he_state_output.chunks_exact_mut(BYTESIZE))
    {
        he_sbox_eval_by_patched_wwlp_cbs(
            &input_byte,
            &mut output_byte,
            fourier_bsk,
            auto_keys,
            ss_key,
            ggsw_base_log,
            ggsw_level,
            log_lut_count,
        );
    }
}


pub fn known_rotate_keyed_lut<Scalar, AccCont, OutputCont>(
    input_cleartext: [u8; BLOCKSIZE_IN_BYTE],
    vec_keyed_sbox: &Vec<GlweCiphertextList<AccCont>>,
    lwe_state_output: &mut LweCiphertextList<OutputCont>,
) where
    Scalar: UnsignedTorus + CastFrom<usize>,
    AccCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let polynomial_size = vec_keyed_sbox.get(0).unwrap().polynomial_size();

    let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
    for i in 0..BLOCKSIZE_IN_BIT {
        let mut lwe_bit = lwe_state_output.get_mut(i);
        let byte_idx = i / BYTESIZE;
        let bit_idx = i % BYTESIZE;

        let keyed_acc_list = vec_keyed_sbox.get(byte_idx).unwrap();

        let acc_idx = bit_idx / num_par_lut;
        let lut_idx = bit_idx % num_par_lut;
        let deg = lut_idx * (1 << BYTESIZE) + input_cleartext[byte_idx] as usize;

        extract_lwe_sample_from_glwe_ciphertext(
            &keyed_acc_list.get(acc_idx),
            &mut lwe_bit,
            MonomialDegree(deg),
        );
    }
}


pub fn known_rotate_keyed_lut_for_half_cbs<Scalar, AccCont, OutputCont>(
    input_cleartext: [u8; BLOCKSIZE_IN_BYTE],
    vec_keyed_sbox_glev: &Vec<Vec<GlweCiphertextList<AccCont>>>,
    lev_state_output: &mut Vec<LweCiphertextList<OutputCont>>,
) where
    Scalar: UnsignedTorus + CastFrom<usize>,
    AccCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let polynomial_size = vec_keyed_sbox_glev.get(0).unwrap().get(0).unwrap().polynomial_size();

    for (level_minus_one, vec_keyed_sbox_acc_list) in vec_keyed_sbox_glev.iter().enumerate() {
        let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
        for i in 0..BLOCKSIZE_IN_BIT {
            let mut lwe_bit = lev_state_output.get_mut(i).unwrap().get_mut(level_minus_one);
            let byte_idx = i / BYTESIZE;
            let bit_idx = i % BYTESIZE;

            let keyed_acc_list = vec_keyed_sbox_acc_list.get(byte_idx).unwrap();

            let acc_idx = bit_idx / num_par_lut;
            let lut_idx = bit_idx % num_par_lut;
            let deg = lut_idx * (1 << BYTESIZE) + input_cleartext[byte_idx] as usize;

            extract_lwe_sample_from_glwe_ciphertext(
                &keyed_acc_list.get(acc_idx),
                &mut lwe_bit,
                MonomialDegree(deg),
            );
        }
    }
}


pub fn convert_lev_state_to_ggsw<Scalar, InputCont, OutputCont>(
    lev_state: &Vec<LweCiphertextList<InputCont>>,
    ggsw_out: &mut GgswCiphertextList<OutputCont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: FourierGgswCiphertextListView,
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let lev_lwe_size = lev_state.get(0).unwrap().lwe_size();
    let glwe_size = ggsw_out.glwe_size();
    let polynomial_size = ggsw_out.polynomial_size();
    let ciphertext_modulus = ggsw_out.ciphertext_modulus();

    let ggsw_level = ggsw_out.decomposition_level_count();

    assert_eq!(lev_state.get(0).unwrap().lwe_ciphertext_count().0, ggsw_level.0);
    assert_eq!(lev_state.get(0).unwrap().ciphertext_modulus(), ciphertext_modulus);
    assert_eq!(lev_lwe_size.to_lwe_dimension().0, glwe_size.to_glwe_dimension().0 * polynomial_size.0);

    let state_size = lev_state.len();
    assert_eq!(ggsw_out.ggsw_ciphertext_count().0, state_size);

    let mut buf_lwe = LweCiphertext::new(
        Scalar::ZERO,
        lev_lwe_size,
        ciphertext_modulus,
    );
    for (lev, mut ggsw) in lev_state.iter().zip(ggsw_out.iter_mut()) {
        let mut glev = GlweCiphertextList::new(
            Scalar::ZERO,
            glwe_size,
            polynomial_size,
            GlweCiphertextCount(ggsw_level.0),
            ciphertext_modulus,
        );

        for (lwe, mut glwe) in lev.iter().zip(glev.iter_mut()) {
            lwe_preprocessing(&lwe, &mut buf_lwe, polynomial_size);
            convert_lwe_to_glwe_const(&buf_lwe, &mut glwe);
            trace_assign(&mut glwe, auto_keys);
        }

        switch_scheme(&glev, &mut ggsw, ss_key);
    }
}


pub fn blind_rotate_keyed_sboxes<Scalar, InputCont, AccCont, OutputCont>(
    ggsw_list: &GgswCiphertextList<InputCont>,
    vec_keyed_sbox_acc: &Vec<GlweCiphertextList<AccCont>>,
    vec_keyed_sbox_mult_by_2_acc: &Vec<GlweCiphertextList<AccCont>>,
    vec_keyed_sbox_mult_by_3_acc: &Vec<GlweCiphertextList<AccCont>>,
    lwe_state_output: &mut LweCiphertextList<OutputCont>,
    lwe_state_output_mult_by_2: &mut LweCiphertextList<OutputCont>,
    lwe_state_output_mult_by_3: &mut LweCiphertextList<OutputCont>,
) where
    Scalar: UnsignedTorus + CastFrom<usize>,
    InputCont: Container<Element=Scalar>,
    AccCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let glwe_size = ggsw_list.glwe_size();
    let polynomial_size = ggsw_list.polynomial_size();
    let ciphertext_modulus = lwe_state_output.ciphertext_modulus();
    let cbs_base_log = ggsw_list.decomposition_base_log();
    let cbs_level = ggsw_list.decomposition_level_count();

    let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
    for (((ggsw_chunk,
        (keyed_sbox_acc_list, mut output_chunk)),
        (keyed_sbox_mult_by_2_acc_list, mut output_mult_by_2_chunk)),
        (keyed_sbox_mult_by_3_acc_list, mut output_mult_by_3_chunk),
    )
        in ggsw_list.chunks_exact(BYTESIZE)
        .zip(
            vec_keyed_sbox_acc.iter()
            .zip(lwe_state_output.chunks_exact_mut(BYTESIZE))
        )
        .zip(
            vec_keyed_sbox_mult_by_2_acc.iter()
            .zip(lwe_state_output_mult_by_2.chunks_exact_mut(BYTESIZE))
        )
        .zip(
            vec_keyed_sbox_mult_by_3_acc.iter()
            .zip(lwe_state_output_mult_by_3.chunks_exact_mut(BYTESIZE))
        )
    {
        let mut fourier_ggsw_chunk = FourierGgswCiphertextList::new(
            vec![c64::default();
            BYTESIZE * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * cbs_level.0
            ],
            BYTESIZE,
            glwe_size,
            polynomial_size,
            cbs_base_log,
            cbs_level,
        );

        for (mut fourier_ggsw, ggsw) in fourier_ggsw_chunk.as_mut_view().into_ggsw_iter().zip(ggsw_chunk.iter()) {
            convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
        }

        for (acc_idx, ((keyed_acc, keyed_acc_mult_by_2), keyed_acc_mult_by_3))
            in keyed_sbox_acc_list.iter()
            .zip(keyed_sbox_mult_by_2_acc_list.iter())
            .zip(keyed_sbox_mult_by_3_acc_list.iter())
            .enumerate()
        {
            let mut acc = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
            acc.as_mut().clone_from_slice(keyed_acc.as_ref());

            let mut acc_mult_by_2 = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
            acc_mult_by_2.as_mut().clone_from_slice(keyed_acc_mult_by_2.as_ref());

            let mut acc_mult_by_3 = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
            acc_mult_by_3.as_mut().clone_from_slice(keyed_acc_mult_by_3.as_ref());

            for (i, fourier_ggsw_bit) in fourier_ggsw_chunk.as_view().into_ggsw_iter().enumerate() {
                let mut buf = acc.clone();
                glwe_ciphertext_monic_monomial_div_assign(&mut buf, MonomialDegree(1 << i));
                glwe_ciphertext_sub_assign(&mut buf, &acc);
                add_external_product_assign(&mut acc, &fourier_ggsw_bit, &buf);

                let mut buf2 = acc_mult_by_2.clone();
                glwe_ciphertext_monic_monomial_div_assign(&mut buf2, MonomialDegree(1 << i));
                glwe_ciphertext_sub_assign(&mut buf2, &acc_mult_by_2);
                add_external_product_assign(&mut acc_mult_by_2, &fourier_ggsw_bit, &buf2);

                let mut buf3 = acc_mult_by_3.clone();
                glwe_ciphertext_monic_monomial_div_assign(&mut buf3, MonomialDegree(1 << i));
                glwe_ciphertext_sub_assign(&mut buf3, &acc_mult_by_3);
                add_external_product_assign(&mut acc_mult_by_3, &fourier_ggsw_bit, &buf3);
            }

            for i in 0..num_par_lut {
                let bit_idx = acc_idx * num_par_lut + i;
                let mut lwe_out = output_chunk.get_mut(bit_idx);
                extract_lwe_sample_from_glwe_ciphertext(&acc, &mut lwe_out, MonomialDegree(i * (1 << BYTESIZE)));

                let mut lwe_out2 = output_mult_by_2_chunk.get_mut(bit_idx);
                extract_lwe_sample_from_glwe_ciphertext(&acc_mult_by_2, &mut lwe_out2, MonomialDegree(i * (1 << BYTESIZE)));

                let mut lwe_out3 = output_mult_by_3_chunk.get_mut(bit_idx);
                extract_lwe_sample_from_glwe_ciphertext(&acc_mult_by_3, &mut lwe_out3, MonomialDegree(i * (1 << BYTESIZE)));
            }
        }
    }
}


pub fn he_sub_bytes_8_to_24_by_patched_wwlp_cbs<Scalar, InputCont, OutputCont>(
    he_state_input: &LweCiphertextList<InputCont>,
    he_state_output: &mut LweCiphertextList<OutputCont>,
    he_state_output_mult_by_2: &mut LweCiphertextList<OutputCont>,
    he_state_output_mult_by_3: &mut LweCiphertextList<OutputCont>,
    fourier_bsk: FourierLweBootstrapKeyView,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: FourierGgswCiphertextListView,
    ggsw_base_log: DecompositionBaseLog,
    ggsw_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
) where
    Scalar: UnsignedTorus + CastInto<usize> + CastFrom<usize>,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    for (((input, mut output), mut output_mult_by_2), mut output_mult_by_3)
        in he_state_input.chunks_exact(BYTESIZE)
        .zip(he_state_output.chunks_exact_mut(BYTESIZE))
        .zip(he_state_output_mult_by_2.chunks_exact_mut(BYTESIZE))
        .zip(he_state_output_mult_by_3.chunks_exact_mut(BYTESIZE))
    {
        he_sbox_8_to_24_eval_by_patched_wwlp_cbs(
            &input,
            &mut output,
            &mut output_mult_by_2,
            &mut output_mult_by_3,
            fourier_bsk,
            auto_keys,
            ss_key,
            ggsw_base_log,
            ggsw_level,
            log_lut_count,
        );
    }
}

fn get_he_state_byte<Scalar, Cont>(
    he_state: &LweCiphertextList<Cont>,
    row: usize,
    col: usize,
) -> LweCiphertextListView<Scalar>
where
    Scalar: UnsignedInteger,
    Cont: Container<Element=Scalar>,
{
    let byte_idx = 4 * col + row;
    he_state.get_sub((8*byte_idx)..(8*byte_idx + 8))
}

fn get_he_state_byte_mut<Scalar, Cont>(
    he_state: &mut LweCiphertextList<Cont>,
    row: usize,
    col: usize,
) -> LweCiphertextListMutView<Scalar>
where
    Scalar: UnsignedInteger,
    Cont: ContainerMut<Element=Scalar>,
{
    let byte_idx = 4 * col + row;
    he_state.get_sub_mut((8*byte_idx)..(8*byte_idx+8))
}

pub fn he_shift_rows<Scalar, Cont>(he_state: &mut LweCiphertextList<Cont>)
where
    Scalar: UnsignedInteger,
    Cont: ContainerMut<Element=Scalar>,
{
    let mut buf = LweCiphertextList::new(
        Scalar::ZERO,
        he_state.lwe_size(),
        LweCiphertextCount(BLOCKSIZE_IN_BIT),
        he_state.ciphertext_modulus(),
    );
    buf.as_mut().clone_from_slice(he_state.as_ref());

    for row in 1..4 {
        for col in 0..4 {
            let mut dst = get_he_state_byte_mut(he_state, row, col);
            let src = get_he_state_byte(&buf, row, (row + col) % 4);
            dst.as_mut().clone_from_slice(src.as_ref());
        }
    }
}

pub fn lev_shift_rows<Scalar, Cont>(lev_state: &mut Vec<LweCiphertextList<Cont>>)
where
    Scalar: UnsignedInteger,
    Cont: ContainerMut<Element=Scalar>,
{
    let lwe_size = lev_state.get(0).unwrap().lwe_size();
    let lwe_ciphertext_count = lev_state.get(0).unwrap().lwe_ciphertext_count();
    let ciphertext_modulus = lev_state.get(0).unwrap().ciphertext_modulus();

    let mut buf = vec![LweCiphertextList::new(
        Scalar::ZERO,
        lwe_size,
        lwe_ciphertext_count,
        ciphertext_modulus,
    ); BLOCKSIZE_IN_BIT];
    for (buf_bit, lev_bit) in buf.iter_mut().zip(lev_state.iter()) {
        buf_bit.as_mut().clone_from_slice(lev_bit.as_ref());
    }

    for row in 0..NUM_ROWS {
        for col in 0..NUM_COLUMNS {
            for bit_idx in 0..BYTESIZE {
                let byte_idx = 4 * col + row;
                let idx = BYTESIZE * byte_idx + bit_idx;

                let shift_byte_idx = 4 * ((col + row) % 4) + row;
                let shift_idx = BYTESIZE * shift_byte_idx + bit_idx;

                lev_state.get_mut(idx).unwrap().as_mut()
                    .clone_from_slice(buf.get(shift_idx).unwrap().as_ref());
            }
        }
    }
}

pub fn he_mix_columns<Scalar, Cont>(he_state: &mut LweCiphertextList<Cont>)
where
    Scalar: UnsignedInteger,
    Cont: ContainerMut<Element=Scalar>,
{
    let mut buf = LweCiphertextList::new(Scalar::ZERO, he_state.lwe_size(), LweCiphertextCount(BLOCKSIZE_IN_BIT), he_state.ciphertext_modulus());
    buf.as_mut().clone_from_slice(he_state.as_ref());

    for row in 0..4 {
        for col in 0..4 {
            let mut tmp = LweCiphertextList::new(Scalar::ZERO, buf.lwe_size(), LweCiphertextCount(BYTESIZE), buf.ciphertext_modulus());

            tmp.as_mut().clone_from_slice(get_he_state_byte(&buf, row, col).as_ref());
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 1) % 4, col),
            );
            he_mult_by_two_assign(&mut tmp);
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 1) % 4, col),
            );
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 2) % 4, col),
            );
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 3) % 4, col),
            );

            get_he_state_byte_mut(he_state, row, col).as_mut().clone_from_slice(tmp.as_ref());
        }
    }
}

pub fn he_mix_columns_precomp<Scalar, Cont, ContMut>(
    he_state: &mut LweCiphertextList<ContMut>,
    he_state_mult_by_2: &LweCiphertextList<Cont>,
    he_state_mult_by_3: &LweCiphertextList<Cont>,
) where
    Scalar: UnsignedInteger,
    Cont: Container<Element=Scalar>,
    ContMut: ContainerMut<Element=Scalar>,
{
    let mut buf = LweCiphertextList::new(Scalar::ZERO, he_state.lwe_size(), LweCiphertextCount(BLOCKSIZE_IN_BIT), he_state.ciphertext_modulus());
    buf.as_mut().clone_from_slice(he_state.as_ref());

    for row in 0..NUM_ROWS {
        for col in 0..NUM_COLUMNS {
            let mut tmp = LweCiphertextList::new(Scalar::ZERO, he_state.lwe_size(), LweCiphertextCount(BYTESIZE), he_state.ciphertext_modulus());

            tmp.as_mut().clone_from_slice(get_he_state_byte(&he_state_mult_by_2, row, col).as_ref());
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&he_state_mult_by_3, (row + 1) % 4, col),
            );
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 2) % 4, col),
            );
            lwe_ciphertext_list_add_assign(
                &mut tmp,
                get_he_state_byte(&buf, (row + 3) % 4, col),
            );

            get_he_state_byte_mut(he_state, row, col).as_mut().clone_from_slice(tmp.as_ref());
        }
    }
}

pub fn lev_mix_columns_precomp<Scalar, Cont, ContMut>(
    lev_state: &mut Vec<LweCiphertextList<ContMut>>,
    lev_state_mult_by_2: &Vec<LweCiphertextList<Cont>>,
    lev_state_mult_by_3: &Vec<LweCiphertextList<Cont>>,
) where
    Scalar: UnsignedInteger,
    Cont: Container<Element=Scalar>,
    ContMut: ContainerMut<Element=Scalar>,
{
    let lwe_size = lev_state.get(0).unwrap().lwe_size();
    let lwe_ciphertext_count = lev_state.get(0).unwrap().lwe_ciphertext_count();
    let ciphertext_modulus = lev_state.get(0).unwrap().ciphertext_modulus();

    let mut buf = vec![LweCiphertextList::new(
        Scalar::ZERO,
        lwe_size,
        lwe_ciphertext_count,
        ciphertext_modulus,
    ); BLOCKSIZE_IN_BIT];
    for (buf_bit, lev_bit) in buf.iter_mut().zip(lev_state.iter()) {
        buf_bit.as_mut().clone_from_slice(lev_bit.as_ref());
    }

    for row in 0..NUM_ROWS {
        for col in 0..NUM_COLUMNS {
            for bit_idx in 0..BYTESIZE {
                let mut tmp = LweCiphertextList::new(
                    Scalar::ZERO,
                    lwe_size,
                    lwe_ciphertext_count,
                    ciphertext_modulus,
                );

                let byte_idx = 4 * col + row;
                let idx = BYTESIZE * byte_idx + bit_idx;
                tmp.as_mut().clone_from_slice(lev_state_mult_by_2.get(idx).unwrap().as_ref());

                let byte_idx = 4 * col + (row + 1) % 4;
                let idx = BYTESIZE * byte_idx + bit_idx;
                lwe_ciphertext_list_add_assign(
                    &mut tmp,
                    lev_state_mult_by_3.get(idx).unwrap().as_view(),
                );

                let byte_idx = 4 * col + (row + 2) % 4;
                let idx = BYTESIZE * byte_idx + bit_idx;
                lwe_ciphertext_list_add_assign(
                    &mut tmp,
                    buf.get(idx).unwrap().as_view(),
                );

                let byte_idx = 4 * col + (row + 3) % 4;
                let idx = BYTESIZE * byte_idx + bit_idx;
                lwe_ciphertext_list_add_assign(
                    &mut tmp,
                    buf.get(idx).unwrap().as_view(),
                );

                let byte_idx = 4 * col + row;
                let idx = BYTESIZE * byte_idx + bit_idx;
                lev_state.get_mut(idx).unwrap().as_mut().clone_from_slice(tmp.as_ref());
            }
        }
    }
}

fn he_mult_by_two<Scalar, Cont>(he_byte: &LweCiphertextList<Cont>) -> LweCiphertextListOwned<Scalar>
where
    Scalar: UnsignedInteger,
    Cont: Container<Element=Scalar>,
{
    debug_assert!(he_byte.entity_count() == BYTESIZE);

    // 2 * (a7, a6, …, a0)
    // = (a6, a5, …, a0, 0) + (0, 0, 0, a7, a7, 0, a7, a7)
    // = (a6, a5, a4, a3 + a7, a2 + a7, a1, a0 + a7, a7)
    let mut output = LweCiphertextList::new(
        Scalar::ZERO,
        he_byte.lwe_size(),
        LweCiphertextCount(BYTESIZE),
        he_byte.ciphertext_modulus(),
    );

    for i in 1..BYTESIZE {
        output.get_mut(i).as_mut().clone_from_slice(he_byte.get(i-1).as_ref());
    }

    let he_msb = he_byte.get(BYTESIZE-1);
    for i in [0, 1, 3, 4] {
        lwe_ciphertext_add_assign(&mut output.get_mut(i), &he_msb);
    }

    output
}

fn he_mult_by_two_assign<Scalar, Cont>(he_byte: &mut LweCiphertextList<Cont>)
where
    Scalar: UnsignedInteger,
    Cont: ContainerMut<Element=Scalar>,
{
    let buf = he_mult_by_two(&he_byte);
    he_byte.as_mut().clone_from_slice(buf.as_ref());
}

fn he_sbox_eval_by_patched_wwlp_cbs<Scalar, InCont, OutCont>(
    input: &LweCiphertextList<InCont>,
    output: &mut LweCiphertextList<OutCont>,
    fourier_bsk: FourierLweBootstrapKeyView,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: FourierGgswCiphertextListView,
    ggsw_base_log: DecompositionBaseLog,
    ggsw_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
) where
    Scalar: UnsignedTorus + CastInto<usize> + CastFrom<usize>,
    InCont: Container<Element=Scalar>,
    OutCont: ContainerMut<Element=Scalar>,
{
    let glwe_size = fourier_bsk.glwe_size();
    let polynomial_size = fourier_bsk.polynomial_size();
    let ciphertext_modulus = output.ciphertext_modulus();

    let mut vec_glev = vec![
        GlweCiphertextList::new(
            Scalar::ZERO,
            glwe_size,
            polynomial_size,
            GlweCiphertextCount(ggsw_level.0),
            ciphertext_modulus,
        ); BYTESIZE];
    for (input_bit, glev) in input.iter().zip(vec_glev.iter_mut()) {
        let glev_mut_view = GlweCiphertextListMutView::from_container(
            glev.as_mut(),
            glwe_size,
            polynomial_size,
            ciphertext_modulus,
        );

        lwe_msb_bit_to_glev_by_trace_with_preprocessing(
            input_bit.as_view(),
            glev_mut_view,
            fourier_bsk,
            auto_keys,
            ggsw_base_log,
            ggsw_level,
            log_lut_count,
        );
    }

    let mut ggsw_bit_list = GgswCiphertextList::new(
        Scalar::ZERO,
        glwe_size,
        polynomial_size,
        ggsw_base_log,
        ggsw_level,
        GgswCiphertextCount(vec_glev.len()),
        ciphertext_modulus,
    );
    for (mut ggsw, glev) in ggsw_bit_list.iter_mut().zip(vec_glev.iter()) {
        switch_scheme(&glev, &mut ggsw, ss_key);
    }

    let mut fourier_ggsw_bit_list = FourierGgswCiphertextList::new(
        vec![c64::default();
        BYTESIZE * polynomial_size.to_fourier_polynomial_size().0
            * glwe_size.0
            * glwe_size.0
            * ggsw_level.0
        ],
        BYTESIZE,
        glwe_size,
        polynomial_size,
        ggsw_base_log,
        ggsw_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_bit_list.as_mut_view().into_ggsw_iter().zip(ggsw_bit_list.iter()) {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }

    let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
    let num_accumulator = if BYTESIZE % num_par_lut == 0 {
        BYTESIZE / num_par_lut
    } else {
        BYTESIZE / num_par_lut + 1
    };

    for acc_idx in 0..num_accumulator {
        let accumulator = (0..polynomial_size.0).map(|i| {
            let lut_idx = acc_idx * num_par_lut + i / (1 << BYTESIZE);
            (((AES128_SBOX[i % (1 << BYTESIZE)] & (1 << lut_idx)) as usize) << ((Scalar::BITS - 1) - lut_idx)).cast_into()
        }).collect::<Vec<Scalar>>();
        let accumulator_plaintext = PlaintextList::from_container(accumulator);
        let mut accumulator = allocate_and_trivially_encrypt_new_glwe_ciphertext(glwe_size, &accumulator_plaintext, ciphertext_modulus);

        for (i, fourier_ggsw_bit) in fourier_ggsw_bit_list.as_view().into_ggsw_iter().into_iter().enumerate() {
            let mut buf = accumulator.clone();
            glwe_ciphertext_monic_monomial_div_assign(&mut buf, MonomialDegree(1 << i));
            glwe_ciphertext_sub_assign(&mut buf, &accumulator);
            add_external_product_assign(&mut accumulator, &fourier_ggsw_bit, &buf);
        }

        for i in 0..num_par_lut {
            let bit_idx = acc_idx * num_par_lut + i;
            let mut lwe_out = output.get_mut(bit_idx);
            extract_lwe_sample_from_glwe_ciphertext(&accumulator, &mut lwe_out, MonomialDegree(i * (1 << BYTESIZE)));
        }
    }
}


pub fn he_sbox_8_to_24_eval_by_patched_wwlp_cbs<Scalar, InputCont, OutputCont>(
    input: &LweCiphertextList<InputCont>,
    output: &mut LweCiphertextList<OutputCont>,
    output_mult_by_2: &mut LweCiphertextList<OutputCont>,
    output_mult_by_3: &mut LweCiphertextList<OutputCont>,
    fourier_bsk: FourierLweBootstrapKeyView,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    ss_key: FourierGgswCiphertextListView,
    ggsw_base_log: DecompositionBaseLog,
    ggsw_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
) where
    Scalar: UnsignedTorus + CastInto<usize> + CastFrom<usize>,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let glwe_size = fourier_bsk.glwe_size();
    let polynomial_size = fourier_bsk.polynomial_size();
    let ciphertext_modulus = output.ciphertext_modulus();

    let mut vec_glev = vec![
        GlweCiphertextList::new(
            Scalar::ZERO,
            glwe_size,
            polynomial_size,
            GlweCiphertextCount(ggsw_level.0),
            ciphertext_modulus,
        ); BYTESIZE];
    for (input_bit, glev) in input.iter().zip(vec_glev.iter_mut()) {
        let glev_mut_view = GlweCiphertextListMutView::from_container(
            glev.as_mut(),
            glwe_size,
            polynomial_size,
            ciphertext_modulus,
        );

        lwe_msb_bit_to_glev_by_trace_with_preprocessing(
            input_bit.as_view(),
            glev_mut_view,
            fourier_bsk,
            auto_keys,
            ggsw_base_log,
            ggsw_level,
            log_lut_count,
        );
    }

    let mut ggsw_bit_list = GgswCiphertextList::new(
        Scalar::ZERO,
        glwe_size,
        polynomial_size,
        ggsw_base_log,
        ggsw_level,
        GgswCiphertextCount(vec_glev.len()),
        ciphertext_modulus,
    );
    for (mut ggsw, glev) in ggsw_bit_list.iter_mut().zip(vec_glev.iter()) {
        switch_scheme(&glev, &mut ggsw, ss_key);
    }

    let mut fourier_ggsw_bit_list = FourierGgswCiphertextList::new(
        vec![c64::default();
            BYTESIZE * polynomial_size.to_fourier_polynomial_size().0
                * glwe_size.0
                * glwe_size.0
                * ggsw_level.0
        ],
        BYTESIZE,
        glwe_size,
        polynomial_size,
        ggsw_base_log,
        ggsw_level,
    );
    for (mut fourier_ggsw, ggsw) in fourier_ggsw_bit_list.as_mut_view().into_ggsw_iter().zip(ggsw_bit_list.iter()) {
        convert_standard_ggsw_ciphertext_to_fourier(&ggsw, &mut fourier_ggsw);
    }

    evaluate_8_to_8_lut(
        fourier_ggsw_bit_list.as_mut_view(),
        output,
        &AES128_SBOX,
        Scalar::BITS - 1,
        ciphertext_modulus,
    );

    evaluate_8_to_8_lut(
        fourier_ggsw_bit_list.as_mut_view(),
        output_mult_by_2,
        &AES128_SBOX_MULT_BY_2,
        Scalar::BITS - 1,
        ciphertext_modulus,
    );

    evaluate_8_to_8_lut(
        fourier_ggsw_bit_list.as_mut_view(),
        output_mult_by_3,
        &AES128_SBOX_MULT_BY_3,
        Scalar::BITS - 1,
        ciphertext_modulus,
    );
}


fn evaluate_8_to_8_lut<Scalar, ContMut>(
    fourier_ggsw_bit_list: FourierGgswCiphertextListMutView,
    output: &mut LweCiphertextList<ContMut>,
    lut: &[u8; 256],
    log_scale: usize,
    ciphertext_modulus: CiphertextModulus<Scalar>,
) where
    Scalar: UnsignedTorus + CastFrom<usize>,
    ContMut: ContainerMut<Element=Scalar>,
{
    let glwe_size = fourier_ggsw_bit_list.glwe_size();
    let polynomial_size = fourier_ggsw_bit_list.polynomial_size();
    let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
    let num_accumulator = if BYTESIZE % num_par_lut == 0 {
        BYTESIZE / num_par_lut
    } else {
        BYTESIZE / num_par_lut + 1
    };

    for acc_idx in 0..num_accumulator {
        let accumulator = (0..polynomial_size.0).map(|i| {
            let lut_idx = acc_idx * num_par_lut + i / (1 << BYTESIZE);
            (((lut[i % (1 << BYTESIZE)] & (1 << lut_idx)) as usize) << (log_scale - lut_idx)).cast_into()
        }).collect::<Vec<Scalar>>();
        let accumulator_plaintext = PlaintextList::from_container(accumulator);
        let mut accumulator = allocate_and_trivially_encrypt_new_glwe_ciphertext(glwe_size, &accumulator_plaintext, ciphertext_modulus);

        for (i, fourier_ggsw_bit) in fourier_ggsw_bit_list.as_view().into_ggsw_iter().into_iter().enumerate() {
            let mut buf = accumulator.clone();
            glwe_ciphertext_monic_monomial_div_assign(&mut buf, MonomialDegree(1 << i));
            glwe_ciphertext_sub_assign(&mut buf, &accumulator);
            add_external_product_assign(&mut accumulator, &fourier_ggsw_bit, &buf);
        }

        for i in 0..num_par_lut {
            let bit_idx = acc_idx * num_par_lut + i;
            let mut lwe_out = output.get_mut(bit_idx);
            extract_lwe_sample_from_glwe_ciphertext(&accumulator, &mut lwe_out, MonomialDegree(i * (1 << BYTESIZE)));
        }
    }
}


pub fn generate_vec_keyed_lut_accumulator<Scalar, KeyCont, G>(
    keyed_lut_list: [[u8; 1 << BYTESIZE]; BLOCKSIZE_IN_BYTE],
    log_scale: usize,
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    glwe_modular_std_dev: impl DispersionParameter,
    ciphertext_modulus: CiphertextModulus<Scalar>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> Vec::<GlweCiphertextListOwned<Scalar>>
where
    Scalar: UnsignedTorus + CastFrom<usize>,
    KeyCont: Container<Element=Scalar>,
    G: ByteRandomGenerator,
{
    let glwe_size = glwe_secret_key.glwe_dimension().to_glwe_size();
    let polynomial_size = glwe_secret_key.polynomial_size();

    let mut vec_keyed_lut_acc = Vec::<GlweCiphertextListOwned<Scalar>>::with_capacity(BLOCKSIZE_IN_BIT);

    let num_par_lut = polynomial_size.0 / (1 << BYTESIZE);
    let num_acc_per_byte = if BYTESIZE % num_par_lut == 0 {
        BYTESIZE / num_par_lut
    } else {
        BYTESIZE / num_par_lut + 1
    };

    for byte_idx in 0..BLOCKSIZE_IN_BYTE {
        let mut keyed_lut_acc_list = GlweCiphertextList::new(Scalar::ZERO, glwe_size, polynomial_size, GlweCiphertextCount(num_acc_per_byte), ciphertext_modulus);

        let keyed_lut = keyed_lut_list[byte_idx];
        for (acc_idx, mut keyed_lut_acc) in keyed_lut_acc_list.iter_mut().enumerate() {
            let acc = (0..polynomial_size.0).map(|i| {
                let lut_idx = acc_idx * num_par_lut + i / (1 << BYTESIZE);
                (((keyed_lut[i % (1 << BYTESIZE)] & (1 << lut_idx)) as usize) << (log_scale - lut_idx)).cast_into()
            }).collect::<Vec<Scalar>>();
            let acc = PlaintextList::from_container(acc);
            encrypt_glwe_ciphertext(glwe_secret_key, &mut keyed_lut_acc, &acc, glwe_modular_std_dev, encryption_generator);
        }
        vec_keyed_lut_acc.push(keyed_lut_acc_list);
    }

    vec_keyed_lut_acc
}


pub fn generate_vec_keyed_lut_glev<Scalar, KeyCont, G>(
    keyed_lut_list: [[u8; 1 << BYTESIZE]; BLOCKSIZE_IN_BYTE],
    glev_base_log: DecompositionBaseLog,
    glev_level: DecompositionLevelCount,
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    glwe_modular_std_dev: impl DispersionParameter,
    ciphertext_modulus: CiphertextModulus<Scalar>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> Vec::<Vec<GlweCiphertextListOwned<Scalar>>>
where
    Scalar: UnsignedTorus + CastFrom<usize>,
    KeyCont: Container<Element=Scalar>,
    G: ByteRandomGenerator,
{
    let mut vec_keyed_lut_glev = Vec::<Vec::<GlweCiphertextListOwned<Scalar>>>::with_capacity(glev_level.0);

    for level_minus_one in 0..glev_level.0 {
        let level = level_minus_one + 1;
        let log_scale = Scalar::BITS - level * glev_base_log.0;

        let vec_keyed_lut_acc_list = generate_vec_keyed_lut_accumulator(
            keyed_lut_list,
            log_scale,
            glwe_secret_key,
            glwe_modular_std_dev,
            ciphertext_modulus,
            encryption_generator,
        );

        vec_keyed_lut_glev.push(vec_keyed_lut_acc_list);
    }

    vec_keyed_lut_glev
}


pub fn get_he_state_error<Scalar, StateCont, SkCont>(
    he_state: &LweCiphertextList<StateCont>,
    plain_state: StateByteMat,
    lwe_sk: &LweSecretKey<SkCont>,
) -> (Vec::<usize>, Scalar)
where
    Scalar: UnsignedInteger + CastFrom<u8>,
    StateCont: Container<Element=Scalar>,
    SkCont: Container<Element=Scalar>,
{
    let plain_state = byte_mat_to_bit_array(plain_state);
    let mut vec_err = Vec::<usize>::with_capacity(BLOCKSIZE_IN_BIT);
    let mut max_err = Scalar::ZERO;

    for (correct_val, he_bit) in plain_state.iter().zip(he_state.iter()) {
        let correct_val = Scalar::cast_from(*correct_val);
        let (_decoded, bit_err, abs_err) = get_val_and_bit_and_abs_err(
            lwe_sk,
            &he_bit,
            correct_val,
            Scalar::ONE << (Scalar::BITS - 1),
        );
        vec_err.push(bit_err as usize);
        max_err = std::cmp::max(max_err, abs_err);
    }

    (vec_err, max_err)
}


pub fn get_lev_int_state_error<Scalar, LevCont, KeyCont>(
    lev_state: &Vec<LweCiphertextList<LevCont>>,
    int_state: [u8; BLOCKSIZE_IN_BIT],
    lev_base_log: DecompositionBaseLog,
    lwe_sk: &LweSecretKey<KeyCont>,
) -> (Vec::<Scalar>, Scalar)
where
    Scalar: UnsignedInteger + CastFrom<u8>,
    LevCont: Container<Element=Scalar>,
    KeyCont: Container<Element=Scalar>,
{
    let mut vec_err = Vec::<Scalar>::with_capacity(BLOCKSIZE_IN_BIT);
    let mut max_err = Scalar::ZERO;

    for (lev_bit, correct_val) in lev_state.iter().zip(int_state.iter()) {
        let mut cur_max_err = Scalar::ZERO;
        for (level_minus_one, lwe_bit) in lev_bit.iter().enumerate() {
            let level = level_minus_one + 1;
            let log_scale = Scalar::BITS - level * lev_base_log.0;

            let (_, err) = get_val_and_abs_err(lwe_sk, &lwe_bit, (*correct_val).cast_into(), Scalar::ONE << log_scale);
            cur_max_err = std::cmp::max(cur_max_err, err);
        }

        vec_err.push(cur_max_err);
        max_err = std::cmp::max(max_err, cur_max_err);
    }

    (vec_err, max_err)
}


pub fn get_he_state_and_error<Scalar, StateCont, SkCont>(
    he_state: &LweCiphertextList<StateCont>,
    plain_state: StateByteMat,
    lwe_sk: &LweSecretKey<SkCont>,
) -> (Vec::<Scalar>, Vec::<usize>, Scalar)
where
    Scalar: UnsignedInteger + CastFrom<u8>,
    StateCont: Container<Element=Scalar>,
    SkCont: Container<Element=Scalar>,
{
    let plain_state = byte_mat_to_bit_array(plain_state);
    let mut vec_out = Vec::<Scalar>::with_capacity(BLOCKSIZE_IN_BIT);
    let mut vec_err = Vec::<usize>::with_capacity(BLOCKSIZE_IN_BIT);
    let mut max_err = Scalar::ZERO;

    for (correct_val, he_bit) in plain_state.iter().zip(he_state.iter()) {
        let correct_val = Scalar::cast_from(*correct_val);
        let (decoded, bit_err, abs_err) = get_val_and_bit_and_abs_err(
            lwe_sk,
            &he_bit,
            correct_val,
            Scalar::ONE << (Scalar::BITS - 1),
        );
        vec_out.push(decoded);
        vec_err.push(bit_err as usize);
        max_err = std::cmp::max(max_err, abs_err);
    }

    (vec_out, vec_err, max_err)
}
