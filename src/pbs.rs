use tfhe::core_crypto::{
    algorithms::polynomial_algorithms::*,
    fft_impl::{
        common::fast_pbs_modulus_switch,
        fft64::{
            crypto::bootstrap::FourierLweBootstrapKeyView,
            math::fft::FftView,
        },
    },
    prelude::*
};
use aligned_vec::CACHELINE_ALIGN;
use dyn_stack::{PodStack, ReborrowMut};
use crate::utils::*;

pub fn generate_accumulator<Scalar, F>(
    polynomial_size: PolynomialSize,
    glwe_size: GlweSize,
    message_modulus: usize,
    ciphertext_modulus: CiphertextModulus<Scalar>,
    delta: Scalar,
    f: F,
) -> GlweCiphertextOwned<Scalar>
where
    Scalar: UnsignedTorus + CastFrom<usize>,
    F: Fn(Scalar) -> Scalar,
{
    // N/(p/2) = size of each block, to correct noise from the input we introduce the
    // notion of box, which manages redundancy to yield a denoised value
    // for several noisy values around a true input value.
    let box_size = polynomial_size.0 / message_modulus;

    // Create the accumulator
    let mut accumulator_scalar = vec![Scalar::ZERO; polynomial_size.0];

    // Fill each box with the encoded denoised value
    for i in 0..message_modulus {
        let index = i * box_size;
        accumulator_scalar[index..index + box_size]
            .iter_mut()
            .for_each(|a| *a = f(Scalar::cast_from(i)) * delta);
    }

    let half_box_size = box_size / 2;

    if ciphertext_modulus.is_compatible_with_native_modulus() {
        // Negate the first half_box_size coefficients to manage negacyclicity and rotate
        for a_i in accumulator_scalar[0..half_box_size].iter_mut() {
            *a_i = (*a_i).wrapping_neg();
        }
    } else {
        let modulus: Scalar = ciphertext_modulus.get_custom_modulus().cast_into();
        for a_i in accumulator_scalar[0..half_box_size].iter_mut() {
            *a_i = (*a_i).wrapping_neg_custom_mod(modulus);
        }
    }

    // Rotate the accumulator
    accumulator_scalar.rotate_left(half_box_size);

    let accumulator_plaintext = PlaintextList::from_container(accumulator_scalar);

    allocate_and_trivially_encrypt_new_glwe_ciphertext(
        glwe_size,
        &accumulator_plaintext,
        ciphertext_modulus,
    )
}

pub fn gen_blind_rotate_local_assign<Scalar: UnsignedTorus + CastInto<usize>>(
    bsk: FourierLweBootstrapKeyView<'_>,
    mut lut: GlweCiphertextMutView<'_, Scalar>,
    mod_switch_offset: ModulusSwitchOffset,
    log_lut_count: LutCountLog,
    lwe: &[Scalar],
    fft: FftView<'_>,
    mut stack: PodStack<'_>,
) {
    let (lwe_body, lwe_mask) = lwe.split_last().unwrap();

    let lut_poly_size = lut.polynomial_size();
    let ciphertext_modulus = lut.ciphertext_modulus();
    assert!(ciphertext_modulus.is_compatible_with_native_modulus());
    let monomial_degree = MonomialDegree(fast_pbs_modulus_switch(
        *lwe_body,
        lut_poly_size,
        mod_switch_offset,
        log_lut_count,
    ));

    lut.as_mut_polynomial_list()
        .iter_mut()
        .for_each(|mut poly| {
            let (mut tmp_poly, _) = stack
                .rb_mut()
                .make_aligned_raw(poly.as_ref().len(), CACHELINE_ALIGN);

            let mut tmp_poly = Polynomial::from_container(&mut *tmp_poly);
            tmp_poly.as_mut().copy_from_slice(poly.as_ref());
            polynomial_wrapping_monic_monomial_div(&mut poly, &tmp_poly, monomial_degree);
        });

    // We initialize the ct_0 used for the successive cmuxes
    let mut ct0 = lut;
    let (mut ct1, mut stack) = stack.make_aligned_raw(ct0.as_ref().len(), CACHELINE_ALIGN);
    let mut ct1 =
        GlweCiphertextMutView::from_container(&mut *ct1, lut_poly_size, ciphertext_modulus);

    for (lwe_mask_element, bootstrap_key_ggsw) in izip!(lwe_mask.iter(), bsk.into_ggsw_iter())
    {
        if *lwe_mask_element != Scalar::ZERO {
            let monomial_degree = MonomialDegree(fast_pbs_modulus_switch(
                *lwe_mask_element,
                lut_poly_size,
                mod_switch_offset,
                log_lut_count,
            ));

            // we effectively inline the body of cmux here, merging the initial subtraction
            // operation with the monic polynomial multiplication, then performing the external
            // product manually

            // We rotate ct_1 and subtract ct_0 (first step of cmux) by performing
            // ct_1 <- (ct_0 * X^{a_hat}) - ct_0
            for (mut ct1_poly, ct0_poly) in izip!(
                ct1.as_mut_polynomial_list().iter_mut(),
                ct0.as_polynomial_list().iter(),
            ) {
                polynomial_wrapping_monic_monomial_mul_and_subtract(
                    &mut ct1_poly,
                    &ct0_poly,
                    monomial_degree,
                );
            }

            // as_mut_view is required to keep borrow rules consistent
            // second step of cmux
            tfhe::core_crypto::fft_impl::fft64::crypto::ggsw::add_external_product_assign(
                ct0.as_mut_view(),
                bootstrap_key_ggsw,
                ct1.as_view(),
                fft,
                stack.rb_mut(),
            );
        }
    }

    if !ciphertext_modulus.is_native_modulus() {
        // When we convert back from the fourier domain, integer values will contain up to 53
        // MSBs with information. In our representation of power of 2 moduli < native modulus we
        // fill the MSBs and leave the LSBs empty, this usage of the signed decomposer allows to
        // round while keeping the data in the MSBs
        let signed_decomposer = SignedDecomposer::new(
            DecompositionBaseLog(ciphertext_modulus.get_custom_modulus().ilog2() as usize),
            DecompositionLevelCount(1),
        );
        ct0.as_mut()
            .iter_mut()
            .for_each(|x| *x = signed_decomposer.closest_representable(*x));
    }
}

pub fn lwe_msb_bit_refresh<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
    refresh_bsk: FourierLweBootstrapKeyView,
) where
    Scalar: UnsignedTorus + CastInto<usize>,
    InputCont: Container<Element = Scalar>,
    OutputCont: ContainerMut<Element = Scalar>,
{
    assert_eq!(input.lwe_size(), refresh_bsk.input_lwe_dimension().to_lwe_size());
    assert_eq!(output.lwe_size(), refresh_bsk.output_lwe_dimension().to_lwe_size());
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());

    let glwe_size = refresh_bsk.glwe_size();
    let polynomial_size = refresh_bsk.polynomial_size();
    let half_box_size = polynomial_size.0 / 2;
    let ciphertext_modulus = input.ciphertext_modulus();

    let mut accumulator = (0..polynomial_size.0).map(|_| {
        (Scalar::ONE).wrapping_neg() << (Scalar::BITS - 2)
    }).collect::<Vec<Scalar>>();

    for a_i in accumulator[0..half_box_size].iter_mut() {
        *a_i = (*a_i).wrapping_neg();
    }
    accumulator.rotate_left(half_box_size);

    let accumulator_plaintext = PlaintextList::from_container(accumulator);
    let accumulator = allocate_and_trivially_encrypt_new_glwe_ciphertext(
        glwe_size,
        &accumulator_plaintext,
        ciphertext_modulus,
    );

    programmable_bootstrap_lwe_ciphertext(input, output, &accumulator, &refresh_bsk);
    lwe_ciphertext_plaintext_add_assign(output, Plaintext(Scalar::ONE << (Scalar::BITS - 2)));
}

pub fn lwe_msb_bit_to_lev<Scalar, InputCont, OutputCont>(
    lwe: &LweCiphertext<InputCont>,
    lev: &mut LweCiphertextList<OutputCont>,
    fourier_bsk: FourierLweBootstrapKeyView,
    lev_base_log: DecompositionBaseLog,
    lev_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
) where
    Scalar: UnsignedTorus + CastInto<usize>,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(lwe.lwe_size(), fourier_bsk.input_lwe_dimension().to_lwe_size());
    assert_eq!(lev.entity_count(), lev_level.0);
    assert_eq!(lwe.ciphertext_modulus(), lev.ciphertext_modulus());

    let glwe_size = fourier_bsk.glwe_size();
    let polynomial_size = fourier_bsk.polynomial_size();
    let half_box_size = polynomial_size.0 / 2;
    let ciphertext_modulus = lwe.ciphertext_modulus();

    let lut_count = 1 << log_lut_count.0;
    for (acc_idx, mut lev_chunk) in lev.chunks_mut(lut_count).enumerate() {
        let mut accumulator = (0..polynomial_size.0).map(|i| {
            let k = i % lut_count;
            let log_scale = Scalar::BITS - (acc_idx * lut_count + k + 1) * lev_base_log.0;
            (Scalar::ONE).wrapping_neg() << (log_scale - 1)
        }).collect::<Vec<Scalar>>();

        for a_i in accumulator[0..half_box_size].iter_mut() {
            *a_i = (*a_i).wrapping_neg();
        }
        accumulator.rotate_left(half_box_size);

        let accumulator_plaintext = PlaintextList::from_container(accumulator);
        let accumulator = allocate_and_trivially_encrypt_new_glwe_ciphertext(
            glwe_size,
            &accumulator_plaintext,
            ciphertext_modulus,
        );

        let mut buffers = ComputationBuffers::new();
        let fft = Fft::new(polynomial_size);
        let fft = fft.as_view();

        buffers.resize(
            programmable_bootstrap_lwe_ciphertext_mem_optimized_requirement::<Scalar>(
                glwe_size,
                polynomial_size,
                fft,
            )
            .unwrap()
            .unaligned_bytes_required(),
        );
        let stack = buffers.stack();

        let (mut local_accumulator_data, stack) = stack.collect_aligned(CACHELINE_ALIGN, accumulator.as_ref().iter().copied());
        let mut local_accumulator = GlweCiphertextMutView::from_container(
            &mut *local_accumulator_data,
            polynomial_size,
            ciphertext_modulus,
        );

        gen_blind_rotate_local_assign(
            fourier_bsk,
            local_accumulator.as_mut_view(),
            ModulusSwitchOffset(0),
            log_lut_count,
            lwe.as_ref(),
            fft,
            stack,
        );

        for (k, mut lwe_out) in lev_chunk.iter_mut().enumerate() {
            let cur_level = acc_idx * lut_count + k + 1;
            let log_scale = Scalar::BITS - cur_level * lev_base_log.0;

            let mut buf = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
            glwe_ciphertext_clone_from(&mut buf, &local_accumulator);
            glwe_ciphertext_monic_monomial_div_assign(&mut buf, MonomialDegree(k));
            glwe_ciphertext_plaintext_add_assign(&mut buf, Plaintext(Scalar::ONE << (log_scale - 1)));

            extract_lwe_sample_from_glwe_ciphertext(&buf, &mut lwe_out, MonomialDegree(0));
        }
    }
}
