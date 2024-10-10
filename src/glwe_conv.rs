use std::collections::HashMap;
use aligned_vec::ABox;
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::c64,
    algorithms::slice_algorithms::slice_wrapping_opposite_assign,
};
use crate::{
    automorphism::*, keyswitch_glwe_ciphertext, mod_switch::*, utils::*, FourierGlweKeyswitchKey
};

pub fn convert_lwe_to_glwe_const<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    let lwe_dimension = input.lwe_size().to_lwe_dimension().0;
    let glwe_dimension = output.glwe_size().to_glwe_dimension().0;
    let polynomial_size = output.polynomial_size().0;

    assert_eq!(lwe_dimension, glwe_dimension * polynomial_size);
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());

    let (lwe_mask, lwe_body) = input.get_mask_and_body();
    let (mut glwe_mask, mut glwe_body) = output.get_mut_mask_and_body();

    // Set body
    *glwe_body.as_mut().get_mut(0).unwrap() = *lwe_body.data;

    // Set mask
    let lwe_mask = lwe_mask.as_ref();
    let glwe_mask = glwe_mask.as_mut();
    for (glwe_poly, lwe_poly) in glwe_mask.chunks_exact_mut(polynomial_size)
        .zip(lwe_mask.chunks_exact(polynomial_size))
    {
        glwe_poly.clone_from_slice(lwe_poly);
        glwe_poly.reverse();
        slice_wrapping_opposite_assign(&mut glwe_poly[0..(polynomial_size - 1)]);
        glwe_poly.rotate_left(polynomial_size - 1);
    }
}


pub fn convert_lwe_to_glwe_by_trace_with_preprocessing<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "only native ciphertext modulus is supported"
    );

    let lwe_size = input.lwe_size();
    let lwe_dimension = lwe_size.to_lwe_dimension();
    let glwe_size = output.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let polynomial_size = output.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    assert_eq!(lwe_dimension.0, glwe_dimension.0 * polynomial_size.0);

    // Pre-processing
    let mut buf = LweCiphertext::new(Scalar::ZERO, lwe_size, ciphertext_modulus);
    lwe_preprocessing(input, &mut buf, polynomial_size);

    // LWEtoGLWEConst
    convert_lwe_to_glwe_const(&buf, output);

    // Clear coefficients except the constant
    trace_assign(output, auto_keys);
}


pub fn convert_lwe_to_glwe_by_trace_with_preprocessing_high_prec<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
    glwe_ksk_to_large: &FourierGlweKeyswitchKey<ABox<[c64]>>,
    glwe_ksk_from_large: &FourierGlweKeyswitchKey<ABox<[c64]>>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "only native ciphertext modulus is supported"
    );
    assert_eq!(glwe_ksk_to_large.input_glwe_size(), glwe_ksk_from_large.output_glwe_size());
    assert_eq!(glwe_ksk_to_large.output_glwe_size(), glwe_ksk_from_large.input_glwe_size());

    let lwe_size = input.lwe_size();
    let lwe_dimension = lwe_size.to_lwe_dimension();
    let glwe_size = output.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let large_glwe_size = glwe_ksk_to_large.output_glwe_size();
    let polynomial_size = output.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    assert_eq!(lwe_dimension.0, glwe_dimension.0 * polynomial_size.0);

    // LWEtoGLWEConst
    let mut buf = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
    convert_lwe_to_glwe_const(&input, &mut buf);

    // GLWE KS to Large
    let mut buf_large = GlweCiphertext::new(Scalar::ZERO, large_glwe_size, polynomial_size, ciphertext_modulus);
    keyswitch_glwe_ciphertext(glwe_ksk_to_large, &buf, &mut buf_large);

    // Pre-processing
    glwe_preprocessing_assign(&mut buf_large);

    // Clear coefficients except the constant
    trace_assign(&mut buf_large, auto_keys);

    // GLWE KS from Large
    keyswitch_glwe_ciphertext(glwe_ksk_from_large, &buf_large, output);
}


pub fn convert_lwes_to_glwe_by_trace_with_preprocessing<Scalar, InputCont, OutputCont>(
    input: &LweCiphertextList<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "only native ciphertext modulus is supported"
    );

    let lwe_size = input.lwe_size();
    let lwe_dimension = lwe_size.to_lwe_dimension();
    let glwe_size = output.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let polynomial_size = output.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    assert_eq!(lwe_dimension.0, glwe_dimension.0 * polynomial_size.0);

    let lwe_count = input.lwe_ciphertext_count().0;
    let mut buf = LweCiphertext::new(Scalar::ZERO, lwe_size, ciphertext_modulus);
    let mut input_glwes = GlweCiphertextList::new(Scalar::ZERO, glwe_size, polynomial_size, GlweCiphertextCount(lwe_count), ciphertext_modulus);

    for (input_lwe, mut input_glwe) in input.iter().zip(input_glwes.iter_mut()) {
        lwe_preprocessing(&input_lwe, &mut buf, polynomial_size);
        convert_lwe_to_glwe_const(&buf, &mut input_glwe);
    }

    let mut buf = pack_lwes(&input_glwes, auto_keys);
    trace_partial_assign(&mut buf, auto_keys, lwe_count);
    glwe_ciphertext_clone_from(output, &buf);
}


fn pack_lwes<Scalar, Cont>(
    input: &GlweCiphertextList<Cont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
) -> GlweCiphertextOwned<Scalar> where
    Scalar: UnsignedTorus,
    Cont: Container<Element=Scalar>,
{
    let glwe_size = input.glwe_size();
    let polynomial_size = input.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    let mut output = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);

    let lwe_count = input.glwe_ciphertext_count().0;
    if lwe_count == 1 {
        glwe_ciphertext_clone_from(&mut output, &input.get(0));
    } else {
        assert_eq!(lwe_count % 2, 0);

        let half_lwe_count = lwe_count / 2;
        let mut input_even = GlweCiphertextList::new(Scalar::ZERO, glwe_size, polynomial_size, GlweCiphertextCount(half_lwe_count), ciphertext_modulus);
        let mut input_odd = GlweCiphertextList::new(Scalar::ZERO, glwe_size, polynomial_size, GlweCiphertextCount(half_lwe_count), ciphertext_modulus);

        for (i, (mut lwe_even, mut lwe_odd)) in input_even.iter_mut().zip(input_odd.iter_mut()).enumerate() {
            glwe_ciphertext_clone_from(&mut lwe_even, &input.get(2*i));
            glwe_ciphertext_clone_from(&mut lwe_odd, &input.get(2*i+1));
        }

        let output_even = pack_lwes(&input_even, auto_keys);
        let output_odd = pack_lwes(&input_odd, auto_keys);

        let mut buf = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
        glwe_ciphertext_sub_assign(&mut buf, &output_odd);
        glwe_ciphertext_monic_monomial_mul_assign(&mut buf, MonomialDegree(polynomial_size.0 / lwe_count));
        glwe_ciphertext_add_assign(&mut buf, &output_even);
        let auto_key = auto_keys.get(&(lwe_count + 1)).unwrap();
        auto_key.auto(&mut output, &buf);

        glwe_ciphertext_clone_from(&mut buf, &output_odd);
        glwe_ciphertext_monic_monomial_mul_assign(&mut buf, MonomialDegree(polynomial_size.0 / lwe_count));
        glwe_ciphertext_add_assign(&mut buf, &output_even);

        glwe_ciphertext_add_assign(&mut output, &buf);
    }

    output
}

