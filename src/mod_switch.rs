use tfhe::core_crypto::prelude::*;

pub fn lwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );
    assert!(
        output.ciphertext_modulus().is_non_native_power_of_two(),
        "output ciphertext modulus is not non-native power-of-two"
    );

    let output_ciphertext_modulus = output.ciphertext_modulus();
    let divisor = output_ciphertext_modulus.get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut())
    {
        *dst = *src - (*src) % divisor;
    }
}

pub fn lwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_non_native_power_of_two(),
        "input ciphertext modulus is not non-native power-of-two"
    );
    assert!(
        output.ciphertext_modulus().is_native_modulus(),
        "output ciphertext modulus is not native"
    );

    let input_ciphertext_modulus = input.ciphertext_modulus();
    let scaling_factor = input_ciphertext_modulus.get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut()) {
        *dst = (*src).wrapping_div(scaling_factor);
    }
}

pub fn lwe_preprocessing_assign<Scalar, ContMut>(
    input: &mut LweCiphertext<ContMut>,
    polynomial_size: PolynomialSize,
) where
    Scalar: UnsignedInteger,
    ContMut: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );
    assert!(
        Scalar::BITS > polynomial_size.0.ilog2() as usize
    );

    let log_small_q = Scalar::BITS - polynomial_size.0.ilog2() as usize;
    let small_ciphertext_modulus = CiphertextModulus::<Scalar>::try_new_power_of_2(log_small_q).unwrap();

    let mut buf = LweCiphertext::new(Scalar::ZERO, input.lwe_size(), small_ciphertext_modulus);

    lwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two(&input, &mut buf);
    lwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native(&buf, input);
}

pub fn lwe_preprocessing<Scalar, InputCont, OutputCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
    polynomial_size: PolynomialSize,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );

    output.as_mut().clone_from_slice(input.as_ref());
    lwe_preprocessing_assign(output, polynomial_size);
}

pub fn glwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two<Scalar, InputCont, OutputCont>(
    input: &GlweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );
    assert!(
        output.ciphertext_modulus().is_non_native_power_of_two(),
        "output ciphertext modulus is not non-native power-of-two"
    );
    assert_eq!(
        input.polynomial_size(),
        output.polynomial_size(),
    );
    assert_eq!(
        input.glwe_size(),
        output.glwe_size(),
    );

    let output_ciphertext_modulus = output.ciphertext_modulus();
    let divisor = output_ciphertext_modulus.get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut())
    {
        *dst = *src - (*src) % divisor;
    }
}

pub fn glwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native<Scalar, InputCont, OutputCont>(
    input: &GlweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_non_native_power_of_two(),
        "input ciphertext modulus is not non-native power-of-two"
    );
    assert!(
        output.ciphertext_modulus().is_native_modulus(),
        "output ciphertext modulus is not native"
    );
    assert_eq!(
        input.polynomial_size(),
        output.polynomial_size(),
    );
    assert_eq!(
        input.glwe_size(),
        output.glwe_size(),
    );

    let input_ciphertext_modulus = input.ciphertext_modulus();
    let scaling_factor = input_ciphertext_modulus.get_power_of_two_scaling_to_native_torus();
    for (src, dst) in input.as_ref().iter().zip(output.as_mut().iter_mut()) {
        *dst = (*src).wrapping_div(scaling_factor);
    }
}

pub fn glwe_preprocessing_assign<Scalar, ContMut>(
    input: &mut GlweCiphertext<ContMut>,
) where
    Scalar: UnsignedInteger,
    ContMut: ContainerMut<Element=Scalar>,
{
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );

    let polynomial_size = input.polynomial_size();

    assert!(
        Scalar::BITS > polynomial_size.0.ilog2() as usize
    );

    let log_small_q = Scalar::BITS - polynomial_size.0.ilog2() as usize;
    let small_ciphertext_modulus = CiphertextModulus::<Scalar>::try_new_power_of_2(log_small_q).unwrap();
    let glwe_size = input.glwe_size();

    let mut buf = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, small_ciphertext_modulus);

    glwe_ciphertext_mod_switch_from_native_to_non_native_power_of_two(&input, &mut buf);
    glwe_ciphertext_mod_raise_from_non_native_power_of_two_to_native(&buf, input);
}

pub fn glwe_preprocessing<Scalar, InputCont, OutputCont>(
    input: &GlweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    Scalar: UnsignedInteger,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    assert!(
        input.ciphertext_modulus().is_native_modulus(),
        "input ciphertext modulus is not native"
    );
    assert_eq!(input.polynomial_size(), output.polynomial_size());
    assert_eq!(input.glwe_size(), output.glwe_size());

    output.as_mut().clone_from_slice(input.as_ref());
    glwe_preprocessing_assign(output);
}
