use tfhe::core_crypto::{
    prelude::*,
    algorithms::polynomial_algorithms::*,
};

/* -------- Error Tracking -------- */
pub fn get_val_and_bit_err<Scalar, KeyCont, C>(
    lwe_secret_key: &LweSecretKey<KeyCont>,
    lwe_ctxt: &LweCiphertext<C>,
    correct_val: Scalar,
    delta: Scalar,
) -> (Scalar, u32)
where
    Scalar: UnsignedInteger,
    KeyCont: Container<Element=Scalar>,
    C: Container<Element=Scalar>,
{
    let scaling = lwe_ctxt.ciphertext_modulus().get_power_of_two_scaling_to_native_torus();

    let decrypted = decrypt_lwe_ciphertext(&lwe_secret_key, &lwe_ctxt).0;
    let decrypted = decrypted * scaling;

    let abs_err = {
        let correct_val = correct_val * delta * scaling;
        let d0 = decrypted.wrapping_sub(correct_val);
        let d1 = correct_val.wrapping_sub(decrypted);
        std::cmp::min(d0, d1) / scaling
    };
    let bit_err = if abs_err != Scalar::ZERO {Scalar::BITS as u32 - abs_err.leading_zeros()} else {0};
    let rounding = (decrypted & (delta >> 1)) << 1;
    let decoded = (decrypted.wrapping_add(rounding)) / delta;

    return (decoded, bit_err);
}

pub fn get_val_and_abs_err<Scalar, KeyCont, C>(
    lwe_secret_key: &LweSecretKey<KeyCont>,
    lwe_ctxt: &LweCiphertext<C>,
    correct_val: Scalar,
    delta: Scalar,
) -> (Scalar, Scalar)
where
    Scalar: UnsignedInteger,
    KeyCont: Container<Element=Scalar>,
    C: Container<Element=Scalar>,
{
    let scaling = lwe_ctxt.ciphertext_modulus().get_power_of_two_scaling_to_native_torus();

    let decrypted = decrypt_lwe_ciphertext(&lwe_secret_key, &lwe_ctxt).0;
    let decrypted = decrypted * scaling;

    let abs_err = {
        let correct_val = correct_val * delta * scaling;
        let d0 = decrypted.wrapping_sub(correct_val);
        let d1 = correct_val.wrapping_sub(decrypted);
        std::cmp::min(d0, d1) / scaling
    };
    let rounding = (decrypted & (delta >> 1)) << 1;
    let decoded = (decrypted.wrapping_add(rounding)) / delta;

    return (decoded, abs_err);
}

pub fn get_val_and_bit_and_abs_err<Scalar, KeyCont, C>(
    lwe_secret_key: &LweSecretKey<KeyCont>,
    lwe_ctxt: &LweCiphertext<C>,
    correct_val: Scalar,
    delta: Scalar,
) -> (Scalar, u32, Scalar)
where
    Scalar: UnsignedInteger,
    KeyCont: Container<Element=Scalar>,
    C: Container<Element=Scalar>,
{
    let scaling = lwe_ctxt.ciphertext_modulus().get_power_of_two_scaling_to_native_torus();

    let decrypted = decrypt_lwe_ciphertext(&lwe_secret_key, &lwe_ctxt).0;
    let decrypted = decrypted * scaling;

    let abs_err = {
        let correct_val = correct_val * delta * scaling;
        let d0 = decrypted.wrapping_sub(correct_val);
        let d1 = correct_val.wrapping_sub(decrypted);
        std::cmp::min(d0, d1) / scaling
    };
    let bit_err = if abs_err != Scalar::ZERO {Scalar::BITS as u32 - abs_err.leading_zeros()} else {0};
    let rounding = (decrypted & (delta >> 1)) << 1;
    let decoded = (decrypted.wrapping_add(rounding)) / delta;

    return (decoded, bit_err, abs_err);
}

pub fn get_glwe_max_err<Scalar, KeyCont, InputCont, PtCont>(
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    input: &GlweCiphertext<InputCont>,
    correct_val_list: &PlaintextList<PtCont>,
) -> Scalar
where
    Scalar: UnsignedTorus,
    KeyCont: Container<Element=Scalar>,
    InputCont: Container<Element=Scalar>,
    PtCont: Container<Element=Scalar>,
{
    assert!(input.ciphertext_modulus().is_compatible_with_native_modulus());
    assert_eq!(glwe_secret_key.glwe_dimension(), input.glwe_size().to_glwe_dimension());
    assert_eq!(glwe_secret_key.polynomial_size(), input.polynomial_size());
    let polynomial_size = input.polynomial_size().0;

    let mut dec = PlaintextList::new(Scalar::ZERO, PlaintextCount(polynomial_size));
    decrypt_glwe_ciphertext(&glwe_secret_key, &input, &mut dec);

    let scaling = input.ciphertext_modulus().get_power_of_two_scaling_to_native_torus();

    let mut max_err = Scalar::ZERO;
    for (decrypted, correct_val) in dec.iter()
        .zip(correct_val_list.iter())
    {
        let decrypted = *decrypted.0 * scaling;
        let correct_val = *correct_val.0 * scaling;

        let abs_err = {
            let d0 = decrypted.wrapping_sub(correct_val);
            let d1 = correct_val.wrapping_sub(decrypted);
            std::cmp::min(d0, d1) / scaling
        };
        max_err = std::cmp::max(max_err, abs_err);
    }

    max_err
}

pub fn get_glwe_l2_err<Scalar, KeyCont, InputCont, PtCont>(
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    input: &GlweCiphertext<InputCont>,
    correct_val_list: &PlaintextList<PtCont>,
) -> f64
where
    Scalar: UnsignedTorus + CastInto<u128>,
    KeyCont: Container<Element=Scalar>,
    InputCont: Container<Element=Scalar>,
    PtCont: Container<Element=Scalar>,
{
    assert!(input.ciphertext_modulus().is_compatible_with_native_modulus());
    assert_eq!(glwe_secret_key.glwe_dimension(), input.glwe_size().to_glwe_dimension());
    assert_eq!(glwe_secret_key.polynomial_size(), input.polynomial_size());
    let polynomial_size = input.polynomial_size().0;

    let mut dec = PlaintextList::new(Scalar::ZERO, PlaintextCount(polynomial_size));
    decrypt_glwe_ciphertext(&glwe_secret_key, &input, &mut dec);

    let scaling = input.ciphertext_modulus().get_power_of_two_scaling_to_native_torus();

    let mut avg_err_square = u128::ZERO;
    for (decrypted, correct_val) in dec.iter()
        .zip(correct_val_list.iter())
    {
        let decrypted = *decrypted.0 * scaling;
        let correct_val = *correct_val.0 * scaling;

        let abs_err = {
            let d0 = decrypted.wrapping_sub(correct_val);
            let d1 = correct_val.wrapping_sub(decrypted);
            std::cmp::min(d0, d1) / scaling
        };
        let abs_err: u128 = Scalar::cast_into(abs_err);

        avg_err_square += abs_err * abs_err;
    }

    (avg_err_square as f64).sqrt()
}

pub fn get_max_err_ggsw_bit<Scalar: UnsignedTorus>(
    glwe_secret_key: &GlweSecretKeyOwned<Scalar>,
    ggsw_in: GgswCiphertextView<Scalar>,
    correct_val: Scalar,
) -> Scalar {
    assert!(correct_val == Scalar::ZERO || correct_val == Scalar::ONE);

    let decomp_base_log = ggsw_in.decomposition_base_log().0;
    let decomp_level = ggsw_in.decomposition_level_count().0;

    let glwe_size = ggsw_in.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let polynomial_size = ggsw_in.polynomial_size().0;

    let glwe_list = ggsw_in.as_glwe_list();
    let glwe_sk_poly_list = glwe_secret_key.as_polynomial_list();

    let mut max_err = Scalar::ZERO;
    for row in 0..(glwe_size.0 * decomp_level) {
        let level = row / glwe_size.0 + 1;
        let log_scale = Scalar::BITS - level * decomp_base_log;
        let delta = Scalar::ONE << log_scale;

        let glwe = glwe_list.get(row);
        let mut dec_pt = PlaintextList::new(Scalar::ZERO, PlaintextCount(ggsw_in.polynomial_size().0));
        decrypt_glwe_ciphertext(glwe_secret_key, &glwe, &mut dec_pt);

        if correct_val == Scalar::ZERO {
            for i in 0..polynomial_size {
                let decrypted = *dec_pt.get(i).0;
                let abs_err = std::cmp::min(decrypted, decrypted.wrapping_neg());
                max_err = std::cmp::max(max_err, abs_err);
            }
        } else {
            let col = row % glwe_size.0;
            if col < glwe_dimension.0 {
                let glwe_sk_poly = glwe_sk_poly_list.get(col);
                for i in 0..polynomial_size {
                    let decrypted = *dec_pt.get(i).0;
                    let abs_err = {
                        let glwe_sk_val = glwe_sk_poly.as_ref().get(i).unwrap();
                        let correct_val = (*glwe_sk_val).wrapping_neg() << log_scale;
                        let d0 = decrypted.wrapping_sub(correct_val);
                        let d1 = correct_val.wrapping_sub(decrypted);
                        std::cmp::min(d0, d1)
                    };
                    max_err = std::cmp::max(max_err, abs_err);
                }
            } else {
                for i in 0..polynomial_size {
                    let decrypted = *dec_pt.get(i).0;
                    let abs_err = {
                        let correct_val = if i == 0 {delta} else {Scalar::ZERO};
                        let d0 = decrypted.wrapping_sub(correct_val);
                        let d1 = correct_val.wrapping_sub(decrypted);
                        std::cmp::min(d0, d1)
                    };
                    max_err = std::cmp::max(max_err, abs_err);
                }
            }
        }
    }

    max_err
}

pub fn get_max_err_ggsw_int<Scalar, InputCont, KeyCont>(
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    ggsw_in: &GgswCiphertext<InputCont>,
    correct_val: Scalar,
) -> Scalar
where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    KeyCont: Container<Element=Scalar>,
{
    let decomp_base_log = ggsw_in.decomposition_base_log().0;
    let decomp_level = ggsw_in.decomposition_level_count().0;

    let glwe_size = ggsw_in.glwe_size();
    let glwe_dimension = glwe_size.to_glwe_dimension();
    let polynomial_size = ggsw_in.polynomial_size().0;

    let glwe_list = ggsw_in.as_glwe_list();
    let glwe_sk_poly_list = glwe_secret_key.as_polynomial_list();

    let mut max_err = Scalar::ZERO;
    for row in 0..(glwe_size.0 * decomp_level) {
        let level = row / glwe_size.0 + 1;
        let log_scale = Scalar::BITS - level * decomp_base_log;
        let delta = Scalar::ONE << log_scale;

        let glwe = glwe_list.get(row);
        let mut dec_pt = PlaintextList::new(Scalar::ZERO, PlaintextCount(ggsw_in.polynomial_size().0));
        decrypt_glwe_ciphertext(glwe_secret_key, &glwe, &mut dec_pt);

        if correct_val == Scalar::ZERO {
            for i in 0..polynomial_size {
                let decrypted = *dec_pt.get(i).0;
                let abs_err: Scalar = std::cmp::min(decrypted, decrypted.wrapping_neg());
                max_err = std::cmp::max(max_err, abs_err);
            }
        } else {
            let col = row % glwe_size.0;
            if col < glwe_dimension.0 {
                let glwe_sk_poly = glwe_sk_poly_list.get(col);
                for i in 0..polynomial_size {
                    let decrypted = *dec_pt.get(i).0;
                    let abs_err = {
                        let glwe_sk_val = glwe_sk_poly.as_ref().get(i).unwrap();
                        let correct_val = (*glwe_sk_val * correct_val).wrapping_neg() << log_scale;
                        let d0 = decrypted.wrapping_sub(correct_val);
                        let d1 = correct_val.wrapping_sub(decrypted);
                        std::cmp::min(d0, d1)
                    };
                    max_err = std::cmp::max(max_err, abs_err);
                }
            } else {
                for i in 0..polynomial_size {
                    let decrypted = *dec_pt.get(i).0;
                    let abs_err = {
                        let correct_val = if i == 0 {correct_val * delta} else {Scalar::ZERO};
                        let d0 = decrypted.wrapping_sub(correct_val);
                        let d1 = correct_val.wrapping_sub(decrypted);
                        std::cmp::min(d0, d1)
                    };
                    max_err = std::cmp::max(max_err, abs_err);
                }
            }
        }
    }
    max_err
}

/* -------- LWE List -------- */
pub fn lwe_ciphertext_list_add_assign<Scalar, LhsContMut, RhsCont>(
    lhs: &mut LweCiphertextList<LhsContMut>,
    rhs: LweCiphertextList<RhsCont>,
) where
    Scalar: UnsignedInteger,
    LhsContMut: ContainerMut<Element=Scalar>,
    RhsCont: Container<Element=Scalar>,
{
    for (mut lwe_lhs, lwe_rhs) in lhs.iter_mut().zip(rhs.iter()) {
        lwe_ciphertext_add_assign(&mut lwe_lhs, &lwe_rhs);
    }
}

/* -------- GLWE -------- */
pub fn glwe_ciphertext_add_monic_mul_assign<Scalar, LhsCont, RhsCont>(
    lhs: &mut GlweCiphertext<LhsCont>,
    rhs: &GlweCiphertext<RhsCont>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    LhsCont: ContainerMut<Element=Scalar>,
    RhsCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(lhs.glwe_size(), rhs.glwe_size());
    assert_eq!(lhs.polynomial_size(), rhs.polynomial_size());
    assert_eq!(lhs.ciphertext_modulus(), rhs.ciphertext_modulus());

    for (mut lhs_poly, rhs_poly) in lhs.as_mut_polynomial_list().iter_mut().zip(rhs.as_polynomial_list().iter()) {
        let mut buf = Polynomial::new(Scalar::ZERO, rhs.polynomial_size());
        polynomial_wrapping_monic_monomial_mul_and_subtract(&mut buf, &rhs_poly, monomial_degree);
        polynomial_wrapping_sub_assign(&mut lhs_poly, &buf);
    }
}

pub fn glwe_ciphertext_monic_monomial_mul<Scalar, LhsCont, RhsCont>(
    glwe_lhs: &mut GlweCiphertext<LhsCont>,
    glwe_rhs: &GlweCiphertext<RhsCont>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    LhsCont: ContainerMut<Element=Scalar>,
    RhsCont: Container<Element=Scalar>,
{
    for (mut poly_lhs, poly_rhs) in glwe_lhs.as_mut_polynomial_list().iter_mut()
        .zip(glwe_rhs.as_polynomial_list().iter())
    {
        polynomial_wrapping_monic_monomial_mul(&mut poly_lhs, &poly_rhs, monomial_degree);
    }
}

pub fn glwe_ciphertext_monic_monomial_mul_assign<Scalar, ContMut>(
    glwe: &mut GlweCiphertext<ContMut>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    ContMut: ContainerMut<Element=Scalar>,
{
    for mut poly in glwe.as_mut_polynomial_list().iter_mut() {
        polynomial_wrapping_monic_monomial_mul_assign(&mut poly, monomial_degree);
    }
}

pub fn glwe_ciphertext_monic_monomial_div<Scalar, LhsCont, RhsCont>(
    glwe_lhs: &mut GlweCiphertext<LhsCont>,
    glwe_rhs: &GlweCiphertext<RhsCont>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    LhsCont: ContainerMut<Element=Scalar>,
    RhsCont: Container<Element=Scalar>,
{
    for (mut poly_lhs, poly_rhs) in glwe_lhs.as_mut_polynomial_list().iter_mut()
        .zip(glwe_rhs.as_polynomial_list().iter())
    {
        polynomial_wrapping_monic_monomial_div(&mut poly_lhs, &poly_rhs, monomial_degree);
    }
}

pub fn glwe_ciphertext_monic_monomial_div_assign<Scalar, ContMut>(
    glwe: &mut GlweCiphertext<ContMut>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    ContMut: ContainerMut<Element=Scalar>,
{
    for mut poly in glwe.as_mut_polynomial_list().iter_mut() {
        polynomial_wrapping_monic_monomial_div_assign(&mut poly, monomial_degree);
    }
}

pub fn glwe_ciphertext_clone_from<Scalar, OutputCont, InputCont>(
    dst: &mut GlweCiphertext<OutputCont>,
    src: &GlweCiphertext<InputCont>,
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    debug_assert!(dst.glwe_size() == src.glwe_size());
    debug_assert!(dst.polynomial_size() == src.polynomial_size());
    dst.as_mut().clone_from_slice(src.as_ref());
}

pub fn encode_bits_into_glwe_ciphertext<Scalar, G>(
    glwe_secret_key: &GlweSecretKeyOwned<Scalar>,
    bit_list: &[Scalar],
    ggsw_bit_decomp_base_log: DecompositionBaseLog,
    ggsw_bit_decomp_level_count: DecompositionLevelCount,
    noise_parameters: impl DispersionParameter,
    generator: &mut EncryptionRandomGenerator<G>,
    ciphertext_modulus: CiphertextModulus<Scalar>,
) -> Vec<GlweCiphertextListOwned<Scalar>>
where
    Scalar: UnsignedTorus,
    G: ByteRandomGenerator,
{
    let glwe_size = glwe_secret_key.glwe_dimension().to_glwe_size();
    let polynomial_size = glwe_secret_key.polynomial_size();
    let num_glwe_list = bit_list.len() / polynomial_size.0;
    let num_glwe_list = if bit_list.len() % polynomial_size.0 == 0 {num_glwe_list} else {num_glwe_list + 1};

    let mut vec_glwe_list = vec![GlweCiphertextList::new(
        Scalar::ZERO,
        glwe_size,
        polynomial_size,
        GlweCiphertextCount(ggsw_bit_decomp_level_count.0),
        ciphertext_modulus,
    ); num_glwe_list];

    for (idx, glwe_list) in vec_glwe_list.iter_mut().enumerate() {
        for (k, mut glwe) in glwe_list.iter_mut().enumerate() {
            let log_scale = Scalar::BITS - ggsw_bit_decomp_base_log.0 * (k + 1) - log2(polynomial_size.0);
            let pt = PlaintextList::from_container(
                (0..polynomial_size.0).map(|i| {
                    let bit_idx = idx * polynomial_size.0 + i;
                    if bit_idx < bit_list.len() {
                        bit_list[bit_idx] << log_scale
                    } else {
                        Scalar::ZERO
                    }
                }).collect::<Vec<Scalar>>()
            );

            encrypt_glwe_ciphertext(&glwe_secret_key, &mut glwe, &pt, noise_parameters, generator);
        }
    }

    vec_glwe_list
}

/* -------- Automorphism -------- */
#[inline]
pub const fn log2(input: usize) -> usize {
    core::mem::size_of::<usize>() * 8 - (input.leading_zeros() as usize) - 1
}

/// Evaluate f(x) on x^k, where k is odd
pub(crate) fn eval_x_k<Scalar>(poly: PolynomialView<'_, Scalar>, k: usize) -> PolynomialOwned<Scalar>
where
    Scalar: UnsignedTorus,
{
    let mut out = PolynomialOwned::new(Scalar::ZERO, poly.polynomial_size());
    eval_x_k_in_memory(&mut out, poly, k);
    out
}

/// Evaluate f(x) on x^k, where k is odd
pub(crate) fn eval_x_k_in_memory<Scalar>(out: &mut PolynomialOwned<Scalar>, poly: PolynomialView<'_, Scalar>, k: usize)
where
    Scalar: UnsignedTorus,
{
    assert_eq!(k % 2, 1);
    assert!(poly.polynomial_size().0.is_power_of_two());
    *out.as_mut().get_mut(0).unwrap() = *poly.as_ref().get(0).unwrap();
    for i in 1..poly.polynomial_size().0 {
        // i-th term becomes ik-th term, but reduced by n
        let j = i * k % poly.polynomial_size().0;
        let sign = if ((i * k) / poly.polynomial_size().0) % 2 == 0
        { Scalar::ONE } else { Scalar::MAX };
        let c = *poly.as_ref().get(i).unwrap();
        *out.as_mut().get_mut(j).unwrap() = sign.wrapping_mul(c);
    }
}

/* -------- Polynomial Algorithm -------- */
/// Multiply (mod $(X^{N}+1)$), the input polynomial with a monic monomial of a given degree i.e.
/// $X^{degree}$, then subtract the input from the result and assign to the output.
///
/// output = input * X^degree - input
///
/// # Note
/// (Code from tfhe::core_crypto::algorithms::polynomial_algorithms)
///
/// Computations wrap around (similar to computing modulo $2^{n\_{bits}}$) when exceeding the
/// unsigned integer capacity.
pub(crate) fn polynomial_wrapping_monic_monomial_mul_and_subtract<Scalar, OutputCont, InputCont>(
    output: &mut Polynomial<OutputCont>,
    input: &Polynomial<InputCont>,
    monomial_degree: MonomialDegree,
) where
    Scalar: UnsignedInteger,
    OutputCont: ContainerMut<Element = Scalar>,
    InputCont: Container<Element = Scalar>,
{
    /// performs the operation: dst = -src - src_orig, with wrapping arithmetic
    fn copy_with_neg_and_subtract<Scalar: UnsignedInteger>(
        dst: &mut [Scalar],
        src: &[Scalar],
        src_orig: &[Scalar],
    ) {
        for ((dst, src), src_orig) in dst.iter_mut().zip(src).zip(src_orig) {
            *dst = src.wrapping_neg().wrapping_sub(*src_orig);
        }
    }

    /// performs the operation: dst = src - src_orig, with wrapping arithmetic
    fn copy_without_neg_and_subtract<Scalar: UnsignedInteger>(
        dst: &mut [Scalar],
        src: &[Scalar],
        src_orig: &[Scalar],
    ) {
        for ((dst, src), src_orig) in dst.iter_mut().zip(src).zip(src_orig) {
            *dst = src.wrapping_sub(*src_orig);
        }
    }

    assert!(
        output.polynomial_size() == input.polynomial_size(),
        "Output polynomial size {:?} is not the same as input polynomial size {:?}.",
        output.polynomial_size(),
        input.polynomial_size(),
    );

    let polynomial_size = output.polynomial_size().0;
    let remaining_degree = monomial_degree.0 % polynomial_size;

    let full_cycles_count = monomial_degree.0 / polynomial_size;
    if full_cycles_count % 2 == 0 {
        copy_with_neg_and_subtract(
            &mut output[..remaining_degree],
            &input[polynomial_size - remaining_degree..],
            &input[..remaining_degree],
        );
        copy_without_neg_and_subtract(
            &mut output[remaining_degree..],
            &input[..polynomial_size - remaining_degree],
            &input[remaining_degree..],
        );
    } else {
        copy_without_neg_and_subtract(
            &mut output[..remaining_degree],
            &input[polynomial_size - remaining_degree..],
            &input[..remaining_degree],
        );
        copy_with_neg_and_subtract(
            &mut output[remaining_degree..],
            &input[..polynomial_size - remaining_degree],
            &input[remaining_degree..],
        );
    }
}

/* -------- Macro -------- */
// https://docs.rs/itertools/0.7.8/src/itertools/lib.rs.html#247-269
#[allow(unused_macros)]
macro_rules! izip {
    (@ __closure @ ($a:expr)) => { |a| (a,) };
    (@ __closure @ ($a:expr, $b:expr)) => { |(a, b)| (a, b) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr)) => { |((a, b), c)| (a, b, c) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr)) => { |(((a, b), c), d)| (a, b, c, d) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr)) => { |((((a, b), c), d), e)| (a, b, c, d, e) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr)) => { |(((((a, b), c), d), e), f)| (a, b, c, d, e, f) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr)) => { |((((((a, b), c), d), e), f), g)| (a, b, c, d, e, f, e) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr)) => { |(((((((a, b), c), d), e), f), g), h)| (a, b, c, d, e, f, g, h) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr)) => { |((((((((a, b), c), d), e), f), g), h), i)| (a, b, c, d, e, f, g, h, i) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr)) => { |(((((((((a, b), c), d), e), f), g), h), i), j)| (a, b, c, d, e, f, g, h, i, j) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr)) => { |((((((((((a, b), c), d), e), f), g), h), i), j), k)| (a, b, c, d, e, f, g, h, i, j, k) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr)) => { |(((((((((((a, b), c), d), e), f), g), h), i), j), k), l)| (a, b, c, d, e, f, g, h, i, j, k, l) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr)) => { |((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m)| (a, b, c, d, e, f, g, h, i, j, k, l, m) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr, $n:expr)) => { |(((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m), n)| (a, b, c, d, e, f, g, h, i, j, k, l, m, n) };
    (@ __closure @ ($a:expr, $b:expr, $c:expr, $d:expr, $e: expr, $f:expr, $g:expr, $h:expr, $i: expr, $j: expr, $k: expr, $l: expr, $m:expr, $n:expr, $o:expr)) => { |((((((((((((((a, b), c), d), e), f), g), h), i), j), k), l), m), n), o)| (a, b, c, d, e, f, g, h, i, j, k, l, m, n, o) };

    ( $first:expr $(,)?) => {
        {
            #[allow(unused_imports)]
            use $crate::core_crypto::commons::utils::ZipChecked;
            ::core::iter::IntoIterator::into_iter($first)
        }
    };
    ( $first:expr, $($rest:expr),+ $(,)?) => {
        {
            #[allow(unused_imports)]
            use tfhe::core_crypto::commons::utils::ZipChecked;
            ::core::iter::IntoIterator::into_iter($first)
                $(.zip_checked($rest))*
                .map($crate::utils::izip!(@ __closure @ ($first, $($rest),*)))
        }
    };
}

#[allow(unused_imports)]
pub(crate) use izip;
