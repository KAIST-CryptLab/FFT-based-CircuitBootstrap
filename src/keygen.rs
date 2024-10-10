use tfhe::core_crypto::prelude::*;

use crate::{allocate_and_generate_new_glwe_keyswitch_key, convert_standard_glwe_keyswitch_key_to_fourier, FftType, FourierGlweKeyswitchKey, FourierGlweKeyswitchKeyOwned};

pub fn keygen_pbs<Scalar: UnsignedTorus, G: ByteRandomGenerator>(
    lwe_dimension: LweDimension,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    lwe_modular_std_dev: impl DispersionParameter,
    glwe_modular_std_dev: impl DispersionParameter,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    secret_generator: &mut SecretRandomGenerator<G>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> (
    LweSecretKey<Vec<Scalar>>,
    GlweSecretKey<Vec<Scalar>>,
    LweSecretKey<Vec<Scalar>>,
    FourierLweBootstrapKeyOwned,
    LweKeyswitchKey<Vec<Scalar>>,
) {
    let small_lwe_secret_key: LweSecretKey<Vec<Scalar>> = LweSecretKey::generate_new_binary(lwe_dimension, secret_generator);
    let glwe_secret_key: GlweSecretKey<Vec<Scalar>> = GlweSecretKey::generate_new_binary(glwe_dimension, polynomial_size, secret_generator);
    let large_lwe_secret_key: LweSecretKey<Vec<Scalar>> = glwe_secret_key.clone().into_lwe_secret_key();

    let lwe_secret_key = large_lwe_secret_key;
    let lwe_secret_key_after_ks = small_lwe_secret_key;

    let bootstrap_key = allocate_and_generate_new_lwe_bootstrap_key(
        &lwe_secret_key_after_ks,
        &glwe_secret_key,
        pbs_base_log,
        pbs_level,
        glwe_modular_std_dev,
        CiphertextModulus::<Scalar>::new_native(),
        encryption_generator,
    );

    let mut fourier_bsk = FourierLweBootstrapKey::new(
        bootstrap_key.input_lwe_dimension(),
        bootstrap_key.glwe_size(),
        bootstrap_key.polynomial_size(),
        bootstrap_key.decomposition_base_log(),
        bootstrap_key.decomposition_level_count(),
    );
    convert_standard_lwe_bootstrap_key_to_fourier(&bootstrap_key, &mut fourier_bsk);
    drop(bootstrap_key);

    let ksk = allocate_and_generate_new_lwe_keyswitch_key(
        &lwe_secret_key,
        &lwe_secret_key_after_ks,
        ks_base_log,
        ks_level,
        lwe_modular_std_dev,
        CiphertextModulus::<Scalar>::new_native(),
        encryption_generator,
    );

    (lwe_secret_key, glwe_secret_key, lwe_secret_key_after_ks, fourier_bsk, ksk)
}

pub fn keygen_pbs_without_ksk<Scalar: UnsignedTorus, G: ByteRandomGenerator>(
    lwe_dimension: LweDimension,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    glwe_modular_std_dev: impl DispersionParameter,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    secret_generator: &mut SecretRandomGenerator<G>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> (
    LweSecretKey<Vec<Scalar>>,
    GlweSecretKey<Vec<Scalar>>,
    LweSecretKey<Vec<Scalar>>,
    FourierLweBootstrapKeyOwned,
) {
    let small_lwe_secret_key: LweSecretKey<Vec<Scalar>> = LweSecretKey::generate_new_binary(lwe_dimension, secret_generator);
    let glwe_secret_key: GlweSecretKey<Vec<Scalar>> = GlweSecretKey::generate_new_binary(glwe_dimension, polynomial_size, secret_generator);
    let large_lwe_secret_key: LweSecretKey<Vec<Scalar>> = glwe_secret_key.clone().into_lwe_secret_key();

    let lwe_secret_key = large_lwe_secret_key;
    let lwe_secret_key_after_ks = small_lwe_secret_key;

    let bootstrap_key = allocate_and_generate_new_lwe_bootstrap_key(
        &lwe_secret_key_after_ks,
        &glwe_secret_key,
        pbs_base_log,
        pbs_level,
        glwe_modular_std_dev,
        CiphertextModulus::<Scalar>::new_native(),
        encryption_generator,
    );

    let mut fourier_bsk = FourierLweBootstrapKey::new(
        bootstrap_key.input_lwe_dimension(),
        bootstrap_key.glwe_size(),
        bootstrap_key.polynomial_size(),
        bootstrap_key.decomposition_base_log(),
        bootstrap_key.decomposition_level_count(),
    );
    convert_standard_lwe_bootstrap_key_to_fourier(&bootstrap_key, &mut fourier_bsk);
    drop(bootstrap_key);

    (lwe_secret_key, glwe_secret_key, lwe_secret_key_after_ks, fourier_bsk)
}

pub fn keygen_pbs_with_glwe_ds<Scalar: UnsignedTorus, G: ByteRandomGenerator>(
    lwe_dimension: LweDimension,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    lwe_modular_std_dev: impl DispersionParameter,
    glwe_modular_std_dev: impl DispersionParameter,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    glwe_ds_base_log: DecompositionBaseLog,
    glwe_ds_level: DecompositionLevelCount,
    common_polynomial_size: PolynomialSize,
    fft_type: FftType,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
    secret_generator: &mut SecretRandomGenerator<G>,
    encryption_generator: &mut EncryptionRandomGenerator<G>,
) -> (
    LweSecretKey<Vec<Scalar>>,
    GlweSecretKey<Vec<Scalar>>,
    LweSecretKey<Vec<Scalar>>,
    FourierLweBootstrapKeyOwned,
    FourierGlweKeyswitchKeyOwned,
) {
    assert_eq!(lwe_dimension.0 % common_polynomial_size.0, 0);
    assert_eq!((glwe_dimension.0 * polynomial_size.0) % common_polynomial_size.0, 0);

    let small_lwe_secret_key: LweSecretKey<Vec<Scalar>> = LweSecretKey::generate_new_binary(lwe_dimension, secret_generator);
    let glwe_secret_key: GlweSecretKey<Vec<Scalar>> = GlweSecretKey::generate_new_binary(glwe_dimension, polynomial_size, secret_generator);
    let large_lwe_secret_key: LweSecretKey<Vec<Scalar>> = glwe_secret_key.clone().into_lwe_secret_key();

    let lwe_secret_key = large_lwe_secret_key;
    let lwe_secret_key_after_ks = small_lwe_secret_key;

    let bootstrap_key = allocate_and_generate_new_lwe_bootstrap_key(
        &lwe_secret_key_after_ks,
        &glwe_secret_key,
        pbs_base_log,
        pbs_level,
        glwe_modular_std_dev,
        ciphertext_modulus,
        encryption_generator,
    );

    let mut fourier_bsk = FourierLweBootstrapKey::new(
        bootstrap_key.input_lwe_dimension(),
        bootstrap_key.glwe_size(),
        bootstrap_key.polynomial_size(),
        bootstrap_key.decomposition_base_log(),
        bootstrap_key.decomposition_level_count(),
    );
    convert_standard_lwe_bootstrap_key_to_fourier(&bootstrap_key, &mut fourier_bsk);
    drop(bootstrap_key);

    let lwe_secret_key_view = GlweSecretKey::from_container(lwe_secret_key.as_ref(), common_polynomial_size);
    let lwe_secret_key_after_ks_view = GlweSecretKey::from_container(lwe_secret_key_after_ks.as_ref(), common_polynomial_size);
    let glwe_ksk = allocate_and_generate_new_glwe_keyswitch_key(
        &lwe_secret_key_view,
        &lwe_secret_key_after_ks_view,
        glwe_ds_base_log,
        glwe_ds_level,
        lwe_modular_std_dev,
        ciphertext_modulus,
        encryption_generator,
    );

    let mut fourier_glwe_ksk = FourierGlweKeyswitchKey::new(
        lwe_secret_key_view.glwe_dimension().to_glwe_size(),
        lwe_secret_key_after_ks_view.glwe_dimension().to_glwe_size(),
        common_polynomial_size,
        glwe_ds_base_log,
        glwe_ds_level,
        fft_type,
    );
    convert_standard_glwe_keyswitch_key_to_fourier(&glwe_ksk, &mut fourier_glwe_ksk);

    (lwe_secret_key, glwe_secret_key, lwe_secret_key_after_ks, fourier_bsk, fourier_glwe_ksk)
}