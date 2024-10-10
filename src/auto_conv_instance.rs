use crate::auto_conv_params::*;
use crate::FftType;
use tfhe::core_crypto::prelude::*;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref PKSK_PARAM_2_2: PkskConvParam<u64> = PkskConvParam::new(
        LweDimension(742), // lwe_dimension
        StandardDev(0.000007069849454709433), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(24), // pksk_base_log
        DecompositionLevelCount(1), // pksk_level
        DecompositionBaseLog(3), // ks_base_log
        DecompositionLevelCount(5), // ks_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_2_2_VANILLA: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(12), // auto_base_log
        DecompositionLevelCount(3), // auto_level
        FftType::Vanilla, // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_2_2_LEV_3: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(13), // auto_base_log
        DecompositionLevelCount(3), // auto_level
        FftType::Split(42), // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref PKSK_PARAM_3_3: PkskConvParam<u64> = PkskConvParam::new(
        LweDimension(864), // lwe_dimension
        StandardDev(0.000000757998020150446), // lwe_modular_std_dev
        PolynomialSize(8192), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(30), // pksk_base_log
        DecompositionLevelCount(1), // pksk_level
        DecompositionBaseLog(3), // ks_base_log
        DecompositionLevelCount(6), // ks_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_3_3_LEV_4: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(8192), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(12), // auto_base_log
        DecompositionLevelCount(4), // auto_level
        FftType::Split(43), // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_3_3_LEV_5: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(8192), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(10), // auto_base_log
        DecompositionLevelCount(5), // auto_level
        FftType::Split(41), // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref PKSK_PARAM_4_4: PkskConvParam<u64> = PkskConvParam::new(
        LweDimension(996), // lwe_dimension
        StandardDev(0.00000006767666038309478), // lwe_modular_std_dev
        PolynomialSize(32768), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(32), // pksk_base_log
        DecompositionLevelCount(1), // pksk_level
        DecompositionBaseLog(3), // ks_base_log
        DecompositionLevelCount(7), // ks_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_4_4_LEV_3: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(32768), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(15), // auto_base_log
        DecompositionLevelCount(3), // auto_level
        FftType::Split16, // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref AUTO_PARAM_4_4_LEV_4: AutoConvParam<u64> = AutoConvParam::new(
        PolynomialSize(32768), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.0000000000000000002168404344971009), // glwe_modular_std_dev
        DecompositionBaseLog(13), // auto_base_log
        DecompositionLevelCount(4), // auto_level
        FftType::Split16, // fft_type
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );
}
