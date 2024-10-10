use crate::wwlp_cbs_params::*;
use crate::FftType;
use tfhe::core_crypto::prelude::*;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref CBS_WOPBS_2_2: CBSParam<u64> = CBSParam::new(
        LweDimension(769), // lwe_dimension
        StandardDev(0.0000043131554647504185), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(15), // pbs_base_log
        DecompositionLevelCount(2), // pbs_level
        DecompositionBaseLog(6), // ks_base_log
        DecompositionLevelCount(2), // ks_level
        DecompositionBaseLog(15), // pfks_base_log
        DecompositionLevelCount(2), // pfks_level
        DecompositionBaseLog(5), // cbs_base_log
        DecompositionLevelCount(3), // cbs_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref CBS_WOPBS_3_3: CBSParam<u64> = CBSParam::new(
        LweDimension(873), // lwe_dimension
        StandardDev(0.0000006428797112843789), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(9), // pbs_base_log
        DecompositionLevelCount(4), // pbs_level
        DecompositionBaseLog(10), // ks_base_log
        DecompositionLevelCount(1), // ks_level
        DecompositionBaseLog(9), // pfks_base_log
        DecompositionLevelCount(4), // pfks_level
        DecompositionBaseLog(6), // cbs_base_log
        DecompositionLevelCount(3), // cbs_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref CBS_WOPBS_4_4: CBSParam<u64> = CBSParam::new(
        LweDimension(953), // lwe_dimension
        StandardDev(0.0000001486733969411098), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(9), // pbs_base_log
        DecompositionLevelCount(4), // pbs_level
        DecompositionBaseLog(11), // ks_base_log
        DecompositionLevelCount(1), // ks_level
        DecompositionBaseLog(9), // pfks_base_log
        DecompositionLevelCount(4), // pfks_level
        DecompositionBaseLog(4), // cbs_base_log
        DecompositionLevelCount(6), // cbs_level
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref WWLP_CBS_WOPBS_2_2: WWLpCBSParam<u64> = WWLpCBSParam::new(
        LweDimension(769), // lwe_dimension
        StandardDev(0.0000043131554647504185), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        DecompositionBaseLog(15), // pbs_base_log
        DecompositionLevelCount(2), // pbs_level
        DecompositionBaseLog(23), // refresh_base_log
        DecompositionLevelCount(1), // refresh_level
        DecompositionBaseLog(6), // ks_base_log
        DecompositionLevelCount(2), // ks_level
        DecompositionBaseLog(7), // auto_base_log
        DecompositionLevelCount(7), // auto_level
        FftType::Split(37), // fft_type_auto
        DecompositionBaseLog(8), // ss_base_log
        DecompositionLevelCount(6), // ss_level
        DecompositionBaseLog(5), // cbs_base_log
        DecompositionLevelCount(3), // cbs_level
        LutCountLog(2), // log_lut_count
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref HIGHPREC_WWLP_CBS_WOPBS_3_3: HighPrecWWLpCBSParam<u64> = HighPrecWWLpCBSParam::new(
        LweDimension(873), // lwe_dimension
        StandardDev(0.0000006428797112843789), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        GlweDimension(2), // large_glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        StandardDev(0.0000000000000000002168404344971009), // large_glwe_modular_std_dev
        DecompositionBaseLog(9), // pbs_base_log
        DecompositionLevelCount(4), // pbs_level
        DecompositionBaseLog(23), // refresh_base_log
        DecompositionLevelCount(1), // refresh_level
        DecompositionBaseLog(10), // ks_base_log
        DecompositionLevelCount(1), // ks_level
        DecompositionBaseLog(15), // glwe_ds_to_large_base_log
        DecompositionLevelCount(3), // glwe_ds_to_large_level
        FftType::Split(44), // fft_type_to_large
        DecompositionBaseLog(6), // auto_base_log
        DecompositionLevelCount(10), // auto_level
        FftType::Split(36), // fft_type_auto
        DecompositionBaseLog(5), // glwe_ds_from_large_base_log
        DecompositionLevelCount(10), // glwe_ds_from_large_level
        FftType::Split(35), // fft_type_from_large
        DecompositionBaseLog(6), // ss_base_log
        DecompositionLevelCount(9), // ss_level
        DecompositionBaseLog(6), // cbs_base_log
        DecompositionLevelCount(3), // cbs_level
        LutCountLog(2), // log_lut_count
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );

    pub static ref HIGHPREC_WWLP_CBS_WOPBS_4_4: HighPrecWWLpCBSParam<u64> = HighPrecWWLpCBSParam::new(
        LweDimension(953), // lwe_dimension
        StandardDev(0.0000001486733969411098), // lwe_modular_std_dev
        PolynomialSize(2048), // polynomial_size
        GlweDimension(1), // glwe_dimension
        GlweDimension(2), // large_glwe_dimension
        StandardDev(0.00000000000000029403601535432533), // glwe_modular_std_dev
        StandardDev(0.0000000000000000002168404344971009), // large_glwe_modular_std_dev
        DecompositionBaseLog(9), // pbs_base_log
        DecompositionLevelCount(4), // pbs_level
        DecompositionBaseLog(23), // refresh_base_log
        DecompositionLevelCount(1), // refresh_level
        DecompositionBaseLog(11), // ks_base_log
        DecompositionLevelCount(1), // ks_level
        DecompositionBaseLog(15), // glwe_ds_to_large_base_log
        DecompositionLevelCount(3), // glwe_ds_to_large_level
        FftType::Split(44), // fft_type_to_large
        // DecompositionBaseLog(6), // auto_base_log
        // DecompositionLevelCount(10), // auto_level
        // FftType::Split(36), // fft_type_auto
        DecompositionBaseLog(12), // auto_base_log
        DecompositionLevelCount(4), // auto_level
        FftType::Split(42), // fft_type_auto
        // DecompositionBaseLog(5), // glwe_ds_from_large_base_log
        // DecompositionLevelCount(10), // glwe_ds_from_large_level
        // FftType::Split(35), // fft_type_from_large
        DecompositionBaseLog(10), // glwe_ds_from_large_base_log
        DecompositionLevelCount(4), // glwe_ds_from_large_level
        FftType::Split(40), // fft_type_from_large
        // DecompositionBaseLog(6), // ss_base_log
        // DecompositionLevelCount(9), // ss_level
        DecompositionBaseLog(10), // ss_base_log
        DecompositionLevelCount(4), // ss_level
        // DecompositionBaseLog(4), // cbs_base_log
        // DecompositionLevelCount(6), // cbs_level
        DecompositionBaseLog(3), // cbs_base_log
        DecompositionLevelCount(8), // cbs_level
        LutCountLog(3), // log_lut_count
        CiphertextModulus::<u64>::new_native(), // ciphertext_modulus
    );
}
