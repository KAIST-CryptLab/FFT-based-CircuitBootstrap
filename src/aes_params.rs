use tfhe::core_crypto::prelude::*;
use crate::FftType;

#[derive(Clone, Copy)]
pub struct AesParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    glwe_ds_base_log: DecompositionBaseLog,
    glwe_ds_level: DecompositionLevelCount,
    common_polynomial_size: PolynomialSize,
    fft_type_ds: FftType,
    auto_base_log: DecompositionBaseLog,
    auto_level: DecompositionLevelCount,
    fft_type_auto: FftType,
    ss_base_log: DecompositionBaseLog,
    ss_level: DecompositionLevelCount,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
}

impl<Scalar: UnsignedInteger> AesParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        pbs_base_log: DecompositionBaseLog,
        pbs_level: DecompositionLevelCount,
        glwe_ds_base_log: DecompositionBaseLog,
        glwe_ds_level: DecompositionLevelCount,
        common_polynomial_size: PolynomialSize,
        fft_type_ds: FftType,
        auto_base_log: DecompositionBaseLog,
        auto_level: DecompositionLevelCount,
        fft_type_auto: FftType,
        ss_base_log: DecompositionBaseLog,
        ss_level: DecompositionLevelCount,
        cbs_base_log: DecompositionBaseLog,
        cbs_level: DecompositionLevelCount,
        log_lut_count: LutCountLog,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
    ) -> Self {
        assert!(lwe_dimension.0 % common_polynomial_size.0 == 0);
        assert!(polynomial_size.0 % common_polynomial_size.0 == 0);

        AesParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            glwe_ds_base_log,
            glwe_ds_level,
            common_polynomial_size,
            fft_type_ds,
            auto_base_log,
            auto_level,
            fft_type_auto,
            ss_base_log,
            ss_level,
            cbs_base_log,
            cbs_level,
            log_lut_count,
            ciphertext_modulus,
        }
    }

    pub fn lwe_dimension(&self) -> LweDimension {
        self.lwe_dimension
    }

    pub fn lwe_modular_std_dev(&self) -> StandardDev {
        self.lwe_modular_std_dev
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn glwe_dimension(&self) -> GlweDimension {
        self.glwe_dimension
    }

    pub fn glwe_modular_std_dev(&self) -> StandardDev {
        self.glwe_modular_std_dev
    }

    pub fn pbs_base_log(&self) -> DecompositionBaseLog {
        self.pbs_base_log
    }

    pub fn pbs_level(&self) -> DecompositionLevelCount {
        self.pbs_level
    }

    pub fn glwe_ds_base_log(&self) -> DecompositionBaseLog {
        self.glwe_ds_base_log
    }

    pub fn glwe_ds_level(&self) -> DecompositionLevelCount {
        self.glwe_ds_level
    }

    pub fn common_polynomial_size(&self) -> PolynomialSize {
        self.common_polynomial_size
    }

    pub fn fft_type_ds(&self) -> FftType {
        self.fft_type_ds
    }

    pub fn auto_base_log(&self) -> DecompositionBaseLog {
        self.auto_base_log
    }

    pub fn auto_level(&self) -> DecompositionLevelCount {
        self.auto_level
    }

    pub fn fft_type_auto(&self) -> FftType {
        self.fft_type_auto
    }

    pub fn ss_base_log(&self) -> DecompositionBaseLog {
        self.ss_base_log
    }

    pub fn ss_level(&self) -> DecompositionLevelCount {
        self.ss_level
    }

    pub fn cbs_base_log(&self) -> DecompositionBaseLog {
        self.cbs_base_log
    }

    pub fn cbs_level(&self) -> DecompositionLevelCount {
        self.cbs_level
    }

    pub fn log_lut_count(&self) -> LutCountLog {
        self.log_lut_count
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus::<Scalar> {
        self.ciphertext_modulus
    }
}

#[derive(Clone, Copy)]
pub struct AesHalfCBSParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    glwe_ds_base_log: DecompositionBaseLog,
    glwe_ds_level: DecompositionLevelCount,
    common_polynomial_size: PolynomialSize,
    fft_type_ds: FftType,
    auto_base_log: DecompositionBaseLog,
    auto_level: DecompositionLevelCount,
    fft_type_auto: FftType,
    ss_base_log: DecompositionBaseLog,
    ss_level: DecompositionLevelCount,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
    half_cbs_auto_base_log: DecompositionBaseLog,
    half_cbs_auto_level: DecompositionLevelCount,
    half_cbs_fft_type_auto: FftType,
    half_cbs_ss_base_log: DecompositionBaseLog,
    half_cbs_ss_level: DecompositionLevelCount,
    half_cbs_base_log: DecompositionBaseLog,
    half_cbs_level: DecompositionLevelCount,
}

impl<Scalar: UnsignedInteger> AesHalfCBSParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        pbs_base_log: DecompositionBaseLog,
        pbs_level: DecompositionLevelCount,
        glwe_ds_base_log: DecompositionBaseLog,
        glwe_ds_level: DecompositionLevelCount,
        common_polynomial_size: PolynomialSize,
        fft_type_ds: FftType,
        auto_base_log: DecompositionBaseLog,
        auto_level: DecompositionLevelCount,
        fft_type_auto: FftType,
        ss_base_log: DecompositionBaseLog,
        ss_level: DecompositionLevelCount,
        cbs_base_log: DecompositionBaseLog,
        cbs_level: DecompositionLevelCount,
        log_lut_count: LutCountLog,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
        half_cbs_auto_base_log: DecompositionBaseLog,
        half_cbs_auto_level: DecompositionLevelCount,
        half_cbs_fft_type_auto: FftType,
        half_cbs_ss_base_log: DecompositionBaseLog,
        half_cbs_ss_level: DecompositionLevelCount,
        half_cbs_base_log: DecompositionBaseLog,
        half_cbs_level: DecompositionLevelCount,
    ) -> Self {
        assert!(lwe_dimension.0 % common_polynomial_size.0 == 0);
        assert!(polynomial_size.0 % common_polynomial_size.0 == 0);

        AesHalfCBSParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            glwe_ds_base_log,
            glwe_ds_level,
            common_polynomial_size,
            fft_type_ds,
            auto_base_log,
            auto_level,
            fft_type_auto,
            ss_base_log,
            ss_level,
            cbs_base_log,
            cbs_level,
            log_lut_count,
            ciphertext_modulus,
            half_cbs_auto_base_log,
            half_cbs_auto_level,
            half_cbs_fft_type_auto,
            half_cbs_ss_base_log,
            half_cbs_ss_level,
            half_cbs_base_log,
            half_cbs_level,
        }
    }

    pub fn lwe_dimension(&self) -> LweDimension {
        self.lwe_dimension
    }

    pub fn lwe_modular_std_dev(&self) -> StandardDev {
        self.lwe_modular_std_dev
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn glwe_dimension(&self) -> GlweDimension {
        self.glwe_dimension
    }

    pub fn glwe_modular_std_dev(&self) -> StandardDev {
        self.glwe_modular_std_dev
    }

    pub fn pbs_base_log(&self) -> DecompositionBaseLog {
        self.pbs_base_log
    }

    pub fn pbs_level(&self) -> DecompositionLevelCount {
        self.pbs_level
    }

    pub fn glwe_ds_base_log(&self) -> DecompositionBaseLog {
        self.glwe_ds_base_log
    }

    pub fn glwe_ds_level(&self) -> DecompositionLevelCount {
        self.glwe_ds_level
    }

    pub fn common_polynomial_size(&self) -> PolynomialSize {
        self.common_polynomial_size
    }

    pub fn fft_type_ds(&self) -> FftType {
        self.fft_type_ds
    }

    pub fn auto_base_log(&self) -> DecompositionBaseLog {
        self.auto_base_log
    }

    pub fn auto_level(&self) -> DecompositionLevelCount {
        self.auto_level
    }

    pub fn fft_type_auto(&self) -> FftType {
        self.fft_type_auto
    }

    pub fn ss_base_log(&self) -> DecompositionBaseLog {
        self.ss_base_log
    }

    pub fn ss_level(&self) -> DecompositionLevelCount {
        self.ss_level
    }

    pub fn cbs_base_log(&self) -> DecompositionBaseLog {
        self.cbs_base_log
    }

    pub fn cbs_level(&self) -> DecompositionLevelCount {
        self.cbs_level
    }

    pub fn log_lut_count(&self) -> LutCountLog {
        self.log_lut_count
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus::<Scalar> {
        self.ciphertext_modulus
    }

    pub fn half_cbs_auto_base_log(&self) -> DecompositionBaseLog {
        self.half_cbs_auto_base_log
    }

    pub fn half_cbs_auto_level(&self) -> DecompositionLevelCount {
        self.half_cbs_auto_level
    }

    pub fn half_cbs_fft_type_auto(&self) -> FftType {
        self.half_cbs_fft_type_auto
    }

    pub fn half_cbs_ss_base_log(&self) -> DecompositionBaseLog {
        self.half_cbs_ss_base_log
    }

    pub fn half_cbs_ss_level(&self) -> DecompositionLevelCount {
        self.half_cbs_ss_level
    }

    pub fn half_cbs_base_log(&self) -> DecompositionBaseLog {
        self.half_cbs_base_log
    }

    pub fn half_cbs_level(&self) -> DecompositionLevelCount {
        self.half_cbs_level
    }
}