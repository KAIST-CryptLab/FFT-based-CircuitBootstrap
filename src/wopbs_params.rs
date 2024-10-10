use tfhe::core_crypto::prelude::*;
use crate::FftType;

#[derive(Clone, Copy)]
pub struct WopbsParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    pfks_base_log: DecompositionBaseLog,
    pfks_level: DecompositionLevelCount,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
    message_size: usize,
}

#[derive(Clone, Copy)]
pub struct ImprovedWopbsParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    auto_base_log: DecompositionBaseLog,
    auto_level: DecompositionLevelCount,
    fft_type_auto: FftType,
    ss_base_log: DecompositionBaseLog,
    ss_level: DecompositionLevelCount,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
    message_size: usize,
}

#[derive(Clone, Copy)]
pub struct HighPrecImprovedWopbsParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    large_glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    large_glwe_modular_std_dev: StandardDev,
    pbs_base_log: DecompositionBaseLog,
    pbs_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    glwe_ds_to_large_base_log: DecompositionBaseLog,
    glwe_ds_to_large_level: DecompositionLevelCount,
    fft_type_to_large: FftType,
    auto_base_log: DecompositionBaseLog,
    auto_level: DecompositionLevelCount,
    fft_type_auto: FftType,
    glwe_ds_from_large_base_log: DecompositionBaseLog,
    glwe_ds_from_large_level: DecompositionLevelCount,
    fft_type_from_large: FftType,
    ss_base_log: DecompositionBaseLog,
    ss_level: DecompositionLevelCount,
    cbs_base_log: DecompositionBaseLog,
    cbs_level: DecompositionLevelCount,
    log_lut_count: LutCountLog,
    ciphertext_modulus: CiphertextModulus::<Scalar>,
    message_size: usize,
}

impl<Scalar: UnsignedInteger> WopbsParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        pbs_base_log: DecompositionBaseLog,
        pbs_level: DecompositionLevelCount,
        ks_base_log: DecompositionBaseLog,
        ks_level: DecompositionLevelCount,
        pfks_base_log: DecompositionBaseLog,
        pfks_level: DecompositionLevelCount,
        cbs_base_log: DecompositionBaseLog,
        cbs_level: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
        message_size: usize,
    ) -> Self {
        WopbsParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            ks_base_log,
            ks_level,
            pfks_base_log,
            pfks_level,
            cbs_base_log,
            cbs_level,
            ciphertext_modulus,
            message_size,
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

    pub fn glwe_dimension(&self) -> GlweDimension  {
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

    pub fn ks_base_log(&self) -> DecompositionBaseLog {
        self.ks_base_log
    }

    pub fn ks_level(&self) -> DecompositionLevelCount {
        self.ks_level
    }

    pub fn pfks_base_log(&self) -> DecompositionBaseLog {
        self.pfks_base_log
    }

    pub fn pfks_level(&self) -> DecompositionLevelCount {
        self.pfks_level
    }

    pub fn cbs_base_log(&self) -> DecompositionBaseLog {
        self.cbs_base_log
    }

    pub fn cbs_level(&self) -> DecompositionLevelCount {
        self.cbs_level
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus::<Scalar> {
        self.ciphertext_modulus
    }

    pub fn message_size(&self) -> usize {
        self.message_size
    }
}

impl<Scalar: UnsignedInteger> ImprovedWopbsParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        pbs_base_log: DecompositionBaseLog,
        pbs_level: DecompositionLevelCount,
        ks_base_log: DecompositionBaseLog,
        ks_level: DecompositionLevelCount,
        auto_base_log: DecompositionBaseLog,
        auto_level: DecompositionLevelCount,
        fft_type_auto: FftType,
        ss_base_log: DecompositionBaseLog,
        ss_level: DecompositionLevelCount,
        cbs_base_log: DecompositionBaseLog,
        cbs_level: DecompositionLevelCount,
        log_lut_count: LutCountLog,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
        message_size: usize,
    ) -> Self {
        ImprovedWopbsParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            ks_base_log,
            ks_level,
            auto_base_log,
            auto_level,
            fft_type_auto,
            ss_base_log,
            ss_level,
            cbs_base_log,
            cbs_level,
            log_lut_count,
            ciphertext_modulus,
            message_size,
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

    pub fn glwe_dimension(&self) -> GlweDimension  {
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

    pub fn ks_base_log(&self) -> DecompositionBaseLog {
        self.ks_base_log
    }

    pub fn ks_level(&self) -> DecompositionLevelCount {
        self.ks_level
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

    pub fn message_size(&self) -> usize {
        self.message_size
    }

}

impl<Scalar: UnsignedInteger> HighPrecImprovedWopbsParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        large_glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        large_glwe_modular_std_dev: StandardDev,
        pbs_base_log: DecompositionBaseLog,
        pbs_level: DecompositionLevelCount,
        ks_base_log: DecompositionBaseLog,
        ks_level: DecompositionLevelCount,
        glwe_ds_to_large_base_log: DecompositionBaseLog,
        glwe_ds_to_large_level: DecompositionLevelCount,
        fft_type_to_large: FftType,
        auto_base_log: DecompositionBaseLog,
        auto_level: DecompositionLevelCount,
        fft_type_auto: FftType,
        glwe_ds_from_large_base_log: DecompositionBaseLog,
        glwe_ds_from_large_level: DecompositionLevelCount,
        fft_type_from_large: FftType,
        ss_base_log: DecompositionBaseLog,
        ss_level: DecompositionLevelCount,
        cbs_base_log: DecompositionBaseLog,
        cbs_level: DecompositionLevelCount,
        log_lut_count: LutCountLog,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
        message_size: usize,
    ) -> Self {
        HighPrecImprovedWopbsParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            large_glwe_dimension,
            glwe_modular_std_dev,
            large_glwe_modular_std_dev,
            pbs_base_log,
            pbs_level,
            ks_base_log,
            ks_level,
            glwe_ds_to_large_base_log,
            glwe_ds_to_large_level,
            fft_type_to_large,
            auto_base_log,
            auto_level,
            fft_type_auto,
            glwe_ds_from_large_base_log,
            glwe_ds_from_large_level,
            fft_type_from_large,
            ss_base_log,
            ss_level,
            cbs_base_log,
            cbs_level,
            log_lut_count,
            ciphertext_modulus,
            message_size,
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

    pub fn glwe_dimension(&self) -> GlweDimension  {
        self.glwe_dimension
    }

    pub fn glwe_modular_std_dev(&self) -> StandardDev {
        self.glwe_modular_std_dev
    }

    pub fn large_glwe_dimension(&self) -> GlweDimension {
        self.large_glwe_dimension
    }

    pub fn large_glwe_modular_std_dev(&self) -> StandardDev {
        self.large_glwe_modular_std_dev
    }

    pub fn pbs_base_log(&self) -> DecompositionBaseLog {
        self.pbs_base_log
    }

    pub fn pbs_level(&self) -> DecompositionLevelCount {
        self.pbs_level
    }

    pub fn glwe_ds_to_large_base_log(&self) -> DecompositionBaseLog {
        self.glwe_ds_to_large_base_log
    }

    pub fn glwe_ds_to_large_level(&self) -> DecompositionLevelCount {
        self.glwe_ds_to_large_level
    }

    pub fn fft_type_to_large(&self) -> FftType {
        self.fft_type_to_large
    }

    pub fn glwe_ds_from_large_base_log(&self) -> DecompositionBaseLog {
        self.glwe_ds_from_large_base_log
    }

    pub fn glwe_ds_from_large_level(&self) -> DecompositionLevelCount {
        self.glwe_ds_from_large_level
    }

    pub fn fft_type_from_large(&self) -> FftType {
        self.fft_type_from_large
    }

    pub fn ks_base_log(&self) -> DecompositionBaseLog {
        self.ks_base_log
    }

    pub fn ks_level(&self) -> DecompositionLevelCount {
        self.ks_level
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

    pub fn message_size(&self) -> usize {
        self.message_size
    }
}
