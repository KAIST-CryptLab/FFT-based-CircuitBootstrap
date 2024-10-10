use tfhe::core_crypto::prelude::*;
use crate::FftType;

#[derive(Clone, Copy)]
pub struct PkskConvParam<Scalar: UnsignedInteger> {
    lwe_dimension: LweDimension,
    lwe_modular_std_dev: StandardDev,
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    pksk_base_log: DecompositionBaseLog,
    pksk_level: DecompositionLevelCount,
    ks_base_log: DecompositionBaseLog,
    ks_level: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus<Scalar>,
}

#[derive(Clone, Copy)]
pub struct AutoConvParam<Scalar: UnsignedInteger> {
    polynomial_size: PolynomialSize,
    glwe_dimension: GlweDimension,
    glwe_modular_std_dev: StandardDev,
    auto_base_log: DecompositionBaseLog,
    auto_level: DecompositionLevelCount,
    fft_type: FftType,
    ciphertext_modulus: CiphertextModulus<Scalar>,
}

impl<Scalar: UnsignedInteger> PkskConvParam<Scalar> {
    pub fn new(
        lwe_dimension: LweDimension,
        lwe_modular_std_dev: StandardDev,
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        pksk_base_log: DecompositionBaseLog,
        pksk_level: DecompositionLevelCount,
        ks_base_log: DecompositionBaseLog,
        ks_level: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus<Scalar>,
    ) -> Self {
        PkskConvParam {
            lwe_dimension,
            lwe_modular_std_dev,
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            pksk_base_log,
            pksk_level,
            ks_base_log,
            ks_level,
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

    pub fn pksk_base_log(&self) -> DecompositionBaseLog {
        self.pksk_base_log
    }

    pub fn pksk_level(&self) -> DecompositionLevelCount {
        self.pksk_level
    }

    pub fn ks_base_log(&self) -> DecompositionBaseLog {
        self.ks_base_log
    }

    pub fn ks_level(&self) -> DecompositionLevelCount {
        self.ks_level
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus<Scalar> {
        self.ciphertext_modulus
    }
}

impl<Scalar: UnsignedInteger> AutoConvParam<Scalar> {
    pub fn new(
        polynomial_size: PolynomialSize,
        glwe_dimension: GlweDimension,
        glwe_modular_std_dev: StandardDev,
        auto_base_log: DecompositionBaseLog,
        auto_level: DecompositionLevelCount,
        fft_type: FftType,
        ciphertext_modulus: CiphertextModulus<Scalar>,
    ) -> Self {
        AutoConvParam {
            polynomial_size,
            glwe_dimension,
            glwe_modular_std_dev,
            auto_base_log,
            auto_level,
            fft_type,
            ciphertext_modulus,
        }
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

    pub fn auto_base_log(&self) -> DecompositionBaseLog {
        self.auto_base_log
    }

    pub fn auto_level(&self) -> DecompositionLevelCount {
        self.auto_level
    }

    pub fn fft_type(&self) -> FftType {
        self.fft_type
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus<Scalar> {
        self.ciphertext_modulus
    }
}
