use std::collections::HashMap;
use aligned_vec::ABox;
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::c64,
};
use crate::{utils::*, glwe_keyswitch::*, fourier_glwe_keyswitch::*};

// The following codes generalize rlweExpand
// from https://github.com/KULeuven-COSIC/SortingHat
// to automorphism on arbitrary GLWE dimension
pub struct AutomorphKey<C: Container<Element=c64>> {
    ksk: FourierGlweKeyswitchKey<C>,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    glwe_dimension: GlweDimension,
    polynomial_size: PolynomialSize,
    auto_k: usize,
}

impl AutomorphKey<ABox<[c64]>> {
    pub fn allocate(
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        glwe_dimension: GlweDimension,
        polynomial_size: PolynomialSize,
        auto_k: usize,
        fft_type: FftType,
    ) -> Self {
        let glwe_size = glwe_dimension.to_glwe_size();
        let ksk = FourierGlweKeyswitchKey::new(glwe_size, glwe_size, polynomial_size, decomp_base_log, decomp_level_count, fft_type);
        AutomorphKey {
            ksk: ksk,
            decomp_base_log,
            decomp_level_count,
            glwe_dimension,
            polynomial_size,
            auto_k: auto_k,
        }
    }

    pub fn decomposition_base_log(&self) -> DecompositionBaseLog {
        self.decomp_base_log
    }

    pub fn decomposition_level_count(&self) -> DecompositionLevelCount {
        self.decomp_level_count
    }

    pub fn glwe_dimension(&self) -> GlweDimension {
        self.glwe_dimension
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    /// Fill this object with the appropriate key switching key
    /// that is used for the automorphism operation
    /// where after_key is {S_i(X)} and before_key is computed as {S_i(X^k)}.
    fn fill_with_automorph_key<Scalar: UnsignedTorus + Sync + Send, G: ByteRandomGenerator>(
        &mut self,
        before_key: &mut GlweSecretKeyOwned<Scalar>,
        after_key: &GlweSecretKeyOwned<Scalar>,
        k: usize,
        noise_parameters: impl DispersionParameter,
        generator: &mut EncryptionRandomGenerator<G>,
    ) {
        debug_assert!(self.glwe_dimension == before_key.glwe_dimension());
        debug_assert!(self.glwe_dimension == after_key.glwe_dimension());
        debug_assert!(self.polynomial_size == before_key.polynomial_size());
        debug_assert!(self.polynomial_size == after_key.polynomial_size());

        let mut before_poly_list = PolynomialList::new(
            Scalar::ZERO,
            self.polynomial_size,
            PolynomialCount(self.glwe_dimension.0),
        );
        for (mut before_poly, after_poly) in before_poly_list.iter_mut()
            .zip(after_key.as_polynomial_list().iter())
        {
            let out = eval_x_k(after_poly.as_view(), k);
            before_poly.as_mut().clone_from_slice(out.as_ref());
        }
        *before_key = GlweSecretKey::from_container(before_poly_list.into_container(), self.polynomial_size);

        self.fill_with_keyswitch_key(before_key, after_key, noise_parameters, generator);
        self.auto_k = k;
    }

    /// Fill this object with the appropriate keyswitching key
    /// that transforms ciphertexts under before_key to ciphertexts under after_key.
    fn fill_with_keyswitch_key<Scalar: UnsignedTorus + Sync + Send, G: ByteRandomGenerator>(
        &mut self,
        before_key: &GlweSecretKeyOwned<Scalar>,
        after_key: &GlweSecretKeyOwned<Scalar>,
        noise_parameters: impl DispersionParameter,
        generator: &mut EncryptionRandomGenerator<G>
    ) {
        debug_assert!(self.glwe_dimension == before_key.glwe_dimension());
        debug_assert!(self.glwe_dimension == after_key.glwe_dimension());
        debug_assert!(self.polynomial_size == before_key.polynomial_size());
        debug_assert!(self.polynomial_size == after_key.polynomial_size());

        let decomp_level_count = self.decomp_level_count;
        let decomp_base_log = self.decomp_base_log;
        let ciphertext_modulus = CiphertextModulus::new_native();

        let standard_ksk = allocate_and_generate_new_glwe_keyswitch_key(
            before_key,
            after_key,
            decomp_base_log,
            decomp_level_count,
            noise_parameters,
            ciphertext_modulus,
            generator,
        );
        convert_standard_glwe_keyswitch_key_to_fourier(&standard_ksk, &mut self.ksk);
    }

    fn keyswitch_ciphertext<Scalar, InputCont, OutputCont>(
        &self,
        after: &mut GlweCiphertext<OutputCont>,
        before: &GlweCiphertext<InputCont>,
    ) where
        Scalar: UnsignedTorus + Sync + Send,
        InputCont: Container<Element=Scalar>,
        OutputCont: ContainerMut<Element=Scalar>,
    {
        keyswitch_glwe_ciphertext(&self.ksk, before, after);
    }

    pub fn auto<Scalar, InputCont, OutputCont>(
        &self,
        after: &mut GlweCiphertext<OutputCont>,
        before: &GlweCiphertext<InputCont>,
    ) where
        Scalar: UnsignedTorus + Sync + Send,
        InputCont: Container<Element=Scalar>,
        OutputCont: ContainerMut<Element=Scalar>,
    {
        let mut before_power = GlweCiphertextOwned::new(Scalar::ZERO, before.glwe_size(), before.polynomial_size(), before.ciphertext_modulus());
        for (mut poly_power, poly) in before_power.as_mut_polynomial_list().iter_mut().zip(before.as_polynomial_list().iter()) {
            poly_power.as_mut().clone_from_slice(eval_x_k(poly, self.auto_k).as_ref());
        }

        self.keyswitch_ciphertext(after, &before_power);
    }
}

pub fn gen_all_auto_keys<Scalar, G>(
    decomp_base_log: DecompositionBaseLog,
    decomp_level: DecompositionLevelCount,
    fft_type: FftType,
    glwe_secret_key: &GlweSecretKeyOwned<Scalar>,
    noise_parameters: impl DispersionParameter,
    generator: &mut EncryptionRandomGenerator<G>,
) -> HashMap<usize, AutomorphKey<ABox<[c64]>>>
where
    Scalar: UnsignedTorus + Sync + Send,
    G: ByteRandomGenerator,
{
    let glwe_dimension = glwe_secret_key.glwe_dimension();
    let polynomial_size = glwe_secret_key.polynomial_size();

    let mut hm = HashMap::new();
    for i in 1..=(polynomial_size.0).ilog2() as usize {
        let k = polynomial_size.0 / (1 << (i - 1)) + 1;
        let mut glwe_ksk = AutomorphKey::allocate(decomp_base_log, decomp_level, glwe_dimension, polynomial_size, i, fft_type);
        let mut before_key = glwe_secret_key.clone();

        glwe_ksk.fill_with_automorph_key(&mut before_key, &glwe_secret_key, k, noise_parameters, generator);
        hm.insert(k, glwe_ksk);
    }

    hm
}

pub fn trace<Scalar, Cont>(
    glwe_in: &GlweCiphertext<Cont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
) -> GlweCiphertextOwned<Scalar>
where
    Scalar: UnsignedTorus + Sync + Send,
    Cont: Container<Element=Scalar>,
{
    let mut out = GlweCiphertext::new(Scalar::ZERO, glwe_in.glwe_size(), glwe_in.polynomial_size(), glwe_in.ciphertext_modulus());
    glwe_ciphertext_clone_from(&mut out, glwe_in);
    trace_assign(&mut out, auto_keys);

    out
}

pub fn trace_assign<Scalar, ContMut>(
    glwe_in: &mut GlweCiphertext<ContMut>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
) where
    Scalar: UnsignedTorus + Sync + Send,
    ContMut: ContainerMut<Element=Scalar>,
{
    trace_partial_assign(glwe_in, auto_keys, 1);
}

pub fn trace_partial_assign<Scalar, Cont>(
    input: &mut GlweCiphertext<Cont>,
    auto_keys: &HashMap<usize, AutomorphKey<ABox<[c64]>>>,
    n: usize,
) where
    Scalar: UnsignedTorus,
    Cont: ContainerMut<Element=Scalar>,
{
    let glwe_size = input.glwe_size();
    let polynomial_size = input.polynomial_size();
    let ciphertext_modulus = input.ciphertext_modulus();

    assert!(polynomial_size.0 % n == 0);

    let mut buf = GlweCiphertextOwned::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
    let mut out: GlweCiphertext<Vec<Scalar>> = GlweCiphertext::new(Scalar::ZERO, glwe_size, polynomial_size, ciphertext_modulus);
    glwe_ciphertext_clone_from(&mut out, input);

    let log_polynomial_size = polynomial_size.0.ilog2() as usize;
    let log_n = n.ilog2() as usize;
    for i in 1..=(log_polynomial_size - log_n) {
        let k = polynomial_size.0 / (1 << (i - 1)) + 1;
        let auto_key = auto_keys.get(&k).unwrap();
        auto_key.auto(&mut buf, &out);
        glwe_ciphertext_add_assign(&mut out, &buf);
    }

    glwe_ciphertext_clone_from(input, &out);
}
