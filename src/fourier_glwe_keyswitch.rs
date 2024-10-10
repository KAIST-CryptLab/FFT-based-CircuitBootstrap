use aligned_vec::{avec, ABox};
use dyn_stack::ReborrowMut;
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::c64,
};

use crate::{
    convert_lwe_to_glwe_const, fourier_glev_ciphertext::*, fourier_glwe_ciphertext::*, fourier_poly_mult_and_add, glev_ciphertext::*, GlweKeyswitchKey
};

#[derive(Debug, Clone, Copy)]
pub enum FftType {
    Vanilla,
    Split(usize),
    Split16,
}

impl FftType {
    pub fn num_split(&self) -> usize {
        match self {
            FftType::Vanilla => 1,
            FftType::Split(_) => 2,
            FftType::Split16 => 4,
        }
    }

    pub fn split_base_log(&self) -> usize {
        match self {
            FftType::Vanilla => 64,
            FftType::Split(b) => *b,
            FftType::Split16 => 16,
        }
    }
}

pub struct FourierGlweKeyswitchKey<C: Container<Element = c64>>
{
    data: C,
    input_glwe_size: GlweSize,
    output_glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    fft_type: FftType,
}

impl<C: Container<Element = c64>> FourierGlweKeyswitchKey<C> {
    pub fn from_container(
        container: C,
        input_glwe_size: GlweSize,
        output_glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        fft_type: FftType,
    ) -> Self {
        let fourier_glev_elem_count = input_glwe_size.to_glwe_dimension().0
            * output_glwe_size.0
            * polynomial_size.to_fourier_polynomial_size().0
            * decomp_level_count.0;
        assert_eq!(
            container.container_len(),
            fourier_glev_elem_count * fft_type.num_split(),
        );

        Self {
            data: container,
            input_glwe_size: input_glwe_size,
            output_glwe_size: output_glwe_size,
            polynomial_size: polynomial_size,
            decomp_base_log: decomp_base_log,
            decomp_level_count: decomp_level_count,
            fft_type,
        }
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn input_glwe_size(&self) -> GlweSize {
        self.input_glwe_size
    }

    pub fn output_glwe_size(&self) -> GlweSize {
        self.output_glwe_size
    }

    pub fn decomp_base_log(&self) -> DecompositionBaseLog {
        self.decomp_base_log
    }

    pub fn decomp_level_count(&self) -> DecompositionLevelCount {
        self.decomp_level_count
    }

    pub fn fft_type(&self) -> FftType {
        self.fft_type
    }

    pub fn as_fourier_glev_ciphertext_list(&self) -> FourierGlevCiphertextListView {
        FourierGlevCiphertextList::from_container(
            self.data.as_ref(),
            self.output_glwe_size,
            self.polynomial_size,
            self.decomp_base_log,
            self.decomp_level_count,
        )
    }
}

impl<C: ContainerMut<Element = c64>> FourierGlweKeyswitchKey<C> {
    pub fn as_mut_fourier_glev_ciphertext_list(&mut self) -> FourierGlevCiphertextListMutView {
        FourierGlevCiphertextList::from_container(
            self.data.as_mut(),
            self.output_glwe_size,
            self.polynomial_size,
            self.decomp_base_log,
            self.decomp_level_count,
        )
    }
}

pub type FourierGlweKeyswitchKeyOwned = FourierGlweKeyswitchKey<ABox<[c64]>>;

impl FourierGlweKeyswitchKeyOwned {
    pub fn new(
        input_glwe_size: GlweSize,
        output_glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        fft_type: FftType,
    ) -> Self {
        let count = input_glwe_size.to_glwe_dimension().0
            * output_glwe_size.0
            * polynomial_size.to_fourier_polynomial_size().0
            * decomp_level_count.0
            * fft_type.num_split();

        Self {
            data: avec![c64::default(); count].into_boxed_slice(),
            input_glwe_size: input_glwe_size,
            output_glwe_size: output_glwe_size,
            polynomial_size: polynomial_size,
            decomp_base_log: decomp_base_log,
            decomp_level_count: decomp_level_count,
            fft_type: fft_type,
        }
    }
}

pub fn convert_standard_glwe_keyswitch_key_to_fourier<Scalar, InputCont, OutputCont>(
    input_ksk: &GlweKeyswitchKey<InputCont>,
    output_ksk: &mut FourierGlweKeyswitchKey<OutputCont>,
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=c64>,
{
    assert_eq!(Scalar::BITS, 64, "current fourier GLWE ksk works on q = 2^64");
    assert_eq!(input_ksk.polynomial_size(), output_ksk.polynomial_size());
    assert_eq!(input_ksk.input_glwe_dimension().to_glwe_size(), output_ksk.input_glwe_size());
    assert_eq!(input_ksk.output_glwe_dimension().to_glwe_size(), output_ksk.output_glwe_size());
    assert_eq!(input_ksk.decomp_base_log(), output_ksk.decomp_base_log());
    assert_eq!(input_ksk.decomp_level_count(), output_ksk.decomp_level_count());

    let polynomial_size = output_ksk.polynomial_size();
    let output_glwe_size = output_ksk.output_glwe_size();
    let decomp_base_log = output_ksk.decomp_base_log();
    let decomp_level = output_ksk.decomp_level_count();
    let ciphertext_modulus = input_ksk.ciphertext_modulus();

    let fft_type = output_ksk.fft_type();
    let num_split = fft_type.num_split();
    let split_base_log = fft_type.split_base_log();

    for (input_glev, mut output_split_fourier_glev_list) in input_ksk.as_glev_ciphertext_list().iter()
        .zip(output_ksk.as_mut_fourier_glev_ciphertext_list().chunks_exact_mut(num_split))
    {
        for (k, mut output_split_fourier_glev) in output_split_fourier_glev_list.iter_mut().enumerate() {
            let mut input_split_glev = GlevCiphertext::new(Scalar::ZERO, output_glwe_size, polynomial_size, decomp_base_log, decomp_level, ciphertext_modulus);

            for (src, dst) in input_glev.as_ref().iter()
                .zip(input_split_glev.as_mut().iter_mut())
            {
                match fft_type {
                    FftType::Vanilla => {
                        *dst = *src;
                    }
                    FftType::Split(_) => {
                        let (shift_up_bit, shift_down_bit) = if k == 0 {
                            (Scalar::BITS - split_base_log, Scalar::BITS - split_base_log)
                        } else {
                            (0, split_base_log)
                        };
                        *dst = ((*src) << shift_up_bit) >> shift_down_bit;
                    }
                    FftType::Split16 => {
                        let shift_up_bit = split_base_log * (num_split - (k + 1));
                        let shift_down_bit = split_base_log * (num_split - 1);
                        *dst = ((*src) << shift_up_bit) >> shift_down_bit;
                    }
                }
            }

            convert_standard_glev_ciphertext_to_fourier(&input_split_glev, &mut output_split_fourier_glev);
        }
    }
}

pub fn keyswitch_glwe_ciphertext<Scalar, KSKeyCont, InputCont, OutputCont>(
    glwe_keyswitch_key: &FourierGlweKeyswitchKey<KSKeyCont>,
    input: &GlweCiphertext<InputCont>,
    output: &mut GlweCiphertext<OutputCont>,
) where
    Scalar: UnsignedTorus,
    KSKeyCont: Container<Element=c64>,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(
        glwe_keyswitch_key.input_glwe_size(),
        input.glwe_size(),
    );
    assert_eq!(
        glwe_keyswitch_key.output_glwe_size(),
        output.glwe_size(),
    );
    assert_eq!(
        glwe_keyswitch_key.polynomial_size(),
        input.polynomial_size(),
    );
    assert_eq!(
        glwe_keyswitch_key.polynomial_size(),
        output.polynomial_size(),
    );
    assert_eq!(
        input.ciphertext_modulus(),
        output.ciphertext_modulus(),
    );

    let polynomial_size = glwe_keyswitch_key.polynomial_size();
    let output_glwe_size = glwe_keyswitch_key.output_glwe_size();
    let decomp_base_log = glwe_keyswitch_key.decomp_base_log();
    let decomp_level = glwe_keyswitch_key.decomp_level_count();
    let ciphertext_modulus = input.ciphertext_modulus();

    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();

    let mut buffers = ComputationBuffers::new();
    buffers.resize(
        fft.backward_scratch()
        .unwrap()
        .unaligned_bytes_required(),
    );
    let mut stack = buffers.stack();

    output.as_mut().fill(Scalar::ZERO);
    output.get_mut_body().as_mut().clone_from_slice(input.get_body().as_ref());

    let decomposer = SignedDecomposer::new(decomp_base_log, decomp_level);
    let fft_type = glwe_keyswitch_key.fft_type();
    let num_split = fft_type.num_split();
    let split_base_log = fft_type.split_base_log();

    let mut buffer_fourier_glwe_list = FourierGlweCiphertextList::new(output_glwe_size, polynomial_size, FourierGlweCiphertextCount(num_split));

    let input_mask = input.get_mask();
    for (input_mask_poly, fourier_glev_split_list) in input_mask.as_polynomial_list().iter()
        .zip(glwe_keyswitch_key.as_fourier_glev_ciphertext_list().chunks_exact(num_split))
    {
        let mut input_decomp_poly_list = PolynomialList::new(
            Scalar::ZERO,
            polynomial_size,
            PolynomialCount(decomp_level.0),
        );

        for (i, val) in input_mask_poly.iter().enumerate() {
            let decomposition_iter = decomposer.decompose(*val);

            for (j, decomp_val) in decomposition_iter.into_iter().enumerate() {
                *input_decomp_poly_list.get_mut(j).as_mut().get_mut(i).unwrap() = decomp_val.value();
            }
        }

        let mut fourier_input_decomp_poly_list = FourierPolynomialList {
            data: avec![
                c64::default();
                polynomial_size.to_fourier_polynomial_size().0
                    * decomp_level.0
            ].into_boxed_slice(),
            polynomial_size: polynomial_size,
        };

        for (decomp_poly, mut fourier_decomp_poly) in input_decomp_poly_list.iter()
            .zip(fourier_input_decomp_poly_list.iter_mut())
        {
            fft.forward_as_integer(
                fourier_decomp_poly.as_mut_view(),
                decomp_poly.as_view(),
                stack.rb_mut(),
            );
        }

        for (mut buffer_fourier_glwe, fourier_glev_split) in buffer_fourier_glwe_list.iter_mut()
            .zip(fourier_glev_split_list.iter())
        {
            for (fourier_decomp_poly, fourier_glwe) in fourier_input_decomp_poly_list.iter_mut()
                .zip(fourier_glev_split.as_fourier_glwe_ciphertext_list().iter().rev())
            {
                for (mut buffer_poly, fourier_poly) in buffer_fourier_glwe.as_mut_fourier_polynomial_list().iter_mut()
                    .zip(fourier_glwe.as_fourier_polynomial_list().iter())
                {
                    fourier_poly_mult_and_add(&mut buffer_poly, &fourier_decomp_poly, &fourier_poly);
                }
            }
        }
    }

    let mut buffer_glwe_list = GlweCiphertextList::new(Scalar::ZERO, output_glwe_size, polynomial_size, GlweCiphertextCount(num_split), ciphertext_modulus);
    for (k, (mut buffer_glwe, buffer_fourier_glwe)) in buffer_glwe_list.iter_mut()
        .zip(buffer_fourier_glwe_list.iter())
        .enumerate()
    {
        for (mut buffer_poly, buffer_fourier_poly) in buffer_glwe.as_mut_polynomial_list().iter_mut()
            .zip(buffer_fourier_glwe.as_fourier_polynomial_list().iter())
        {
            fft.backward_as_torus(buffer_poly.as_mut_view(), buffer_fourier_poly.as_view(), stack.rb_mut());
        }

        let log_scaling = match fft_type {
            FftType::Vanilla => 0,
            FftType::Split(_) => if k == 0 {0} else {split_base_log},
            FftType::Split16 => k * split_base_log,
        };
        glwe_ciphertext_cleartext_mul_assign(&mut buffer_glwe, Cleartext(Scalar::ONE << log_scaling));
        glwe_ciphertext_add_assign(output, &buffer_glwe);
    }
}

pub fn keyswitch_lwe_ciphertext_by_glwe_keyswitch<Scalar, InputCont, OutputCont, KSKeyCont>(
    input: &LweCiphertext<InputCont>,
    output: &mut LweCiphertext<OutputCont>,
    glwe_keyswitch_key: &FourierGlweKeyswitchKey<KSKeyCont>
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element = Scalar>,
    OutputCont: ContainerMut<Element = Scalar>,
    KSKeyCont: Container<Element = c64>,
{
    assert_eq!(input.ciphertext_modulus(), output.ciphertext_modulus());
    let ciphertext_modulus = input.ciphertext_modulus();

    let polynomial_size = glwe_keyswitch_key.polynomial_size();

    let input_lwe_dimension = input.lwe_size().to_lwe_dimension();
    let output_lwe_dimension = output.lwe_size().to_lwe_dimension();

    assert_eq!(input_lwe_dimension.0 % polynomial_size.0, 0);
    assert_eq!(output_lwe_dimension.0 % polynomial_size.0, 0);

    let input_glwe_dimension = GlweDimension(input_lwe_dimension.0 / polynomial_size.0);
    let input_glwe_size = input_glwe_dimension.to_glwe_size();
    let output_glwe_dimension = GlweDimension(output_lwe_dimension.0 / polynomial_size.0);
    let output_glwe_size = output_glwe_dimension.to_glwe_size();

    assert_eq!(glwe_keyswitch_key.input_glwe_size(), input_glwe_size);
    assert_eq!(glwe_keyswitch_key.output_glwe_size(), output_glwe_size);

    let mut input_buf = GlweCiphertext::new(Scalar::ZERO, input_glwe_size, polynomial_size, ciphertext_modulus);
    let mut output_buf = GlweCiphertext::new(Scalar::ZERO, output_glwe_size, polynomial_size, ciphertext_modulus);

    convert_lwe_to_glwe_const(&input, &mut input_buf);
    keyswitch_glwe_ciphertext(glwe_keyswitch_key, &input_buf, &mut output_buf);
    extract_lwe_sample_from_glwe_ciphertext(&output_buf, output, MonomialDegree(0));
}
