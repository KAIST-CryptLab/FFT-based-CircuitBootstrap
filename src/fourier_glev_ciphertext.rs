use aligned_vec::{avec, AVec};
use dyn_stack::ReborrowMut;
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::c64,
};

use crate::{FourierGlweCiphertextList, FourierGlweCiphertextListMutView, FourierGlweCiphertextListView, GlevCiphertext};

pub struct FourierGlevCiphertext<C: Container<Element = c64>> {
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomposition_base_log: DecompositionBaseLog,
    decomposition_level_count: DecompositionLevelCount,
}

impl<C: Container<Element = c64>> AsRef<[c64]> for FourierGlevCiphertext<C> {
    fn as_ref(&self) -> &[c64] {
        self.data.as_ref()
    }
}

impl<C: ContainerMut<Element = c64>> AsMut<[c64]> for FourierGlevCiphertext<C> {
    fn as_mut(&mut self) -> &mut [c64] {
        self.data.as_mut()
    }
}

pub fn fourier_glev_ciphertext_size(glwe_size: GlweSize, polynomial_size: PolynomialSize, decomposition_level_count: DecompositionLevelCount) -> usize {
    glwe_size.0
        * polynomial_size.to_fourier_polynomial_size().0
        * decomposition_level_count.0
}

/// A [`FourierGlevCiphertext`] owning the memory for its own storage.
pub type FourierGlevCiphertextOwned = FourierGlevCiphertext<AVec<c64>>;
/// A [`FourierGlevCiphertext`] immutably borrowing memory for its own storage.
pub type FourierGlevCiphertextView<'data> = FourierGlevCiphertext<&'data [c64]>;
/// A [`FourierGlevCiphertext`] mutably borrowing memory for its own storage.
pub type FourierGlevCiphertextMutView<'data> = FourierGlevCiphertext<&'data mut [c64]>;

impl<C: Container<Element = c64>> FourierGlevCiphertext<C>
{
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomposition_base_log: DecompositionBaseLog,
        decomposition_level_count: DecompositionLevelCount,
    ) -> FourierGlevCiphertext<C> {
        let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
        assert_eq!(
            container.container_len(),
            glwe_size.0 * fourier_poly_size * decomposition_level_count.0
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
            decomposition_base_log: decomposition_base_log,
            decomposition_level_count: decomposition_level_count,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn decomposition_base_log(&self) -> DecompositionBaseLog {
        self.decomposition_base_log
    }

    pub fn decomposition_level_count(&self) -> DecompositionLevelCount {
        self.decomposition_level_count
    }

    pub fn as_fourier_glwe_ciphertext_list(&self) -> FourierGlweCiphertextListView {
        FourierGlweCiphertextList::from_container(
            self.data.as_ref(),
            self.glwe_size,
            self.polynomial_size,
        )
    }
}

impl<C: ContainerMut<Element = c64>> FourierGlevCiphertext<C> {
    pub fn as_mut_fourier_glwe_ciphertext_list(&mut self) -> FourierGlweCiphertextListMutView {
        FourierGlweCiphertextList::from_container(
            self.data.as_mut(),
            self.glwe_size,
            self.polynomial_size,
        )
    }
}

impl FourierGlevCiphertextOwned {
    pub fn new(
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomposition_base_log: DecompositionBaseLog,
        decomposition_level_count: DecompositionLevelCount,
    ) -> FourierGlevCiphertextOwned {
        Self::from_container(
            avec![
                c64::default();
                glwe_size.0
                    * polynomial_size.to_fourier_polynomial_size().0
                    * decomposition_level_count.0
            ],
            glwe_size,
            polynomial_size,
            decomposition_base_log,
            decomposition_level_count,
        )
    }
}

pub fn convert_standard_glev_ciphertext_to_fourier<Scalar, InputCont, OutputCont>(
    standard: &GlevCiphertext<InputCont>,
    fourier: &mut FourierGlevCiphertext<OutputCont>,
) where
    Scalar: UnsignedTorus,
    InputCont: Container<Element = Scalar>,
    OutputCont: ContainerMut<Element = c64>,
{
    assert_eq!(standard.glwe_size(), fourier.glwe_size());
    assert_eq!(standard.polynomial_size(), fourier.polynomial_size());
    assert_eq!(standard.decomposition_base_log(), fourier.decomposition_base_log());
    assert_eq!(standard.decomposition_level_count(), fourier.decomposition_level_count());

    let polynomial_size = standard.polynomial_size();
    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();

    let mut buffers = ComputationBuffers::new();
    buffers.resize(
        fft.forward_scratch()
        .unwrap()
        .unaligned_bytes_required(),
    );
    let mut stack = buffers.stack();

    for (glwe, mut fourier_glwe) in standard.as_glwe_ciphertext_list().iter()
        .zip(fourier.as_mut_fourier_glwe_ciphertext_list().iter_mut())
    {
        for (poly, mut fourier_poly) in glwe.as_polynomial_list().iter()
            .zip(fourier_glwe.as_mut_fourier_polynomial_list().iter_mut())
        {
            fft.forward_as_torus(fourier_poly.as_mut_view(), poly.as_view(), stack.rb_mut());
        }
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`FourierGlevCiphertext`] entities.
#[derive(Clone, Copy)]
pub struct FourierGlevCiphertextCreationMetadata(
    pub GlweSize,
    pub PolynomialSize,
    pub DecompositionBaseLog,
    pub DecompositionLevelCount,
);

impl<C: Container<Element = c64>> CreateFrom<C> for FourierGlevCiphertext<C> {
    type Metadata = FourierGlevCiphertextCreationMetadata;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let FourierGlevCiphertextCreationMetadata(glwe_size, polynomial_size, decomposition_base_log, decomposition_level_count) = meta;
        Self::from_container(from, glwe_size, polynomial_size, decomposition_base_log, decomposition_level_count)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct FourierGlevCiphertextCount(pub usize);



pub struct FourierGlevCiphertextList<C: Container<Element = c64>> {
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomposition_base_log: DecompositionBaseLog,
    decomposition_level_count: DecompositionLevelCount,
}

impl<C: Container<Element = c64>> AsRef<[c64]> for FourierGlevCiphertextList<C> {
    fn as_ref(&self) -> &[c64] {
        self.data.as_ref()
    }
}

impl<C: ContainerMut<Element = c64>> AsMut<[c64]> for FourierGlevCiphertextList<C> {
    fn as_mut(&mut self) -> &mut [c64] {
        self.data.as_mut()
    }
}


/// A [`FourierGlevCiphertextList`] owning the memory for its own storage.
pub type FourierGlevCiphertextListOwned = FourierGlevCiphertextList<AVec<c64>>;
/// A [`FourierGlevCiphertext`] immutably borrowing memory for its own storage.
pub type FourierGlevCiphertextListView<'data> = FourierGlevCiphertextList<&'data [c64]>;
/// A [`FourierGlevCiphertext`] mutably borrowing memory for its own storage.
pub type FourierGlevCiphertextListMutView<'data> = FourierGlevCiphertextList<&'data mut [c64]>;

impl<C: Container<Element = c64>> FourierGlevCiphertextList<C> {
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomposition_base_log: DecompositionBaseLog,
        decomposition_level_count: DecompositionLevelCount,
    ) -> FourierGlevCiphertextList<C> {
        let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
        assert_eq!(
            container.container_len() % (glwe_size.0 * fourier_poly_size),
            0,
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
            decomposition_base_log: decomposition_base_log,
            decomposition_level_count: decomposition_level_count,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn decomposition_base_log(&self) -> DecompositionBaseLog {
        self.decomposition_base_log
    }

    pub fn decomposition_level_count(&self) -> DecompositionLevelCount {
        self.decomposition_level_count
    }

    pub fn fourier_glev_ciphertext_count(&self) -> FourierGlevCiphertextCount {
        let fourier_poly_size = self.polynomial_size.to_fourier_polynomial_size().0;
        let count = self.data.container_len() / (self.glwe_size.0 * fourier_poly_size * self.decomposition_level_count.0);

        FourierGlevCiphertextCount(count)
    }
}

impl FourierGlevCiphertextListOwned {
    pub fn new(
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomposition_base_log: DecompositionBaseLog,
        decomposition_level_count: DecompositionLevelCount,
        fourier_glev_ciphertext_count: FourierGlevCiphertextCount,
    ) -> FourierGlevCiphertextListOwned {
        Self::from_container(
            avec![
                c64::default();
                fourier_glev_ciphertext_count.0
                    * glwe_size.0
                    * polynomial_size.to_fourier_polynomial_size().0
                    * decomposition_level_count.0
            ],
            glwe_size,
            polynomial_size,
            decomposition_base_log,
            decomposition_level_count,
        )
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`FourierGlevCiphertextList`] entities.
#[derive(Clone, Copy)]
pub struct FourierGlevCiphertextListCreationMetadata(
    pub GlweSize,
    pub PolynomialSize,
    pub DecompositionBaseLog,
    pub DecompositionLevelCount,
);

impl<C: Container<Element = c64>> CreateFrom<C>
    for FourierGlevCiphertextList<C>
{
    type Metadata = FourierGlevCiphertextListCreationMetadata;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let FourierGlevCiphertextListCreationMetadata(glwe_size, polynomial_size, decomp_base_log, decomp_level_count) =
            meta;
        Self::from_container(from, glwe_size, polynomial_size, decomp_base_log, decomp_level_count)
    }
}

impl<C: Container<Element = c64>> ContiguousEntityContainer
    for FourierGlevCiphertextList<C>
{
    type Element = C::Element;

    type EntityViewMetadata = FourierGlevCiphertextCreationMetadata;

    type EntityView<'this> = FourierGlevCiphertextView<'this>
    where
        Self: 'this;

    type SelfViewMetadata = FourierGlevCiphertextListCreationMetadata;

    type SelfView<'this> = FourierGlevCiphertextListView<'this>
    where
        Self: 'this;

    fn get_entity_view_creation_metadata(&self) -> Self::EntityViewMetadata {
        FourierGlevCiphertextCreationMetadata(self.glwe_size(), self.polynomial_size(), self.decomposition_base_log(), self.decomposition_level_count())
    }

    fn get_entity_view_pod_size(&self) -> usize {
        fourier_glev_ciphertext_size(self.glwe_size(), self.polynomial_size(), self.decomposition_level_count())
    }

    fn  get_self_view_creation_metadata(&self) -> Self::SelfViewMetadata {
        FourierGlevCiphertextListCreationMetadata(
            self.glwe_size(),
            self.polynomial_size(),
            self.decomposition_base_log(),
            self.decomposition_level_count(),
        )
    }
}

impl<C: ContainerMut<Element = c64>> ContiguousEntityContainerMut
    for FourierGlevCiphertextList<C>
{
    type EntityMutView<'this> = FourierGlevCiphertextMutView<'this>
    where
        Self: 'this;

    type SelfMutView<'this> = FourierGlevCiphertextListMutView<'this>
    where
        Self: 'this;
}
