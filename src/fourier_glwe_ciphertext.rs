use aligned_vec::{avec, AVec};
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::c64,
};

// Extension of tfhe::core_crypto::fft_impl::fft64::math::FourierPolynomialList
pub struct FourierPolynomialList<C: Container<Element = c64>> {
    pub data: C,
    pub polynomial_size: PolynomialSize,
}

pub type FourierPolynomialListView<'data> = FourierPolynomialList<&'data [c64]>;
pub type FourierPolynomialListMutView<'data> = FourierPolynomialList<&'data mut [c64]>;

impl<C: Container<Element = c64>> FourierPolynomialList<C> {
    pub fn polynomial_count(&self) -> PolynomialCount {
        PolynomialCount(
            self.data.container_len() / self.polynomial_size.to_fourier_polynomial_size().0
        )
    }

    pub fn iter(
        &self
    ) -> impl DoubleEndedIterator<Item = FourierPolynomial<&'_ [c64]>> {
        assert_eq!(
            self.data.container_len() % self.polynomial_size.to_fourier_polynomial_size().0,
            0,
        );
        self.data
            .as_ref()
            .chunks_exact(self.polynomial_size.to_fourier_polynomial_size().0)
            .map(move |slice| FourierPolynomial { data: slice })
    }
}

impl<C: ContainerMut<Element = c64>> FourierPolynomialList<C> {
    pub fn iter_mut(
        &mut self,
    ) -> impl DoubleEndedIterator<Item = FourierPolynomial<&'_ mut [c64]>> {
        assert_eq!(
            self.data.container_len() % self.polynomial_size.to_fourier_polynomial_size().0,
            0,
        );
        self.data
            .as_mut()
            .chunks_exact_mut(self.polynomial_size.to_fourier_polynomial_size().0)
            .map(move |slice| FourierPolynomial { data: slice })
    }
}

pub struct FourierGlweCiphertext<C: Container<Element = c64>> {
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
}

/// A [`FourierGlweCiphertext`] owning the memory for its own storage.
pub type FourierGlweCiphertextOwned = FourierGlweCiphertext<AVec<c64>>;
/// A [`FourierGlweCiphertext`] immutably borrowing memory for its own storage.
pub type FourierGlweCiphertextView<'data> = FourierGlweCiphertext<&'data [c64]>;
/// A [`FourierGlevCiphertext`] mutably borrowing memory for its own storage.
pub type FourierGlweCiphertextMutView<'data> = FourierGlweCiphertext<&'data mut [c64]>;

impl<C: Container<Element = c64>> AsRef<[c64]> for FourierGlweCiphertext<C> {
    fn as_ref(&self) -> &[c64] {
        self.data.as_ref()
    }
}

impl<C: ContainerMut<Element = c64>> AsMut<[c64]> for FourierGlweCiphertext<C> {
    fn as_mut(&mut self) -> &mut [c64] {
        self.data.as_mut()
    }
}

pub fn fourier_glwe_ciphertext_size(glwe_size: GlweSize, polynomial_size: PolynomialSize) -> usize {
    glwe_size.0 * polynomial_size.to_fourier_polynomial_size().0
}

impl<C: Container<Element = c64>> FourierGlweCiphertext<C>
{
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
    ) -> FourierGlweCiphertext<C> {
        let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
        assert_eq!(
            container.container_len(),
            glwe_size.0 * fourier_poly_size
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn as_fourier_polynomial_list(&self) -> FourierPolynomialListView {
        FourierPolynomialList {
            data: self.data.as_ref(),
            polynomial_size: self.polynomial_size,
        }
    }
}

impl<C: ContainerMut<Element = c64>> FourierGlweCiphertext<C>
{
    pub fn as_mut_fourier_polynomial_list(&mut self) -> FourierPolynomialListMutView {
        FourierPolynomialList {
            data: self.data.as_mut(),
            polynomial_size: self.polynomial_size,
        }
    }
}

impl FourierGlweCiphertextOwned {
    pub fn new(
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
    ) -> FourierGlweCiphertextOwned {
        let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
        Self::from_container(
            avec![
                c64::default();
                glwe_size.0  * fourier_poly_size
            ],
            glwe_size,
            polynomial_size,
        )
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`FourierGlweCiphertext`] entities.
#[derive(Clone, Copy)]
pub struct FourierGlweCiphertextCreationMetadata(
    pub GlweSize,
    pub PolynomialSize,
);

impl<C: Container<Element = c64>> CreateFrom<C> for FourierGlweCiphertext<C> {
    type Metadata = FourierGlweCiphertextCreationMetadata;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let FourierGlweCiphertextCreationMetadata(glwe_size, polynomial_size) = meta;
        Self::from_container(from, glwe_size, polynomial_size)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct FourierGlweCiphertextCount(pub usize);


pub struct FourierGlweCiphertextList<C: Container<Element = c64>>
{
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
}

impl<C: Container<Element = c64>> AsRef<[c64]> for FourierGlweCiphertextList<C> {
    fn as_ref(&self) -> &[c64] {
        self.data.as_ref()
    }
}

impl<C: ContainerMut<Element = c64>> AsMut<[c64]> for FourierGlweCiphertextList<C> {
    fn as_mut(&mut self) -> &mut [c64] {
        self.data.as_mut()
    }
}


/// A [`FourierGlweCiphertextList`] owning the memory for its own storage.
pub type FourierGlweCiphertextListOwned = FourierGlweCiphertextList<AVec<c64>>;
/// A [`FourierGlweCiphertext`] immutably borrowing memory for its own storage.
pub type FourierGlweCiphertextListView<'data> = FourierGlweCiphertextList<&'data [c64]>;
/// A [`FourierGlweCiphertext`] mutably borrowing memory for its own storage.
pub type FourierGlweCiphertextListMutView<'data> = FourierGlweCiphertextList<&'data mut [c64]>;

impl<C: Container<Element = c64>> FourierGlweCiphertextList<C> {
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
    ) -> FourierGlweCiphertextList<C> {
        let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
        assert_eq!(
            container.container_len() % (glwe_size.0 * fourier_poly_size),
            0,
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn fourier_glwe_ciphertext_count(&self) -> FourierGlweCiphertextCount {
        let fourier_poly_size = self.polynomial_size.to_fourier_polynomial_size().0;
        let count = self.data.container_len() / (self.glwe_size.0 * fourier_poly_size);

        FourierGlweCiphertextCount(count)
    }
}

impl FourierGlweCiphertextListOwned {
    pub fn new(
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        fourier_glwe_ciphertext_count: FourierGlweCiphertextCount,
    ) -> FourierGlweCiphertextListOwned {
        Self::from_container(
            avec![
                c64::default();
                fourier_glwe_ciphertext_count.0
                    * glwe_size.0
                    * polynomial_size.to_fourier_polynomial_size().0
            ],
            glwe_size,
            polynomial_size,
        )
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`FourierGlweCiphertextList`] entities.
#[derive(Clone, Copy)]
pub struct FourierGlweCiphertextListCreationMetadata(
    pub GlweSize,
    pub PolynomialSize,
    pub FourierGlweCiphertextCount,
);

impl<C: Container<Element = c64>> CreateFrom<C>
    for FourierGlweCiphertextList<C>
{
    type Metadata = FourierGlweCiphertextListCreationMetadata;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let FourierGlweCiphertextListCreationMetadata(glwe_size, polynomial_size, _fourier_glwe_ciphertext_count) = meta;
        Self::from_container(from, glwe_size, polynomial_size)
    }
}

impl<C: Container<Element = c64>> ContiguousEntityContainer
    for FourierGlweCiphertextList<C>
{
    type Element = C::Element;

    type EntityViewMetadata = FourierGlweCiphertextCreationMetadata;

    type EntityView<'this> = FourierGlweCiphertextView<'this>
    where
        C: 'this;

    type SelfViewMetadata = FourierGlweCiphertextListCreationMetadata;

    type SelfView<'this> = FourierGlweCiphertextListView<'this>
    where
        C: 'this;

    fn get_entity_view_creation_metadata(&self) -> Self::EntityViewMetadata {
        FourierGlweCiphertextCreationMetadata(self.glwe_size(), self.polynomial_size())
    }

    fn get_entity_view_pod_size(&self) -> usize {
        fourier_glwe_ciphertext_size(self.glwe_size(), self.polynomial_size())
    }

    fn get_self_view_creation_metadata(&self) -> Self::SelfViewMetadata {
        FourierGlweCiphertextListCreationMetadata(
            self.glwe_size(),
            self.polynomial_size(),
            self.fourier_glwe_ciphertext_count(),
        )
    }
}

impl<C: ContainerMut<Element = c64>> ContiguousEntityContainerMut
    for FourierGlweCiphertextList<C>
{
    type EntityMutView<'this> = FourierGlweCiphertextMutView<'this>
    where
        C: 'this;

    type SelfMutView<'this> = FourierGlweCiphertextListMutView<'this>
    where
        C: 'this;
}
