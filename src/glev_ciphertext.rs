use tfhe::core_crypto::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct GlevCiphertext<C: Container>
where
    C::Element: UnsignedInteger,
{
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus::<C::Element>,
}

/// A [`GlevCiphertext`] owning the memory for its own storage.
pub type GlevCiphertextOwned<Scalar> = GlevCiphertext<Vec<Scalar>>;
/// A [`GlevCiphertext`] immutably borrowing memory for its own storage.
pub type GlevCiphertextView<'data, Scalar> = GlevCiphertext<&'data [Scalar]>;
/// A [`GlevCiphertext`] mutably borrowing memory for its own storage.
pub type GlevCiphertextMutView<'data, Scalar> = GlevCiphertext<&'data mut [Scalar]>;

impl<T: UnsignedInteger, C: Container<Element=T>> AsRef<[T]> for GlevCiphertext<C> {
    fn as_ref(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<T: UnsignedInteger, C: ContainerMut<Element=T>> AsMut<[T]> for GlevCiphertext<C> {
    fn as_mut(&mut self) -> &mut [T] {
        self.data.as_mut()
    }
}

impl<Scalar: UnsignedInteger, C: Container<Element=Scalar>> GlevCiphertext<C>
{
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
    ) -> GlevCiphertext<C> {
        assert_eq!(
            container.container_len(),
            glwe_size.0 * polynomial_size.0 * decomp_level_count.0,
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
            decomp_base_log: decomp_base_log,
            decomp_level_count: decomp_level_count,
            ciphertext_modulus: ciphertext_modulus,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn decomposition_base_log(&self) -> DecompositionBaseLog {
        self.decomp_base_log
    }

    pub fn decomposition_level_count(&self) -> DecompositionLevelCount {
        self.decomp_level_count
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus::<Scalar> {
        self.ciphertext_modulus
    }

    pub fn into_container(self) -> C {
        self.data
    }

    pub fn as_glwe_ciphertext_list(&self) -> GlweCiphertextList<&'_ [Scalar]> {
        GlweCiphertextList::from_container(
            self.data.as_ref(),
            self.glwe_size,
            self.polynomial_size,
            self.ciphertext_modulus,
        )
    }
}

impl<Scalar: UnsignedInteger, C: ContainerMut<Element=Scalar>> GlevCiphertext<C>
{
    pub fn as_mut_glwe_ciphertext_list(&mut self) -> GlweCiphertextList<&'_ mut [Scalar]>  {
        GlweCiphertextList::from_container(
            self.data.as_mut(),
            self.glwe_size,
            self.polynomial_size,
            self.ciphertext_modulus,
        )
    }
}

impl<Scalar: UnsignedInteger> GlevCiphertextOwned<Scalar> {
    pub fn new(
        fill_with: Scalar,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
    ) -> GlevCiphertextOwned<Scalar> {
        Self::from_container(
            vec![fill_with; glwe_size.0 * polynomial_size.0 * decomp_level_count.0],
            glwe_size,
            polynomial_size,
            decomp_base_log,
            decomp_level_count,
            ciphertext_modulus,
        )
    }
}

pub fn glev_ciphertext_size(
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomp_level_count: DecompositionLevelCount,
) -> usize {
    glwe_size.0 * polynomial_size.0 * decomp_level_count.0
}

pub fn encrypt_glev_ciphertext<Scalar, KeyCont, InputCont, OutputCont, Gen>(
    glwe_secret_key: &GlweSecretKey<KeyCont>,
    output: &mut GlevCiphertext<OutputCont>,
    pt: &PlaintextList<InputCont>,
    noise_parameters: impl DispersionParameter,
    generator: &mut EncryptionRandomGenerator<Gen>,
) where
    Scalar: UnsignedTorus,
    KeyCont: Container<Element=Scalar>,
    InputCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
    Gen: ByteRandomGenerator,
{
    assert_eq!(
        output.glwe_size(),
        glwe_secret_key.glwe_dimension().to_glwe_size(),
    );
    assert_eq!(
        output.polynomial_size(),
        glwe_secret_key.polynomial_size(),
    );
    assert_eq!(
        pt.plaintext_count().0,
        output.polynomial_size().0,
    );

    let polynomial_size = output.polynomial_size().0;
    let decomp_base_log = output.decomposition_base_log().0;

    let mut glev = output.as_mut_glwe_ciphertext_list();
    for (k, mut glwe) in glev.iter_mut().enumerate() {
        let level = k + 1;
        let log_scale = Scalar::BITS - level * decomp_base_log;

        let scaled_pt = PlaintextList::from_container((0..polynomial_size).map(|i| {
            *pt.get(i).0 << log_scale
        }).collect::<Vec<Scalar>>());
        encrypt_glwe_ciphertext(glwe_secret_key, &mut glwe, &scaled_pt, noise_parameters, generator);
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`GlevCiphertext`] entities.
#[derive(Clone, Copy)]
pub struct GlevCiphertextCreationMetadata<Scalar: UnsignedInteger>(
    pub GlweSize,
    pub PolynomialSize,
    pub DecompositionBaseLog,
    pub DecompositionLevelCount,
    pub CiphertextModulus<Scalar>,
);

impl<Scalar: UnsignedInteger, C: Container<Element = Scalar>> CreateFrom<C> for GlevCiphertext<C> {
    type Metadata = GlevCiphertextCreationMetadata<Scalar>;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let GlevCiphertextCreationMetadata(glwe_size, polynomial_size, decomposition_base_log, decomposition_level_count, ciphertext_modulus) = meta;
        Self::from_container(from, glwe_size, polynomial_size, decomposition_base_log, decomposition_level_count, ciphertext_modulus)
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Debug)]
pub struct GlevCiphertextCount(pub usize);

pub struct GlevCiphertextList<C: Container>
where
    C::Element: UnsignedInteger,
{
    data: C,
    glwe_size: GlweSize,
    polynomial_size: PolynomialSize,
    decomp_base_log: DecompositionBaseLog,
    decomp_level_count: DecompositionLevelCount,
    ciphertext_modulus: CiphertextModulus::<C::Element>,
}

/// A [`GlevCiphertext`] owning the memory for its own storage.
pub type GlevCiphertextListOwned<Scalar> = GlevCiphertextList<Vec<Scalar>>;
/// A [`GlevCiphertext`] immutably borrowing memory for its own storage.
pub type GlevCiphertextListView<'data, Scalar> = GlevCiphertextList<&'data [Scalar]>;
/// A [`GlevCiphertext`] mutably borrowing memory for its own storage.
pub type GlevCiphertextListMutView<'data, Scalar> = GlevCiphertextList<&'data mut [Scalar]>;


impl<T: UnsignedInteger, C: Container<Element= T>> AsRef<[T]> for GlevCiphertextList<C> {
    fn as_ref(&self) -> &[T] {
        self.data.as_ref()
    }
}

impl<T: UnsignedInteger, C: ContainerMut<Element= T>> AsMut<[T]> for GlevCiphertextList<C> {
    fn as_mut(&mut self) -> &mut [T] {
        self.data.as_mut()
    }
}

impl<Scalar: UnsignedInteger, C: Container<Element = Scalar>> GlevCiphertextList<C> {
    pub fn from_container(
        container: C,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
    ) -> GlevCiphertextList<C> {
        let glev_size = glwe_size.0 * polynomial_size.0 * decomp_level_count.0;
        assert_eq!(
            container.container_len() % glev_size, 0
        );

        Self {
            data: container,
            glwe_size: glwe_size,
            polynomial_size: polynomial_size,
            decomp_base_log: decomp_base_log,
            decomp_level_count: decomp_level_count,
            ciphertext_modulus: ciphertext_modulus,
        }
    }

    pub fn glwe_size(&self) -> GlweSize {
        self.glwe_size
    }

    pub fn polynomial_size(&self) -> PolynomialSize {
        self.polynomial_size
    }

    pub fn decomposition_base_log(&self) -> DecompositionBaseLog {
        self.decomp_base_log
    }

    pub fn decomposition_level_count(&self) -> DecompositionLevelCount {
        self.decomp_level_count
    }

    pub fn ciphertext_modulus(&self) -> CiphertextModulus::<Scalar> {
        self.ciphertext_modulus
    }

    pub fn into_container(self) -> C {
        self.data
    }
}

impl<Scalar: UnsignedInteger> GlevCiphertextListOwned<Scalar> {
    pub fn new(
        fill_with: Scalar,
        glwe_size: GlweSize,
        polynomial_size: PolynomialSize,
        decomp_base_log: DecompositionBaseLog,
        decomp_level_count: DecompositionLevelCount,
        glev_ciphertext_count: GlevCiphertextCount,
        ciphertext_modulus: CiphertextModulus::<Scalar>,
    ) -> GlevCiphertextListOwned<Scalar> {
        Self::from_container(
            vec![fill_with;
                glev_ciphertext_count.0
                    * glwe_size.0
                    * polynomial_size.0
                    * decomp_level_count.0
            ],
            glwe_size,
            polynomial_size,
            decomp_base_log,
            decomp_level_count,
            ciphertext_modulus,
        )
    }
}

/// Metadata used in the [`CreateFrom`] implementation to create [`GlevCiphertextList`] entities.
#[derive(Clone, Copy)]
pub struct GlevCiphertextListCreationMetadata<Scalar: UnsignedInteger>(
    pub GlweSize,
    pub PolynomialSize,
    pub DecompositionBaseLog,
    pub DecompositionLevelCount,
    pub CiphertextModulus<Scalar>,
);

impl<Scalar: UnsignedInteger, C: Container<Element = Scalar>> CreateFrom<C>
    for GlevCiphertextList<C>
{
    type Metadata = GlevCiphertextListCreationMetadata<Scalar>;

    #[inline]
    fn create_from(from: C, meta: Self::Metadata) -> Self {
        let GlevCiphertextListCreationMetadata(glwe_size, polynomial_size, decomp_base_log, decomp_level_count, ciphertext_modulus) =
            meta;
        Self::from_container(from, glwe_size, polynomial_size, decomp_base_log, decomp_level_count, ciphertext_modulus)
    }
}

impl<Scalar: UnsignedInteger, C: Container<Element = Scalar>> ContiguousEntityContainer
    for GlevCiphertextList<C>
{
    type Element = C::Element;

    type EntityViewMetadata = GlevCiphertextCreationMetadata<Self::Element>;

    type EntityView<'this> = GlevCiphertextView<'this, Self::Element>
    where
        Self: 'this;

    type SelfViewMetadata = GlevCiphertextListCreationMetadata<Self::Element>;

    type SelfView<'this> = GlevCiphertextListView<'this, Self::Element>
    where
        Self: 'this;

    fn get_entity_view_creation_metadata(&self) -> Self::EntityViewMetadata {
        GlevCiphertextCreationMetadata(self.glwe_size(), self.polynomial_size(), self.decomposition_base_log(), self.decomposition_level_count(), self.ciphertext_modulus())
    }

    fn get_entity_view_pod_size(&self) -> usize {
        glev_ciphertext_size(self.glwe_size(), self.polynomial_size(), self.decomposition_level_count())
    }

    fn get_self_view_creation_metadata(&self) -> Self::SelfViewMetadata {
        GlevCiphertextListCreationMetadata(
            self.glwe_size(),
            self.polynomial_size(),
            self.decomposition_base_log(),
            self.decomposition_level_count(),
            self.ciphertext_modulus(),
        )
    }
}

impl<Scalar: UnsignedInteger, C: ContainerMut<Element = Scalar>> ContiguousEntityContainerMut
    for GlevCiphertextList<C>
{
    type EntityMutView<'this> = GlevCiphertextMutView<'this, Self::Element>
    where
        Self: 'this;

    type SelfMutView<'this> = GlevCiphertextListMutView<'this, Self::Element>
    where
        Self: 'this;
}
