use aligned_vec::CACHELINE_ALIGN;
use dyn_stack::{ReborrowMut, SizeOverflow, StackReq};
use tfhe::core_crypto::{
    prelude::*,
    fft_impl::fft64::{
        c64,
        math::{
            fft::FftView,
            polynomial::{FourierPolynomialView, FourierPolynomialMutView},
        },
    },
};
use crate::izip;

pub fn fourier_poly_mult_and_backward<Scalar, LhsCont, RhsCont, OutputCont>(
    output: &mut Polynomial<OutputCont>,
    lhs_int: &FourierPolynomial<LhsCont>,
    rhs_torus: &FourierPolynomial<RhsCont>,
) where
    Scalar: UnsignedTorus,
    LhsCont: Container<Element=c64>,
    RhsCont: Container<Element=c64>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(lhs_int.polynomial_size(), rhs_torus.polynomial_size());
    assert_eq!(lhs_int.polynomial_size(), output.polynomial_size());

    let lhs_int = lhs_int.as_view();
    let rhs_torus = rhs_torus.as_view();

    let polynomial_size = lhs_int.polynomial_size();
    let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;

    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();

    let mut buffers = ComputationBuffers::new();
    buffers.resize(
        fourier_poly_mult_scratch(fft)
        .unwrap()
        .unaligned_bytes_required(),
    );

    let stack = buffers.stack();
    let (mut output_buffer, substack0) = stack.make_aligned_raw::<c64>(
        fourier_poly_size,
        CACHELINE_ALIGN,
    );
    let output_buffer = &mut *output_buffer;

    update_with_fmadd(
        output_buffer,
        lhs_int.data,
        rhs_torus.data,
        true,
        fourier_poly_size,
    );

    let fourier_output = FourierPolynomialView { data: output_buffer };
    fft.backward_as_torus(output.as_mut_view(), fourier_output, substack0);
}

pub fn fourier_poly_mult<LhsCont, RhsCont, OutputCont>(
    output: &mut FourierPolynomial<OutputCont>,
    lhs_int: &FourierPolynomial<LhsCont>,
    rhs_torus: &FourierPolynomial<RhsCont>,
) where
    LhsCont: Container<Element=c64>,
    RhsCont: Container<Element=c64>,
    OutputCont: ContainerMut<Element=c64>,
{
    assert_eq!(lhs_int.polynomial_size(), rhs_torus.polynomial_size());
    assert_eq!(lhs_int.polynomial_size(), output.polynomial_size());

    let lhs_int = lhs_int.as_view();
    let rhs_torus = rhs_torus.as_view();
    let output = output.as_mut_view();

    let polynomial_size = lhs_int.polynomial_size();
    let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;

    update_with_fmadd(
        &mut *output.data,
        lhs_int.data,
        rhs_torus.data,
        true,
        fourier_poly_size,
    );
}

pub fn fourier_poly_mult_and_add<LhsCont, RhsCont, OutputCont>(
    output: &mut FourierPolynomial<OutputCont>,
    lhs_int: &FourierPolynomial<LhsCont>,
    rhs_torus: &FourierPolynomial<RhsCont>,
) where
    LhsCont: Container<Element=c64>,
    RhsCont: Container<Element=c64>,
    OutputCont: ContainerMut<Element=c64>,
{
    assert_eq!(lhs_int.polynomial_size(), rhs_torus.polynomial_size());
    assert_eq!(lhs_int.polynomial_size(), output.polynomial_size());

    let lhs_int = lhs_int.as_view();
    let rhs_torus = rhs_torus.as_view();
    let output = output.as_mut_view();

    let polynomial_size = lhs_int.polynomial_size();
    let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;

    update_with_fmadd(
        &mut *output.data,
        lhs_int.data,
        rhs_torus.data,
        false,
        fourier_poly_size,
    );
}

pub fn polynomial_mul_by_fft<Scalar, LhsCont, RhsCont, OutputCont>(
    output: &mut Polynomial<OutputCont>,
    lhs_int: &Polynomial<LhsCont>,
    rhs_torus: &Polynomial<RhsCont>,
) where
    Scalar: UnsignedTorus,
    LhsCont: Container<Element=Scalar>,
    RhsCont: Container<Element=Scalar>,
    OutputCont: ContainerMut<Element=Scalar>,
{
    assert_eq!(lhs_int.polynomial_size(), rhs_torus.polynomial_size());
    assert_eq!(lhs_int.polynomial_size(), output.polynomial_size());

    output.as_mut().fill(Scalar::ZERO);
    let polynomial_size = lhs_int.polynomial_size();

    let fft = Fft::new(polynomial_size);
    let fft = fft.as_view();

    let mut buffers = ComputationBuffers::new();
    buffers.resize(
        polynomial_mul_by_fft_scratch(polynomial_size, fft)
        .unwrap()
        .unaligned_bytes_required(),
    );
    let mut stack = buffers.stack();

    let fourier_poly_size = polynomial_size.to_fourier_polynomial_size().0;
    let align = CACHELINE_ALIGN;

    let (mut fourier_lhs, mut substack0) = stack.rb_mut().make_aligned_raw::<c64>(fourier_poly_size, align);
    let (mut fourier_rhs, mut substack1) = substack0.rb_mut().make_aligned_raw::<c64>(fourier_poly_size, align);
    let (mut fourier_out, mut substack2) = substack1.rb_mut().make_aligned_raw::<c64>(fourier_poly_size, align);
    let fourier_out = &mut *fourier_out;

    let fourier_lhs = fft
        .forward_as_integer(
            FourierPolynomialMutView { data: &mut fourier_lhs },
            lhs_int.as_view(),
            substack2.rb_mut(),
        )
        .data;
    let fourier_rhs = fft
        .forward_as_torus(
            FourierPolynomialMutView { data: &mut fourier_rhs },
            rhs_torus.as_view(),
            substack2.rb_mut(),
        )
        .data;

    update_with_fmadd(
        fourier_out,
        fourier_lhs,
        fourier_rhs,
        true,
        fourier_poly_size,
    );

    let fourier_out = FourierPolynomialView { data: fourier_out };
    fft.backward_as_torus(output.as_mut_view(), fourier_out, substack2.rb_mut());
}

pub fn fourier_poly_mult_scratch(
    fft: FftView<'_>,
) -> Result<StackReq, SizeOverflow> {
    let align = CACHELINE_ALIGN;
    let fourier_polynomial_size = fft.polynomial_size().to_fourier_polynomial_size().0;
    let fourier_scratch = StackReq::try_new_aligned::<c64>(fourier_polynomial_size, align)?;

    let substack0 = fft.backward_scratch()?;
    substack0.try_and(fourier_scratch)
}

pub fn polynomial_mul_by_fft_scratch(
    polynomial_size: PolynomialSize,
    fft: FftView<'_>,
) -> Result<StackReq, SizeOverflow> {
    let align = CACHELINE_ALIGN;
    let fourier_polynomial_size = polynomial_size.to_fourier_polynomial_size().0;
    let fourier_scratch = StackReq::try_new_aligned::<c64>(fourier_polynomial_size, align)?;

    let substack2 = StackReq::try_any_of([
        fft.forward_scratch()?,
        fft.backward_scratch()?,
    ])?;
    let substack1 = substack2.try_and(fourier_scratch)?;
    let substack0 = substack1.try_and(fourier_scratch)?;
    substack0.try_and(fourier_scratch)
}

// From tfhe::core_crypto::fft64::crypto::ggsw
#[cfg_attr(__profiling, inline(never))]
pub(crate) fn update_with_fmadd(
    output_fft_buffer: &mut [c64],
    lhs_polynomial_list: &[c64],
    fourier: &[c64],
    is_output_uninit: bool,
    fourier_poly_size: usize,
) {
    struct Impl<'a> {
        output_fft_buffer: &'a mut [c64],
        lhs_polynomial_list: &'a [c64],
        fourier: &'a [c64],
        is_output_uninit: bool,
        fourier_poly_size: usize,
    }

    impl pulp::WithSimd for Impl<'_> {
        type Output = ();

        #[inline(always)]
        fn with_simd<S: pulp::Simd>(self, simd: S) -> Self::Output {
            // Introducing a function boundary here means that the slices
            // get `noalias` markers, possibly allowing better optimizations from LLVM.
            //
            // see:
            // https://github.com/rust-lang/rust/blob/56e1aaadb31542b32953292001be2312810e88fd/library/core/src/slice/mod.rs#L960-L966
            #[inline(always)]
            fn implementation<S: pulp::Simd>(
                simd: S,
                output_fft_buffer: &mut [c64],
                lhs_polynomial_list: &[c64],
                fourier: &[c64],
                is_output_uninit: bool,
                fourier_poly_size: usize,
            ) {
                let rhs = S::c64s_as_simd(fourier).0;

                if is_output_uninit {
                    for (output_fourier, ggsw_poly) in izip!(
                        output_fft_buffer.into_chunks(fourier_poly_size),
                        lhs_polynomial_list.into_chunks(fourier_poly_size)
                    ) {
                        let out = S::c64s_as_mut_simd(output_fourier).0;
                        let lhs = S::c64s_as_simd(ggsw_poly).0;

                        for (out, &lhs, &rhs) in izip!(out, lhs, rhs) {
                            *out = simd.c64s_mul(lhs, rhs);
                        }
                    }
                } else {
                    for (output_fourier, ggsw_poly) in izip!(
                        output_fft_buffer.into_chunks(fourier_poly_size),
                        lhs_polynomial_list.into_chunks(fourier_poly_size)
                    ) {
                        let out = S::c64s_as_mut_simd(output_fourier).0;
                        let lhs = S::c64s_as_simd(ggsw_poly).0;

                        for (out, &lhs, &rhs) in izip!(out, lhs, rhs) {
                            *out = simd.c64s_mul_add_e(lhs, rhs, *out);
                        }
                    }
                }
            }

            implementation(
                simd,
                self.output_fft_buffer,
                self.lhs_polynomial_list,
                self.fourier,
                self.is_output_uninit,
                self.fourier_poly_size,
            );
        }
    }

    pulp::Arch::new().dispatch(Impl {
        output_fft_buffer,
        lhs_polynomial_list,
        fourier,
        is_output_uninit,
        fourier_poly_size,
    });
}
