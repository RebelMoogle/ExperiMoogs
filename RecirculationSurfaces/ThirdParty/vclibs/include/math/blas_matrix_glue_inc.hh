//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================

/** \file

    This file is included by blas_matrix.hh.

    DO NOT INCLUDE IT DIRECTLY!

    \arg Map BLAS functions defined for vectors VC::math::blas::vec
    and matrix types, e.g., VC::math::blas::ge_mat to reference types
    VC::math::blas::vector_reference_t, VC::math::blas::matrix_reference_t.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/// \sa ::VC::math::blas::swap() \ingroup vc_blas
template <typename T,typename VX,typename VY>
void swap(const vector_reference_t<T,VX>& _x,
          const vector_reference_t<T,VY>& _y) {
  swap(_x.vector(),_x.data(),_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::scal() \ingroup vc_blas
template <typename T,typename V>
void scal(const T& _alpha,const vector_reference_t<T,V>& _x) {
  scal(_alpha,_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::copy() \ingroup vc_blas
  template <typename T,typename VX,typename VY>
void copy(const vector_const_reference_t<T,VX>& _x,
          const vector_reference_t<T,VY>& _y) {
  copy(_x.vector(),_x.data(),_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::copy_scal() \ingroup vc_blas
template <typename T,typename VX,typename VY>
void copy_scal(const T& _alpha,
          const vector_const_reference_t<T,VX>& _x,
          const vector_reference_t<T,VY>& _y) {
  copy_scal(_alpha,_x.vector(),_x.data(),_y.vector(),_y.data());
}


/// \sa ::VC::math::blas::axpy() \ingroup vc_blas
template <typename T,typename VX,typename VY>
void axpy(const T& _alpha,
          const vector_const_reference_t<T,VX>& _x,
          const vector_reference_t<T,VY>& _y) {
  axpy(_alpha,_x.vector(),_x.data(),_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::dot() \ingroup vc_blas
template <typename T,typename VX,typename VY>
T dot(const vector_const_reference_t<T,VX>& _x,
      const vector_const_reference_t<T,VY>& _y) {
  return dot(_x.vector(),_x.data(),_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::ddot() \ingroup vc_blas
template <typename T,typename VX,typename VY>
double ddot(const vector_const_reference_t<T,VX>& _x,
            const vector_const_reference_t<T,VY>& _y) {
  return ddot(_x.vector(),_x.data(),_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::nrm2() \ingroup vc_blas
template <typename T,typename V>
T nrm2(const vector_const_reference_t<T,V>& _x) {
  return nrm2(_x.vector(),_x.data());
}
/// \sa ::VC::math::blas::asum() \ingroup vc_blas
template <typename T,typename V>
T asum(const vector_const_reference_t<T,V>& _x) {
  return asum(_x.vector(),_x.data());
}

//-----------------------------------------------------------------------------

/// \sa ::VC::math::blas::gemv() \ingroup vc_blas
template <typename T,
          TransposeFlag TRANS,int M,int N,int LD,int NX,int INCX,int NY,int INCY>
void mv(const T& _alpha,
        const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
        const vector_const_reference_t<T,vec<NX,INCX> >& _x,
        const T& _beta,
        const vector_reference_t<T,vec<NY,INCY> >& _y) {
  mv(_alpha,_a.matrix(),_a.data(),_x.vector(),_x.data(),_beta,_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::gbmv() \ingroup vc_blas
template <typename T,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD,
          int NX,int INCX,int NY,int INCY>
void mv(const T& _alpha,
        const matrix_const_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> >& _a,
        const vector_const_reference_t<T,vec<NX,INCX> >& _x,
        const T& _beta,
        const vector_reference_t<T,vec<NY,INCY> >& _y) {
  mv(_alpha,_a.matrix(),_a.data(),_x.vector(),_x.data(),_beta,_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::symv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int LD,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const matrix_const_reference_t<T,sy_mat<UPLO,N,LD> >& _a,
        const vector_const_reference_t<T,vec<NV,INCX> >& _x,
        const T& _beta,
        const vector_reference_t<T,vec<NV,INCY> >& _y) {
  mv(_alpha,_a.matrix(),_a.data(),_x.vector(),_x.data(),_beta,_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::sbmv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a,
        const vector_const_reference_t<T,vec<NV,INCX> >& _x,
        const T& _beta,
        const vector_reference_t<T,vec<NV,INCY> >& _y) {
  mv(_alpha,_a.matrix(),_a.data(),_x.vector(),_x.data(),_beta,_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::spmv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const matrix_const_reference_t<T,sp_mat<UPLO,N> >& _a,
        const vector_const_reference_t<T,vec<NV,INCX> >& _x,
        const T& _beta,
        const vector_reference_t<T,vec<NV,INCY> >& _y) {
  mv(_alpha,_a.matrix(),_a.data(),_x.vector(),_x.data(),_beta,_y.vector(),_y.data());
}

/// \sa ::VC::math::blas::trmv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD,int NV,int INCX>
void mv(const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,N,LD> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  mv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::tbmv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int K,int LD,int NV,int INCX>
void mv(const matrix_const_reference_t<T,tb_mat<UPLO,TRANS,DIAG,N,K,LD> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  mv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::tpmv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int NV,int INCX>
void mv(const matrix_const_reference_t<T,tp_mat<UPLO,TRANS,DIAG,N> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  mv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::trsv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD,int NV,int INCX>
void sv(const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,N,LD> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  sv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::tbsv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int K,int LD,int NV,int INCX>
void sv(const matrix_const_reference_t<T,tb_mat<UPLO,TRANS,DIAG,N,K,LD> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  sv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::tpsv() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int NV,int INCX>
void sv(const matrix_const_reference_t<T,tp_mat<UPLO,TRANS,DIAG,N> >& _a,
        const vector_reference_t<T,vec<NV,INCX> >& _x) {
  sv(_a.matrix(),_a.data(),_x.vector(),_x.data());
}

/// \sa ::VC::math::blas::syr() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int LD,int NV,int INCX>
void r(const T& _alpha,
       const vector_const_reference_t<T,vec<NV,INCX> >& _x,
       const matrix_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  r(_alpha,_x.vector(),_x.data(),_a.matrix(),_a.data());
}

/// \sa ::VC::math::blas::spr() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int NV,int INCX>
void r(const T& _alpha,
       const vector_const_reference_t<T,vec<NV,INCX> >& _x,
       const matrix_reference_t<T,sp_mat<UPLO,N> >& _a) {
  r(_alpha,_x.vector(),_x.data(),_a.matrix(),_a.data());
}

/// \sa ::VC::math::blas::syr2() \ingroup vc_blas
  template <typename T,UpperLowerFlag UPLO,int N,int LD,int NX,int NY,int INCX,int INCY>
void r2(const T& _alpha,
        const vector_const_reference_t<T,vec<NX,INCX> >& _x,
        const vector_const_reference_t<T,vec<NY,INCY> >& _y,
        const matrix_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  r2(_alpha,_x.vector(),_x.data(),_y.vector(),_y.data(),_a.matrix(),_a.data());
}
/// \sa ::VC::math::blas::spr2() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int N,int NX,int NY,int INCX,int INCY>
void r2(const T& _alpha,
        const vector_const_reference_t<T,vec<NX,INCX> >& _x,
        const vector_const_reference_t<T,vec<NY,INCY> >& _y,
        const matrix_reference_t<T,sp_mat<UPLO,N> >& _a) {
  syr2(_alpha,_x.vector(),_x.data(),_y.vector(),_y.data(),_a.matrix(),_a.data());
}

//-----------------------------------------------------------------------------

/// set all elememts to _a \ingroup vc_blas
template <typename T,int N,int INC>
void ld_all(const T& _a,const vector_reference_t<T,vec<N,INC> >& _x) {
  ld_all(_a,_x.vector(),_x.data());
}

/// set all but _i-th element zero \ingroup vc_blas
template <typename T,int N,int INC>
void ld_unit(int _i,const T& _a,const vector_reference_t<T,vec<N,INC> >& _x) {
  ld_unit(_i,_a,_x.vector(),_x.data());
}

/// load random values in [0,1] \sa ::VC::math::blas::ld_rand() \ingroup vc_blas
  template <typename T,int N,int INC>
void ld_rand(const vector_reference_t<T,vec<N,INC> >& _x) {
  ld_rand(_x.vector(),_x.data());
}

/// add _alpha to all elements \ingroup vc_blas
template <typename T,typename V>
void adds(const T& _alpha,const vector_reference_t<T,V>& _x) {
  adds(_alpha,_x.vector(),_x.data());
}

//-----------------------------------------------------------------------------

/// \sa ::VC::math::blas::gemm() \ingroup vc_blas
template <typename T,
          TransposeFlag TRANSA,TransposeFlag TRANSB,
          int MA,int NA,int MB,int NB,int MC,int NC,
          int LDA,int LDB,int LDC>
void mm(const T& _alpha,
        const matrix_const_reference_t<T,ge_mat<TRANSA,MA,NA,LDA> >& _a,
        const matrix_const_reference_t<T,ge_mat<TRANSB,MB,NB,LDB> >& _b,
        const T& _beta,
        const matrix_reference_t<T,ge_mat<NoT,MC,NC,LDC> >& _c) {
  mm(_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data(),_beta,_c.matrix(),_c.data());
}

/// \sa ::VC::math::blas::symm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,
          int NA,int MB,int NB,int MC,int NC,
          int LDA,int LDB,int LDC>
void mm(SideFlag _side,
        const T& _alpha,
        const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
        const matrix_const_reference_t<T, ge_mat<NoT,MB,NB,LDB> >& _b,
        const T& _beta,
        const matrix_reference_t<T,ge_mat<NoT,MC,NC,LDC> >& _c) {
  mm(_side,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data(),
     _beta,_c.matrix(),_c.data());
}

/// (SideFlag = Left) \sa ::VC::math::blas::symm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,
          int NA,int MB,int NB,int MC,int NC,
          int LDA,int LDB,int LDC>
void mm(const T& _alpha,
        const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
        const matrix_const_reference_t<T, ge_mat<NoT,MB,NB,LDB> >& _b,
        const T& _beta,
        const matrix_reference_t<T,ge_mat<NoT,MC,NC,LDC> >& _c) {
  mm(Left,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data(),
     _beta,_c.matrix(),_c.data());
}

/// \sa ::VC::math::blas::syrk() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,
          int MA,int NA,int NC,int LDA,int LDC>
void rk(const T& _alpha,
        const matrix_const_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
        const T& _beta,
        const matrix_reference_t<T,sy_mat<UPLO,NC,LDC> >& _c) {
  rk(_alpha,_a.matrix(),_a.data(),_beta,_c.matrix(),_c.data());
}

/// \sa ::VC::math::blas::syr2k() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,
          int MA,int NA,int MB,int NB,int NC,int LDA,int LDB,int LDC>
void r2k(TransposeFlag _trans,
         const T& _alpha,
         const matrix_const_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
         const matrix_const_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b,
         const T& _beta,
         const matrix_reference_t<T,sy_mat<UPLO,NC,LDC> >& _c) {
  r2k(_trans,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data(),
      _beta,_c.matrix(),_c.data());
}

/// \sa ::VC::math::blas::trmm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int LDA,int LDB>
void mm(SideFlag _side,
        const T& _alpha,
        const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,NA,LDA> >& _a,
        const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  mm(_side,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data());
}

/// \sa ::VC::math::blas::trmm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int LDA,int LDB>
void mm(const T& _alpha,
        const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,NA,LDA> >& _a,
        const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  mm(Left,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data());
}

/// \sa ::VC::math::blas::trsm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int LDA,int LDB>
void sm(SideFlag _side,
        const T& _alpha,
        const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,NA,LDA> >& _a,
        const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  sm(_side,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data());
}

/// \sa ::VC::math::blas::trsm() \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int LDA,int LDB>
void sm(const T& _alpha,
        const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,NA,LDA> >& _a,
        const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  sm(Left,_alpha,_a.matrix(),_a.data(),_b.matrix(),_b.data());
}

//=============================================================================

/// load zeros \sa ::VC::math::blas::ld_zero() \ingroup vc_blas
template <typename T,typename A>
void ld_zero(const matrix_reference_t<T,A>& _a) {
  ld_zero(_a.matrix(),_a.data());
}

/// load constant \sa ::VC::math::blas::ld_all() \ingroup vc_blas
template <typename T,typename A>
void ld_all(const T& _alpha,const matrix_reference_t<T,A>& _a) {
  if (_alpha==T(0))
    ld_zero(_a.matrix(),_a.data());
  else
    ld_all(_alpha,_a.matrix(),_a.data());
}

/// load identity \sa ::VC::math::blas::ld_eye() \ingroup vc_blas
template <typename T,typename A>
void ld_eye(const matrix_reference_t<T,A>& _a,T _x=T(1)) {
  ld_eye(_a.matrix(),_a.data(),_x);
}

# ifndef DOXYGEN_SKIP
template <typename T>
struct _drand48_op {
  T operator()(T&) { return T(drand48()); }
};
#endif

/// load random values in [0,1] \sa ::VC::math::blas::ld_rand() \ingroup vc_blas
template <typename T,typename A>
void ld_rand(const matrix_reference_t<T,A>& _a) {
  map_f(_drand48_op<T>(),_a);
}

/// load diagonal matrix  \sa ::VC::math::blas::ld_diag() \ingroup vc_blas
template <typename T,typename V,typename A>
void ld_diag(const T& _alpha,
             const vector_const_reference_t<T,V>& _x,
             const matrix_reference_t<T,A>& _a) {
  ld_diag(_alpha,_x.vector(),_x.data(),_a.matrix(),_a.data());
}

/// scale matrix \sa ::VC::math::blas::mscal() \ingroup vc_blas
template <typename T,typename A>
void mscal(const T& _alpha,const matrix_reference_t<T,A>& _a) {
  mscal(_alpha,_a.matrix(),_a.data());
}
/// scale matrix \sa ::VC::math::blas::mscal() \ingroup vc_blas
template <typename T,typename A>
void adds(const T& _alpha,const matrix_reference_t<T,A>& _a) {
  adds(_alpha,_a.matrix(),_a.data());
}
/// elementwise multiplication ::VC::math::blas::emul() \ingroup vc_blas
template <typename T,typename A,typename B>
void emul(const T& _alpha,
           const vector_const_reference_t<T,B>& _b,
           const vector_reference_t<T,A>& _a) {
  emul(_alpha,_b.vector(),_b.data(),_a.vector(),_a.data());
}
/// elementwise multiplication ::VC::math::blas::emul() \ingroup vc_blas
template <typename T,typename A,typename B>
void emul(const T& _alpha,
           const matrix_const_reference_t<T,B>& _b,
           const matrix_reference_t<T,A>& _a) {
  emul(_alpha,_b.matrix(),_b.data(),_a.matrix(),_a.data());
}
/// elementwise multiplication ::VC::math::blas::ediv() \ingroup vc_blas
template <typename T,typename A,typename B>
void ediv(const T& _alpha,
           const matrix_const_reference_t<T,B>& _b,
           const matrix_reference_t<T,A>& _a) {
  ediv(_alpha,_b.matrix(),_b.data(),_a.matrix(),_a.data());
}
/// elementwise multiplication ::VC::math::blas::ediv() \ingroup vc_blas
template <typename T,typename A,typename B>
void ediv(const T& _alpha,
           const vector_const_reference_t<T,B>& _b,
           const vector_reference_t<T,A>& _a) {
  ediv(_alpha,_b.vector(),_b.data(),_a.vector(),_a.data());
}
/// add matrices _a+=_alpha*_b ::VC::math::blas::madd() \ingroup vc_blas
template <typename T,typename A,typename B>
void madd(const T& _alpha,
          const matrix_const_reference_t<T,B>& _b,
          const matrix_reference_t<T,A>& _a) {
  madd(_alpha,_b.matrix(),_b.data(),_a.matrix(),_a.data());
}
/// right-multiplication by diagonal matrix ::VC::math::blas::mscal_cols() \ingroup vc_blas
template <typename T,typename A,typename V>
void mscal_cols(const T& _alpha,
                const vector_const_reference_t<T,V>& _d,
                const matrix_reference_t<T,A>& _a) {
  mscal_cols(_alpha,_d.vector(),_d.data(),_a.matrix(),_a.data());
}
///  \ingroup vc_blas
template <typename T,typename A,typename B>
void copy(const matrix_const_reference_t<T,B>& _b,
          const matrix_reference_t<T,A>& _a) {
  copy(_b.matrix(),_b.data(),_a.matrix(),_a.data());
}

/// add transposed: A=(A+A')*_alpha \ingroup vc_blas
template <typename T,typename A>
void add_trans(const T& _alpha,
               const matrix_reference_t<T,A>& _a) {
  add_trans(_alpha,_a.matrix(),_a.data());
}

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
