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

    \arg Map BLAS functions to vector and matrix types in VC::math::blas
    (VC::math::blas::vec,VC::math::blas::ge_mat,...).

    \arg Specializations may apply!

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/// \sa ::VC::math::blas::swap() \ingroup vc_blas_mat
  template <typename T,typename VX,typename VY>
void swap(const VX& _vx,T* _x,const VY& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_vy.n());
  assert(_x!=_y);
  swap(_vx.n(),_x,_vx.inc(),_y,_vy.inc());
}

/// \sa ::VC::math::blas::scal() \ingroup vc_blas_mat
template <typename T,typename V>
void scal(const T& _alpha,const V& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  scal(_vx.n(),_alpha,_x,_vx.inc());
}

/// \sa ::VC::math::blas::copy() \ingroup vc_blas_mat
template <typename T,typename VX,typename VY>
void copy(const VX& _vx,const T* _x,const VY& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_vy.n());
  assert(_x!=_y);
  copy(_vx.n(),_x,_vx.inc(),_y,_vy.inc());
}

/// \sa ::VC::math::blas::copy()+::VC::math::blas::scal(): y=a*x \ingroup vc_blas_mat
template <typename T,typename VX,typename VY>
void copy_scal(const T& _alpha,const VX& _vx,const T* _x,const VY& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_x!=_y);
  copy(_vx.n(),_x,_vx.inc(),_y,_vy.inc());
  if (_alpha!=T(1))
    scal(_vy.n(),_alpha,_y,_vy.inc());
}


/// \sa ::VC::math::blas::axpy() \ingroup vc_blas_mat
template <typename T,typename VX,typename VY>
void axpy(const T& _alpha,const VX& _vx,const T* _x,const VY& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_vy.n());
  assert(_x!=_y);
  axpy(_vx.n(),_alpha,_x,_vx.inc(),_y,_vy.inc());
}

/// \sa ::VC::math::blas::dot() \ingroup vc_blas_mat
template <typename T,typename VX,typename VY>
T dot(const VX& _vx,const T* _x,const VY& _vy,const T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_vy.n());
  return dot(_vx.n(),_x,_vx.inc(),_y,_vy.inc());
}

/// \sa ::VC::math::blas::dot() \ingroup vc_blas_mat
template <typename T,typename VX,typename VY>
double ddot(const VX& _vx,const T* _x,const VY& _vy,const T* _y) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_vy.n());
  return ddot(_vx.n(),_x,_vx.inc(),_y,_vy.inc());
}

/// \sa ::VC::math::blas::nrm2() \ingroup vc_blas_mat
template <typename T,typename V>
T nrm2(const V& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  return nrm2(_vx.n(),_x,_vx.inc());
}
/// \sa ::VC::math::blas::asum() \ingroup vc_blas_mat
template <typename T,typename V>
T asum(const V& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  return asum(_vx.n(),_x,_vx.inc());
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

/*
  Q: Why do we need template parameter NV? - Should be equal to N?
  A: Yes, it should be equal to N but we might mix fixed and VarInt
     dimensions: then this would not compile! ("no matching function for call").
     Consequently, dimensions are only check at runtime!
 */

# define _CHECK_DIM_XY()                              \
  assert(_ma.n()== (TRANS==NoT ? _vx.n() : _vy.n())); \
  assert(_ma.m()== (TRANS==NoT ? _vy.n() : _vx.n())); \
  assert(_x!=_y)

# define _CHECK_DIM_SXY()                             \
  assert(_ma.n()==_vx.n());                           \
  assert(_ma.n()==_vy.n());                           \
  assert(_x!=_y)

# define _CHECK_DIM_X() assert(_ma.n()== _vx.n())
# define _CHECK_DIM_Y() assert(_ma.n()== _vy.n())

/// \sa ::VC::math::blas::gemv() \ingroup vc_blas_mat
template <typename T,
          TransposeFlag TRANS,int M,int N,int LD,int NX,int INCX,int NY,int INCY>
void mv(const T& _alpha,
        const ge_mat<TRANS,M,N,LD>& _ma,const T* _a,
        const vec<NX,INCX>& _vx,const T* _x,
        const T& _beta,
        const vec<NY,INCY>& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_XY();
  gemv(_ma.trans(),_ma.m(),_ma.n(),_alpha,_a,_ma.ld(),
       _x,_vx.inc(),_beta,_y,_vy.inc());
}

  // specialize for M,N=const, M,N,INCX,INCY=const -- benchmark

/// \sa ::VC::math::blas::gbmv() \ingroup vc_blas_mat
template <typename T,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD,
          int NX,int INCX,int NY,int INCY>
void mv(const T& _alpha,
        const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,const T* _a,
        const vec<NX,INCX>& _vx,const T* _x,
        const T& _beta,
        const vec<NY,INCY>& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_XY();
  gbmv(_ma.trans(),_ma.m(),_ma.n(),_ma.kl(),_ma.ku(),_alpha,_a,_ma.ld(),
       _x,_vx.inc(),_beta,_y,_vy.inc());
}

/// \sa ::VC::math::blas::symv() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const sy_mat<UPLO,N,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,const T* _x,
        const T& _beta,
        const vec<NV,INCY>& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_SXY();
  symv(_ma.uplo(),_ma.n(),_alpha,_a,_ma.ld(),
       _x,_vx.inc(),_beta,_y,_vy.inc());
}

  // specialize for M,N=const, M,N,INCX,INCY=const -- benchmark

/// \sa ::VC::math::blas::sbmv() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const sb_mat<UPLO,N,K,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,const T* _x,
        const T& _beta,
        const vec<NV,INCY>& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_SXY();
  sbmv(_ma.uplo(),_ma.n(),_ma.k(),_alpha,_a,_ma.ld(),
       _x,_vx.inc(),_beta,_y,_vy.inc());
}

  // specialize for K=0,M,N=const, K=0,M,N,INCX,INCY=const -- benchmark

/// \sa ::VC::math::blas::spmv() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int NV,int INCX,int INCY>
void mv(const T& _alpha,
        const sp_mat<UPLO,N>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,const T* _x,
        const T& _beta,
        const vec<NV,INCY>& _vy,T* _y) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_SXY();
  spmv(_ma.uplo(),_ma.n(),_alpha,_a,
       _x,_vx.inc(),_beta,_y,_vy.inc());
}

/// \sa ::VC::math::blas::trmv() \ingroup vc_blas_mat_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD,int NV,int INCX>
void mv(const tr_mat<UPLO,TRANS,DIAG,N,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  trmv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_a,_ma.ld(),_x,_vx.inc());
}

  // have this prototype also for all others; specialize

/// \sa ::VC::math::blas::tbmv() \ingroup vc_blas_mat_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int K,int LD,int NV,int INCX>
void mv(const tb_mat<UPLO,TRANS,DIAG,N,K,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  tbmv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_ma.k(),_a,_ma.ld(),_x,_vx.inc());
}

/// \sa ::VC::math::blas::tpmv() \ingroup vc_blas_mat_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int NV,int INCX>
void mv(const tp_mat<UPLO,TRANS,DIAG,N>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  tpmv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_a,_x,_vx.inc());
}

  // extend these to others...? not here!


/// \sa ::VC::math::blas::trsv() \ingroup vc_blas_mat_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int NV,int LD,int INCX>
void sv(const tr_mat<UPLO,TRANS,DIAG,N,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  trsv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_a,_ma.ld(),_x,_vx.inc());
}

  // have this prototype also for all others; specialize

/// \sa ::VC::math::blas::tbsv() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int K,int LD,int NV,int INCX>
void sv(const tb_mat<UPLO,TRANS,DIAG,N,K,LD>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  tbsv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_ma.k(),_a,_ma.ld(),_x,_vx.inc());
}

/// \sa ::VC::math::blas::tpsv() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int NV,int INCX>
void sv(const tp_mat<UPLO,TRANS,DIAG,N>& _ma,const T* _a,
        const vec<NV,INCX>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  tpsv(_ma.uplo(),_ma.trans(),_ma.diag(),_ma.n(),_a,_x,_vx.inc());
}


/// \sa ::VC::math::blas::syr() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD,int NV,int INCX>
void r(const T& _alpha,
       const vec<NV,INCX>& _vx,const T* _x,
       const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  syr(_ma.uplo(),_ma.n(),_alpha,_x,_vx.inc(),_a,_ma.ld());
}

/// \sa ::VC::math::blas::spr() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int NV,int INCX>
void r(const T& _alpha,
       const vec<NV,INCX>& _vx,const T* _x,
       const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X();
  spr(_ma.uplo(),_ma.n(),_alpha,_x,_vx.inc(),_a);
}

/// \sa ::VC::math::blas::syr2() \ingroup vc_blas_mat
  template <typename T,UpperLowerFlag UPLO,int N,int LD,int NX,int NY,int INCX,int INCY>
void r2(const T& _alpha,
        const vec<NX,INCX>& _vx,const T* _x,
        const vec<NY,INCY>& _vy,const T* _y,
        const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X(); _CHECK_DIM_Y();
  syr2(_ma.uplo(),_ma.n(),_alpha,_x,_vx.inc(),_y,_vy.inc(),_a,_ma.ld());
}

/// \sa ::VC::math::blas::spr2() \ingroup vc_blas_mat
  template <typename T,UpperLowerFlag UPLO,int N,int NX,int NY,int INCX,int INCY>
void r2(const T& _alpha,
        const vec<NX,INCX>& _vx,const T* _x,
        const vec<NY,INCY>& _vy,const T* _y,
        const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  _CHECK_DIM_X(); _CHECK_DIM_Y();
  spr2(_ma.uplo(),_ma.n(),_alpha,_x,_vx.inc(),_y,_vy.inc(),_a);
}


# undef _CHECK_DIM_Y
# undef _CHECK_DIM_X
# undef _CHECK_DIM_XY

//-----------------------------------------------------------------------------

/// \sa ::VC::math::blas::gemm() \ingroup vc_blas_mat
template <typename T,
          TransposeFlag TRANSA,TransposeFlag TRANSB,
          int MA,int NA,int MB,int NB,int MC,int NC,
          int LDA,int LDB,int LDC>
void mm(const T& _alpha,
        const ge_mat<TRANSA,MA,NA,LDA>& _ma,const T* _a,
        const ge_mat<TRANSB,MB,NB,LDB>& _mb,const T* _b,
        const T& _beta,
        const ge_mat<NoT,MC,NC,LDC>& _mc,T* _c) {
  _VC_DBG_BLAS_SCOPE();

  int m,n,k;

  // see BLAS specs: we have to evaluate TRANSA, TRANSB
  if (_ma.trans()==NoT) {
    m=_ma.m();
    k=_ma.n();
  }
  else {
    m=_ma.n();
    k=_ma.m();
  }

  n=_mb.trans()==NoT ? _mb.n() : _mb.m();

  assert(m==_mc.m() && "dimensions don't match");
  assert(n==_mc.n() && "dimensions don't match");
  assert(k==(_mb.trans()==NoT ? _mb.m() : _mb.n()) && "dimensions don't match");

  assert(_c!=_a && _c!=_b);

  gemm(_ma.trans(),_mb.trans(), m,n,k,
       _alpha,_a,_ma.ld(), _b,_mb.ld(), _beta, _c,_mc.ld());
}

/// \sa ::VC::math::blas::symm() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,
          int NA,int MB,int NB,int MC,int NC,
          int LDA,int LDB,int LDC>
void mm(SideFlag _side,
        const T& _alpha,
        const sy_mat<UPLO,NA,LDA>& _ma,const T* _a,
        const ge_mat<NoT,MB,NB,LDB>& _mb,const T* _b,
        const T& _beta,
        const ge_mat<NoT,MC,NC,LDC>& _mc,T* _c) {
  _VC_DBG_BLAS_SCOPE();

# ifndef NDEBUG
  if (_side==Left) {
    assert(_ma.n()==_mb.m() && "dimensions don't match");
    assert(_mc.m()==_ma.n() && "dimensions don't match");
    assert(_mc.n()==_mb.n() && "dimensions don't match");
  }
  else {
    assert(_ma.n()==_mb.n() && "dimensions don't match");
    assert(_mc.m()==_mb.m() && "dimensions don't match");
    assert(_mc.n()==_ma.n() && "dimensions don't match");
  }

  assert(_c!=_a && _c!=_b);
#endif

  symm(_side,_ma.uplo(),_mc.m(),_mc.n(),
       _alpha,_a,_ma.ld(), _b,_mb.ld(), _beta,_c,_mc.ld());
}

/// \sa ::VC::math::blas::syrk() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,
          int MA,int NA,int NC,int LDA,int LDC>
void rk(const T& _alpha,
        const ge_mat<TRANS,MA,NA,LDA>& _ma,const T* _a,
        const T& _beta,
        const sy_mat<UPLO,NC,LDC>& _mc,T* _c) {
  _VC_DBG_BLAS_SCOPE();

  assert(_mc.n()==(_ma.trans()==NoT ? _ma.m() : _ma.n()) && "dimensions don't match");
  assert(_c!=_a);

  syrk(_mc.uplo(),_ma.trans(),_mc.n(),_ma.trans()==NoT ? _ma.n() : _ma.m(),
       _alpha,_a,_ma.ld(),_beta,_c,_mc.ld());
}

/// \sa ::VC::math::blas::syr2k() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,
          int MA,int NA,int MB,int NB,int NC,int LDA,int LDB,int LDC>
void r2k(TransposeFlag _trans,
         const T& _alpha,
         const ge_mat<NoT,MA,NA,LDA>& _ma,const T* _a,
         const ge_mat<NoT,MB,NB,LDB>& _mb,const T* _b,
         const T& _beta,
         const sy_mat<UPLO,NC,LDC>& _mc,T* _c) {
  _VC_DBG_BLAS_SCOPE();

# ifndef NDEBUG
  if (_trans==NoT) {
    assert(_mc.n()==_ma.m() && "dimensions don't match");
    assert(_mc.n()==_mb.m() && "dimensions don't match");
    assert(_ma.n()==_mb.n() && "dimensions don't match");
  }
  else {
    assert(_mc.n()==_ma.n() && "dimensions don't match");
    assert(_mc.n()==_mb.n() && "dimensions don't match");
    assert(_ma.m()==_mb.m() && "dimensions don't match");
  }
  assert(_c!=_a && _c!=_b);
# endif

  syr2k(_mc.uplo(),_trans,_mc.n(),_trans==NoT ? _ma.n() : _ma.m(),
        _alpha,_a,_ma.ld(),_b,_mb.ld(),_beta,_c,_mc.ld());
}

/// \sa ::VC::math::blas::trmm() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int LDA,int LDB>
void mm(SideFlag _side,
        const T& _alpha,
        const tr_mat<UPLO,TRANS,DIAG,NA,LDA>& _ma,const T* _a,
        const ge_mat<NoT,MB,NB,LDB>& _mb,T* _b) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==(_side==Left ? _mb.m() : _mb.n()) && "dimensions don't match");
  assert(_b!=_a);

  trmm(_side,_alpha,_a,_ma.ld(),_b,_mb.ld());
}

/// \sa ::VC::math::blas::trsm() \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int MB,int NB,int MX,int NX,int LDA,int LDB>
void sm(SideFlag _side,
        const T& _alpha,
        const tr_mat<UPLO,TRANS,DIAG,NA,LDA>& _ma,const T* _a,
        const ge_mat<NoT,MB,NB,LDB>& _mb,T* _b) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==(_side==Left ? _mb.m() : _mb.n()) && "dimensions don't match");
  assert(_b!=_a);

  trsm(_side,_alpha,_a,_ma.ld(),_b,_mb.ld());
}

//=============================================================================

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void ld_zero(const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<_ma.m();++i)
      _a[k+i]=T(0);
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
void ld_zero(const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int m=_ma.kl()+_ma.ku()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]=T(0);
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD>
void ld_zero(const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]=T(0);
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
  }
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
void ld_zero(const sb_mat<UPLO,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]=T(0);
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N>
void ld_zero(const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]=T(0);
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int LD>
void ld_zero(const tr_mat<UPLO,TRANS,NoU,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]=T(0);
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
  }
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int K,int LD>
void ld_zero(const tb_mat<UPLO,TRANS,NoU,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]=T(0);
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N>
void ld_zero(const tp_mat<UPLO,TRANS,NoU,N>& _ma,T* _a) {
  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]=T(0);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,int N,int INC>
void ld_all(const T& _alpha,const vec<N,INC>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  int n=_vx.n(), inc=_vx.inc();
  for (int i=0;i<n;++i)
    _x[i*inc]=_alpha;
}

/// set all elements zero and _i-th element _alpha  \ingroup vc_blas_mat
template <typename T,int N,int INC>
void ld_unit(int _i,const T& _alpha,const vec<N,INC>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  assert(0<=_i && _i<_vx.n());

  int n=_vx.n(), inc=_vx.inc();
  for (int i=0;i<n;++i)
    _x[i*inc]=(i==_i) ? _alpha : T(0);

}

/// set random values in [0,1] (`drand48()`) \ingroup vc_blas_mat
template <typename T,int N,int INC>
void ld_rand(const vec<N,INC>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  int n=_vx.n(), inc=_vx.inc();
  for (int i=0;i<n;++i)
    _x[i*inc]=T(drand48());
}


/// set all elements to zero \ingroup vc_blas_mat
template <typename T,int N,int INC>
void adds(const T& _alpha,const vec<N,INC>& _vx,T* _x) {
  _VC_DBG_BLAS_SCOPE();
  int n=_vx.n(), inc=_vx.inc();
  for (int i=0;i<n;++i)
    _x[i*inc]+=_alpha;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// set all elements to _alpha \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void ld_all(const T& _alpha,const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<_ma.m();++i)

      _a[k+i]=_alpha;
}

/// set all elements to _alpha \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD>
void ld_all(const T& _alpha,const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]=_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]=_alpha;
  }
}

/// set all elements to _alpha \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N>
void ld_all(const T& _alpha,const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]=_alpha;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// load identity matrix \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void ld_eye(const ge_mat<TRANS,M,N,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<_ma.m();++i)
      _a[k+i]=(i==j) ? _x : T(0);
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
void ld_eye(const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld(), ku=_ma.ku();
  int m=_ma.kl()+ku+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]=(i==ku) ? _x : T(0);
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD>
void ld_eye(const sy_mat<UPLO,N,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0) {
      _a[k+m0]=_x;
      for (int i=m0+1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    int m=0,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m) {
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
      _a[k+m]=_x;
    }
  }
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
void ld_eye(const sb_mat<UPLO,N,K,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;

  if (UPLO==Lower) {
    for (int j=0;j<_ma.n();++j,k+=ld) {
      _a[k]=_x;
      for (int i=1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    for (int j=0;j<_ma.n();++j,k+=ld) {
      for (int i=0;i<m-1;++i)
        _a[k+i]=T(0);
      _a[k+m-1]=_x;
    }
  }

}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N>
void ld_eye(const sp_mat<UPLO,N>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  ld_zeros(_ma,_a);
  int n=_ma.n(), i=0;
  int m=n;

  for (int j=0;j<n;++j,--m) {
    _a[i]=_x;
    i+=m;
  }
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,
          UpperLowerFlag UPLO, TransposeFlag TRANS,
          int N,int LD>
void ld_eye(const tr_mat<UPLO,TRANS,NoU,N,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0) {
      _a[k+m0]=_x;
      for (int i=m0+1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    int m=0,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m) {
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
      _a[k+m]=_x;
    }
  }
}
/// load identity matrix \ingroup vc_blas_mat
template <typename T,
          UpperLowerFlag UPLO, TransposeFlag TRANS,
          int N,int LD>
void ld_eye(const tr_mat<UPLO,TRANS,UnitDiag,N,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();
  assert((_x==1.0) && "matrix has unit diagonal by definition");

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0) {
      // do not modify diagonal
      for (int i=m0+1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    int m=0,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m) {
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
      // do not modify diagonal
    }
  }
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO, TransposeFlag TRANS,int N,int K,int LD>
void ld_eye(const tb_mat<UPLO,TRANS,NoU,N,K,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;

  if (UPLO==Lower) {
    for (int j=0;j<_ma.n();++j,k+=ld) {
      _a[k]=_x;
      for (int i=1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    for (int j=0;j<_ma.n();++j,k+=ld) {
      for (int i=0;i<m-1;++i)
        _a[k+i]=T(0);
      _a[k+m-1]=_x;
    }
  }
}
/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO, TransposeFlag TRANS,int N,int K,int LD>
void ld_eye(const tb_mat<UPLO,TRANS,UnitDiag,N,K,LD>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();
  assert((_x==1.0) && "matrix has unit diagonal by definition");

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;

  if (UPLO==Lower) {
    for (int j=0;j<_ma.n();++j,k+=ld) {
      // do not modify diagonal
      for (int i=1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    for (int j=0;j<_ma.n();++j,k+=ld) {
      for (int i=0;i<m-1;++i)
        _a[k+i]=T(0);
      // do not modify diagonal
    }
  }
}

/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO, TransposeFlag TRANS,int N>
void ld_eye(const tp_mat<UPLO,TRANS,NoU,N>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();

  ld_zeros(_ma,_a);
  int n=_ma.n(), i=0;
  int m=n;

  for (int j=0;j<n;++j,--m) {
    _a[i]=_x;
    i+=m;
  }
}
/// load identity matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO, TransposeFlag TRANS,int N>
void ld_eye(const tp_mat<UPLO,TRANS,UnitDiag,N>& _ma,T* _a,const T& _x=T(1)) {
  _VC_DBG_BLAS_SCOPE();
  assert((_x==1.0) && "matrix has unit diagonal by definition");
  for (int i=0;i<_ma.size();++i) _a[i]=T(0);
}

//-----------------------------------------------------------------------------

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,TransposeFlag TRANS,int M,int N,int LD>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_ma.m());

  int ld=_ma.ld(), inc=_vx.inc();
  int k=0,ii=0;
  for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc)
    for (int i=0;i<_ma.m();++i)
      _a[k+i]=(i==j) ? _alpha*_x[ii] : T(0);
}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
void ld_diag(const T& _alpha,const vec<NX,INCX>& _vx,const T* _x,
             const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_ma.m());

  int ld=_ma.ld(), ku=_ma.ku(), inc=_vx.inc();
  int m=_ma.kl()+ku+1;
  int k=0,ii=0;
  for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc)
    for (int i=0;i<m;++i)
      _a[k+i]=(i==ku) ? _alpha*_x[ii] : T(0);
}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,UpperLowerFlag UPLO,int N,int LD>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_ma.m());

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld(), inc=_vx.inc();
    int k=0,ii=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0,ii+=inc) {
      _a[k+m0]=_alpha*_x[ii];
      for (int i=m0+1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    int m=0,ld=_ma.ld(), inc=_vx.inc();
    int k=0,ii=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m,ii+=inc) {
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
      _a[k+m]=_alpha*_x[ii];
    }
  }
}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,UpperLowerFlag UPLO,int N,int K,int LD>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const sb_mat<UPLO,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_ma.m());

  int ld=_ma.ld(), inc=_vx.inc();
  int m=_ma.k()+1;
  int k=0,ii=0;

  if (UPLO==Lower) {
    for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc) {
      _a[k]=_alpha*_x[ii];
      for (int i=1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc) {
      for (int i=0;i<m-1;++i)
        _a[k+i]=T(0);
      _a[k+m-1]=_x[ii];
    }
  }
}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,UpperLowerFlag UPLO,int N>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_vx.n()==_ma.m());

  ld_zeros(_ma,_a);
  int n=_ma.n(), i=0, inc=_vx.inc();
  int m=n,ii=0;

  for (int j=0;j<n;++j,--m,ii+=inc) {
    _a[i]=_alpha*_x[ii];
    i+=m;
  }
}

/// load (scaled) identity matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,
          UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const tr_mat<UPLO,TRANS,DIAG,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  static_assert(DIAG==NoU,"cannot modify constant unit diagonal");
  assert(_vx.n()==_ma.m());

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld(), inc=_vx.inc();
    int k=0,ii=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0,ii+=inc) {
      _a[k+m0]=_alpha*_x[ii];
      for (int i=m0+1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    int m=0,ld=_ma.ld(), inc=_vx.inc();
    int k=0,ii=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m,ii+=inc) {
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
      _a[k+m]=_alpha*_x[ii];
    }
  }
}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,
          UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int K,int LD>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const tb_mat<UPLO,TRANS,DIAG,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  static_assert(DIAG==NoU,"cannot modify constant unit diagonal");
  assert(_vx.n()==_ma.m());

  int ld=_ma.ld(), inc=_vx.inc();
  int m=_ma.k()+1;
  int k=0,ii=0;

  if (UPLO==Lower) {
    for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc) {
      _a[k]=_alpha*_x[ii];
      for (int i=1;i<m;++i)
        _a[k+i]=T(0);
    }
  }
  else {
    assert(UPLO==Upper);
    for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc) {
      for (int i=0;i<m-1;++i)
        _a[k+i]=T(0);
      _a[k+m-1]=_alpha*_x[ii];
    }
  }

}

/// load (scaled) diagonal matrix \ingroup vc_blas_mat
template <typename T,int NX,int INCX,UpperLowerFlag UPLO,
          TransposeFlag TRANS,DiagonalFlag DIAG,int N>
void ld_diag(const T& _alpha,
             const vec<NX,INCX>& _vx,const T* _x,
             const tp_mat<UPLO,TRANS,DIAG,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  static_assert(DIAG==NoU,"cannot modify constant unit diagonal");
  assert(_vx.n()==_ma.m());

  ld_zeros(_ma,_a);
  int n=_ma.n(), i=0, inc=_vx.inc();
  int m=n,ii=0;

  for (int j=0;j<n;++j,--m,ii+=inc) {
    _a[i]=_alpha*_x[ii];
    i+=m;
  }
}

//-----------------------------------------------------------------------------

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void mscal(const T& _alpha,const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int ld=_ma.ld();
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<_ma.m();++i)
      _a[k+i]*=_alpha;
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
void mscal(const T& _alpha,const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int ld=_ma.ld();
  int m=_ma.kl()+_ma.ku()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]*=_alpha;
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD>
void mscal(const T& _alpha,const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]*=_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]*=_alpha;
  }
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
void mscal(const T& _alpha,const sb_mat<UPLO,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]*=_alpha;
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N>
void mscal(const T& _alpha,const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]*=_alpha;
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int LD>
void mscal(const T& _alpha,const tr_mat<UPLO,TRANS,NoU,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]*=_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]*=_alpha;
  }
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int K,int LD>
void mscal(const T& _alpha,const tb_mat<UPLO,TRANS,NoU,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]*=_alpha;
}

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N>
void mscal(const T& _alpha,const tp_mat<UPLO,TRANS,NoU,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (_alpha==T(1)) return;

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]*=_alpha;
}

//-----------------------------------------------------------------------------

/// add scalar to all elements  \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void adds(const T& _alpha,const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld();
  int k=0;
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<_ma.m();++i)
      _a[k+i]+=_alpha;
}

/// add scalar to all elements  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N,int LD>
void adds(const T& _alpha,const sy_mat<UPLO,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]+=_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]+=_alpha;
  }
}

/// add scalar to all elements  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,int N>
void adds(const T& _alpha,const sp_mat<UPLO,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]+=_alpha;
}

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

template <typename T>
struct _EmulOp {
  _EmulOp(const T& _alpha) : alpha(_alpha) {}
  T alpha;
  void operator()(T& _a,const T& _b) { _a*=_b*alpha; }
};
template <typename T>
struct _EdivOp {
  _EdivOp(const T& _alpha) : alpha(_alpha) {}
  T alpha;
  void operator()(T& _a,const T& _b) { _a/=_b*alpha; }
};

# endif

/// element-wise multiplication _a*=_alpha*_b (note: few exceptions to MA==MB) \ingroup vc_blas_mat
template <typename T,typename MA,typename MB>
void emul(const T& _alpha,
           const MB& _mb,const T* _b,const MA& _ma,T* _a) {
  __assert_no_unit_diag<MA>();
  Foreach2<T,_EmulOp<T>,MA,MB>(_ma,_a,_mb,(T*) _b,_EmulOp<T>(_alpha));
}
/// element-wise division _a/=_alpha*_b (note: few exceptions to MA==MB) \ingroup vc_blas_mat
template <typename T,typename MA,typename MB>
void ediv(const T& _alpha,
           const MB& _mb,const T* _b,const MA& _ma,T* _a) {
  // require assert on unit diagonal
  Foreach2<T,_EdivOp<T>,MA,MB>(_ma,_a,_mb,(T*) _b,_EdivOp<T>(_alpha));
}

//-----------------------------------------------------------------------------

//
// NOTE: We could simplify madd() by using Foreach2 (single implementation).
//       See above (emul(), ediv()).
//

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANSB,
          int MA,int NA,int LDA,int MB,int NB,int LDB>
void madd(const T& _alpha,
          const ge_mat<TRANSB,MB,NB,LDB>& _mb,const T* _b,
          const ge_mat<NoT,MA,NA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  if (TRANSB==NoT) {
    assert(_ma.m()==_mb.m());
    assert(_ma.n()==_mb.n());

    int lda=_ma.ld(), ldb=_mb.ld();
    int ka=0,kb=0;
    for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb)
      for (int i=0;i<_ma.m();++i)
        _a[ka+i]+=_alpha*_b[kb+i];
  }
  else {
    assert(_a!=_b && "cannot madd(A,trans(A))");
    assert(_ma.n()==_mb.m());
    assert(_ma.m()==_mb.n());

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<_ma.m();++i)
        _a[j*_ma.ld()+i]+=_alpha*_b[i*_mb.ld()+j];
  }
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,
          int MA,int NA,int KLA,int KUA,int LDA,
          int MB,int NB,int KLB,int KUB,int LDB>
void madd(const T& _alpha,
          const gb_mat<TRANS,MB,NB,KLB,KUB,LDB>& _mb,const T* _b,
          const gb_mat<TRANS,MA,NA,KLA,KUA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.m()==_mb.m());
  assert(_ma.n()==_mb.n());
  assert(_ma.kl()==_mb.kl());
  assert(_ma.ku()==_mb.ku());

  int lda=_ma.ld(), ldb=_mb.ld();
  int m=_ma.kl()+_ma.ku()+1;
  int ka=0,kb=0;
  for (int j=0;j<_ma.n();++j,ka+=lda,kb=ldb)
    for (int i=0;i<m;++i)
      _a[ka+i]+=_alpha*_b[kb+i];
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
  template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NB,int LDB>
void madd(const T& _alpha,
          const sy_mat<UPLO,NB,LDB>& _mb,const T* _b,
          const sy_mat<UPLO,NA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_mb.n());

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),lda=_ma.ld(), ldb=_mb.ld();
    int ka=0,kb=0;
    for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb,++m0)
      for (int i=m0;i<m;++i)
        _a[ka+i]+=_alpha*_b[kb+i];
  }
  else {
    assert(UPLO==Upper);
    int m=1,lda=_ma.ld(),ldb=_mb.ld();
    int ka=0,kb=0;
    for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb,++m)
      for (int i=0;i<m;++i)
        _a[ka+i]+=_alpha*_b[kb+i];
  }
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,
          int NA,int KA,int LDA,int NB,int KB,int LDB>
void madd(const T& _alpha,
          const sb_mat<UPLO,NB,KB,LDB>& _mb,const T* _b,
          const sb_mat<UPLO,NA,KA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_mb.n());
  assert(_ma.k()==_mb.k());

  int lda=_ma.ld(),ldb=_mb.ld();
  int m=_ma.k()+1;
  int ka=0,kb=0;
  for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb)
    for (int i=0;i<m;++i)
      _a[ka+i]+=_alpha*_b[kb+i];
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
  template <typename T,UpperLowerFlag UPLO,int NA,int NB>
void madd(const T& _alpha,
          const sp_mat<UPLO,NB>& _mb,const T* _b,
          const sp_mat<UPLO,NA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_mb.n());

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]+=_alpha*_b[i];
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,
          int NA,int LDA,int NB,int LDB>
void madd(const T& _alpha,
          const tr_mat<UPLO,TRANS,NoU,NB,LDB>& _mb,const T* _b,
          const tr_mat<UPLO,TRANS,NoU,NA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_mb.n());

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),lda=_ma.ld(),ldb=_mb.ld();
    int ka=0,kb=0;
    for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb,++m0)
      for (int i=m0;i<m;++i)
        _a[ka+i]+=_alpha*_b[kb+i];
  }
  else {
    assert(UPLO==Upper);
    int m=1,lda=_ma.lda(),ldb=_mb.ldb();
    int ka=0,kb=0;
    for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb,++m)
      for (int i=0;i<m;++i)
        _a[ka+i]+=_alpha*_b[kb+i];
  }
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,
          int NA,int KA,int LDA,int NB,int KB,int LDB>
void madd(const T& _alpha,
          const tb_mat<UPLO,TRANS,NoU,NB,KB,LDB>& _mb,const T* _b,
          const tb_mat<UPLO,TRANS,NoU,NA,KA,LDA>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_mb.n());
  assert(_ma.k()==_mb.k());

  int lda=_ma.ld(),ldb=_mb.ld();
  int m=_ma.k()+1;
  int ka=0,kb=0;
  for (int j=0;j<_ma.n();++j,ka+=lda,kb+=ldb)
    for (int i=0;i<m;++i)
      _a[ka+i]+=_alpha*_b[kb+i];
}

/// addition _a+=_alpha*_b \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int NA,int NB>
void madd(const T& _alpha,
          const tp_mat<UPLO,TRANS,NoU,NB>& _mb,const T* _b,
          const tp_mat<UPLO,TRANS,NoU,NA>& _ma,T* _a) {

  assert(_ma.n()==_mb.n());

  int sz=_ma.n();
  sz*=sz+1;
  sz/=2;
  for (int i=0;i<sz;++i)
    _a[i]+=_alpha*_b[i];
}

//-----------------------------------------------------------------------------

/// right-multiply by diagonal matrix \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N,int LD,int NX,int INCX>
void mscal_cols(const T& _alpha,
                const vec<NX,INCX>& _vd,const T* _d,
                const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int ld=_ma.ld(),inc=_vd.inc();
  int k=0;

  if (TRANS==NoT) {
    assert(_ma.n()==_vd.n());
    int ii=0;
    for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc)
      for (int i=0;i<_ma.m();++i)
        _a[k+i]*=_d[ii]*_alpha;
  }
  else {
    assert(_ma.m()==_vd.n());
    for (int j=0;j<_ma.n();++j,k+=ld)
      for (int i=0;i<_ma.m();++i)
        _a[k+i]*=_d[i*inc]*_alpha;
  }
}


/// right-multiply by diagonal matrix \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,
          int M,int N,int KL,int KU,int LD,int NX,int INCX>
void mscal_cols(const T& _alpha,
                const vec<NX,INCX>& _vd,const T* _d,
                const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(TRANS==NoT && "not implemented");
  assert(_ma.n()==_vd.n());

  int ld=_ma.ld(),inc=_vd.inc();
  int m=_ma.kl()+_ma.ku()+1;
  int k=0,ii=0;
  for (int j=0;j<_ma.n();++j,k+=ld,ii+=inc)
    for (int i=0;i<m;++i)
      _a[k+i]*=_d[ii]*_alpha;
}

  // not implemented: const gb_mat< Transpose ,M,N,KL,KU,LD>

/// multiply by scalar  \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int LD,
          int NX,int INCX>
void mscal_cols(const T& _alpha,
                const vec<NX,INCX>& _vd,const T* _d,
                const tr_mat<UPLO,TRANS,NoU,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();
  assert(_ma.n()==_vd.n());

  int inc=_vd.inc();

  if (UPLO==Lower) {
    int m0=0,m=_ma.n(),ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m0)
      for (int i=m0;i<m;++i)
        _a[k+i]*=_d[((TRANS==NoT) ? j : i)*inc]*_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1,ld=_ma.ld();
    int k=0;
    for (int j=0;j<_ma.n();++j,k+=ld,++m)
      for (int i=0;i<m;++i)
        _a[k+i]*=_d[((TRANS==NoT) ? j : i)*inc]*_alpha;
  }
}

/// right-multiply by diagonal matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,
          int N,int K,int LD,int NX,int INCX>
void mscal_cols(const T& _alpha,
                const vec<NX,INCX>& _vd,const T* _d,
                const tb_mat<UPLO,TRANS,NoU,N,K,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(TRANS==NoT && "not implemented");
  assert(_ma.n()==_vd.n());

  int ld=_ma.ld();
  int m=_ma.k()+1;
  int k=0,inc=_vd.inc();
  for (int j=0;j<_ma.n();++j,k+=ld)
    for (int i=0;i<m;++i)
      _a[k+i]*=_d[j*inc]*_alpha;
}

/// right-multiply by diagonal matrix \ingroup vc_blas_mat
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,int N,int NX,int INCX>
void mscal_cols(const T& _alpha,
                const vec<NX,INCX>& _vd,const T* _d,
                const tp_mat<UPLO,TRANS,NoU,N>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.n()==_vd.n());

  int n=_ma.n(),inc=_vd.inc(),k=0;

  if (UPLO==Lower) {
    int m=_ma.n();
    for (int j=0;j<n;++j,--m)
      for (int i=0;i<m;++i,++k)
        _a[k]*=_d[j*inc]*_alpha;
  }
  else {
    assert(UPLO==Upper);
    int m=1;
    for (int j=0;j<n;++j,++m)
      for (int i=0;i<m;++i,++k)
        _a[k]*=_d[j*inc]*_alpha;
  }
}

//-----------------------------------------------------------------------------

/// add transposed: A=(A+A')*_alpha
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
void add_trans(const T& _alpha,
               const ge_mat<TRANS,M,N,LD>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(_ma.m()==_ma.n());

  int ld=_ma.ld();
  for (int j=0;j<_ma.n();++j)
    for (int i=0;i<=j;++i) {
      T x=_a[j*ld+i], y=_a[i*ld+j];
      _a[j*ld+i]=_a[i*ld+j]=(x+y)*_alpha;
    }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

  //
  // specializations
  //

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS>
void ld_zero(const ge_mat<TRANS,VarInt,VarInt,VarInt>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  int m=_ma.m(), n=_ma.n(), ld=_ma.ld();

  if (ld>m) {
    int k=0;
    for (int j=0;j<n;++j,k+=ld)
      for (int i=0;i<m;++i)
        _a[k+i]=T(0);
  }
  else {
    assert(ld==m);
    int sz=m*n;
    for (int i=0;i<sz;++i)
      _a[i]=T(0);
  }
}

/// set all elements to zero \ingroup vc_blas_mat
template <typename T,TransposeFlag TRANS,int M,int N>
void ld_zero(const ge_mat<TRANS,M,N,M>& _ma,T* _a) {
  _VC_DBG_BLAS_SCOPE();

  assert(M!=VarInt && N!=VarInt && "something went wrong with specialization");
  int sz=_ma.m()*_ma.n();
  for (int i=0;i<sz;++i)
    _a[i]=T(0);
}

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
