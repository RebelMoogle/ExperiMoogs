//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================


#ifndef __VC_MATH_BLAS_HH
#define __VC_MATH_BLAS_HH

#include "lapack_types.hh"
#include "blas_flags.hh"

# ifndef VC_NO_BLAS_REPLACEMENTS
#  include "vcblas.hh"
# endif

namespace VC {
namespace math {

/// BLAS wrappers, helpers, prototypes, matrix package \sa [\ref vc_blas] \ingroup vc_blas
namespace blas {

/** \defgroup vc_blas BLAS: C++ interface
    \ingroup vc_math

    \sa VC::math::blas
 */

 /** \def VC_NO_BLAS_REPLACEMENTS
     \ingroup vc_blas
     Don't use BLAS replacements [\ref vc_vcblas] within BLAS wrappers
     [\ref vc_blas].

     If `VC_NO_BLAS_REPLACEMENTS` is *not* defined (default setting),
     then the C++ BLAS (e.g., VC::math::blas::gemv() for
     `DGEMV`,`SGEMV`) may call replacements (e.g.,
     VC::math::blas::vc_gemv()) for fixed and/or small dimensions and
     let the compiler inline and optimize code.

     \sa [\ref vc_vcblas], [\ref vc_blas]
   */

# ifdef DOXYGEN_SKIP
# ifndef VC_NO_BLAS_REPLACEMENTS
#  define VC_NO_BLAS_REPLACEMENTS
# endif
#  error "doxygen only"
# undef VC_NO_BLAS_REPLACEMENTS
# endif

typedef lapack::INTEGER INTEGER;   //!< integer type used by BLAS

} // namespace blas
} // namespace math
} // namespace VC

#include "blas_prototypes.hh"

//=============================================================================

namespace VC {
namespace math {
namespace blas {

  //
  // C++ wrappers to low-level blas functions
  //


//
// double precision
//
//-----------------------------------------------------------------------------

/** \defgroup vc_blas_1 BLAS level 1 C++ wrapper: vector-vector operations
    \ingroup vc_blas
 */
/** \defgroup vc_blas_2 BLAS level 2 C++ wrapper: matrix-vector operations
    \ingroup vc_blas
 */
/** \defgroup vc_blas_3 BLAS level 3 C++ wrapper: matrix-matrix operations
    \ingroup vc_blas
 */

/// limits for replacing BLAS calls \ingroup vc_vcblas
enum _VC_BLAS_LIMTS {
  _L1_NMAX=512,
  _L2_NMAX=16,
  _L3_NMAX=16
};

  //
  // level 1
  //

/** BLAS swap _x and _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dswap.f">[reference]</a>
*/
inline void
swap(int _n,
     double* _x,int _incx, double* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_swap(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/DSWAP(&n,_x,&incx,_y,&incy);
}

/** BLAS scal _x = _alpha*_x
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dscal.f">[reference]</a>
*/
inline void
scal(int _n,double _alpha, double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_scal(_n,_alpha,_x,(_incx==1 ? 1 : _incx));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx;
  /**/DSCAL(&n,&_alpha,_x,&incx);
}

/** BLAS copy _y = _x
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dcopy.f">[reference]</a>
*/
inline void
copy(int _n,
     const double* _x,int _incx, double* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_copy(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/DCOPY(&n,_x,&incx,_y,&incy);
}

/** BLAS AXPY _y = _alpha*_x+_y.
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/daxpy.f">[reference]</a>
*/
inline void
axpy(int _n,double _alpha,const double* _x,int _incx,double* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_axpy(_n,_alpha,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/DAXPY(&n,&_alpha,_x,&incx,_y,&incy);
}

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
inline double
dot(int _n,
    const double* _x,int _incx, const double* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_dot(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  return /**/DDOT(&n,_x,&incx,_y,&incy);
}

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
inline double
ddot(int _n,
     const double* _x,int _incx, const double* _y,int _incy) {
  return dot(_n,_x,_incx,_y,_incy);
}

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
inline double
ddot(int _n,
     const float* _x,int _incx, const float* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_ddot(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  return /**/DSDOT(&n,_x,&incx,_y,&incy);
}


/** BLAS 2-norm ||_x||_2
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dnrm2.f">[reference]</a>
*/
inline double
nrm2(int _n, const double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_nrm2(_n,_x,(_incx==1 ? 1 : _incx));
  }
# endif
  INTEGER n=_n,incx=_incx;
  return /**/DNRM2(&n,_x,&incx);
}

/** BLAS 1-norm || _x || _1
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dasum.f">[reference]</a>
*/
inline double
asum(int _n, const double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_asum(_n,_x,(_incx==1 ? 1 : _incx));
  }
# endif
  INTEGER n=_n,incx=_incx;
  return /**/DASUM(&n,_x,&incx);
}


//
// level 2
//

/** BLAS DGEMV _y=_alpha*op(_A)*_x+_beta*_y.
    \ingroup vc_blas_2
     <a href="http://www.netlib.org/blas/dgemv.f">[reference]</a>
*/
inline void
gemv(TransposeFlag _trans,int _m,int _n,
     double _alpha,
     const double* _A,int _lda,const double* _x,int _incx,
     double _beta,
     double* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<=_L2_NMAX && _n<_L2_NMAX) {
    vc_gemv(_trans,_m,_n,_alpha,_A,_lda,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char trans=char(_trans);
  INTEGER m=_m,n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/DGEMV(&trans,&m,&n,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSYMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsymv.f">[reference]</a>
*/
inline void
symv(UpperLowerFlag _uplo,int _n,
     double _alpha,
     const double* _A,int _lda,
     const double* _x,int _incx,
     double _beta,
     double* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_symv(_uplo,_n,_alpha,_A,_lda,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/DSYMV(&uplo,&n,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSBMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsbmv.f">[reference]</a>
*/
inline void
sbmv(UpperLowerFlag _uplo,int _n,int _k,
     double _alpha,
     const double* _A,int _lda,
     const double* _x,int _incx,
     double _beta,
     double* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo);
  INTEGER k=_k,n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/DSBMV(&uplo,&n,&k,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSPMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspmv.f">[reference]</a>
*/
inline void
spmv(UpperLowerFlag _uplo,int _n,
     double _alpha,const double* _Ap,
     const double* _x,int _incx,
     double _beta,
     double* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spmv(_uplo,_n,_alpha,_Ap,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/DSPMV(&uplo,&n,&_alpha,_Ap,_x,&incx,&_beta,_y,&incy);
}


/** BLAS DTRMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtrmv.f">[reference]</a>
*/
inline void
trmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const double* _A,int _lda,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_trmv(_uplo,_trans,_diag,_n,_A,_lda,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/DTRMV(&uplo,&trans,&diag,&n,_A,&lda,_x,&incx);
}

/** BLAS DTBMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtbmv.f">[reference]</a>
*/
inline void
tbmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,int _k,
     const double* _A,int _lda,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,k=_k,lda=_lda,incx=_incx;
  /**/DTBMV(&uplo,&trans,&diag,&n,&k,_A,&lda,_x,&incx);
}

/** BLAS DTPMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtpmv.f">[reference]</a>
*/
inline void
tpmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const double* _Ap,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_tpmv(_uplo,_trans,_diag,_n,_Ap,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,incx=_incx;
  /**/DTPMV(&uplo,&trans,&diag,&n,_Ap,_x,&incx);
}



/** BLAS DTRSV _x=op(_A)^(-1)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtrsv.f">[reference]</a>
*/
inline void
trsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const double* _A,int _lda,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_trsv(_uplo,_trans,_diag,_n,_A,_lda,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/DTRSV(&uplo,&trans,&diag,&n,_A,&lda,_x,&incx);
}


/** BLAS DTBSV _x=op(_A)^(-1)*_x
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtbsv.f">[reference]</a>
*/
inline void
tbsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,int _k,
     const double* _A,int _lda,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,k=_k,lda=_lda,incx=_incx;
  /**/DTBSV(&uplo,&trans,&diag,&n,&k,_A,&lda,_x,&incx);
}

/** BLAS DTPSV _A=_x*_x^(T)+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtpsv.f">[reference]</a>
 */
inline void
tpsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const double* _Ap,
     double* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_tpsv(_uplo,_trans,_diag,_n,_Ap,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,incx=_incx;
  /**/DTPSV(&uplo,&trans,&diag,&n,_Ap,_x,&incx);
}



/** BLAS DSYR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsyr.f">[reference]</a>
*/
inline void
syr(UpperLowerFlag _uplo,int _n,
    double _alpha,
    const double* _x,int _incx,
    double* _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_syr(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),_A,_lda);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/DSYR(&uplo,&n,&_alpha,_x,&incx,_A,&lda);
}

/** BLAS DSPR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspr.f">[reference]</a>
*/
inline void
spr(UpperLowerFlag _uplo,int _n,
    double _alpha,
    const double* _x,int _incx,
    double* _Ap) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spr(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),_Ap);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx;
  /**/DSPR(&uplo,&n,&_alpha,_x,&incx,_Ap);
}

/** BLAS DSYR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsyr2.f">[reference]</a>
*/
inline void
syr2(UpperLowerFlag _uplo,int _n,
     double _alpha,
     const double* _x,int _incx,
     const double* _y,int _incy,
     double* _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_syr2(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),
            _y,(_incy==1 ? 1 : _incy),_A,_lda);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/DSYR2(&uplo,&n,&_alpha,_x,&incx,_y,&incy,_A,&lda);
}

/** BLAS DSPR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspr2.f">[reference]</a>
*/
inline void
spr2(UpperLowerFlag _uplo,int _n,
     double _alpha,
     const double* _x,int _incx,
     const double* _y,int _incy,
     double* _Ap) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spr2(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),
            _y,(_incy==1 ? 1 : _incy),_Ap);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/DSPR2(&uplo,&n,&_alpha,_x,&incx,_y,&incy,_Ap);
}

//
// level 3
//

/** BLAS DGEMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dgemm.f">[reference]</a>
*/
inline void
gemm(TransposeFlag _transA,TransposeFlag _transB,int _m,int _n,int _k,
     double _alpha,
     const double* _A,int _lda,
     const double* _B,int _ldb,
     double _beta,
     double* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<_L2_NMAX && _n<_L2_NMAX && _k<_L2_NMAX) {
    vc_gemm(_transA,_transB,_m,_n,_k,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char transA=char(_transA),transB=char(_transB);
  INTEGER m=_m,n=_n,k=_k,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/DGEMM(&transA,&transB,&m,&n,&k,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DSYMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsymm.f">[reference]</a>
*/
inline void
symm(SideFlag _side,UpperLowerFlag _uplo,int _m,int _n,
     double _alpha,
     const double* _A,int _lda,
     const double* _B,int _ldb,
     double _beta,
     double* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<_L2_NMAX && _n<_L2_NMAX) {
    vc_symm(_side,_uplo,_m,_n,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/DSYMM(&side,&uplo,&m,&n,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DSYRK _C=_alpha*_A*_A^T+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsyrk.f">[reference]</a>
*/
inline void
syrk(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
     double _alpha,
     const double* _A,int _lda,
     double _beta,
     double* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _k<_L2_NMAX) {
    vc_syrk(_uplo,_trans,_n,_k,_alpha,_A,_lda,_beta,_C,_ldc);
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans);
  INTEGER k=_k,n=_n,lda=_lda,ldc=_ldc;
  /**/DSYRK(&uplo,&trans,&n,&k,&_alpha,_A,&lda,&_beta,_C,&ldc);
}

/** BLAS DSYR2K _C=_alpha*_A*_B^T+_alpha*_B*_A^T+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsyr2k.f">[reference]</a>
*/
inline void
syr2k(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
      double _alpha,
      const double* _A,int _lda,
      const double* _B,int _ldb,
      double _beta,
      double* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _k<_L2_NMAX) {
    vc_syr2k(_uplo,_trans,_n,_k,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans);
  INTEGER k=_k,n=_n,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/DSYR2K(&uplo,&trans,&n,&k,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DTRMM _B=_alpha*op(_A)*op(_B).
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dtrmm.f">[reference]</a>
*/
inline void
trmm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
     int _m,int _n,
     double _alpha,
     const double* _A,int _lda,
     double* _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _n<_L2_NMAX) {
    vc_trmm(_side,_uplo,_transA,_diag,_m,_n,_alpha,_A,_lda,_B,_ldb);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo),transA=char(_transA),diag=char(_diag);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb;
  /**/DTRMM(&side,&uplo,&transA,&diag,&m,&n,&_alpha,_A,&lda,_B,&ldb);
}

/** blas DTRSM _B=_alpha*op(_A^-1)*op(_B).
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dtrsm.f">[reference]</a>
*/
inline void
trsm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
     int _m,int _n,
     double _alpha,
     const double* _A,int _lda,
     double* _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _n<_L2_NMAX) {
    vc_trsm(_side,_uplo,_transA,_diag,_m,_n,_alpha,_A,_lda,_B,_ldb);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo),transA=char(_transA),diag=char(_diag);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb;
  /**/DTRSM(&side,&uplo,&transA,&diag,&m,&n,&_alpha,_A,&lda,_B,&ldb);
}

//
// single precision
//
//-----------------------------------------------------------------------------

  //
  // level 1
  //

/** BLAS swap _x and _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dswap.f">[reference]</a>
*/
inline void
swap(int _n,
     float* _x,int _incx, float* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_swap(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/SSWAP(&n,_x,&incx,_y,&incy);
}


/** BLAS scal _x = _alpha*_x
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dscal.f">[reference]</a>
*/
inline void
scal(int _n,float _alpha, float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_scal(_n,_alpha,_x,(_incx==1 ? 1 : _incx));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx;
  /**/SSCAL(&n,&_alpha,_x,&incx);
}

/** BLAS copy _y = _x
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dcopy.f">[reference]</a>
*/
inline void
copy(int _n,
     const float* _x,int _incx, float* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_copy(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/SCOPY(&n,_x,&incx,_y,&incy);
}

/** BLAS AXPY _y = _alpha*_x+_y.
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/daxpy.f">[reference]</a>
*/
inline void
axpy(int _n,float _alpha,const float* _x,int _incx,float* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    vc_axpy(_n,_alpha,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
    return;
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/SAXPY(&n,&_alpha,_x,&incx,_y,&incy);
}

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
inline float
dot(int _n,
    const float* _x,int _incx, const float* _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_dot(_n,_x,(_incx==1 ? 1 : _incx), _y,(_incy==1 ? 1 : _incy));
  }
# endif
  INTEGER n=_n,incx=_incx,incy=_incy;
  return /**/SDOT(&n,_x,&incx,_y,&incy);
}

/** BLAS 2-norm ||_x||_2
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dnrm2.f">[reference]</a>
*/
inline float
nrm2(int _n, const float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_nrm2(_n,_x,(_incx==1 ? 1 : _incx));
  }
# endif
  INTEGER n=_n,incx=_incx;
  return /**/SNRM2(&n,_x,&incx);
}

/** BLAS 1-norm || _x || _1
    \ingroup vc_blas_1
    <a href="http://www.netlib.org/blas/dasum.f">[reference]</a>
*/
inline float
asum(int _n, const float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<=_L1_NMAX) {
    return vc_asum(_n,_x,(_incx==1 ? 1 : _incx));
  }
# endif
  INTEGER n=_n,incx=_incx;
  return /**/SASUM(&n,_x,&incx);
}


//
// level 2
//

/** BLAS DGEMV _y=_alpha*op(_A)*_x+_beta*_y.
    \ingroup vc_blas_2
     <a href="http://www.netlib.org/blas/dgemv.f">[reference]</a>
*/
inline void
gemv(TransposeFlag _trans,int _m,int _n,
     float _alpha,
     const float* _A,int _lda,const float* _x,int _incx,
     float _beta,
     float* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<=_L2_NMAX && _n<_L2_NMAX) {
    vc_gemv(_trans,_m,_n,_alpha,_A,_lda,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char trans=char(_trans);
  INTEGER m=_m,n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/SGEMV(&trans,&m,&n,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSYMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsymv.f">[reference]</a>
*/
inline void
symv(UpperLowerFlag _uplo,int _n,
     float _alpha,
     const float* _A,int _lda,
     const float* _x,int _incx,
     float _beta,
     float* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_symv(_uplo,_n,_alpha,_A,_lda,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/SSYMV(&uplo,&n,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSBMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsbmv.f">[reference]</a>
*/
inline void
sbmv(UpperLowerFlag _uplo,int _n,int _k,
     float _alpha,
     const float* _A,int _lda,
     const float* _x,int _incx,
     float _beta,
     float* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo);
  INTEGER k=_k,n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/SSBMV(&uplo,&n,&k,&_alpha,_A,&lda,_x,&incx,&_beta,_y,&incy);
}

/** BLAS DSPMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspmv.f">[reference]</a>
*/
inline void
spmv(UpperLowerFlag _uplo,int _n,
     float _alpha,const float* _Ap,
     const float* _x,int _incx,
     float _beta,
     float* _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spmv(_uplo,_n,_alpha,_Ap,
            _x,(_incx==1 ? 1 : _incx),_beta,_y,(_incy==1 ? 1 : _incy));
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/SSPMV(&uplo,&n,&_alpha,_Ap,_x,&incx,&_beta,_y,&incy);
}


/** BLAS DTRMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtrmv.f">[reference]</a>
*/
inline void
trmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const float* _A,int _lda,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_trmv(_uplo,_trans,_diag,_n,_A,_lda,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/STRMV(&uplo,&trans,&diag,&n,_A,&lda,_x,&incx);
}

/** BLAS DTBMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtbmv.f">[reference]</a>
*/
inline void
tbmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,int _k,
     const float* _A,int _lda,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,k=_k,lda=_lda,incx=_incx;
  /**/STBMV(&uplo,&trans,&diag,&n,&k,_A,&lda,_x,&incx);
}

/** BLAS DTPMV _x=op(_A)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtpmv.f">[reference]</a>
*/
inline void
tpmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const float* _Ap,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_tpmv(_uplo,_trans,_diag,_n,_Ap,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,incx=_incx;
  /**/STPMV(&uplo,&trans,&diag,&n,_Ap,_x,&incx);
}



/** BLAS DTRSV _x=op(_A)^(-1)*_x.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtrsv.f">[reference]</a>
*/
inline void
trsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const float* _A,int _lda,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_trsv(_uplo,_trans,_diag,_n,_A,_lda,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/STRSV(&uplo,&trans,&diag,&n,_A,&lda,_x,&incx);
}


/** BLAS DTBSV _x=op(_A)^(-1)*_x
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtbsv.f">[reference]</a>
*/
inline void
tbsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,int _k,
     const float* _A,int _lda,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,k=_k,lda=_lda,incx=_incx;
  /**/STBSV(&uplo,&trans,&diag,&n,&k,_A,&lda,_x,&incx);
}

/** BLAS DTPSV _A=_x*_x^(T)+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dtpsv.f">[reference]</a>
 */
inline void
tpsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const float* _Ap,
     float* _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_tpsv(_uplo,_trans,_diag,_n,_Ap,_x,(_incx==1 ? 1 : _incx));
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans), diag=char(_diag);
  INTEGER n=_n,incx=_incx;
  /**/STPSV(&uplo,&trans,&diag,&n,_Ap,_x,&incx);
}



/** BLAS DSYR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsyr.f">[reference]</a>
*/
inline void
syr(UpperLowerFlag _uplo,int _n,
    float _alpha,
    const float* _x,int _incx,
    float* _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_syr(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),_A,_lda);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx;
  /**/SSYR(&uplo,&n,&_alpha,_x,&incx,_A,&lda);
}

/** BLAS DSPR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspr.f">[reference]</a>
*/
inline void
spr(UpperLowerFlag _uplo,int _n,
    float _alpha,
    const float* _x,int _incx,
    float* _Ap) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spr(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),_Ap);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx;
  /**/SSPR(&uplo,&n,&_alpha,_x,&incx,_Ap);
}

/** BLAS DSYR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dsyr2.f">[reference]</a>
*/
inline void
syr2(UpperLowerFlag _uplo,int _n,
     float _alpha,
     const float* _x,int _incx,
     const float* _y,int _incy,
     float* _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_syr2(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),
            _y,(_incy==1 ? 1 : _incy),_A,_lda);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,incx=_incx,incy=_incy;
  /**/SSYR2(&uplo,&n,&_alpha,_x,&incx,_y,&incy,_A,&lda);
}

/** BLAS DSPR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_2
    <a href="http://www.netlib.org/blas/dspr2.f">[reference]</a>
*/
inline void
spr2(UpperLowerFlag _uplo,int _n,
     float _alpha,
     const float* _x,int _incx,
     const float* _y,int _incy,
     float* _Ap) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX) {
    vc_spr2(_uplo,_n,_alpha,_x,(_incx==1 ? 1 : _incx),
            _y,(_incy==1 ? 1 : _incy),_Ap);
    return ;
  }
# endif
  char uplo=char(_uplo);
  INTEGER n=_n,incx=_incx,incy=_incy;
  /**/SSPR2(&uplo,&n,&_alpha,_x,&incx,_y,&incy,_Ap);
}

//
// level 3
//

/** BLAS DGEMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dgemm.f">[reference]</a>
*/
inline void
gemm(TransposeFlag _transA,TransposeFlag _transB,int _m,int _n,int _k,
     float _alpha,
     const float* _A,int _lda,
     const float* _B,int _ldb,
     float _beta,
     float* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<_L2_NMAX && _n<_L2_NMAX && _k<_L2_NMAX) {
    vc_gemm(_transA,_transB,_m,_n,_k,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char transA=char(_transA),transB=char(_transB);
  INTEGER m=_m,n=_n,k=_k,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/SGEMM(&transA,&transB,&m,&n,&k,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DSYMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsymm.f">[reference]</a>
*/
inline void
symm(SideFlag _side,UpperLowerFlag _uplo,int _m,int _n,
     float _alpha,
     const float* _A,int _lda,
     const float* _B,int _ldb,
     float _beta,
     float* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_m<_L2_NMAX && _n<_L2_NMAX) {
    vc_symm(_side,_uplo,_m,_n,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/SSYMM(&side,&uplo,&m,&n,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DSYRK _C=_alpha*_A*_A^T+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsyrk.f">[reference]</a>
*/
inline void
syrk(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
     float _alpha,
     const float* _A,int _lda,
     float _beta,
     float* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _k<_L2_NMAX) {
    vc_syrk(_uplo,_trans,_n,_k,_alpha,_A,_lda,_beta,_C,_ldc);
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans);
  INTEGER k=_k,n=_n,lda=_lda,ldc=_ldc;
  /**/SSYRK(&uplo,&trans,&n,&k,&_alpha,_A,&lda,&_beta,_C,&ldc);
}

/** BLAS DSYR2K _C=_alpha*_A*_B^T+_alpha*_B*_A^T+_beta*_C.
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dsyr2k.f">[reference]</a>
*/
inline void
syr2k(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
      float _alpha,
      const float* _A,int _lda,
      const float* _B,int _ldb,
      float _beta,
      float* _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _k<_L2_NMAX) {
    vc_syr2k(_uplo,_trans,_n,_k,_alpha,_A,_lda,_B,_ldb,_beta,_C,_ldc);
    return ;
  }
# endif
  char uplo=char(_uplo),trans=char(_trans);
  INTEGER k=_k,n=_n,lda=_lda,ldb=_ldb,ldc=_ldc;
  /**/SSYR2K(&uplo,&trans,&n,&k,&_alpha,_A,&lda,_B,&ldb,&_beta,_C,&ldc);
}

/** BLAS DTRMM _B=_alpha*op(_A)*op(_B).
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dtrmm.f">[reference]</a>
*/
inline void
trmm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
     int _m,int _n,
     float _alpha,
     const float* _A,int _lda,
     float* _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _n<_L2_NMAX) {
    vc_trmm(_side,_uplo,_transA,_diag,_m,_n,_alpha,_A,_lda,_B,_ldb);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo),transA=char(_transA),diag=char(_diag);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb;
  /**/STRMM(&side,&uplo,&transA,&diag,&m,&n,&_alpha,_A,&lda,_B,&ldb);
}

/** blas DTRSM _B=_alpha*op(_A^-1)*op(_B).
    \ingroup vc_blas_3
    <a href="http://www.netlib.org/blas/dtrsm.f">[reference]</a>
*/
inline void
trsm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
     int _m,int _n,
     float _alpha,
     const float* _A,int _lda,
     float* _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()
# ifndef VC_NO_BLAS_REPLACEMENTS
  if (_n<_L2_NMAX && _n<_L2_NMAX) {
    vc_trsm(_side,_uplo,_transA,_diag,_m,_n,_alpha,_A,_lda,_B,_ldb);
    return ;
  }
# endif
  char side=char(_side),uplo=char(_uplo),transA=char(_transA),diag=char(_diag);
  INTEGER m=_m,n=_n,lda=_lda,ldb=_ldb;
  /**/STRSM(&side,&uplo,&transA,&diag,&m,&n,&_alpha,_A,&lda,_B,&ldb);
}

//=============================================================================



//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

#endif // __VC_MATH_BLAS_HH
