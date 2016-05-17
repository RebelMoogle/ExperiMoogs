//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================


#ifndef __VC_MATH_VCBLAS_HH
#define __VC_MATH_VCBLAS_HH

#include "blas_flags.hh"

#include "../system/platform.hh"

/** \file vcblas.hh Math: replacements for BLAS routines (C++ interface)
 */

//=============================================================================

namespace VC {
namespace math {
namespace blas {

/** \defgroup vc_vcblas BLAS reimplementation
    \ingroup vc_blas_details
    Replacement for BLAS: BLAS calls can be replaced by their counterparts
    (e.g., \c DGEMM replaced by VC::math::blas::vc_gemm()). The replacements
    are using the C++-interface of the blas wrappers [\ref vc_blas].

    The rationale is to use the replacements for small matrices and/or
    fixed dimensions to give the optimizer of the compiler a chance to
    produce inlined code. (Note: for fixed SideFlag, UpperLowerFlag,
    TransposeFlag, DiagonalFlag, most code should be subject to dead
    code elimination.)

    The code is taken from the BLAS reference implementation
    http://www.netlib.org/blas/.

    \sa VC::math::blas, VC_NO_BLAS_REPLACEMENTS
 */


  //
  // level 1
  //

/// inline replacement for VC::math::blas::swap() \ingroup vc_vcblas
template <typename T>
inline void vc_swap(const int _n,T*__RESTRICT _x,const int _incx,T*__RESTRICT _y,const int _incy) {
  _VC_DBG_BLAS_SCOPE()
  for (int i=0;i<_n;++i) {
    T tmp=_x[i*_incx]; _x[i*_incx]=_y[i*_incy]; _y[i*_incx]=tmp;
  }
}

/// inline replacement for VC::math::blas::scal() \ingroup vc_vcblas
template <typename T>
inline void vc_scal(int _n,T _alpha, T*__RESTRICT _x,
int _incx) {
  _VC_DBG_BLAS_SCOPE()
  for (int i=0;i<_n;++i) _x[i*_incx]*=_alpha;
}

/// inline replacement for VC::math::blas::copy() \ingroup vc_vcblas
template <typename T>
inline void vc_copy(int _n, const T*__RESTRICT _x,int _incx, T*__RESTRICT _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
  for (int i=0;i<_n;++i) _y[i*_incy]=_x[i*_incx];
}

/// inline replacement for VC::math::blas::axpy() \ingroup vc_vcblas
template <typename T>
inline void vc_axpy(int _n,T _alpha,const T*__RESTRICT _x,int _incx,T*__RESTRICT _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()
  for (int i=0;i<_n;++i) _y[i*_incy]+=_x[i*_incx]*_alpha;
}

/// inline replacement for VC::math::blas::dot() \ingroup vc_vcblas
template <typename T>
inline T vc_dot(int _n,const T*__RESTRICT _x,int _incx, const T*__RESTRICT _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
  T d=_x[0]*_y[0];
  for (int i=1;i<_n;++i) d+=_x[i*_incx]*_y[i*_incy];
  return d;
}

/// inline replacement for VC::math::blas::ddot() \ingroup vc_vcblas
template <typename T>
inline T vc_ddot(int _n,const T*__RESTRICT _x,int _incx, const T*__RESTRICT _y,int _incy) {
  _VC_DBG_BLAS_SCOPE()
  double d=_x[0]*_y[0];
  for (int i=1;i<_n;++i) d+=_x[i*_incx]*_y[i*_incy];
  return T(d);
}

/// inline replacement for VC::math::blas::nrm2() \ingroup vc_vcblas
template <typename T>
inline T vc_nrm2(int _n, const T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  T d=_x[0]*_x[0];
  for (int i=1;i<_n;++i) d+=_x[i*_incx]*_x[i*_incx];
  return sqrt(d);
}

/// inline replacement for VC::math::blas::asum() \ingroup vc_vcblas
template <typename T>
inline T vc_asum(int _n, const T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()
  T d=fabs(_x[0]);
  for (int i=1;i<_n;++i) d+=fabs(_x[i*_incx]);
  return d;
}

  //
  // level 2
  //

/// inline replacement for VC::math::blas::gemv() \ingroup vc_vcblas
template <typename T>
inline void
vc_gemv(TransposeFlag _trans,int _m,int _n,
        T _alpha,
        const T*__RESTRICT _A,int _lda,const T*__RESTRICT _x,int _incx,
        T _beta,
        T*__RESTRICT _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()

  if (_trans==NoT) { // y = alpha*A*x + beta*y
    for (int i=0;i<_m;++i)
      _y[i*_incy]*=_beta;

    for (int j=0;j<_n;++j) {
      T a=_alpha*_x[j*_incx];
      for (int i=0;i<_m;++i)
        _y[i*_incy]+=a*_A[j*_lda+i]; // as weighted sum of columns
    }
  }
  else { // y = alpha*A'*x + beta*y
    for (int j=0;j<_n;++j)
      _y[j*_incy]*=_beta;

    T accum;
    for (int j=0;j<_n;++j) {
      accum=_A[j*_lda]*_x[0];
      for (int i=1;i<_m;++i)
        accum+=_A[j*_lda+i]*_x[i*_incx];
      _y[j*_incy]+=_alpha*accum;
    }
  }
}

  //
  // not implemented: vc_gbmv()
  //

/// inline replacement for VC::math::blas::symv() \ingroup vc_vcblas
template <typename T>
inline void
vc_symv(UpperLowerFlag _uplo,int _n,
     T _alpha,
     const T*__RESTRICT _A,int _lda,
     const T*__RESTRICT _x,int _incx,
     T _beta,
     T*__RESTRICT _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()

  for (int i=0;i<_n;++i)
    _y[i*_incy]*=_beta;

  if (_uplo==Upper) {
    for (int j=0;j<_n;++j) {
      T a=_alpha*_x[j*_incx];
      T accum=0;
      for (int i=0;i<j;++i) {
        _y[i*_incy]+=a*_A[j*_lda+i]; // as weighted sum of columns
        accum+=_A[j*_lda+i]*_x[i*_incx];
      }
      _y[j*_incy]+=a*_A[j*_lda+j]+_alpha*accum;
    }
  }
  else {
    for (int j=0;j<_n;++j) {
      T a=_alpha*_x[j*_incx];
      T accum=0;
      _y[j*_incy]+=a*_A[j*_lda+j];
      for (int i=j+1;i<_n;++i) {
        _y[i*_incy]+=a*_A[j*_lda+i]; // as weighted sum of columns
        accum+=_A[j*_lda+i]*_x[i*_incx];
      }
      _y[j*_incy]+=_alpha*accum;
    }
  }
}

  //
  // not implemented: vc_sbmv()
  //

/// inline replacement for VC::math::blas::spmv() \ingroup vc_vcblas
template <typename T>
inline void
vc_spmv(UpperLowerFlag _uplo,int _n,
        T _alpha,const T*__RESTRICT _Ap,
        const T*__RESTRICT _x,int _incx,
        T _beta,
        T*__RESTRICT _y, int _incy) {
  _VC_DBG_BLAS_SCOPE()

  for (int i=0;i<_n;++i)
    _y[i*_incy]*=_beta;

  if (_uplo==Upper) {
    int kk=0;
    for (int j=0;j<_n;++j,kk+=j) {
      T a=_alpha*_x[j*_incx];
      T accum=0;
      int k=kk;
      for (int i=0;i<j;++i,++k) {
        _y[i*_incy]+=a*_Ap[k];
        accum+=_Ap[k]*_x[i*_incx];
      }
      _y[j*_incy]+=a*_Ap[kk+j]+_alpha*accum;
    }
  }
  else {
    int kk=0;
    for (int j=0;j<_n;++j,kk+=(_n-j+1)) {
      T a=_alpha*_x[j*_incx];
      T accum=0;
      _y[j*_incy]+=a*_Ap[kk];
      int k=kk+1;
      for (int i=j+1;i<_n;++i,++k) {
        _y[i*_incy]+=a*_Ap[k];
        accum+=_Ap[k]*_x[i*_incx];
      }
      _y[j*_incy]+=_alpha*accum;
    }
  }
}

/// inline replacement for VC::math::blas::trmv() \ingroup vc_vcblas
template <typename T>
inline void
vc_trmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
        const T*__RESTRICT _A,int _lda,
        T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()

    if (_trans==NoT) {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j) {
          T a=_x[j*_incx];
          for (int i=0;i<j;++i)
            _x[i*_incx]+=a*_A[j*_lda+i];
          if (_diag==NoUnitDiag)
            _x[j*_incx]*=_A[j*_lda+j];
        }
      }
      else {
        for (int j=_n-1;j>=0;--j) {
          T a=_x[j*_incx];
          for (int i=_n-1;i>j;--i)
            _x[i*_incx]+=a*_A[j*_lda+i];
          if (_diag==NoUnitDiag)
            _x[j*_incx]*=_A[j*_lda+j];
        }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int j=_n-1;j>=0;--j) {
          T a=_x[j*_incx];
          if (_diag==NoUnitDiag)
            a*=_A[j*_lda+j];
          for (int i=j-1;i>=0;--i)
            a+=_A[j*_lda+i]*_x[i*_incx];
          _x[j*_incx]=a;
        }
      }
      else {
        for (int j=0;j<_n;++j) {
          T a=_x[j*_incx];
          if (_diag==NoUnitDiag)
            a*=_A[j*_lda+j];
          for (int i=j+1;i<_n;++i)
          a+=_A[j*_lda+i]*_x[i*_incx];
          _x[j*_incx]=a;
        }
      }
    }
}

  //
  // not implemented: vc_tbmv()
  //

/// inline replacement for VC::math::blas::tpmv() \ingroup vc_vcblas
template <typename T>
inline void
vc_tpmv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
        const T*__RESTRICT _Ap,
        T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()

  if (_trans==NoT) {
    if (_uplo==Upper) {
      int kk=0;
      for (int j=0;j<_n;++j,kk+=j) {
        /*if (_x[j*_incx]!=T(0))*/ {
          T a=_x[j*_incx];
          int k=kk;
          for (int i=0;i<j;++i,++k)
            _x[i*_incx]+=a*_Ap[k];
          if (_diag==NoUnitDiag)
            _x[j*_incx]*=_Ap[kk+j];
        }
      }
    }
    else {
      int kk=(_n*(_n+1))/2-1;
      for (int j=_n-1;j>=0;--j) {
        T a=_x[j*_incx];
        int k=kk;
        for (int i=_n-1;i>j;--i,--k)
          _x[i*_incx]+=a*_Ap[k];
        if (_diag==NoUnitDiag)
          _x[j*_incx]*=_Ap[kk-_n+j+1];
        kk-=(_n-j);
      }
    }
  }
  else {
    if (_uplo==Upper) {
      int kk=(_n*(_n+1))/2-1;
      for (int j=_n-1;j>=0;--j) {
        T a=_x[j*_incx];
        if (_diag==NoUnitDiag)
          a*=_Ap[kk];
        int k=kk-1;
        for (int i=j-1;i>=0;--i,--k)
          a+=_Ap[k]*_x[i*_incx];
        _x[j*_incx]=a;
        kk-=(j+1);
      }
    }
    else {
      int kk=0;
      for (int j=0;j<_n;++j) {
        T a=_x[j*_incx];
        if (_diag==NoUnitDiag)
          a*=_Ap[kk];
        int k=kk+1;
        for (int i=j+1;i<_n;++i,++k)
          a+=_Ap[k]*_x[i*_incx];
        _x[j*_incx]=a;
        kk+=(_n-j);
      }
    }
  }
}

/// inline replacement for VC::math::blas::trsv() \ingroup vc_vcblas
template <typename T>
inline void
vc_trsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const T*__RESTRICT _A,int _lda,
     T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()

    if (_trans==NoT) {
      if (_uplo==Upper) {
        for (int j=_n-1;j>=0;--j) {
          /*if (_x[j*_incx]!=T(0))*/ {
            if (_diag==NoUnitDiag)
              _x[j*_incx]/=_A[j*_lda+j];
            T a=_x[j*_incx];
            for (int i=j-1;i>=0;--i)
              _x[i*_incx]-=a*_A[j*_lda+i];
          }
        }
      }
      else {
        for (int j=0;j<_n;++j) {
          /*if (_x[j*_incx]!=T(0))*/ {
            if (_diag==NoUnitDiag)
              _x[j*_incx]/=_A[j*_lda+j];
            T a=_x[j*_incx];
            for (int i=j+1;i<_n;++i)
              _x[i*_incx]-=a*_A[j*_lda+i];
          }
        }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j) {
          T a=_x[j*_incx];
          for (int i=0;i<j;++i)
            a-=_A[j*_lda+i]*_x[i*_incx];
          if (_diag==NoUnitDiag)
            a/=_A[j*+_lda+j];
          _x[j*_incx]=a;
        }
      }
      else {
        for (int j=_n-1;j>=0;--j) {
          T a=_x[j*_incx];
          for (int i=_n-1;i>j;--i)
            a-=_A[j*_lda+i]*_x[i*_incx];
          if (_diag==NoUnitDiag)
            a/=_A[j*_lda+j];
          _x[j*_incx]=a;
        }
      }
    }
}

  //
  // not implemented: vc_tbsv()
  //

/// inline replacement for VC::math::blas::tpsv() \ingroup vc_vcblas
template <typename T>
inline void
vc_tpsv(UpperLowerFlag _uplo,TransposeFlag _trans,DiagonalFlag _diag,int _n,
     const T*__RESTRICT _Ap,
     T*__RESTRICT _x,int _incx) {
  _VC_DBG_BLAS_SCOPE()

    if (_trans==NoT) {
      if (_uplo==Upper) {
        int kk=(_n*(_n+1))/2-1;
        for (int j=_n-1;j>=0;--j) {
          /*if (_x[j*_incx]!=T(0))*/ {
            if (_diag==NoUnitDiag)
              _x[j*_incx]/=_Ap[kk];
            T a=_x[j*_incx];
            int k=kk-1;
            for (int i=j-1;i>=0;--i,--k)
              _x[i*_incx]-=a*_Ap[k];
          }
          kk-=(j+1);
        }
      }
      else {
        int kk=0;
        for (int j=0;j<_n;++j)  {
          /*if (_x[j*_incx]!=T(0))*/ {
            if (_diag==NoUnitDiag)
              _x[j*_incx]/=_Ap[kk];
            T a=_x[j*_incx];
            int k=kk+1;
            for (int i=j+1;i<_n;++i,++k)
              _x[i*_incx]-=a*_Ap[k];
          }
          kk+=(_n-j);
        }
      }
    }
    else {
      if (_uplo==Upper) {
        int kk=0;
        for (int j=0;j<_n;++j) {
          T a=_x[j*_incx];
          int k=kk;
          for (int i=0;i<j;++i,++k)
            a-=_Ap[k]*_x[i*_incx];
          if (_diag==NoUnitDiag)
            a/=_Ap[kk+j];
          _x[j*_incx]=a;
          kk+=(j+1);
        }
      }
      else {
        int kk=(_n*(_n+1))/2-1;
        for (int j=_n-1;j>=0;--j) {
          T a=_x[j*_incx];
          int k=kk;
          for (int i=_n-1;i>j;--i,--k)
            a-=_Ap[k]*_x[i*_incx];
          if (_diag==NoUnitDiag)
            a/=_Ap[kk-_n+j+1];
          _x[j*_incx]=a;
          kk-=(_n-j);
        }
      }
    }
}

/// inline replacement for VC::math::blas::syr() \ingroup vc_vcblas
template <typename T>
inline void
vc_syr(UpperLowerFlag _uplo,int _n,
    T _alpha,
    const T*__RESTRICT _x,int _incx,
    T*__RESTRICT _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()

    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {
        T a=_alpha*_x[j*_incx];
        for (int i=0;i<=j;++i)
          _A[j*_lda+i]+=_x[i*_incx]*a;
      }
    }
    else {
      for (int j=0;j<_n;++j)
        /*if (_x[j*_incx]!=T(0))*/ {
          T a=_alpha*_x[j*_incx];
          for (int i=j;i<_n;++i)
            _A[j*_lda+i]+=_x[i*_incx]*a;
        }
    }
}

/// inline replacement for VC::math::blas::spr() \ingroup vc_vcblas
template <typename T>
inline void
vc_spr(UpperLowerFlag _uplo,int _n,
    T _alpha,
    const T*__RESTRICT _x,int _incx,
    T*__RESTRICT _Ap) {
  _VC_DBG_BLAS_SCOPE()

  int kk=0;
  if (_uplo==Upper) {
    for (int j=0;j<_n;++j) {
      /*if (_x[j*_incx]!=T(0))*/ {
        T a=_alpha*_x[j*_incx];
        int k=kk;
        for (int i=0;i<=j;++i,++k)
          _Ap[k]+=_x[i*_incx]*a;
      }
      kk+=(j+1);
    }
  }
  else {
    for (int j=0;j<_n;++j) {
      /*if (_x[j*_incx]!=T(0))*/ {
        T a=_alpha*_x[j*_incx];
        int k=kk;
        for (int  i=j;i<_n;++i,++k)
          _Ap[k]+=_x[i*_incx]*a;
      }
      kk+=(_n-j);
    }
  }
}

/// inline replacement for VC::math::blas::syr2() \ingroup vc_vcblas
template <typename T>
inline void
vc_syr2(UpperLowerFlag _uplo,int _n,
     T _alpha,
     const T*__RESTRICT _x,int _incx,
     const T*__RESTRICT _y,int _incy,
     T*__RESTRICT _A,int _lda) {
  _VC_DBG_BLAS_SCOPE()

  if (_uplo==Upper) {
    for (int j=0;j<_n;++j) {
      /*if (_x[j*_incx]!=T(0) || _y[j*+_incx]!=T(0))*/ {
        T a=_alpha*_y[j*_incy];
        T b=_alpha*_x[j*_incx];
        for (int i=0;i<=j;++i)
          _A[j*_lda+i]+=_x[i*_incx]*a+_y[i*_incy]*b;
      }
    }
  }
  else {
     for (int j=0;j<_n;++j) {
       /*if (_x[j*_incx]!=T(0) || _y[j*+_incx]!=T(0))*/ {
         T a=_alpha*_y[j*_incy];
         T b=_alpha*_x[j*_incx];
         for (int i=j;i<_n;++i)
           _A[j*_lda+i]+=_x[i*_incx]*a+_y[i*_incy]*b;
       }
     }
  }
}

/// inline replacement for VC::math::blas::spr2() \ingroup vc_vcblas
template <typename T>
inline void
vc_spr2(UpperLowerFlag _uplo,int _n,
        T _alpha,
        const T*__RESTRICT _x,int _incx,
        const T*__RESTRICT _y,int _incy,
        T*__RESTRICT _Ap) {
  _VC_DBG_BLAS_SCOPE()

  int kk=0;

  if (_uplo==Upper) {
    for (int j=0;j<_n;++j) {
      /*if (_x[j*_incx]!=T(0) || _y[j*+_incx]!=T(0))*/ {
        T a=_alpha*_y[j*_incy];
        T b=_alpha*_x[j*_incx];
        int k=kk;
        for (int i=0;i<=j;++i,++k)
          _Ap[k]+=_x[i*_incx]*a+_y[i*_incy]*b;
      }
      kk+=(j+1);
    }
  }
  else {
    for (int j=0;j<_n;++j) {
      /*if (_x[j*_incx]!=T(0) || _y[j*+_incx]!=T(0))*/ {
        T a=_alpha*_y[j*_incy];
        T b=_alpha*_x[j*_incx];
        int k=kk;
        for (int i=j;i<_n;++i,++k)
          _Ap[k]+=_x[i*_incx]*a+_y[i*_incy]*b;
      }
      kk+=(_n-j);
    }
  }
}

  //
  // level 3
  //

/// inline replacement for VC::math::blas::gemm() \ingroup vc_vcblas
template <typename T>
inline void
vc_gemm(TransposeFlag _transA,TransposeFlag _transB,int _m,int _n,int _k,
     T _alpha,
     const T*__RESTRICT _A,int _lda,
     const T*__RESTRICT _B,int _ldb,
     T _beta,
     T*__RESTRICT _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()

  // trivial cases
  if (_alpha==T(0)) {
    if (_beta==T(0)) {
      if (_ldc==_m)
        for (int i=0;i<_m*_n;++i) _C[i]=T(0);
      else
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]=T(0);
    }
    else if (_beta!=T(1)) {
      if (_ldc==_m)
        for (int i=0;i<_m*_n;++i) _C[i]*=_beta;
      else
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]*=_beta;
    }
    return;
  }

  if (_transB==NoT) {
    if (_transA==NoT) {

      for (int j=0;j<_n;++j) {
        if (_beta==T(0)) {
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]*=_beta;
        }

        for (int h=0;h<_k;++h)
          /*if (_B[j*_ldb+h]!=T(0))*/ {
          T a=_alpha*_B[j*_ldb+h];
            for (int i=0;i<_m;++i)
              _C[j*_ldc+i]+=a*_A[h*_lda+i];
          }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        for (int i=0;i<_m;++i) {
          T a=T(0);
          for (int h=0;h<_k;++h)
            a+=_A[i*_lda+h]*_B[j*_ldb+h];
          _C[j*_ldc+i]=(_beta==T(0) ? _alpha*a : _alpha*a+_beta*_C[j*_ldc+i]);
        }
      }
    }
  }
  else {
    if (_transA==NoT) {
      for (int j=0;j<_n;++j) {
        if (_beta==T(0)) {
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]*=_beta;
        }
        for (int h=0;h<_k;++h) {
          /*if (_B[h*_ldb+j]!=T(0))*/ {
            T a=_alpha*_B[h*_ldb+j];
            for (int i=0;i<_m;++i)
              _C[j*_ldc+i]+=a*_A[h*_lda+i];
          }
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        for (int i=0;i<_m;++i) {
          T a=T(0);
          for (int h=0;h<_k;++h)
            a+=_A[i*_lda+h]*_B[h*_ldb+j];

          _C[j*_ldc+i]=(_beta==T(0) ? _alpha*a : _alpha*a+_beta*_C[j*_ldc+i]);
        }
      }
    }
  }
}

/// inline replacement for VC::math::blas::symm() \ingroup vc_vcblas
template <typename T>
inline void
vc_symm(SideFlag _side,UpperLowerFlag _uplo,int _m,int _n,
        T _alpha,
        const T*__RESTRICT _A,int _lda,
        const T*__RESTRICT _B,int _ldb,
        T _beta,
        T*__RESTRICT _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()

  // trivial cases
  if (_alpha==T(0)) {
    if (_beta==T(0)) {
      if (_ldc==_m)
        for (int i=0;i<_m*_n;++i) _C[i]=T(0);
      else
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]=T(0);
    }
    else if (_beta!=T(1)) {
      if (_ldc==_m)
        for (int i=0;i<_m*_n;++i) _C[i]*=_beta;
      else
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i)
            _C[j*_ldc+i]*=_beta;
    }
    return;
  }

  if (_side==Left) {
    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {
        for (int i=0;i<_m;++i) {
          T a=_alpha*_B[j*_ldb+i];
          T b=T(0);
          for (int k=0;k<i;++k) {
            _C[j*_ldc+k]+=a*_A[i*_lda+k];
            b+=_B[j*_ldb+k]*_A[i*_lda+k];
          }
          if (_beta==T(0))
            _C[j*_ldc+i]=a*_A[i*_lda+i]+_alpha*b;
          else
            _C[j*_ldc+i]=_beta*_C[j*_ldc+i]+a*_A[i*_lda+i]+_alpha*b;
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        for (int i=_m-1;i>=0;--i) {
          T a=_alpha*_B[j*_ldb+i];
          T b=T(0);
          for (int k=i+1;k<_m;++k) {
            _C[j*_ldc+k]+=a*_A[i*_lda+k];
            b+=_B[j*_ldb+k]*_A[i*_lda+k];
          }
          if (_beta==T(0))
            _C[j*_ldc+i]=a*_A[i*_lda+i]+_alpha*b;
          else
            _C[j*_ldc+i]=_beta*_C[j*_ldc+i]+a*_A[i*_lda+i]+_alpha*b;
        }
      }
    }
  }
  else {
    for (int j=0;j<_n;++j) {
      T a=_alpha*_A[j*_lda+j];
      if (_beta==T(0)) {
        for (int i=0;i<_m;++i)
          _C[j*_ldc+i]=a*_B[j*_ldb+i];
      }
      else {
        for (int i=0;i<_m;++i)
          _C[j*_ldc+i]=_beta*_C[j*_ldc+i]+a*_B[j*_ldb+i];
      }
      for (int k=0;k<j;++k) {
        T a=_alpha*(_uplo==Upper ? _A[j*_lda+k] : _A[k*_lda+j]);

        for (int i=0;i<_m;++i)
          _C[j*_ldc+i]+=a*_B[k*_ldb+i];
      }
      for (int k=j+1;k<_n;++k) {
        T a=_alpha*(_uplo==Upper ? _A[k*_lda+j] : _A[j*_lda+k]);

        for (int i=0;i<_m;++i)
          _C[j*_ldc+i]+=a*_B[k*_ldb+i];
      }
    }
  }
}

/// inline replacement for VC::math::blas::syrk() \ingroup vc_vcblas
template <typename T>
inline void
vc_syrk(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
        T _alpha,
        const T*__RESTRICT _A,int _lda,
        T _beta,
        T*__RESTRICT _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()

  if (_alpha==T(0)) {
    if (_uplo==Upper) {
      if (_beta==T(0)) {
        for (int j=0;j<_n;++j)
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]=T(0);
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]*=_beta;
      }
    }
    else {
      if (_beta==T(0)) {
        for (int j=0;j<_n;++j)
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]=T(0);
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]*=_beta;
      }
    }
    return;
  }

  if (_trans==NoT) {

    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {

        if (_beta==T(0)) {
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]*=_beta;
        }

        for (int h=0;h<_k;++h) {
          /*if (_A[h*_lda+j]!=T(0))*/ {
            T a=_alpha*_A[h*_lda+j];
            for (int i=0;i<=j;++i)
              _C[j*_ldc+i]+=a*_A[h*_lda+i];
          }
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {

        if (_beta==T(0)) {
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]*=_beta;
        }

        for (int h=0;h<_k;++h) {
          /*if (_A[h*_lda+j]!=T(0))*/ {
            T a=_alpha*_A[h*_lda+j];
            for (int i=j;i<_n;++i)
              _C[j*_ldc+i]+=a*_A[h*_lda+i];
          }
        }
      }
    }
  }
  else {
    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {
        for (int i=0;i<=j;++i) {
          T a=T(0);
          for (int h=0;h<_k;++h)
            a+=_A[i*_lda+h]*_A[j*_lda+h];
          if (_beta==T(0))
            _C[j*_ldc+i]=_alpha*a;
          else
            _C[j*_ldc+i]=_alpha*a+_beta*_C[j*_ldc+i];
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        for (int i=j;i<_n;++i) {
          T a=T(0);
          for (int h=0;h<_k;++h)
            a+=_A[i*_lda+h]*_A[j*_lda+h];
          if (_beta==T(0))
            _C[j*_ldc+i]=_alpha*a;
          else
            _C[j*_ldc+i]=_alpha*a+_beta*_C[j*_ldc+i];
        }
      }
    }
  }
}

/// inline replacement for VC::math::blas::syr2k() \ingroup vc_vcblas
template <typename T>
inline void
vc_syr2k(UpperLowerFlag _uplo,TransposeFlag _trans,int _n,int _k,
         T _alpha,
         const T*__RESTRICT _A,int _lda,
         const T*__RESTRICT _B,int _ldb,
         T _beta,
         T*__RESTRICT _C,int _ldc) {
  _VC_DBG_BLAS_SCOPE()

  if (_alpha==T(0)) {
    if (_uplo==Upper) {
      if (_beta==T(0)) {
        for (int j=0;j<_n;++j)
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]=T(0);
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]*=_beta;
      }
    }
    else {
      if (_beta==T(0)) {
        for (int j=0;j<_n;++j)
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]=T(0);
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]*=_beta;
      }
    }
    return;
  }

  if (_trans==NoT) {
    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {
        if (_beta==T(0)) {
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=0;i<=j;++i)
            _C[j*_ldc+i]*=_beta;
        }

        for (int h=0;h<_k;++h) {
          /*if (_A[h*_lda+j]!=T(0) || _B[h*_ldb+j]!=T(0))*/ {
            T a=_alpha*_B[h*_ldb+j];
            T b=_alpha*_A[h*_lda+j];

            for (int i=0;i<=j;++i) {
              _C[j*_ldc+i]+=_A[h*_lda+i]*a+_B[h*_ldb+i]*b;
            }
          }
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        if (_beta==T(0)) {
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]=T(0);
        }
        else if (_beta!=T(1)) {
          for (int i=j;i<_n;++i)
            _C[j*_ldc+i]*=_beta;
        }

        for (int h=0;h<_k;++h) {
          /*if (_A[h*_lda+j]!=T(0) || _B[h*_ldb+j]!=T(0))*/ {
            T a=_alpha*_B[h*_ldb+j];
            T b=_alpha*_A[h*_lda+j];

            for (int i=j;i<_n;++i) {
              _C[j*_ldc+i]+=_A[h*_lda+i]*a+_B[h*_ldb+i]*b;
            }
          }
        }
      }
    }
  }
  else {
    if (_uplo==Upper) {
      for (int j=0;j<_n;++j) {
        for (int i=0;i<=j;++i) {
          T a=T(0);
          T b=T(0);
          for (int h=0;h<_k;++h) {
            a+=_A[i*_lda+h]*_B[j*_ldb+h];
            b+=_B[i*_lda+h]*_A[j*_ldb+h];
          }
          if (_beta==T(0))
            _C[j*_ldc+i]=_alpha*a+_alpha*b;
          else
            _C[j*_ldc+i]=_beta*_C[j*_ldc+i]+_alpha*a+_alpha*b;
        }
      }
    }
    else {
      for (int j=0;j<_n;++j) {
        for (int i=j;i<_n;++i) {
          T a=T(0);
          T b=T(0);
          for (int h=0;h<_k;++h) {
            a+=_A[i*_lda+h]*_B[j*_ldb+h];
            b+=_B[i*_lda+h]*_A[j*_ldb+h];
          }
          if (_beta==T(0))
            _C[j*_ldc+i]=_alpha*a+_alpha*b;
          else
            _C[j*_ldc+i]=_beta*_C[j*_ldc+i]+_alpha*a+_alpha*b;
        }
      }
    }
  }
}

/// inline replacement for VC::math::blas::trmm() \ingroup vc_vcblas
template <typename T>
inline void
vc_trmm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
     int _m,int _n,
     T _alpha,
     const T*__RESTRICT _A,int _lda,
     T*__RESTRICT _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()

  if (_alpha==T(0)) {
    if (_ldb==_m) {
      for (int i=0;i<_m*_n;++i)
        _B[i]=T(0);
    }
    else {
      for (int j=0;j<_n;++j)
        for (int i=0;i<_m;++i)
          _B[j*_ldb+i]=T(0);
    }
    return;
  }

  if (_side==Left) {
    if (_transA==NoT) {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j)
          for (int k=0;k<_m;++k)
            /*if (_B[j*_ldb+k]!=T(0))*/ {
              T a=_alpha*_B[j*_ldb+k];
              for (int i=0;i<k;++i)
                _B[j*_ldb+i]+=a*_A[k*_lda+i];
              if (_diag==NoUnitDiag)
                a*=_A[k*_lda+k];
              _B[j*_ldb+k]=a;
            }
      }
      else {
        for (int j=0;j<_n;++j)
          for (int k=_m-1;k>=0;--k)
            /*if (_B[j*_ldb+k]!=T(0))*/ {
              T a=_alpha*_B[j*_ldb+k];
              _B[j*_ldb+k]=a;
              if (_diag==NoUnitDiag)
                _B[j*_ldb+k]*=_A[k*_lda+k];
              for (int i=k+1;i<_m;++i)
                _B[j*_ldb+i]+=a*_A[k*_lda+i];
            }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j)
          for (int i=_m-1;i>=0;--i) {
            T a=_B[j*_ldb+i];
            if (_diag==NoUnitDiag)
              a*=_A[i*_lda+i];
            for (int k=0;k<i;++k)
              a+=_A[i*_lda+k]*_B[j*_ldb+k];
            _B[j*_ldb+i]=_alpha*a;
          }
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i) {
            T a=_B[j*_ldb+i];
            if (_diag==NoUnitDiag)
              a*=_A[i*_lda+i];
            for (int k=i+1;k<_m;++k)
              a+=_A[i*_lda+k]*_B[j*_ldb+k];
            _B[j*_ldb+i]=_alpha*a;
          }
      }
    }
  }
  else {
    if (_transA==NoT) {
      if (_uplo==Upper) {
        for (int j=_n-1;j>=0;--j) {
          T a=_alpha;
          if (_diag==NoUnitDiag)
            a*=_A[j*_lda+j];
          for (int i=0;i<_m;++i)
            _B[j*_ldb+i]=a*_B[j*_ldb+i];
          for (int k=0;k<j;++k)
            /*if (_A[j*_lda+k]!=T(0))*/ {
              a=_alpha*_A[j*_lda+k];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]+=a*_B[k*_ldb+i];
            }
        }
      }
      else {
        for (int j=0;j<_n;++j) {
          T a=_alpha;
          if (_diag==NoUnitDiag)
            a*=_A[j*_lda+j];
          for (int i=0;i<_m;++i)
            _B[j*_ldb+i]=a*_B[j*_ldb+i];
          for (int k=j+1;k<_n;++k)
            /*if (_A[j*_lda+k]!=T(0))*/ {
              a=_alpha*_A[j*_lda+k];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]+=a*_B[k*_ldb+i];
            }
        }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int k=0;k<_n;++k) {
          for (int j=0;j<k;++j) {
            /*if (_A[k*_lda+j]!=T(0))*/ {
              T a=_alpha*_A[k*_lda+j];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]+=a*_B[k*_ldb+i];
            }
          }
          T b=_alpha;
          if (_diag==NoUnitDiag)
            b*=_A[k*_lda+k];
          if (b!=T(1))
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=b;
        }
      }
      else {
        for (int k=_n-1;k>=0;--k) {
          for (int j=k+1;j<_n;++j) {
            /*if (_A[k*_lda+j]!=T(0))*/ {
              T a=_alpha*_A[k*_lda+j];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]+=a*_B[k*_ldb+i];
            }
          }
          T b=_alpha;
          if (_diag==NoUnitDiag)
            b*=_A[k*_lda+k];
          if (b!=T(1))
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=b;
        }
      }
    }
  }
}

/// inline replacement for VC::math::blas::trsm() \ingroup vc_vcblas
template <typename T>
inline void
vc_trsm(SideFlag _side,UpperLowerFlag _uplo,TransposeFlag _transA,DiagonalFlag _diag,
        int _m,int _n,
        T _alpha,
        const T*__RESTRICT _A,int _lda,
        T*__RESTRICT _B,int _ldb) {
  _VC_DBG_BLAS_SCOPE()

  if (_alpha==T(0)) {
    if (_ldb==_m) {
      for (int i=0;i<_m*_n;++i)
        _B[i]=T(0);
    }
    else {
      for (int j=0;j<_n;++j)
        for (int i=0;i<_m;++i)
          _B[j*_ldb+i]=T(0);
    }
    return;
  }

  if (_side==Left) {
    if (_transA==NoT) {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j) {
          if (_alpha!=T(1)) {
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=_alpha;
          }
          for (int k=_m-1;k>=0;--k)
            if (_B[j*_ldb+k]!=T(0)) {
              if (_diag==NoUnitDiag)
                _B[j*_ldb+k]/=_A[k*_lda+k];
              for (int i=0;i<k;++i)
                _B[j*_ldb+i]-=_B[j*_ldb+k]*_A[k*_lda+i];
            }
        }
      }
      else {
        for (int j=0;j<_n;++j) {
          if (_alpha!=T(1)) {
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=_alpha;
          }
          for (int k=0;k<_m;++k)
            /*if (_B[j*_ldb+k]!=T(0))*/ {
              if (_diag==NoUnitDiag)
                _B[j*_ldb+k]/=_A[k*_lda+k];
              for (int i=k+1;i<_m;++i)
                _B[j*_ldb+i]-=_B[j*_ldb+k]*_A[k*_lda+i];
            }
        }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j)
          for (int i=0;i<_m;++i) {
            T a=_alpha*_B[j*_ldb+i];
            for (int k=0;k<i;++k)
              a-=_A[i*_lda+k]*_B[j*_ldb+k];
            if (_diag==NoUnitDiag)
              a/=_A[i*_lda+i];
            _B[j*_ldb+i]=a;
          }
      }
      else {
        for (int j=0;j<_n;++j)
          for (int i=_m-1;i>=0;--i) {
            T a=_alpha*_B[j*_ldb+i];
            for (int k=i+1;k<_m;++k)
              a-=_A[i*_lda+k]*_B[j*_ldb+k];
            if (_diag==NoUnitDiag)
              a/=_A[i*_lda+i];
            _B[j*_ldb+i]=a;
          }
      }
    }
  }
  else {
    if (_transA==NoT) {
      if (_uplo==Upper) {
        for (int j=0;j<_n;++j) {
          if (_alpha!=T(1))
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=_alpha;
          for (int k=0;k<j;++k)
            /*if (_A[j*_lda+k]!=T(0))*/
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]-=_A[j*_lda+k]*_B[k*_ldb+i];
          if (_diag==NoUnitDiag) {
            T a=T(1)/_A[j*_lda+j];
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=a;
          }
        }
      }
      else {
        for (int j=_n-1;j>=0;--j) {
          if (_alpha!=T(1))
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=_alpha;
          for (int k=j+1;k<_n;++k)
            /*if (_A[j*_lda+k]!=T(0))*/ {
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]-=_A[j*_lda+k]*_B[k*_ldb+i];
          }
          if (_diag==NoUnitDiag) {
            T a=T(1)/_A[j*_lda+j];
            for (int i=0;i<_m;++i)
              _B[j*_ldb+i]*=a;
          }
        }
      }
    }
    else {
      if (_uplo==Upper) {
        for (int k=_n-1;k>=0;--k) {
          if (_diag==NoUnitDiag) {
            T a=T(1)/_A[k*_lda+k];
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=a;
          }
          for (int j=0;j<k;++j)
            /*if (_A[k*_lda+j]!=T(0))*/ {
              T a=_A[k*_lda+j];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]-=a*_B[k*_ldb+i];
            }
          if (_alpha!=T(1))
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=_alpha;
        }
      }
      else {
        for (int k=0;k<_n;++k) {
          if (_diag==NoUnitDiag) {
            T a=T(1)/_A[k*_lda+k];
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=a;
          }
          for (int j=k+1;j<_n;++j)
            /*if (_A[k*_lda+j]!=T(0))*/ {
              T a=_A[k*_lda+j];
              for (int i=0;i<_m;++i)
                _B[j*_ldb+i]-=a*_B[k*_ldb+i];
            }
          if (_alpha!=T(1))
            for (int i=0;i<_m;++i)
              _B[k*_ldb+i]*=_alpha;
        }
      }
    }
  }
}

//=============================================================================



//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

#endif // __VC_MATH_VCBLAS_HH