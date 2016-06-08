#ifndef __VC_MATH_LU_HH
#define __VC_MATH_LU_HH

#include <limits>
#include <algorithm>
#include <type_traits>

#include "lapack.hh"

/** \file lu.hh Math: LU decomposition
 */

namespace VC {
namespace math {

/** LU decomposition of an `MxM` matrix _a with leading dimension `_LDA`.
    \ingroup vc_math
    This is an `inline` implementation of lu().
*/
template <typename T,int M,int LDA=M>
inline bool _lu(T* _a,int* _ipiv) {
  // LU decomposition from Numerical Recipes
# define AELT(i,j) _a[i+j*LDA]

  T   big, tmp;
  int i,j,k,imax;
  T   vv[M];  // scaling

  for (i=0;i<M;++i) {
    big=fabs(AELT(i,0));
    for (j=1;j<M;++j)
      if ((tmp=fabs(AELT(i,j)))>big) big=tmp;
    if (big==T(0))
      return false;
    vv[i]=T(1)/big;
  }

  for (k=0;k<M;++k) {
    big=vv[k]*fabs(AELT(k,k));
    imax=k;
    for (i=k+1;i<M;++i) {
      if ((tmp=vv[i]*fabs(AELT(i,k)))>big) {
        big=tmp;
        imax=i;
      }
    }

    if (k!=imax) {
      for (j=0;j<M;++j)
        std::swap(AELT(imax,j),AELT(k,j));
      vv[imax]=vv[k];
    }
    _ipiv[k]=imax+1; // 1-based indices

    if (AELT(k,k)==T(0)) AELT(k,k)=std::numeric_limits<T>::min();

    for (i=k+1;i<M;++i) {
      tmp=(AELT(i,k)/=AELT(k,k));
      for (j=k+1;j<M;++j)
        AELT(i,j)-=tmp*AELT(k,j);
    }
  }
  // missing: track swaps for sign of determinant

  return true;
#undef AELT
}

//-----------------------------------------------------------------------------

/** Back-substitution from lu() decomposition.
    \ingroup vc_math
    Simplified Inline version of LAPACK `DGETRS` (lapack::getrs()).
    This is an `inline` implementation of lu_subs().
*/
template <typename T,int M,int LDA=M>
inline void _lu_subs(const T* _a,const int* _ipiv,T* _b) {
  int ii=0,i,j,ip;
  T sum;

  for (int i=0;i<M;++i) { // forward
    ip=_ipiv[i]-1; // 1-based indices
    sum=_b[ip];
    _b[ip]=_b[i];
    if (ii!=0)
      for (int j=ii-1;j<i;++j) sum-=_a[i+j*LDA]*_b[j];
    else if (sum!=T(0))
      ii=i+1;
    _b[i]=sum;
  }
  for (i=M-1;i>=0;--i) { // backward
    sum=_b[i];
    for (j=i+1;j<M;++j) sum-=_a[i+j*LDA]*_b[j];
    _b[i]=sum/_a[i+i*LDA];
  }
}

/** Back-substitution from lu() decomposition.
    \ingroup vc_math
    Simplified Inline version of LAPACK `DGETRS` (lapack::getrs()).
    This is an `inline` implementation of lu_subs().
*/
template <typename T,int M,int LDA=M>
inline void _lu_subs(const T* _a,const int* _ipiv,T* _b,int _n) {
  for (int j=0;j<_n;++j)
    _lu_subs<T,M,LDA>(_a,_ipiv,_b+j*M);
}

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

template <typename T,int M,int LDA=M>
struct lu_t {
  lu_t() {}
  bool operator()(T* _a,int* _ipiv) const {
    int  info;
    if (std::is_same<int,lapack::INTEGER>::value) {
      info=lapack::getrf(M,M,_a,M,_ipiv);
    }
    else {
      lapack::INTEGER pi[M];
      info=lapack::getrf(M,M,_a,M,pi);
      for (int i=0;i<M;++i)
        _ipiv[i]=int(pi[i]);
    }
    assert(info>=0);

    return info==0;
  }
};

#  define SPECIALIZE_LU(M)\
  template <typename T>\
  struct lu_t<T,M,M> {\
  lu_t() {}\
  bool operator()(T* _a,int* _ipiv) {\
    return _lu<T,M,M>(_a,_ipiv);\
  }};

SPECIALIZE_LU(2)
SPECIALIZE_LU(3)
SPECIALIZE_LU(4)
SPECIALIZE_LU(5)
SPECIALIZE_LU(6)
SPECIALIZE_LU(7)
SPECIALIZE_LU(8)
SPECIALIZE_LU(9)
SPECIALIZE_LU(10)
SPECIALIZE_LU(11)
SPECIALIZE_LU(12)

#   ifndef NDEBUG

SPECIALIZE_LU(13)
SPECIALIZE_LU(14)
SPECIALIZE_LU(15)
SPECIALIZE_LU(16)
SPECIALIZE_LU(17)
SPECIALIZE_LU(18)
SPECIALIZE_LU(19)
SPECIALIZE_LU(20)
SPECIALIZE_LU(21)
SPECIALIZE_LU(22)
SPECIALIZE_LU(23)
SPECIALIZE_LU(24)
SPECIALIZE_LU(25)
SPECIALIZE_LU(26)
SPECIALIZE_LU(27)
SPECIALIZE_LU(28)
SPECIALIZE_LU(29)
SPECIALIZE_LU(30)

// TODO: require benchmark!!!
#   endif

#  undef SPECIALIZE_LU

# endif // DOXYGEN_SKIP

/** LU decomposition of an `MxM` matrix _a with leading dimension `LDA`.
    \ingroup vc_math

    \tparam T value type
    \tparam M number of rows
    \tparam LDA leading dimension of matrix

    Factorize square matrix `A=L*U`.

    \param[in,out] _a matrix on **input**, LU decomposition on **output**.
    \param _ipiv stores the pivot indices as **1-based** indices (same as LAPACK)
    (`_ipiv[i]==i+1` if column `i` has not been swapped).

    Matrix storage and meaning of the parameters and result is as for
    `DGETRF` (lapack::getrf()) in LAPACK with the following
    differences:

    \arg Only square matrices are supported, however,
    leading dimension applies.
    \arg In case of a singular matrix the factorization/result
    is undefined.

    **The standard implementation calls lapack::getrf(). For small dimensions
    `M` there exist specializations that call _lu().**

    \return `false` if the matrix `_a` is singular
*/
template <typename T,int M,int LDA=M>
bool lu(T* _a,int* _ipiv) {
  lu_t<T,M,LDA> luinst; return luinst(_a,_ipiv);
}

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

template <typename T,int M,int LDA=M>
struct lu_subs_t {
  lu_subs_t(const T* _a,const int* _ipiv,T* _b,int _n=1) {
    int info;
    if (std::is_same<int,lapack::INTEGER>::value)
      info=lapack::getrs(blas::NoT,M,_n,_a,M,_ipiv,_b,M);
    else {
      lapack::INTEGER pi[M];
      for (int i=0;i<M;++i)
        pi[i]=lapack::INTEGER(_ipiv[i]);
      info=lapack::getrs(blas::NoT,M,_n,_a,M,pi,_b,M);
    }
    assert(info==0 && "invalid argument");
    use_nowarn(info);
  }
};
#  define SPECIALIZE_LUSUBS(M)\
  template <typename T>\
  struct lu_subs_t<T,M,M> {\
  lu_subs_t(const T* _a,const int* _ipiv,T* _b) { _lu_subs(_a,_ipiv,_b); }\
    lu_subs_t(const T* _a,const int* _ipiv,T* _b,int _n) {\
      _lu_subs<T,M,M>(_a,_ipiv,_b,_n);                    \
  }};

SPECIALIZE_LUSUBS(2)
SPECIALIZE_LUSUBS(3)
SPECIALIZE_LUSUBS(4)
SPECIALIZE_LUSUBS(5)
SPECIALIZE_LUSUBS(6)

#   ifndef NDEBUG
SPECIALIZE_LUSUBS(7)
SPECIALIZE_LUSUBS(8)
SPECIALIZE_LUSUBS(9)
SPECIALIZE_LUSUBS(10)
SPECIALIZE_LUSUBS(11)
SPECIALIZE_LUSUBS(12)
// TODO: require benchmark!!!
#   endif

#  undef SPECIALIZE_LU

# endif // DOXYGEN_SKIP

/** Back-substitution from lu() decomposition.
    \ingroup vc_math

    \tparam T value type
    \tparam M number of rows
    \tparam LDA leading dimension of matrix

    Solve `A*x=b` from factorization `A=L*U`.

    \param _a LU decomposition as on exit from lu()
    \param _ipiv pivot indices as on exit from lu()
    \param _b[in,out] on **input** right-hand-side `M`-vector, **on output*
    solution `x`

    **The standard implementation calls lapack::getrs(). For small dimensions
    `M` there exist specializations that call _lu().**
*/
template <typename T,int M,int LDA=M>
void lu_subs(const T* _a,const int* _ipiv,T* _b) {
  lu_subs_t<T,M,LDA>(_a,_ipiv,_b);
}

/** Back-substitution from lu() decomposition.
    \ingroup vc_math

    \tparam T value type
    \tparam M number of rows
    \tparam LDA leading dimension of matrix

    Solve `A*x=b` from factorization `A=L*U`.

    \param _a LU decomposition as on exit from lu()
    \param _ipiv pivot indices as on exit from lu()
    \param _b[in,out] on **input** right-hand-side `M x _n` matrix, **on output*
    solution `x`
    \param _n number of columns on right-hand-side

    **The standard implementation calls lapack::getrs(). For small dimensions
    `M` there exist specializations that call _lu().**
*/
template <typename T,int M,int LDA=M>
void lu_subs(const T* _a,const int* _ipiv,T* _b,int _n) {
  lu_subs_t<T,M,LDA>(_a,_ipiv,_b,_n);
}

//-----------------------------------------------------------------------------
} // namespace VC
} // namespace math
//-----------------------------------------------------------------------------
#endif // __VC_MATH_LU_HH
