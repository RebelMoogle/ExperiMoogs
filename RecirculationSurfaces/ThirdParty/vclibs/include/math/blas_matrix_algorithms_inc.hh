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

    \arg Matrix algorithms (for matrix_reference_t), mostly calling LAPACK.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

#include "../system/platform.hh" // alloca
#include "lapack.hh" // really? INCLUDE ALL THIS HERE?

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/// sum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T sum(const vector_const_reference_t<T,vec<N,INC> >& _x) {
  T s=_x(0);
  for (int j=1;j<_x.n();++j)
    s+=_x(j);
  return s;
}
/// sum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T sum(const vector_reference_t<T,vec<N,INC> >& _x) {
  return sum(_x.const_ref());
}

/// mean of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T mean(const vector_const_reference_t<T,vec<N,INC> >& _x) {
  return sum(_x)/T(N);
}
/// mean of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T mean(const vector_reference_t<T,vec<N,INC> >& _x) {
  return mean(_x.const_ref());
}

/// minimum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T min(const vector_const_reference_t<T,vec<N,INC> >& _x) {
  T s=_x(0);
  for (int j=1;j<_x.n();++j)
    s=std::min(_x(j),s);
  return s;
}
/// minimum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T min(const vector_reference_t<T,vec<N,INC> >& _x) {
  return min(_x.const_ref());
}

/// maximum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T max(const vector_const_reference_t<T,vec<N,INC> >& _x) {
  T s=_x(0);
  for (int j=1;j<_x.n();++j)
    s=std::max(_x(j),s);
  return s;
}
/// maximum of entries \ingroup vc_blas_algorithms
template <typename T,int N,int INC>
T max(const vector_reference_t<T,vec<N,INC> >& _x) {
  return max(_x.const_ref());
}

//-----------------------------------------------------------------------------

/// load column sums to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void sum(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.n());
  for (int j=0;j<_a.n();++j)
    _x(j)=sum(column(_a,j));
}
/// load column sums to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void sum(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.m());
  for (int i=0;i<_a.m();++i)
    _x(i)=sum(row(_a,i));
}
/// load column sums to `x` \ingroup vc_blas_algorithms
template <typename T,TransposeFlag TRANS,int M,int N,int K,int LD,int INC>
void sum(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  return sum(_a.const_ref(),_x);
}

/// load column means to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void mean(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.n());
  for (int j=0;j<_a.n();++j)
    _x(j)=mean(column(_a,j));
}
/// load column means to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void mean(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.m());
  for (int i=0;i<_a.m();++i)
    _x(i)=mean(row(_a,i));
}
/// load column means to `x` \ingroup vc_blas_algorithms
template <typename T,TransposeFlag TRANS,int M,int N,int LD,int K,int INC>
void mean(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  return mean(_a.const_ref(),_x);
}


/// load column minima to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void min(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.n());
  for (int j=0;j<_a.n();++j)
    _x(j)=min(column(_a,j));
}
/// load column minima to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void min(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.m());
  for (int i=0;i<_a.m();++i)
    _x(i)=min(row(_a,i));
}
/// load column minima to `x` \ingroup vc_blas_algorithms
template <typename T,TransposeFlag TRANS,int M,int N,int LD,int K,int INC>
void min(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  return min(_a.const_ref(),_x);
}


/// load column maxima to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void max(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.n());
  for (int j=0;j<_a.n();++j)
    _x(j)=max(column(_a,j));
}
/// load column maxima to `x` \ingroup vc_blas_algorithms
template <typename T,int M,int N,int LD,int K,int INC>
void max(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  assert(_x.n()==_a.m());
  for (int i=0;i<_a.m();++i)
    _x(i)=max(row(_a,i));
}
/// load column maxima to `x` \ingroup vc_blas_algorithms
template <typename T,TransposeFlag TRANS,int M,int N,int LD,int K,int INC>
void max(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
         const vector_reference_t<T,vec<K,INC> >& _x) {
  return max(_a.const_ref(),_x);
}

//-----------------------------------------------------------------------------

/** LU factorization.
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, contents are overwritten by LU on output
    (ususally, _a should be a quadratic matrix)
    \param[out] _ipiv store permutation from pivoting (array of length
    <tt>min(_a.m(),_a.n()</tt>)
    \return `true` on success, `false` if `_a` was singular
    \sa VC::math::lapack::getrf()
 */
template <typename T,int M,int N,int LD>
bool lu(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
        lapack::INTEGER* _ipiv) {
  int info=lapack::getrf(_a.m(),_a.n(),_a.data(),_a.ld(),_ipiv);
  assert(info>=0 && "invalid argument");
  return info==0;
}

/** LU factorization.
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, contents are overwritten by LU on output
    (ususally, _a should be a quadratic matrix)
    \param[out] _ipiv store permutation from pivoting (array of length
    <tt>min(_a.m(),_a.n()</tt>)
    \return `true` on success, `false` if `_a` was singular
    \sa VC::math::lapack::gbtrf(), lu
 */
template <typename T,int M,int N,int KL,int KU,int LD>
bool lu(const matrix_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >& _a,
        lapack::INTEGER* _ipiv) {
  assert(_a.ld()>=2*_a.kl()+_a.ku()+1 && "must preallocate factors, see gbtrf()");
  int info=lapack::gbtrf(_a.m(),_a.n(),_a.kl(),_a.ku(),_a.data(),
                         _a.ld(),_ipiv);
  assert(info>=0 && "invalid argument");
  return info==0;
}

/** Solve Ax=B (or A'*X=B) by backsubstitution from LU (lu()).
    \ingroup vc_blas_algorithms
    \param[in] _a LU decomposed matrix as result of lu()
    \param[in] _ipiv permutation as result of lu()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::getrs(), lu()
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int MB,int NB,int LDB>
void lu_subs(const matrix_const_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
             lapack::INTEGER* _ipiv,
             const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  assert(_a.m()==_a.n());
  assert(_b.m()==_a.m());
  int info=lapack::getrs(TRANS,_a.n(),_b.n(),_a.data(),_a.ld(),_ipiv,
                 _b.data(),_b.ld());
  assert(info==0 && "invalid argument");
  use_nowarn(info);
}
/** Solve Ax=b (or A'*X=b) by backsubstitution from LU (lu()).
    \ingroup vc_blas_algorithms
    \param[in] _a LU decomposed matrix as result of lu()
    \param[in] _ipiv permutation as result of lu()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::getrs(), lu()
 */

template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int NB,int INCB>
void lu_subs(const matrix_const_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
             lapack::INTEGER* _ipiv,
             const vector_reference_t<T,vec<NB,INCB> >& _b) {
  lu_subs(_a,_ipiv,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B (or A'*X=B) by backsubstitution from LU (lu()).
    \ingroup vc_blas_algorithms
    \param[in] _a LU decomposed matrix as result of lu()
    \param[in] _ipiv permutation as result of lu()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::gbtrs(), lu()
 */

template <typename T,TransposeFlag TRANS,int MA,int NA,int KL,int KU,int LDA,int MB,int NB,int LDB>
void lu_subs(const matrix_const_reference_t<T,gb_mat<TRANS,MA,NA,KL,KU,LDA> >& _a,
             lapack::INTEGER* _ipiv,
             const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  assert(_a.m()==_a.n());
  assert(_b.m()==_a.m());
  int info=lapack::gbtrs(TRANS,_a.n(),_a.kl(),_a.ku(),_b.n(),_a.data(),_a.ld(),_ipiv,
                 _b.data(),_b.ld());
  assert(info==0 && "invalid argument");
  use_nowarn(info);
}
/** Solve Ax=b (or A'*X=b) by backsubstitution from LU (lu()).
    \ingroup vc_blas_algorithms
    \param[in] _a LU decomposed matrix as result of lu()
    \param[in] _ipiv permutation as result of lu()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::gbtrs(), lu()
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int KL,int KU,int LDA,int NB,int INCB>
void lu_subs(const matrix_const_reference_t<T,gb_mat<TRANS,MA,NA,KL,KU,LDA> >& _a,
             lapack::INTEGER* _ipiv,
             const vector_reference_t<T,vec<NB,INCB> >& _b) {
  lu_subs(_a,_ipiv,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B (or A'*X=B) by LU factorization (lu()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, overwritten by LU factorization on output (as by lu())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was singular
    \sa lu(), lu_subs
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int MB,int NB,int LDB>
bool lu_solve(const matrix_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
              const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  lapack::INTEGER* ipiv;
// # if VC_ALLOCA_LIMIT>0
//   if (_a.m()*sizeof(lapack::INTEGER)<=VC_ALLOCA_LIMIT)
//     ipiv=(lapack::INTEGER*) alloca(_a.m()*sizeof(lapack::INTEGER));
//   else
// #endif
    _ws.alloc(&ipiv,_a.m());

  if (!lu(_a,ipiv))
    return false;
  else
    lu_subs(_a.const_ref(),ipiv,_b);
  return true;
}
/** Solve Ax=b (or A'*X=b) by LU factorization (lu()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, overwritten by LU factorization on output (as by lu())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was singular
    \sa lu(), lu_subs()
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int NB,int INCB>
bool lu_solve(const matrix_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
              const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return lu_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B (or A'*X=B) by LU factorization (lu()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, overwritten by LU factorization on output (as by lu())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was singular
    \sa lu(), lu_subs()
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int KL,int KU,int MB,int NB,int LDB>
bool lu_solve(const matrix_reference_t<T,gb_mat<TRANS,MA,NA,KL,KU,LDA> >& _a,
              const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  lapack::INTEGER* ipiv;
// # if VC_ALLOCA_LIMIT>0
//   if (_a.m()*sizeof(lapack::INTEGER)<=VC_ALLOCA_LIMIT)
//     ipiv=(lapack::INTEGER*) alloca(_a.m()*sizeof(lapack::INTEGER));
//   else
// #endif
    _ws.alloc(&ipiv,_a.m());

  if (!lu(_a,ipiv))
    return false;
  else
    lu_subs(_a,ipiv,_b);
  return true;
}

/** Solve Ax=B (or A'*X=B) by LU factorization (lu()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, overwritten by LU factorization on output (as by lu())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was singular
    \sa lu(), lu_subs()
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA,int KL,int KU,int NB,int INCB>
bool lu_solve(const matrix_reference_t<T,gb_mat<TRANS,MA,NA,KL,KU,LDA> >& _a,
              const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return lu_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

  //
  // TODO lu_solve tridiagonal KL=KU=1 --> copy diagonals, gtrsv, (difference: _a unmodified)
  //

/** Compute inverse of _a doing LU factorization.
    \ingroup vc_blas_algorithms
    \param[in,out] _a inverse on output
    \sa lu(), lapack::getri(), temporary workspace is allocated from stack
 */
template <typename T,TransposeFlag TRANS,int MA,int NA,int LDA>
bool lu_invert(const matrix_reference_t<T,ge_mat<TRANS,MA,NA,LDA> >& _a,
               lapack::Workspace& _ws=lapack::Workspace::instance()) {
  static int lwork=0;
  if (lwork==0) {
    T work;
    int info=lapack::getri(_a.n(),(T*) 0,_a.ld(),(lapack::INTEGER*) 0,&work,-1);
    assert(info==0 && "workspace query failed");
    use_nowarn(info);
    lwork=int(work);
  }
  lapack::INTEGER* ipiv;
  T* work;
  _ws.alloc(&work,lwork, &ipiv,_a.m());

  if (!lu(_a,ipiv))
    return false;
  else {
    int info=lapack::getri(_a.n(),_a.data(),_a.ld(),ipiv,work,lwork);
    assert(info==0);
    use_nowarn(info);
  }
  return true;
}

  //
  // TODO: generalization / switch to tuned algorithms for small matrices
  //

//-----------------------------------------------------------------------------

/** Cholesky factorization A=LL' (UPLO==Lower) or A=R'R (UPLO==Upper).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, contents are overwritten by LL' on output
    \return `true` on success, `false` if `_a` was not positive definite
    \sa VC::math::lapack::potrf()
 */
template <typename T,UpperLowerFlag UPLO,int N,int LD>
bool chol(const matrix_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  int info=lapack::potrf(_a.uplo(),_a.n(),_a.data(),_a.ld());
  assert(info>=0 && "invalid argument");
  return info==0;
}
/** Cholesky factorization A=LL' (UPLO==Lower) or A=R'R (UPLO==Upper).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, contents are overwritten by LL' on output
    \return `true` on success, `false` if `_a` was not positive definite
    \sa VC::math::lapack::pbtrf()
 */
template <typename T,UpperLowerFlag UPLO,int N,int KD,int LD>
bool chol(const matrix_reference_t<T,sb_mat<UPLO,N,KD,LD> >& _a) {
  int info=lapack::pbtrf(_a.uplo(),_a.n(),_a.kd(),_a.data(),_a.ld());
  assert(info>=0 && "invalid argument");
  return info==0;
}
/** Cholesky factorization A=LL' (UPLO==Lower) or A=R'R (UPLO==Upper).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, contents are overwritten by LL' on output
    \return `true` on success, `false` if `_a` was not positive definiteg
    \sa VC::math::lapack::pbtrf()
 */
template <typename T,UpperLowerFlag UPLO,int N>
bool chol(const matrix_reference_t<T,sp_mat<UPLO,N> >& _a) {
  int info=lapack::pptrf(_a.uplo(),_a.n(),_a.data());
  assert(info>=0 && "invalid argument");
  return info==0;
}

/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::potrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int MB,int NB,int LDB>
void chol_subs(const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
               const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  assert(_b.m()==_a.n());
  int info=lapack::potrs(_a.uplo(),_a.n(),_b.n(),_a.data(),_a.ld(),_b.data(),_b.ld());
  assert(info==0 && "invalid argument");
  use_nowarn(info);
}
/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::potrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int MB,int NB,int INCB>
void chol_subs(const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
               const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_subs(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::pbtrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int K,int LDA,int MB,int NB,int LDB>
void chol_subs(const matrix_const_reference_t<T,sb_mat<UPLO,NA,K,LDA> >& _a,
               const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  assert(_b.m()==_a.n());
  int info=lapack::pbtrs(_a.uplo(),_a.n(),_a.k(),_b.n(),_a.data(),_a.ld(),_b.data(),_b.ld());
  assert(info==0 && "invalid argument");
  use_nowarn(info);
}
/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::pbtrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int K,int LDA,int MB,int NB,int INCB>
void chol_subs(const matrix_const_reference_t<T,sb_mat<UPLO,NA,K,LDA> >& _a,
               const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_subs(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::pptrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int MB,int NB,int LDB>
void chol_subs(const matrix_const_reference_t<T,sp_mat<UPLO,NA> >& _a,
               const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  assert(_b.m()==_a.n());
  int info=lapack::pptrs(_a.uplo(),_a.n(),_b.n(),_a.data(),_b.data(),_b.ld());
  assert(info==0 && "invalid argument");
  use_nowarn(info);
}
/** Solve Ax=B (or A'*X=B) by backsubstitution from Cholesky facorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in] _a decomposed matrix (LL' or R'R) as result of chol()
    \param[in,out] _b right-hand-side on input, solution on output
    \sa VC::math::lapack::pptrs(), chol()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int MB,int NB,int INCB>
void chol_subs(const matrix_const_reference_t<T,sp_mat<UPLO,NA> >& _a,
               const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_subs(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}


/** Solve Ax=B by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int MB,int NB,int LDB>
bool chol_solve(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
                const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  if (!chol(_a))
    return false;
  else
    chol_subs(_a.const_ref(),_b);
  return true;
}
/** Solve Ax=b by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return true on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NB,int INCB>
bool chol_solve(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
                const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int K,int LDA,int MB,int NB,int LDB>
bool chol_solve(const matrix_reference_t<T,sb_mat<UPLO,NA,K,LDA> >& _a,
                const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  if (!chol(_a))
    return false;
  else
    chol_subs(_a.const_ref(),_b);
  return true;
}
/** Solve Ax=b by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int K,int LDA,int NB,int INCB>
bool chol_solve(const matrix_reference_t<T,sb_mat<UPLO,NA,K,LDA> >& _a,
                const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

/** Solve Ax=B by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return true on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int MB,int NB,int LDB>
bool chol_solve(const matrix_reference_t<T,sp_mat<UPLO,NA> >& _a,
                const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b) {
  if (!chol(_a))
    return false;
  else
    chol_subs(_a.const_ref(),_b);
  return true;
}
/** Solve Ax=b by Cholesky factorization (chol()).
    \ingroup vc_blas_algorithms
    \param[in,out] _a input SPD matrix, overwritten by Cholesky factorization on output (as by chol())
    \param[in,out] _b right-hand-side on input, solution on output (as for lu_subs())
    \return `true` on success, `false` if `_a` was not positive definite
    \sa chol(), chol_subs()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int NB,int INCB>
bool chol_solve(const matrix_reference_t<T,sp_mat<UPLO,NA> >& _a,
                const vector_reference_t<T,vec<NB,INCB> >& _b) {
  return chol_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.inc()));
}

  //
  // TODO chol_solve tridiagonal K=1 --> copy diagonals, pttrs
  //

/** Compute inverse of SPD _a doing Cholesky factorization.
    \ingroup vc_blas_algorithms
    \param[in,out] _a inverse on output
    \sa chol(), lapack::potri()
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA>
bool chol_invert(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a) {
  if (!chol(_a))
    return false;
  else {
    int info=lapack::potri(_a.uplo(),_a.n(),_a.data(),_a.ld());
    assert(info==0);
    use_nowarn(info);
  }
  return true;
}

/** Compute inverse of SPD _a doing Cholesky factorization.
    \ingroup vc_blas_algorithms
    \param[in,out] _a inverse on output
    \sa chol(), lapack::pptri()
 */
template <typename T,UpperLowerFlag UPLO,int NA>
bool chol_invert(const matrix_reference_t<T,sp_mat<UPLO,NA> >& _a) {
  if (!chol(_a))
    return false;
  else {
    int info=lapack::pptri(_a.uplo(),_a.n(),_a.data());
    assert(info==0);
    use_nowarn(info);
  }
  return true;
}

  //
  // TODO: generalization / switch to tuned algorithms for small matrices
  //

//-----------------------------------------------------------------------------

/** Solve least-squares problem \f$||Ax-b||^2\rightarrow\min\f$.
    \ingroup vc_blas_algorithms
    Uses complete orthogonal factorization, yields least-norm solution for
    underdetermined system, see lapack::gelsy().
    \param[in,out] _a matrix A (will be destroyed)
    \param[in,out] _b right hand side on input, solution on output
    \param _rcond inverse condition number to determine effective rank
    \param _work ...
    \param _lwork temporary storage as for lapack::gelsy(),
    default uses lapack::Workspace::instance()
    \return effective rank of _a w.r.t. _rcond (==_a.n() on success)
 */
template <typename T,int MA,int NA,int LDA,int MB,int NB,int LDB>
int lsq_solve(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
              const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b,
              const T& _rcond,
              T* _work=0,int _lwork=0,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_b.m()==_a.m());
  assert(_rcond>=T(0));
  INTEGER* jpvt; //=(INTEGER*) allocXa(_a.n()*sizeof(INTEGER));
  int rank=-1,lwork=_lwork;
  T* work=_work;
  if (work==0) {
    lwork=lapack::gelsy_lwork(_a.m(),_a.n(),_b.n());
    _ws.alloc(&work,lwork,&jpvt,_a.n());
  }
  else
    _ws.alloc(&jpvt,_a.n());

  for (int i=0;i<_a.n();++i)
    jpvt[i]=0; // treat all columns as "free"

  int info=lapack::gelsy(_a.m(),_a.n(),_b.n(),_a.data(),_a.ld(),
                         _b.data(),_b.ld(),jpvt,_rcond,rank,work,lwork);
  assert(info==0);
  use_nowarn(info);

  return rank;
}

/** Solve least-squares problem \f$||Ax-b||^2\rightarrow\min\f$.
    \ingroup vc_blas_algorithms
    Uses complete orthogonal factorization, yields least-norm solution for
    underdetermined system, see lapack::gelsy().
    \param[in,out] _a matrix A (will be destroyed)
    \param[in,out] _b right hand side on input, solution on output
    \param _rcond inverse condition number to determine effective rank
    \param _work ...
    \param _lwork temporary storage as for lapack::gelsy(),
    default uses lapack::Workspace::instance()
    \return effective rank of _a w.r.t. _rcond (==_a.n() on success)
 */
template <typename T,int MA,int NA,int LDA,int NB,int INCB>
int lsq_solve(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
              const vector_reference_t<T,vec<NB,INCB> >& _b,
              const T& _rcond,
              T* _work=0,int _lwork=0,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  return lsq_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.n()),
                   _rcond,_work,_lwork,_ws);
}


//-----------------------------------------------------------------------------

/** Solve least-squares problem \f$||Ax-b||^2\rightarrow\min\f$.
    \ingroup vc_blas_algorithms
    Uses SVD, yields least-norm solution for underdetermined system,
    see lapack::gelsd().
    \param[in,out] _a matrix A (will be destroyed)
    \param[in,out] _b right hand side on input, solution on output
    \param[out] _s singular values of _a in descreasing order
    \param _rcond inverse condition number to determine effective rank
    \param _work ...
    \param _iwork ...
    \param _lwork temporary storage as for lapack::gelsd(),
    default uses lapack::Workspace::instance()
    \return effective rank of _a w.r.t. _rcond (==_a.n() on success)
 */
template <typename T,int MA,int NA,int LDA,int MB,int NB,int LDB,int NS>
int svd_solve(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
              const matrix_reference_t<T,ge_mat<NoT,MB,NB,LDB> >& _b,
              const vector_reference_t<T,vec<NS,1> >& _s,
              const T& _rcond,
              T* _work=0,int _lwork=0,INTEGER* _iwork=0,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_b.m()==_a.m());
  assert(_s.n()>=std::min(_a.m(),_a.n()));
  assert(_rcond>=T(0));
  INTEGER* iwork=_iwork;
  int rank=-1,lwork=_lwork;
  T* work=_work;
  if (work==0 || iwork==0) {
    lwork=lapack::gelsd_lwork(_a.m(),_a.n(),_b.n());
    int liwork=lapack::gelsd_liwork(_a.m(),_a.n(),_b.n());
    _ws.alloc(&work,lwork,&iwork,liwork);
  }
  int info=lapack::gelsd(_a.m(),_a.n(),_b.n(),_a.data(),_a.ld(),
                         _b.data(),_b.ld(),_s,_rcond,rank,work,lwork,iwork);
  assert(info==0);
  use_nowarn(info);

  return rank;
}

/** Solve least-squares problem \f$||Ax-b||^2\rightarrow\min\f$.
    \ingroup vc_blas_algorithms
    Uses SVD, yields least-norm solution for underdetermined system,
    see lapack::gelsd().
    \param[in,out] _a matrix A (will be destroyed)
    \param[in,out] _b right hand side on input, solution on output
    \param[out] _s singular values of _a in descreasing order
    \param _rcond inverse condition number to determine effective rank
    \param _work ...
    \param _iwork ...
    \param _lwork temporary storage as for lapack::gelsd(),
    default uses lapack::Workspace::instance()
    \return effective rank of _a w.r.t. `_rcond` (`==_a.n()` on success)
 */
  template <typename T,int MA,int NA,int LDA,int NB,int INCB,int NS>
int svd_solve(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
              const vector_reference_t<T,vec<NB,INCB> >& _b,
              const vector_reference_t<T,vec<NS,1> >& _s,
              const T& _rcond,
              T* _work=0,int _lwork=0,INTEGER* _iwork=0,
              lapack::Workspace& _ws=lapack::Workspace::instance()) {
  return svd_solve(_a,mat<GE>::ref(_b.data(),_b.n(),1,_b.n()),
                   _s,_rcond,_work,_lwork,_iwork,_ws);
  }

//-----------------------------------------------------------------------------

/** Compute eigenvalues and eigenvectors of symmetric matrix.
    \ingroup vc_blas_algorithms
    \param[in,out] _a input matrix, on output _a holds orthonormal eigenvectors
    if `_compute_evecs==true`, otherwise the input triangle is destroyed
    \param[out] _w eigenvalues in ascending order
    \param _compute_evecs compute only eigenvalues if `false`
    \param _work ...
    \param _lwork temporary storage as for lapack::syev(),
    default uses lapack::Workspace::instance()
    \return `true` on success, `false` on failure (failed to converge)
    \sa lapack::syev
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NW>
bool eig(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
         const vector_reference_t<T,vec<NW,1> >& _w,
         bool _compute_evecs,
         T* _work=0,int _lwork=0,
         lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_w.n()==_a.n());
  int lwork=_lwork;
  T* work=_work;
  if (work==0) {
    lwork=lapack::syev_lwork(_a.n());
    _ws.alloc(&work,lwork);
  }
  int info=lapack::syev(_compute_evecs ? 'V' : 'N',UPLO,
                        _a.n(),_a.data(),_a.ld(),
                        _w.data(),work,lwork);
  assert(info>=0 && "invalid value");

  return info==0;
}

/** Compute eigenvalues and eigenvectors of symmetric matrix.
    \ingroup vc_blas_algorithms
    \arg Faster than eig() for large matrices. Slow for small matrices.
    \arg Does not provide all options
    of lapack::syevr() to find \a some eigenvalues.
    \param[in,out] _a input matrix, respective triangle is destroyed
    \param[out] _w eigenvalues in ascending order
    \param[in,out] _v eigenvectors on output
    \param _work ...
    \param _lwork ...
    \param _iwork ...
    \param _liwork temporary storage as for lapack::syev(),
    default uses lapack::Workspace::instance()
    \return `true` on success, `false` on failure (failed to converge)
    \sa lapack::syevr
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NW,int MV,int NV,int LDV>
bool eig(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
         const vector_reference_t<T,vec<NW,1> >& _w,
         const matrix_reference_t<T,ge_mat<NoT,MV,NV,LDV> >& _v,
         T* _work=0,int _lwork=0,lapack::INTEGER* _iwork=0,int _liwork=0,
         lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_w.n()==_a.n());
  assert(_v.m()==_v.n());
  assert(_v.n()==_a.n());

  int lwork=_lwork, liwork=_liwork,m;
  T* work=_work;
  INTEGER* iwork=_iwork,*isuppz=0; //=(INTEGER*) allocXa(sizeof(lapack::INTEGER)*2*_a.n());
  if (work==0 || _iwork==0) {
    lwork=lapack::syevr_lwork(_a.n(),liwork);
    _ws.alloc(&work,lwork, &iwork,liwork, &isuppz,2*_a.n());
  }
  else {
    _ws.alloc(&isuppz,2*_a.n());
  }
  assert(isuppz!=0);
  int info=lapack::syevr('V','A',UPLO,_a.n(),_a.data(),_a.ld(),
                         T(0),T(0),0,0,T(0),m,_w.data(),_v.data(),_v.ld(),
                         isuppz,work,lwork,iwork,liwork);

  assert(info>=0 && "invalid value");
  assert(m==_a.n() || info>0);

  return info==0;
}
/** Compute eigenvalues of symmetric matrix.
    \ingroup vc_blas_algorithms
    \arg Faster than eig() for large matrices, slow for small matrices.
    \arg Does not provide all options
    of lapack::syevr() to find \a some eigenvalues.
    \param[in,out] _a input matrix, respective triangle is destroyed
    \param[out] _w eigenvalues in ascending order
    \param _work ...
    \param _lwork ...
    \param _iwork ...
    \param _liwork temporary storage as for lapack::syev(),
    default uses lapack::Workspace::instance()
    \return `true` on success, `false` on failure (failed to converge)
    \sa lapack::syevr
 */
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NW,int MV,int NV,int LDV>
bool eig(const matrix_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
         const vector_reference_t<T,vec<NW,1> >& _w,
         T* _work=0,int _lwork=0,lapack::INTEGER* _iwork=0,int _liwork=0,
         lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_w.n()==_a.n());

  int lwork=_lwork, liwork=_liwork,m;
  T* work=_work;
  lapack::INTEGER* iwork=_iwork,*isuppz;
  //=(lapack::INTEGER*) allocXa(sizeof(lapack::INTEGER)*2*_a.n());
  if (work==0 || _iwork==0) {
    lwork=lapack::syevr_lwork(_a.n(),liwork);
    _ws.alloc(&work,lwork,&iwork,liwork,&isuppz,2*_a.n());
  }
  else
    _ws.alloc(&isuppz,2*_a.n());

  int info=lapack::syevr('N','A',UPLO,_a.n(),_a.data(),_a.ld(),
                         T(0),T(0),0,0,T(0),m,_w.data(),0,0,
                         isuppz,work,lwork,iwork,liwork);

  assert(info>=0 && "invalid value");
  assert(m==_a.n() || info>0);

  return info==0;
}

//-----------------------------------------------------------------------------

/** Compute SVD A=U*S*V'.
    \ingroup vc_blas_algorithms
    Computes the SVD.
    \arg An "economy size" version is computed for <tt>_u.n()==k</tt> and/or
    <tt>_vt.m()=k ,k=min(_a.m(),_a.n())</tt>.
    \arg The matrix _vt contains the \a transposed V'.
    \arg The input matrix _a is destroyed.
    \param[in,out] _a m x n input, will be destroyed
    \param[out] _u m x m or m x min(m,n) orthogonal matrix U (output)
    \param[out] _s singular values as on diaglional of S, descending (output)
    \param[out] _vt n x n or min(m,n) x n orthogonal matrix V' (output)
    \param _work ...
    \param _lwork temporary storage as for lapack::syev(),
    default uses lapack::Workspace::instance()
    \return `true` on success or `false` (failed to converge)
    \sa lapack::gesvd()
 */
template <typename T,
          int MA,int NA,int LDA,
          int NS,int MU,int NU,int LDU, int MVT,int NVT,int LDVT>
bool svd(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
         const matrix_reference_t<T,ge_mat<NoT,MU,NU,LDU> >& _u,
         const vector_reference_t<T,vec<NS,1> >& _s,
         const matrix_reference_t<T,ge_mat<NoT,MVT,NVT,LDVT> >& _vt,
         T* _work=0,int _lwork=0,
         lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_s.n()>=std::min(_a.m(),_a.n()));
  int k=std::min(_a.m(),_a.n());
  use_nowarn(k);
  assert(_u.m()==_u.n() || _u.n()==k);
  assert(_vt.m()==_vt.n() || _vt.m()==k);
  assert(_u.m()==_a.m());
  assert(_vt.n()==_a.n());

  int lwork=_lwork;
  T* work=_work;
  if (work==0) {
    lwork=lapack::gesvd_lwork(_a.m(),_a.n());
    _ws.alloc(&work,lwork);
  }
  int info=lapack::gesvd(_u.n()==_a.m() ? 'A' : 'S',
                         _vt.m()==_a.n() ? 'A' : 'S',
                         _a.m(),_a.n(),_a.data(),_a.ld(),
                         _s.data(),_u.data(),_u.ld(),_vt.data(),_vt.ld(),
                         work,lwork);
  assert(info>=0 && "invalid parameters");

  return info==0;
}

/** Compute singular values of _a.
    \ingroup vc_blas_algorithms
    Computes \a complete SVD (versus economy size).
    \param[in,out] _a m x n input, will be destroyed
    \param[out] _s singular values, descending (output)
    \param _work ...
    \param _lwork temporary storage as for lapack::syev(),
    default uses lapack::Workspace::instance()
    \return `true` on success or `false` (failed to converge)
    \sa lapack::gesvd()
 */
template <typename T,
          int MA,int NA,int LDA,int NS>
bool svd(const matrix_reference_t<T,ge_mat<NoT,MA,NA,LDA> >& _a,
         const vector_reference_t<T,vec<NS,1> >& _s,
         T* _work=0,int _lwork=0,
         lapack::Workspace& _ws=lapack::Workspace::instance()) {
  assert(_s.n()>=std::min(_a.m(),_a.n()));

  int lwork=_lwork;
  T* work=_work;
  if (work==0) {
    lwork=lapack::gesvd_lwork(_a.m(),_a.n());
    _ws.alloc(&work,lwork);
  }
  int info=lapack::gesvd('N','N',_a.m(),_a.n(),_a.data(),_a.ld(),
                         _s.data(),0,1,0,1,
                         work,lwork);
  assert(info>=0 && "invalid parameters");

  return info==0;
}

//-----------------------------------------------------------------------------

  // generalized eigenvectors

//-----------------------------------------------------------------------------

  // generalized linar least squares (fit quadric)

//-----------------------------------------------------------------------------

/** Determine effective rank from singular values _s of a matrix.
    \param _s singular values as returned by svd()
    \param _tol tolerance >= 0
    \return effective rank
    \sa rank(), svd()
 */
template <typename T,int N>
int rank_s(const vector_const_reference_t<T,vec<N,1> >& _s,T _tol) {
  assert(_tol>T(0));

  int k=_s.n()-1;
  while (k>=0 && _s(k)<=_tol) --k;

  return k+1;
}

/** Determine effective rank of _a.
    Copies _a to temporary storage (_a is not overwritten), applies svd()
    and counts singular values which are > _tol.
    \param[in] _a input matrix
    \param _tol tolerance (default: max(m,n)*sigma_max*epsilon())
    \return effective rank
    \sa rank_s(), svd()
 */
template <typename T,typename A>
int rank(const matrix_const_reference_t<T,A>& _a,T _tol=-1) {
  // copy _a -- a,s,work shared the same space; have to allocate all at once
  int lwork=lapack::gesvd_lwork(_a.m(),_a.n());
  T*  work;
  int k=std::min(_a.m(),_a.n());
  lapack::Workspace::instance().alloc(&work,lwork+k+_a.m()*_a.n());
  typename mat<GE>::ref_t<T>::type a=mat<GE>::ref(work+lwork+k,_a.m(),_a.n());

  a=_a;
  // do svd
  typename vec<>::ref_t<T>::type s=vec<>::ref(work+lwork,k);
  bool success=svd(a,s,work,lwork);
  assert(success);
  use_nowarn(success);

  if (_tol<T(0))
    _tol=std::max(a.m(),a.n())*s(0)*std::numeric_limits<T>::epsilon();

  // count effektive rank
  --k;
  while (k>=0 && s(k)<=_tol) --k;

  return k+1;
}

/** Determine condition number of _a w.r.t. L_2 norm.
    Copies _a to temporary storage (_a is not overwritten) and applies svd().
    \param[in] _a input matrix
    \return condition number w.r.t. L_2 norm
    \sa svd()
 */
template <typename T,typename A>
T cond2(const matrix_const_reference_t<T,A>& _a) {
  // copy _a -- a,s,work shared the same space; have to allocate all at once
  int lwork=lapack::gesvd_lwork(_a.m(),_a.n());
  T*  work;
  int k=std::min(_a.m(),_a.n());
  lapack::Workspace::instance().alloc(&work,lwork+k+_a.m()*_a.n());
  typename mat<GE>::ref_t<T>::type a=mat<GE>::ref(work+lwork+k,_a.m(),_a.n());

  a=_a;
  // do svd
  typename vec<>::ref_t<T>::type s=vec<>::ref(work+lwork,k);
  bool success=svd(a,s,work,lwork);
  assert(success);
  use_nowarn(success);

  return s(0)/s(k-1);
}

/** Determine L_2 norm of matrix _a.
    Copies _a to temporary storage (_a is not overwritten) and applies svd().
    \param[in] _a input matrix
    \return ||_a||_2
    \sa svd()
 */
template <typename T,typename A>
T norm2(const matrix_const_reference_t<T,A>& _a) {
  // copy _a -- a,s,work shared the same space; have to allocate all at once
  int lwork=lapack::gesvd_lwork(_a.m(),_a.n());
  T*  work;
  int k=std::min(_a.m(),_a.n());
  lapack::Workspace::instance().alloc(&work,lwork+k+_a.m()*_a.n());
  typename mat<GE>::ref_t<T>::type a=mat<GE>::ref(work+lwork+k,_a.m(),_a.n());

  a=_a;
  // do svd
  typename vec<>::ref_t<T>::type s=vec<>::ref(work+lwork,k);
  bool success=svd(a,s,work,lwork);
  assert(success);
  use_nowarn(success);

  return s(0);
}

/// Determine L_2 norm of vector _x.
template <typename T,typename V>
T norm2(const vector_const_reference_t<T,V>& _x) {
  return _x.nrm2();
}

/** Determine Frobenius norm of _a.
    \param[in] _a input matrix
    \return ||_a||_2
    \sa svd()
 */
template <typename T,TransposeFlag TRANS,int N,int M,int LD>
  T normf(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a) {
  T sum(0);
  for (int j=0;j<_a.n();++j)
    for (int i=0;i<_a.m();++i) {
      T d=_a.data()[j*_a.ld()+i];
      sum+=d*d;
    }
  return sqrt(sum);
}

/** Get pseudo-inverse _pinva of _a.
    Copies _a to temporary storage (_a is not overwritten) and applies svd().
    \param[in] _a m x n input matrix
    \param[out] _pinva n x m pseudo-innverse of _a as output
    \param  _tol tolerance (default: max(m,n)*sigma_max*epsilon())
    \return `false` on failure
 */
template <typename T,typename A,int M,int N,int LD>
bool pinv(const matrix_const_reference_t<T,A>& _a,
          const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _pinva,
          T _tol=T(-1)) {
  assert(_a.m()==_pinva.n() && _a.n()==_pinva.m());

  // copy _a -- a,s,u,v,work shared the same space; have to allocate all at once
  int lwork=lapack::gesvd_lwork(_a.m(),_a.n());
  T*  work;
  int m=_a.m(),n=_a.n(),k=std::min(_a.m(),_a.n());
  lapack::Workspace::instance().alloc(&work,lwork+ k + m*k + n*k + m*n);
  typename mat<GE>::ref_t<T>::type a=mat<GE>::ref(work+lwork+k+m*k+n*k,_a.m(),_a.n());

  a=_a;
  // do svd
  typename vec<>::ref_t<T>::type     s=vec<>::ref(work+lwork,k);
  typename mat<GE>::ref_t<T>::type   u=mat<GE>::ref(work+lwork+k,m,k);
  typename mat<GE>::ref_t<T>::type  vt=mat<GE>::ref(work+lwork+k+m*k,k,n);
  bool success=svd(a,u,s,vt,work,lwork);
  if (!success)
    return false; // this should never happen!

  // invert diag(s)
  if (_tol<T(0))
    _tol=std::max(a.m(),a.n())*s(0)*std::numeric_limits<T>::epsilon();

  for (int i=0;i<k;++i)
    if (fabs(s.data()[i])>_tol)
      s.data()[i]=T(1)/s.data()[i];
    else
      s.data()[i]=T(0);

  // get pseudo-inverse
  trans(vt)*=diag(s);
  _pinva=trans(vt)*trans(u);
  //_pinva.mm(T(1),trans(vt).const_ref(),trans(u).const_ref(),T(0));

  return true;
}

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

template <typename T,int N>
int rank_s(const vector_reference_t<T,vec<N,1> >& _s,T _tol) {
  return rank_s(_s.const_ref(),_tol);
}
template <typename T,typename A>
int rank(const matrix_reference_t<T,A>& _a,T _tol=-1) {
  return rank(_a.const_ref(),_tol);
}
template <typename T,typename A>
T cond2(const matrix_reference_t<T,A>& _a) {
  return cond2(_a.const_ref());
}
template <typename T,typename A>
T norm2(const matrix_reference_t<T,A>& _a) {
  return norm2(_a.const_ref());
}
/// Determine L_2 norm of vector _x.
template <typename T,typename V>
T norm2(const vector_reference_t<T,V>& _x) {
  return norm2(_x.const_ref());
}
template <typename T,TransposeFlag TRANS,int N,int M,int LD>
  T normf(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a) {
  return normf(_a.const_ref());
}
template <typename T,typename A,int M,int N,int LD>
bool pinv(const matrix_reference_t<T,A>& _a,
          const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _pinva,
          T _tol=T(-1)) {
  return pinv(_a.const_ref(),_pinva,_tol);
}

# endif

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
