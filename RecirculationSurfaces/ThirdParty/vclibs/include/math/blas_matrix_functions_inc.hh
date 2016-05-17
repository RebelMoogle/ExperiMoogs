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

    \arg Non-BLAS functions in VC::math::blas dealing with vectors
    and matrices.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

#ifdef DOXYGEN_SKIP

  //
  // documentation only
  //

/// get reference to column vector (GE only) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_reference_t<T,V>
column(const matrix_reference_t<T,A>& _a,int _j);
/// get reference to column vector (GE only) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_const_reference_t<T,V>
column(const matrix_const_reference_t<T,A>& _a,int _j);

/// get reference to row vector (GE only) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_reference_t<T,V>
row(const matrix_reference_t<T,A>& _a,int _j);
/// get reference to row vector (GE only) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_const_reference_t<T,V>
row(const matrix_const_reference_t<T,A>& _a,int _j);

/// get reference to diagonal as vector (GE,GB,SY,SB,TR,TB w/ NoU) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_reference_t<T,V>
diag(const matrix_reference_t<T,A>& _a,int _j);
/// get reference to diagonal as vector (GE,GB,SY,SB,TR,TB w/ NoU) \ingroup vc_blas
template <typename T,typename V,typname A>
vector_const_reference_t<T,V>
diag(const matrix_const_reference_t<T,A>& _a,int _j);

/// get reference to subvector _x(_i1:_i2) \ingroup vc_blas
template <typename T,typename V>
vector_reference_t<T,V>
block(const vector_reference_t<T,V>& _x,int _i1,_int _i2);
/// get reference to subvector _x(_i1:_i2) \ingroup vc_blas
template <typename T,typename V>
const_vector_reference_t<T,V>
block(const const_vector_reference_t<T,V>& _x,int _i1,_int _i2);

/// get reference to submatrix _a(_i1:_i2, _j1:_j2) (GE only), see cp_block() \ingroup vc_blas
template <typename T,typename V,typname A,typname B>
matrix_reference_t<T,B>
block(const matrix_reference_t<T,A>& _a,int _i1,_int _i2,int _j1,int _j2);
/// get reference to submatrix _a(_i1:_i2, _j1:_j2) (GE only), see cp_block() \ingroup vc_blas
template <typename T,typename V,typname A,typname B>
matrix_const_reference_t<T,B>
block(const matrix_const_reference_t<T,A>& _a,int _i1,_int _i2,int _j1,int _j2);

/// copy column _j of _a to vector _x \ingroup vc_blas
template <typename T,typename V,typname A>
void cp_column(const matrix_const_reference_t<T,A>& _a,
               const vector_reference_t<T,V>& _x,int _j);
/// copy row _i of _a to vector _x \ingroup vc_blas
template <typename T,typename V,typname A>
void cp_row(const matrix_const_reference_t<T,A>& _a,
            const vector_reference_t<T,V>& _x,int _i);
/// copy diagonal of _a to vector _x \ingroup vc_blas
template <typename T,typename V,typname A>
void cp_diag(const matrix_const_reference_t<T,A>& _a,
             const vector_reference_t<T,V>& _x);

/// copy submatrix _a(_i1:_i2, _j1:_j2) to _b (which is GE) \ingroup vc_blas
template <typename T,typename A,typename B>
void cp_block(const matrix_const_reference_t<T,A>& _a,
              int _i1,int _i2,int _j1,int _j2,
              const matrix_reference_t<T,B>& _b);

/// copy matrix _b to _a (_a is GE, or SY->TP, SP->SY, TR->TP, TP->TR) \ingroup vc_blas
template <typename T,typename A,typename B>
void cp(const matrix_const_reference_t<T,B>& _b,
        const matrix_reference_t<T,A>& _a);

/// copy data _ary element-wise to vector (no checks) \ingroup vc_blas
template <typename S,typename T,typename V>
void cp(const S* _ary,const vector_reference_t<T,V>& _x);
/// copy vector elements to linear array _ary (no checks) \ingroup vc_blas
template <typename S,typename T,typename V>
void cp(const vector_const_reference_t<T,V>& _x,S* _ary);

/// copy data _ary element-wise to matrix (no checks) \ingroup vc_blas
template <typename S,typename T,typename A>
void cp(const S* _ary,const matrix_reference_t<T,A>& _a);
/// copy matrix elements to linear array _ary (no checks) \ingroup vc_blas
template <typename S,typename T,typename A>
void cp(const matrix_const_reference_t<T,A>& _a,S* _ary);

/// expression for printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename V>
Expression<> print(const vector_const_reference_t<T,V>& _x);
/// expression for printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename A>
Expression<> print(const matrix_const_reference_t& _A);
/// expression for pretty printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename V>
Expression<> pp(const vector_const_reference_t<T,V>& _x);
/// expression for pretty printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename A>
Expression<> pp(const matrix_const_reference_t<T,A>& _A);
/// expression for matlab-style printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename V>
Expression<> pmatlab(const vector_const_reference_t<T,V>& _x);
/// expression for matlab-style printing using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename A>
Expression<> pmatlab(const matrix_const_reference_t<T,A>& _A);
/// expression for printing stored data using ostream::operator<<, e.g., cout << p(x) \ingroup vc_blas_io
typename <typename T,typename V>
Expression<> pdata(const matrix_const_reference_t<T,V>& _x);

#else // => !defined(DOXYGEN_SKIP)

//=============================================================================

//
// reference rows and columns of ge_mat
//

  // NoT

template <typename T,int M,int N,int LD>
vector_reference_t<T,vec<M,1> >
column(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _j) {
  assert(0<=_j && _j<_a.n());
  return vector_reference_t<T,vec<M,1> >
    (vec<M,1>(_a.m(),1),_a.data()+_j*_a.matrix().ld());
}
template <typename T,int M,int N,int LD>
vector_const_reference_t<T,vec<M,1> >
column(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _j) {
  assert(0<=_j && _j<_a.n());
  return vector_const_reference_t<T,vec<M,1> >
    (vec<M,1>(_a.m(),1),_a.data()+_j*_a.matrix().ld());
}

template <typename T,int M,int N,int LD>
vector_reference_t<T,vec<N,LD> >
row(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i) {
  assert(0<=_i && _i<_a.m());
  return vector_reference_t<T,vec<N,LD> >
    (vec<N,LD>(_a.n(),_a.matrix().ld()),_a.data()+_i);
}
template <typename T,int M,int N,int LD>
vector_const_reference_t<T,vec<N,LD> >
row(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i) {
  assert(0<=_i && _i<_a.n());
  return vector_const_reference_t<T,vec<N,LD> >
    (vec<N,LD>(_a.n(),_a.matrix().ld()),_a.data()+_i);
}

  // Transposed

template <typename T,int M,int N,int LD>
vector_reference_t<T,vec<N,LD> >
column(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _j) {
  return row(matrix_reference_t<T,ge_mat<NoT,M,N,LD> >
             (ge_mat<NoT,M,N,LD>(_a.m(),_a.n(),_a.matrix().ld()),_a.data()),_j);
}
template <typename T,int M,int N,int LD>
vector_const_reference_t<T,vec<N,LD> >
column(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _j) {
  return row(matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >
             (ge_mat<NoT,M,N,LD>(_a.m(),_a.n(),_a.matrix().ld()),_a.data()),_j);
}

template <typename T,int M,int N,int LD>
vector_reference_t<T,vec<M,1> >
row(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i) {
 return column(matrix_reference_t<T,ge_mat<NoT,M,N,LD> >
               (ge_mat<NoT,M,N,LD>(_a.m(),_a.n(),_a.matrix().ld()),_a.data()),_i);
}
template <typename T,int M,int N,int LD>
vector_reference_t<T,vec<M,1> >
row(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i) {
  return column(matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >
                (ge_mat<NoT,M,N,LD>(_a.m(),_a.n(),_a.matrix().ld()),_a.data()),_i);
}

//-----------------------------------------------------------------------------

//
// reference diagonal of ge_mat
//

template <typename T,TransposeFlag TRANS,int M,int N,int LD>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}

//
// reference diagonal of gb_mat
//

template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data()+_a.ku());
}
template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data()+_a.ku());
}

//
// reference diagonal of sy_mat
//

template <typename T,int N,int LD,UpperLowerFlag UPLO>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}
template <typename T,int N,int LD,UpperLowerFlag UPLO>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}

//
// reference diagonal of sb_mat
//

template <typename T,int N,int K,int LD,UpperLowerFlag UPLO>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),
     _a.data()+(UPLO==Upper ? _a.k() : 0));
}
template <typename T,int N,int K,int LD,UpperLowerFlag UPLO>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),
     _a.data()+(UPLO==Upper ? _a.k() : 0));
}

//
// reference diagonal of tr_mat
//

template <typename T,int N,int LD,TransposeFlag TRANS,UpperLowerFlag UPLO>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,tr_mat<UPLO,TRANS,NoU,N,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}
template <typename T,int N,int LD,TransposeFlag TRANS,UpperLowerFlag UPLO>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,NoU,N,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),_a.data());
}

//
// reference diagonal of tb_mat
//

template <typename T,int N,int K,int LD,UpperLowerFlag UPLO,TransposeFlag TRANS>
vector_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_reference_t<T,tb_mat<UPLO,TRANS,NoU,N,K,LD> >& _a) {
  return vector_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),
     _a.data()+(UPLO==Upper ? _a.k() : 0));
}
template <typename T,int N,int K,int LD,UpperLowerFlag UPLO,TransposeFlag TRANS>
vector_const_reference_t<T,vec<VarInt,VarInt> >
diag(const matrix_const_reference_t<T,tb_mat<UPLO,TRANS,NoU,N,K,LD> >& _a) {
  return vector_const_reference_t<T,vec<VarInt,VarInt> >
    (vec<VarInt,VarInt>(std::min(_a.m(),_a.n()),_a.matrix().ld()+1),
     _a.data()+(UPLO==Upper ? _a.k() : 0));
}

//-----------------------------------------------------------------------------

//
// reference subvector
//

template <typename T,int N,int INC>
vector_reference_t<T,vec<VarInt,INC> >
block(const vector_reference_t<T,vec<N,INC> >& _x,int _i1,int _i2) {
  _VC_DBG_BLAS_SCOPE();
  assert(0<=_i1); assert(_i2<_x.n()); assert(_i1<=_i2);
  int n=_i2-_i1+1;
  return vector_reference_t<T,vec<VarInt,INC> >
    (vec<VarInt,INC>(n,_x.vector().inc()),_x.data()+_i1*_x.vector().inc());
}
template <typename T,int N,int INC>
vector_const_reference_t<T,vec<VarInt,INC> >
block(const vector_const_reference_t<T,vec<N,INC> >& _x,int _i1,int _i2) {
  _VC_DBG_BLAS_SCOPE();
  assert(0<=_i1); assert(_i2<_x.n()); assert(_i1<=_i2);
  int n=_i2-_i1+1;
  return vector_const_reference_t<T,vec<VarInt,INC> >
    (vec<VarInt,INC>(n,_x.vector().inc()),_x.data()+_i1*_x.vector().inc());
}

// using idx_range

template <typename T,int N,int INC,bool END1,bool END2>
vector_reference_t<T,vec<VarInt,INC> >
block(const vector_reference_t<T,vec<N,INC> >& _x,const idx_range<END1,END2>& _r) {
  _VC_DBG_BLAS_SCOPE();
  int n=_x.n();
  return block(_x,_r.eval_i1(n),_r.eval_i2(n));
}
template <typename T,int N,int INC,bool END1,bool END2>
vector_const_reference_t<T,vec<VarInt,INC> >
block(const vector_const_reference_t<T,vec<N,INC> >& _x,const idx_range<END1,END2>& _r) {
  _VC_DBG_BLAS_SCOPE();
  int n=_x.n();
  return block(_x,_r.eval_i1(n),_r.eval_i2(n));
}

//
// reference submatrix in ge_mat
//

template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      int _i1,int _i2,int _j1,int _j2) {
  _VC_DBG_BLAS_SCOPE();
  assert(0<=_i1); assert(_i2<_a.m_rows()); assert(_i1<=_i2);
  assert(0<=_j1); assert(_j2<_a.n_cols()); assert(_j1<=_j2);
  int m=_i2-_i1+1, n=_j2-_j1+1;
  return  matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
    (ge_mat<TRANS,VarInt,VarInt,LD>(m,n,_a.matrix().ld()),
     _a.data()+_a.matrix().ld()*_j1+_i1);
}

template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      int _i1,int _i2,int _j1,int _j2) {
  _VC_DBG_BLAS_SCOPE();
  assert(0<=_i1); assert(_i2<_a.m_rows()); assert(_i1<=_i2);
  assert(0<=_j1); assert(_j2<_a.n_cols()); assert(_j1<=_j2);
  int m=_i2-_i1+1, n=_j2-_j1+1;
  return  matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
    (ge_mat<TRANS,VarInt,VarInt,LD>(m,n,_a.matrix().ld()),
     _a.data()+_a.matrix().ld()*_j1+_i1);
}

// using idx_range

template <typename T,TransposeFlag TRANS,int M,int N,int LD,
          bool END1,bool END2,bool END3,bool END4>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n();
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),_rj.eval_i1(n),_rj.eval_i2(n));
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,
          bool END1,bool END2,bool END3,bool END4>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n();
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),_rj.eval_i1(n),_rj.eval_i2(n));
}

// A(i1:i2,j)
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,int _j) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m();
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),_j,_j);
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,int _j) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m();
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),_j,_j);
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,const idx_end& _j) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n(), j=n-_j.offset-1;
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),j,j);
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_range<END1,END2>& _ri,const idx_end& _j) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n(), j=n-_j.offset-1;
  return block(_a,_ri.eval_i1(m),_ri.eval_i2(m),j,j);
}

// A(i,j1:j2)
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      int _i,const idx_range<END1,END2>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int n=_a.n();
  return block(_a,_i,_i,_rj.eval_i1(n),_rj.eval_i2(n));
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      int _i,const idx_range<END1,END2>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int n=_a.n();
  return block(_a,_i,_i,_rj.eval_i1(n),_rj.eval_i2(n));
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_end& _i,const idx_range<END1,END2>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n(), i=m-_i.offset-1;
  return block(_a,i,i,_rj.eval_i1(n),_rj.eval_i2(n));
}
template <typename T,TransposeFlag TRANS,int M,int N,int LD,bool END1,bool END2>
matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> >
block(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
      const idx_end& _i,const idx_range<END1,END2>& _rj) {
  _VC_DBG_BLAS_SCOPE();
  int m=_a.m(), n=_a.n(), i=m-_i.offset-1;
  return block(_a,i,i,_rj.eval_i1(n),_rj.eval_i2(n));
}

//-----------------------------------------------------------------------------

template <typename S,typename T,typename V>
void cp(const S* _ary,const vector_reference_t<T,V>& _x) {
  for (int i=0;i<_x.n();++i)
    _x.data()[i*_x.inc()]=T(*_ary++);
}
template <typename S,typename T,typename V>
void cp(const vector_const_reference_t<T,V>& _x,S* _ary) {
  for (int i=0;i<_x.n();++i)
    *_ary++=_x.data()[i*_x.inc()];
}
  // matrix versions are defined in blas_matrix_foreach_inc.hh


//
// general functions for copying columns, rows, diagonals
//

template <typename T,int N,int INC,typename A>
void cp_column(const matrix_const_reference_t<T,A>& _a,
               const vector_reference_t<T,vec<N,INC> >& _x,int _j) {
  _VC_DBG_BLAS_SCOPE();
  if (_a.trans()!=NoT)
    cp_row(_a,_x,_j);
  else {
    assert(_x.n()==_a.m());
    for (int i=0;i<_a.m();++i)
      _at(_x,i)=_at(_a,i,_j);
  }
}

template <typename T,int N,int INC,typename A>
void cp_row(const matrix_const_reference_t<T,A>& _a,
            const vector_reference_t<T,vec<N,INC> >& _x,int _i) {
  _VC_DBG_BLAS_SCOPE();
  if (_a.trans()!=NoT)
    cp_column(_a,_x,_i);
  else {
    assert(_x.n()==_a.n());
    for (int j=0;j<_a.n();++j)
      _at(_x,j)=_at(_a,_i,j);
  }
}

template <typename T,int N,int INC,typename A>
void cp_diag(const matrix_const_reference_t<T,A>& _a,
             const vector_reference_t<T,vec<N,INC> >& _x) {
  _VC_DBG_BLAS_SCOPE();
  int k=std::min(_a.m(),_a.n());
  assert(_x.n()==k);
  for (int i=0;i<k;++i)
    _at(_x,i)=_at(_a,i,i);
}

//-----------------------------------------------------------------------------

//
// copy columns and rows of ge_mat
//

template <typename T,int M,int N,int LD>
void cp_column(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
               const vector_reference_t<T,vec<M,1> >& _x,int _i) {
  copy(column(_a,_i),_x);
}
template <typename T,int M,int N,int LD>
void cp_row(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,
            const vector_reference_t<T,vec<N,LD> >& _x,int _i) {
  copy(row(_a,_i),_x);
}

template <typename T,int M,int N,int LD>
void cp_column(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
               const vector_reference_t<T,vec<N,LD> >& _x,int _i) {
  copy(column(_a,_i),_x);
}
template <typename T,int M,int N,int LD>
void cp_row(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,
            const vector_reference_t<T,vec<M,1> >& _x,int _i) {
  copy(row(_a,_i),_x);
}

//-----------------------------------------------------------------------------

//
// copy diagonal of ge_mat
//

template <typename T,TransposeFlag TRANS,int M,int N,int NV,int INC,int LD>
void cp_diag(const matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a,
             const vector_reference_t<T,vec<NV,INC> >& _x,int /*_i*/) {
  copy(diag(_a),_x); // BLAS vector copy
}

//-----------------------------------------------------------------------------

//
// copy submatrix of ge_mat
//

//
// general function for copying a block to a ge_mat
//

template <typename T,TransposeFlag TRANS,int M,int N,int LD,typename A>
void cp_block(const matrix_const_reference_t<T,A>& _a,
              int _i1,int _i2,int _j1,int _j2,
              const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _b) {
  _VC_DBG_BLAS_SCOPE();

  //
  // slow "general" function based on _at()
  //

  assert(0<=_i1); assert(_i2<_a.m_rows()); assert(_i1<=_i2);
  assert(0<=_j1); assert(_j2<_a.n_cols()); assert(_j1<=_j2);
  int m=_i2-_i1+1, n=_j2-_j1+1;

  assert(_b.m_rows()==m);
  assert(_b.n_cols()==n);

  for (int j=0;j<n;++j)
    for (int i=0;i<m;++i)
      _at(_b,i,j)=_at(_a,_i1+i,_j1+j);
}
template <typename T,TransposeFlag TRANSA,TransposeFlag TRANSB,
          int MA,int NA,int MB,int NB,int LDA,int LDB>
void cp_block(const matrix_const_reference_t<T,ge_mat<TRANSA,MA,NA,LDA> >& _a,
              int _i1,int _i2,int _j1,int _j2,
              const matrix_reference_t<T,ge_mat<TRANSB,MB,NB,LDB> >& _b) {
  _VC_DBG_BLAS_SCOPE();

  cp(block(_a,_i1,_i2,_j1,_j2),_b); // respects TRANSA, TRANSC
}

//-----------------------------------------------------------------------------

//
// copy matrix to ge_mat
//

/// copy matrix to ge_mat \ingroup vc_blas
template <typename T,
          TransposeFlag TRANSA,int MA,int NA,int LDA,
          TransposeFlag TRANSB,int MB,int NB,int LDB>
void cp(const matrix_const_reference_t<T,ge_mat<TRANSA,MB,NB,LDB> >& _b,
        const matrix_reference_t<T,ge_mat<TRANSB,MA,NA,LDA> >& _a) {
  _VC_DBG_BLAS_SCOPE();

  T* a=_a.data();
  const T* b=_b.data();

  assert(a!=b);

  if (TRANSA==TRANSB) {
    assert(_a.m()==_b.m());
    assert(_a.n()==_b.n());

    int m=_a.m(), n=_a.n();
    int lda=_a.ld(), ldb=_b.ld();

    for (int j=0;j<n;++j)
      for (int i=0;i<m;++i)
        a[j*lda+i]=b[j*ldb+i];
  }
  else {
    assert(_a.m()==_b.n());
    assert(_a.n()==_b.m());

    int m=_a.m(), n=_a.n();
    int lda=_a.ld(), ldb=_b.ld();

    for (int j=0;j<n;++j)
      for (int i=0;i<m;++i)
        a[j*lda+i]=b[i*ldb+j];
  }
}

/// copy matrix to ge_mat \ingroup vc_blas
template <typename T,TransposeFlag TRANS,int M,int N,int LD,typename A>
void cp(const matrix_const_reference_t<T,A>& _b,
        const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a) {
  _VC_DBG_BLAS_SCOPE();

  //
  // general function based on cp_column()
  //

  assert(_a.data()!=_b.data());

  for (int j=0;j<_a.n();++j)
    cp_column(_b,column(_a,j),j);
}

//
// SY -> SP
//

/// copy matrix from sy_mat to sp_mat \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int NA,int NB,int LDA>
void cp(const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
        const matrix_reference_t<T,sp_mat<UPLO,NB> >& _b) {

  assert(_a.data()!=_b.data());
  assert(_a.n()==_b.n());

  int n=_a.n(),k=0,ld=_a.ld();
  const T* a=_a.data();
  T* b=_b.data();

  if (UPLO==Upper) {
    for (int j=0;j<n;++j)
      for (int i=0;i<=j;++i)
        b[k++]=a[j*ld+i];
  }
  else {
    for (int j=0;j<n;++j)
      for (int i=j;i<n;++i)
        b[k++]=a[j*ld+i];
  }
}

//
// SY -> SP
//

/// copy matrix from sp_mat to sy_mat \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int NA,int NB,int LDB>
void cp(const matrix_const_reference_t<T,sp_mat<UPLO,NA> >& _a,
        const matrix_reference_t<T,sy_mat<UPLO,NB,LDB> >& _b) {

  assert(_a.data()!=_b.data());
  assert(_a.n()==_b.n());

  int n=_a.n(),k=0,ld=_a.ld();
  const T* a=_a.data();
  T* b=_b.data();

  if (UPLO==Upper) {
    for (int j=0;j<n;++j)
      for (int i=0;i<=j;++i)
        b[j*ld+i]=a[k++];
  }
  else {
    for (int j=0;j<n;++j)
      for (int i=j;i<n;++i)
        b[j*ld+i]=a[k++];
  }
}

//-----------------------------------------------------------------------------


// TODO: HERE IS A LOT MISSING!

//
// SY -> SY
//

/// @todo move/collect this somewhere else; require more of cp()
template <typename T>
struct _Identity {
  T operator()(const T& _b) { return _b;  }
};

/// copy sy_mat \ingroup vc_blas
template <typename T,UpperLowerFlag UPLO,int NA,int LDA,int NB,int LDB>
void cp(const matrix_const_reference_t<T,sy_mat<UPLO,NA,LDA> >& _a,
        const matrix_reference_t<T,sy_mat<UPLO,NB,LDB> >& _b) {
  _Identity<T> id;
  map_f2(id,_b,_a);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//
// helper CP_Block uses
// - block() (reference) for GE matrix (specialization)
// - returns an Expr_cp_block (to call cp_block()) for all other types
//

template <typename T,typename M>
struct CP_Block {

  typedef matrix_reference_t<T,M> ref_t;
  typedef matrix_const_reference_t<T,M> const_ref_t;

  typedef Expression<Type_matrix,Expr_cp_block,T,M> op_t;
  typedef Expression<Type_matrix,Expr_cp_block,T,M> const_op_t;

  template <typename ARG1,typename ARG2,typename ARG3,typename ARG4>
  static op_t
  get(const ref_t& _a,const ARG1& _i1,const ARG2& _i2,
      const ARG1& _j1,const ARG2& _j2) {
    return op_t(_a.const_ref(),_i1,_i2,_j1,_j2);
  }
  template <typename ARG1,typename ARG2>
  static op_t
  get(const ref_t& _a,const ARG1& _i,const ARG2& _j) {
    return op_t(_a.const_ref(),_i,_j);
  }

  template <typename ARG1,typename ARG2,typename ARG3,typename ARG4>
  static op_t
  get(const const_ref_t& _a,const ARG1& _i1,const ARG2& _i2,
      const ARG1& _j1,const ARG2& _j2) {
    return op_t(_a,_i1,_i2,_j1,_j2);
  }
  template <typename ARG1,typename ARG2>
  static op_t
  get(const const_ref_t& _a,const ARG1& _i,const ARG2& _j) {
    return op_t(_a,_i,_j);
  }
};

template <typename T,
          TransposeFlag TRANS,int M,int N,int LD>
struct CP_Block<T,ge_mat<TRANS,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<TRANS,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> > const_ref_t;

  typedef matrix_reference_t<T,ge_mat<TRANS,VarInt,VarInt,LD> > op_t;
  typedef matrix_const_reference_t<T,ge_mat<TRANS,VarInt,VarInt,VarInt> > const_op_t;

  template <typename ARG1,typename ARG2,typename ARG3,typename ARG4>
  static op_t
  get(const ref_t& _a,
      const ARG1& _i1,const ARG2& _i2,const ARG3& _j1,const ARG4& _j2) {
    return block(_a,_i1,_i2,_j1,_j2);
  }
  template <typename ARG1,typename ARG2>
  static op_t
  get(const ref_t& _a,const ARG1& _i,const ARG2& _j) {
    return block(_a,_i,_j);
  }

  template <typename ARG1,typename ARG2,typename ARG3,typename ARG4>
  static op_t
  get(const const_ref_t& _a,
      const ARG1& _i1,const ARG2& _i2,const ARG3& _j1,const ARG4& _j2) {
    return block(_a,_i1,_i2,_j1,_j2);
  }
  template <typename ARG1,typename ARG2>
  static op_t
  get(const const_ref_t& _a,const ARG1& _i,const ARG2& _j) {
    return block(_a,_i,_j);
  }
};

//-----------------------------------------------------------------------------

//
// helper CP_Column uses
// - column() (reference) for GE matrix (specialization)
// - returns an Expr_cp_column (to call cp_column()) for all other types
//

template <typename T,typename M>
struct CP_Column {

  typedef matrix_reference_t<T,M> ref_t;
  typedef matrix_const_reference_t<T,M> const_ref_t;

  typedef Expression<Type_matrix,Expr_cp_column,T,M> op_t;
  typedef Expression<Type_matrix,Expr_cp_column,T,M> const_op_t;

  static op_t
  get(const ref_t& _a,int _j) { return op_t(_a.const_ref(),_j);  }

  static op_t
  get(const const_ref_t& _a,int _j) { return op_t(_a.const_ref(),_j);  }
};

template <typename T,int M,int N,int LD>
struct CP_Column<T,ge_mat<NoT,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<NoT,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<M,1> > op_t;
  typedef vector_const_reference_t<T,vec<M,1> > const_op_t;

  static op_t get(const ref_t& _a,int _j) { return column(_a,_j); }
  static op_t get(const const_ref_t& _a,int _j) { return column(_a,_j); }
};
template <typename T,int M,int N,int LD>
struct CP_Column<T,ge_mat<Transpose,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<Transpose,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<N,LD> > op_t;
  typedef vector_const_reference_t<T,vec<N,LD> > const_op_t;

  static op_t get(const ref_t& _a,int _j) { return column(_a,_j); }
  static op_t get(const const_ref_t& _a,int _j) { return column(_a,_j); }
};

//-----------------------------------------------------------------------------

//
// helper CP_Row uses
// - row() (reference) for GE matrix (specialization)
// - returns an Expr_cp_row (to call cp_row()) for all other types
//

template <typename T,typename M>
struct CP_Row {

  typedef matrix_reference_t<T,M> ref_t;
  typedef matrix_const_reference_t<T,M> const_ref_t;

  typedef Expression<Type_matrix,Expr_cp_row,T,M> op_t;
  typedef Expression<Type_matrix,Expr_cp_row,T,M> const_op_t;

  static op_t
  get(const ref_t& _a,int _i) { return op_t(_a.const_ref(),_i);  }

  static op_t
  get(const const_ref_t& _a,int _i) { return op_t(_a.const_ref(),_i);  }
};

template <typename T,int M,int N,int LD>
struct CP_Row<T,ge_mat<NoT,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<NoT,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<N,LD> > op_t;
  typedef vector_const_reference_t<T,vec<N,LD> > const_op_t;

  static op_t get(const ref_t& _a,int _i) { return row(_a,_i); }
  static op_t get(const const_ref_t& _a,int _i) { return row(_a,_i); }
};
template <typename T,int M,int N,int LD>
struct CP_Row<T,ge_mat<Transpose,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<Transpose,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<M,1> > op_t;
  typedef vector_const_reference_t<T,vec<M,1> > const_op_t;

  static op_t get(const ref_t& _a,int _i) { return row(_a,_i); }
  static op_t get(const const_ref_t& _a,int _i) { return row(_a,_i); }
};

//-----------------------------------------------------------------------------

//
// helper CP_Diag uses
// - diag() (reference) for GE,GB,TR,TB,SY,SB (no unit diagonal) matrix (specialization)
// - returns an Expr_cp_diag (to call cp_diag()) for all other types
//

template <typename T,typename M>
struct CP_Diag {

  typedef matrix_reference_t<T,M> ref_t;
  typedef matrix_const_reference_t<T,M> const_ref_t;

  typedef Expression<Type_matrix,Expr_cp_diag,T,M> op_t;
  typedef Expression<Type_matrix,Expr_cp_diag,T,M> const_op_t;

  static op_t get(const ref_t& _a) { return op_t(_a.const_ref()); }
  static op_t get(const const_ref_t& _a) { return op_t(_a.const_ref()); }
};

template <typename T,
          TransposeFlag TRANS,int M,int N,int LD>
struct CP_Diag<T,ge_mat<TRANS,M,N,LD> > {

  typedef matrix_reference_t<T,ge_mat<TRANS,M,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,ge_mat<TRANS,M,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

template <typename T,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
struct CP_Diag<T,gb_mat<TRANS,M,N,KL,KU,LD> > {

  typedef matrix_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> > ref_t;
  typedef matrix_const_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

template <typename T,
          UpperLowerFlag UPLO,int N,int LD>
struct CP_Diag<T,sy_mat<UPLO,N,LD> > {

  typedef matrix_reference_t<T,sy_mat<UPLO,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,sy_mat<UPLO,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

template <typename T,
          int N,int K,int LD,UpperLowerFlag UPLO>
struct CP_Diag<T,sb_mat<UPLO,N,K,LD> > {

  typedef matrix_reference_t<T,sb_mat<UPLO,N,K,LD> > ref_t;
  typedef matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

template <typename T,
          int N,int LD,TransposeFlag TRANS,UpperLowerFlag UPLO>
struct CP_Diag<T,tr_mat<UPLO,TRANS,NoU,N,LD> > {

  typedef matrix_reference_t<T,tr_mat<UPLO,TRANS,NoU,N,LD> > ref_t;
  typedef matrix_const_reference_t<T,tr_mat<UPLO,TRANS,NoU,N,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

template <typename T,
          int N,int K,int LD,UpperLowerFlag UPLO,TransposeFlag TRANS>
struct CP_Diag<T,tb_mat<UPLO,TRANS,NoU,N,K,LD> > {

  typedef matrix_reference_t<T,tb_mat<UPLO,TRANS,NoU,N,K,LD> > ref_t;
  typedef matrix_const_reference_t<T,tb_mat<UPLO,TRANS,NoU,N,K,LD> > const_ref_t;

  typedef vector_reference_t<T,vec<VarInt,VarInt> > op_t;
  typedef vector_const_reference_t<T,vec<VarInt,VarInt> > const_op_t;

  static op_t get(const ref_t& _a) { return diag(_a); }
  static op_t get(const const_ref_t& _a) { return diag(_a); }
};

///-----------------------------------------------------------------------------
# endif // !defined(DOXYGEN_SKIP)
//-----------------------------------------------------------------------------

# ifdef DOXYGEN_SKIP
/** Get reference to symmetric matrix defined by triangle of quadratic ge_mat.
    \ingroup vc_blas
    Same for matrix_const_reference_t.

    \todo C++0x: matrix_reference methods sy_upper, etc. (as for GE_Matrix)
*/
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,sy_mat<Upper,N,LD> >
sy_upper(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);
/** Get reference to symmetric matrix defined by triangle of quadratic ge_mat.
    \ingroup vc_blas
    Same for matrix_const_reference_t.
*/
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,sy_mat<Lower,N,LD> >
sy_lower(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);

# else

# define DEF_GE_SY(ref,func,uplo) \
template <typename T,TransposeFlag TRANS,int M,int N,int LD>\
ref<T,sy_mat<uplo,N,LD> > func(const ref<T,ge_mat<TRANS,M,N,LD> >& _a) {\
  assert(_a.m()==_a.n());\
  return ref<T,sy_mat<uplo,N,LD> >(sy_mat<uplo,N,LD>(_a.m(),_a.ld()),_a.data()); \
}
DEF_GE_SY(matrix_reference_t,sy_upper,Upper)
DEF_GE_SY(matrix_reference_t,sy_lower,Lower)
DEF_GE_SY(matrix_const_reference_t,sy_upper,Upper)
DEF_GE_SY(matrix_const_reference_t,sy_lower,Lower)

# undef _GE_SY
#endif

# ifdef DOXYGEN_SKIP
/** Get reference to upper triangle of quadratic ge_mat.
    \ingroup vc_blas
    Same for matrix_const_reference_t.
*/
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,tr_mat<Upper,TRANS,NoU,N,LD> >
tr_upper(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);
/// Same as tr_upper() with UnitDiag \ingroup vc_blas
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,tr_mat<Upper,TRANS,UnitDiag,N,LD> >
tr_upper_u(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);
/** Get reference to lower triangle of quadratic ge_mat.
    \ingroup vc_blas
    Same for matrix_const_reference_t.
*/
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,tr_mat<Lower,TRANS,NoU,N,LD> >
tr_lower(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);
/// Same as tr_upper() with UnitDiag \ingroup vc_blas
template <typename T,TransposeFlag TRANS,int M,int N,int LD>
matrix_reference_t<T,tr_mat<Lower,TRANS,UnitDiag,N,LD> >
tr_lower_u(const matrix_reference_t<T,ge_mat<TRANS,M,N,LD> >& _a);
# else

# define DEF_GE_TR(ref,func,uplo,diag)                          \
template <typename T,TransposeFlag TRANS,int M,int N,int LD>\
ref<T,tr_mat<uplo,TRANS,diag,N,LD> > func(const ref<T,ge_mat<TRANS,M,N,LD> >& _a) { \
  assert(_a.m()==_a.n());\
  return ref<T,tr_mat<uplo,TRANS,diag,N,LD> >\
  (tr_mat<uplo,TRANS,diag,N,LD>(_a.m(),_a.ld()),_a.data());     \
}

DEF_GE_TR(matrix_reference_t,tr_upper,Upper,NoU)
DEF_GE_TR(matrix_reference_t,tr_upper_u,Upper,UnitDiag)
DEF_GE_TR(matrix_reference_t,tr_lower,Lower,NoU)
DEF_GE_TR(matrix_reference_t,tr_lower_u,Lower,UnitDiag)
DEF_GE_TR(matrix_const_reference_t,tr_upper,Upper,NoU)
DEF_GE_TR(matrix_const_reference_t,tr_upper_u,Upper,UnitDiag)
DEF_GE_TR(matrix_const_reference_t,tr_lower,Lower,NoU)
DEF_GE_TR(matrix_const_reference_t,tr_lower_u,Lower,UnitDiag)


# undef DEG_GE_TR

#endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

/// print expression: print in "simplest" format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_print,vector_const_reference_t<T,A> >
print(const vector_const_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_print,vector_const_reference_t<T,A> >(_a);
}
/// pretty print expression \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,A> >
pp(const vector_const_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,A> >(_a);
}
/// print in matlab format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_pmatlab,vector_const_reference_t<T,A> >
pmatlab(const vector_const_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_pmatlab,vector_const_reference_t<T,A> >(_a);
}


/// print expression: print in "simplest" format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_print,vector_const_reference_t<T,A> >
print(const vector_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_print,vector_const_reference_t<T,A> >(_a);
}
/// pretty print expression \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,A> >
pp(const vector_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,A> >(_a);
}
/// print in matlab format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_vector,Expr_pmatlab,vector_const_reference_t<T,A> >
pmatlab(const vector_reference_t<T,A>& _a) {
  return Expression<Type_vector,Expr_pmatlab,vector_const_reference_t<T,A> >(_a);
}

/// print expression: print in "simplest" format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_print,matrix_const_reference_t<T,A> >
print(const matrix_const_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_print,matrix_const_reference_t<T,A> >(_a);
}
/// pretty print expression \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >
pp(const matrix_const_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >(_a);
}
/// print in matlab format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_pmatlab,matrix_const_reference_t<T,A> >
pmatlab(const matrix_const_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_pmatlab,matrix_const_reference_t<T,A> >(_a);
}
/// print all stored data values \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_pdata,matrix_const_reference_t<T,A> >
pdata(const matrix_const_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_pdata,matrix_const_reference_t<T,A> >(_a);
}

/// print expression: print in "simplest" format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_print,matrix_const_reference_t<T,A> >
print(const matrix_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_print,matrix_const_reference_t<T,A> >(_a);
}
/// pretty print expression \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >
pp(const matrix_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >(_a);
}
/// print in matlab format \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_pmatlab,matrix_const_reference_t<T,A> >
pmatlab(const matrix_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_pmatlab,matrix_const_reference_t<T,A> >(_a);
}
/// print all stored data values \ingroup vc_blas
template <typename T,typename A>
Expression<Type_matrix,Expr_pdata,matrix_const_reference_t<T,A> >
pdata(const matrix_reference_t<T,A>& _a) {
  return Expression<Type_matrix,Expr_pdata,matrix_const_reference_t<T,A> >(_a);
}

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
