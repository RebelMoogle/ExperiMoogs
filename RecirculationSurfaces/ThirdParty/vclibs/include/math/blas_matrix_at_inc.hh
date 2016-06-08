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

    \arg Addressing elements in matrices, see
    http://www.netlib.org/lapack/lug/node121.html

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

/// get element _x(_i) \ingroup vc_blas
template <typename T,typename V>
T& at(const vector_reference_t<T,V>& _x,int _i);
/// get element _x(_i) \ingroup vc_blas
template <typename T,typename V>
T at(const vector_const_reference_t<T,V>& _x,int _i);

/// get element _a(_i,_j) \ingroup vc_blas
template <typename T,typename M>
T& at(const matrix_reference_t<T,M>& _a,int _i,int _j);
/// get element _a(_i,_j) \ingroup vc_blas
template <typename T,typename M>
T at(const matrix_const_reference_t<T,M>& _a,int _i,int _j);

#endif

//=============================================================================

//
// get element of a vec
//

template <typename T,int N,int INC>
T& _at(const vector_reference_t<T,vec<N,INC> >& _x,int _i) {
  return _x.data()[_i*_x.inc()]; // "assertion-free"
}

template <typename T,int N,int INC>
const T& _at(const vector_const_reference_t<T,vec<N,INC> >& _x,int _i) {
  return _x.data()[_i*_x.inc()]; // "assertion-free"
}

template <typename T,int N,int INC>
T& at(const vector_reference_t<T,vec<N,INC> >& _x,int _i) {
  assert(0<=_i && _i<_x.n());
  return _at(_x,_i);
}

template <typename T,int N,int INC>
const T& at(const vector_const_reference_t<T,vec<N,INC> >& _x,int _i) {
  assert(0<=_i && _i<_x.n());
  return _at(_x,_i);
}

//-----------------------------------------------------------------------------

//
// get element of a ge_mat
//

  // NoT

template <typename T,int M,int N,int LD>
T& _at(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i,int _j) {
  return _a.data()[_j*_a.matrix().ld()+_i]; // "assertion-free"
}

template <typename T,int M,int N,int LD>
const T& _at(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i,int _j) {
  return _a.data()[_j*_a.matrix().ld()+_i]; // "assertion-free"
}


template <typename T,int M,int N,int LD>
T& at(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.m());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,int M,int N,int LD>
const T& at(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.m());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

  // Transpose

template <typename T,int M,int N,int LD>
T& _at(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i,int _j) {
  return  _a.data()[_i*_a.matrix().ld()+_j]; // "assertion-free"
}

template <typename T,int M,int N,int LD>
const T& _at(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i,int _j) {
  return  _a.data()[_i*_a.matrix().ld()+_j]; // "assertion-free"
}

template <typename T,int M,int N,int LD>
T& at(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.m());
  return _at(_a,_i,_j);
}

template <typename T,int M,int N,int LD>
const T& at(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.m());
  return _at(_a,_i,_j);
}


//-----------------------------------------------------------------------------

//
// get element of a gb_mat
//

template <typename T,int M,int N,int KL,int KU,int LD>
T& _at(const matrix_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >& _a,int _i,int _j) {
  int i=_a.ku()+_i-_j, j=_j;
  assert(0<=i && i<=_a.kl()+_a.ku());
  return _a.data()[j*_a.matrix().lb()+i];
}
template <typename T,int M,int N,int KL,int KU,int LD>
T& _at(const matrix_reference_t<T,gb_mat<Transpose,M,N,KL,KU,LD> >& _a,int _j,int _i) {
  int i=_a.ku()+_i-_j, j=_j;
  assert(0<=i && i<=_a.kl()+_a.ku());
  return _a.data()[j*_a.matrix().lb()+i];
}

template <typename T,int M,int N,int KL,int KU,int LD>
T _at(const matrix_const_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >& _a,int _i,int _j) {
  assert(_j-_a.ku()<=_i);
  assert(_i<=_j+_a.kl());
  int i=_a.ku()+_i-_j, j=_j;
  if (i<0 || i>_a.ku()+_a.kl()) return T(0);
  return _a.data()[j*_a.matrix().lb()+i];
}
template <typename T,int M,int N,int KL,int KU,int LD>
T _at(const matrix_const_reference_t<T,gb_mat<Transpose,M,N,KL,KU,LD> >& _a,int _j,int _i) {
  assert(_j-_a.ku()<=_i);
  assert(_i<=_j+_a.kl());
  int i=_a.ku()+_i-_j, j=_j;
  if (i<0 || i>_a.ku()+_a.kl()) return T(0);
  return _a.data()[j*_a.matrix().lb()+i];
}

template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
T& at(const matrix_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.m());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
T at(const matrix_const_reference_t<T,gb_mat<TRANS,M,N,KL,KU,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.m());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//
// get element of a sy_mat
//

template <typename T,int N,int LD>
T& _at(const matrix_reference_t<T,sy_mat<Upper,N,LD> >& _a,int _i,int _j) {
  assert(_i<=_j); // think it's better to forbid this
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,int N,int LD>
const T& _at(const matrix_const_reference_t<T,sy_mat<Upper,N,LD> >& _a,int _i,int _j) {
  if (_i>_j) std::swap(_i,_j);
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,int N,int LD>
T& _at(const matrix_reference_t<T,sy_mat<Lower,N,LD> >& _a,int _i,int _j) {
  assert(_j<=_i); // think it's better to forbid this
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,int N,int LD>
const T& _at(const matrix_const_reference_t<T,sy_mat<Lower,N,LD> >& _a,int _i,int _j) {
  if (_i<_j) std::swap(_i,_j);
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,UpperLowerFlag UPLO,int N,int LD>
T& at(const matrix_reference_t<T,sy_mat<UPLO,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,UpperLowerFlag UPLO,int N,int LD>
const T& at(const matrix_const_reference_t<T,sy_mat<UPLO,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//
// get element of a sb_mat
//

template <typename T,int N,int K,int LD>
T& _at(const matrix_reference_t<T,sb_mat<Upper,N,K,LD> >& _a,int _i,int _j) {
  assert(_i<=_j); // think it's better to forbid this
  int i=_a.k()+_i-_j, j=_j;
  assert(0<=i && i<=_a.k());
  return _a.data()[j*_a.matrix().lb()+i];
}
template <typename T,int N,int K,int LD>
T& _at(const matrix_reference_t<T,sb_mat<Lower,N,K,LD> >& _a,int _i,int _j) {
  assert(_j<=_i); // think it's better to forbid this
  int i=_i-_j, j=_j;
  assert(0<=i && i<=_a.k());
  return _a.data()[j*_a.matrix().lb()+i];
}

template <typename T,int N,int K,int LD>
T _at(const matrix_const_reference_t<T,sb_mat<Upper,N,K,LD> >& _a,int _i,int _j) {
  if (_i>_j) std::swap(_i,_j);
  assert(_i>=_j-_a.k());
  int i=_a.k()+_i-_j, j=_j;
  if (i<0 || i>_a.k()) return T(0);
  return _a.data()[j*_a.matrix().lb()+i];
}
template <typename T,int N,int K,int LD>
T _at(const matrix_const_reference_t<T,sb_mat<Lower,N,K,LD> >& _a,int _i,int _j) {
  if (_i<_j) std::swap(_i,_j);
  asset(_i<=_j+_a.k());
  int i=_i-_j, j=_j;
  if (i<0 || i>_a.k()) return T(0);
  return _a.data()[j*_a.matrix().lb()+i];
}

template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
T& at(const matrix_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
const T& at(const matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//
// get element of a sp_mat
//

template <typename T,int N>
T& _at(const matrix_reference_t<T,sp_mat<Upper,N> >& _a,int _i,int _j) {
  assert(_i<=_j); // think it's better to forbid this
  return _a.data()[_j*(_j+1)/2+_i];
}

template <typename T,int N>
const T& _at(const matrix_const_reference_t<T,sp_mat<Upper,N> >& _a,int _i,int _j) {
  if (_i>_j) std::swap(_i,_j);
  return _a.data()[_j*(_j+1)/2+_i];
}

template <typename T,int N>
T& _at(const matrix_reference_t<T,sp_mat<Lower,N> >& _a,int _i,int _j) {
  assert(_j<=_i); // think it's better to forbid this
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}

template <typename T,int N>
const T& _at(const matrix_const_reference_t<T,sp_mat<Lower,N> >& _a,int _i,int _j) {
  if (_i<_j) std::swap(_i,_j);
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}

template <typename T,UpperLowerFlag UPLO,int N>
T& at(const matrix_reference_t<T,sp_mat<UPLO,N> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

  template <typename T,UpperLowerFlag UPLO,int N>
const T& at(const matrix_const_reference_t<T,sp_mat<UPLO,N> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//
// get element of a tr_mat
//

template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tr_mat<Upper,NoT,DIAG,N,LD> >& _a,int _i,int _j) {
  assert(DIAG==NoU ? (_i<=_j) : _i<_j);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tr_mat<Lower,NoT,DIAG,N,LD> >& _a,int _i,int _j) {
  assert(DIAG==NoU ? (_j<=_i) : _j<_i);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tr_mat<Upper,Transpose,DIAG,N,LD> >& _a,int _j,int _i) {
  assert(DIAG==NoU ? (_i<=_j) : _i<_j);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tr_mat<Lower,Transpose,DIAG,N,LD> >& _a,int _j,int _i) {
  assert(DIAG==NoU ? (_j<=_i) : _j<_i);
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tr_mat<Upper,NoT,DIAG,N,LD> >& _a,int _i,int _j) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tr_mat<Lower,NoT,DIAG,N,LD> >& _a,int _i,int _j) {
  if (_j>_i) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tr_mat<Upper,Transpose,DIAG,N,LD> >& _a,int _j,int _i) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*_a.matrix().ld()+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tr_mat<Lower,Transpose,DIAG,N,LD> >& _a,int _j,int _i) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*_a.matrix().ld()+_i];
}

template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD>
T& at(const matrix_reference_t<T,tr_mat<UPLO,TRANS,DIAG,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD>
T at(const matrix_const_reference_t<T,tr_mat<UPLO,TRANS,DIAG,N,LD> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//
// get element of a tb_mat
//

//-----------------------------------------------------------------------------

//
// get element of a tp_mat
//

template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tp_mat<Upper,NoT,DIAG,N> >& _a,int _i,int _j) {
  assert(DIAG==NoU ? (_i<=_j) : _i<_j);
  return _a.data()[_j*(_j+1)/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tp_mat<Lower,NoT,DIAG,N> >& _a,int _i,int _j) {
  assert(DIAG==NoU ? (_j<=_i) : _j<_i);
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tp_mat<Upper,Transpose,DIAG,N> >& _a,int _j,int _i) {
  assert(DIAG==NoU ? (_i<=_j) : _i<_j);
  return _a.data()[_j*(_j+1)/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T& _at(const matrix_reference_t<T,tp_mat<Lower,Transpose,DIAG,N> >& _a,int _j,int _i) {
  assert(DIAG==NoU ? (_j<=_i) : _j<_i);
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}

template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tp_mat<Upper,NoT,DIAG,N> >& _a,int _i,int _j) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*(_j+1)/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tp_mat<Lower,NoT,DIAG,N> >& _a,int _i,int _j) {
  if (_j>_i) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tp_mat<Upper,Transpose,DIAG,N> >& _a,int _j,int _i) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[_j*(_j+1)/2+_i];
}
template <typename T,DiagonalFlag DIAG,int N,int LD>
T _at(const matrix_const_reference_t<T,tp_mat<Lower,Transpose,DIAG,N> >& _a,int _j,int _i) {
  if (_i>_j) return T(0);
  if (DIAG==UnitDiag && _i==_j) return T(1);
  return _a.data()[(2*_a.n()-1-_j)*_j/2+_i];
}

template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD>
T& at(const matrix_reference_t<T,tp_mat<UPLO,TRANS,DIAG,N> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int N,int LD>
T at(const matrix_const_reference_t<T,tp_mat<UPLO,TRANS,DIAG,N> >& _a,int _i,int _j) {
  assert(0<=_i && _i<_a.n());
  assert(0<=_j && _j<_a.n());
  return _at(_a,_i,_j);
}

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC