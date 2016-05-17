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

    \arg Rules for transposing matrices (matrix_const_reference_t)

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

# ifdef DOXYGEN_SKIP

  /** \brief Yields reference to transposed matrix.
      \ingroup vc_blas


      \arg  Changes the TransposeFlag of a matrix (reference to ge_mat,
      gb_mat, tr_mat, tb_mat, tp_mat; no effect to sy_mat, sb_mat,
      sp_mat).

      \arg \b Important: only the flag is changed. The dimension
      properties m() and n() stay the same! (Matrix-matrix
      multiplication mm() is aware of that.)

      \arg The meaning of this flag depends on the BLAS algorithm that
      the matrix is used with.

      Hints for use:

      \arg Use only with constant references (matrix_const_reference_t).
      (Although trans() is also defined for matrix_reference_t.)

      \arg Use only on arguments to BLAS functions or their wrappers
      (like methods/operators of matrix_reference_t, e.g.,
      matrix_reference_t::mm).

      \sa matrix_const_reference_t, matrix_reference_t
   */
matrix_const_reference_t<>
trans(const matrix_const_reference_t<>& _a);

/** \brief Yields reference to 1xN matrix.
    \ingroup vc_blas
    Transposing a vector (vector_const_reference_t = column vector) yields
    a 1xN ge_mat matrix.
    \sa matrix_const_reference_t<>trans(const matrix_const_reference_t<>& _a)
    \sa matrix_const_reference_t, matrix_reference_t
*/
matrix_const_reference_t<>
trans(const vector_const_reference_t<>& _a);

# else

// vector

template <typename T,int N,int INC>
matrix_const_reference_t<T,ge_mat<NoT,1,N,INC> >
trans(const vector_const_reference_t<T,vec<N,INC> >& _x) {
  return matrix_const_reference_t<T,ge_mat<NoT,1,N,INC> >
    (ge_mat<NoT,1,N,INC>(1,_x.vector().n(),_x.vector().inc()),_x.data());
}

template <typename T,int N,int INC>
matrix_reference_t<T,ge_mat<NoT,1,N,INC> >
trans(const vector_reference_t<T,vec<N,INC> >& _x) {
  return matrix_reference_t<T,ge_mat<NoT,1,N,INC> >
    (ge_mat<NoT,1,N,INC>(1,_x.vector().n(),_x.vector().inc()),_x.data());
}



// GE

template <typename T,int M,int N,int LD>
matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >
trans(const matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >& _a) {
  return matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >
    (ge_mat<Transpose,M,N,LD>(_a.matrix().m(),_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}
template <typename T,int M,int N,int LD>
matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >
trans(const matrix_const_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a) {
  return matrix_const_reference_t<T,ge_mat<NoT,M,N,LD> >
    (ge_mat<NoT,M,N,LD>(_a.matrix().m(),_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}


template <typename T,int M,int N,int LD>
matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >
trans(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _a) {
  return matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >
    (ge_mat<Transpose,M,N,LD>(_a.matrix().m(),_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}
template <typename T,int M,int N,int LD>
matrix_reference_t<T,ge_mat<NoT,M,N,LD> >
trans(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _a) {
  return matrix_reference_t<T,ge_mat<NoT,M,N,LD> >
    (ge_mat<NoT,M,N,LD>(_a.matrix().m(),_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}

// GB

template <typename T,int M,int N,int KL,int KU,int LD>
matrix_const_reference_t<T,gb_mat<Transpose,M,N,KL,KU,LD> >
trans(const matrix_const_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >& _a) {
  return matrix_const_reference_t<T,gb_mat<Transpose,M,N,KL,KU,LD> >
    (_a.matrix().m(),_a.matrix().n(),_a.matrix().kl(),_a.matrix().ku(),
     _a.matrix().ld(),_a.data());
}
template <typename T,int M,int N,int KL,int KU,int LD>
matrix_const_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >
trans(const matrix_const_reference_t<T,gb_mat<Transpose,M,N,KL,KU,LD> >& _a) {
  return matrix_const_reference_t<T,gb_mat<NoT,M,N,KL,KU,LD> >
    (gb_mat<NoT,M,N,KL,KU,LD>
     (_a.matrix().m(),_a.matrix().n(),_a.matrix().kl(),_a.matrix().ku(),
      _a.matrix().ld()),_a.data());
}

// SY,SB,SP -- no effect

template <typename T,UpperLowerFlag UPLO,int N,int LD>
const matrix_const_reference_t<T,sy_mat<UPLO,N,LD> >&
trans(const matrix_const_reference_t<T,sy_mat<UPLO,N,LD> >& _a) {
  return _a;
}
template <typename T,UpperLowerFlag UPLO,int N,int K,int LD>
const matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> >&
  trans(const matrix_const_reference_t<T,sb_mat<UPLO,N,K,LD> >& _a) {
  return _a;
}
template <typename T,UpperLowerFlag UPLO,int N>
const matrix_const_reference_t<T,sp_mat<UPLO,N> >&
trans(const matrix_const_reference_t<T,sp_mat<UPLO,N> >& _a) {
  return _a;
}

// TR

template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int LD>
matrix_const_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >
trans(const matrix_const_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >& _a) {
  return matrix_const_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >
    (tr_mat<UPLO,Transpose,DIAG,N,LD>(_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}
template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int LD>
matrix_const_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >
trans(const matrix_const_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >& _a) {
  return matrix_const_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >
    (tr_mat<UPLO,NoT,DIAG,N,LD>(_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}


template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int LD>
matrix_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >
trans(const matrix_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >& _a) {
  return matrix_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >
    (tr_mat<UPLO,Transpose,DIAG,N,LD>(_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}
template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int LD>
matrix_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >
trans(const matrix_reference_t<T,tr_mat<UPLO,Transpose,DIAG,N,LD> >& _a) {
  return matrix_reference_t<T,tr_mat<UPLO,NoT,DIAG,N,LD> >
    (tr_mat<UPLO,NoT,DIAG,N,LD>(_a.matrix().n(),_a.matrix().ld()),
     _a.data());
}

// TB

template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int K,int LD>
matrix_const_reference_t<T,tb_mat<UPLO,Transpose,DIAG,N,K,LD> >
trans(const matrix_const_reference_t<T,tb_mat<UPLO,NoT,DIAG,N,K,LD> >& _a) {
  return matrix_const_reference_t<T,tb_mat<UPLO,Transpose,DIAG,N,K,LD> >
    (tb_mat<UPLO,Transpose,DIAG,N,K,LD>(_a.matrix().n(),_a.matrix().k(),
                                        _a.matrix().ld()),
     _a.data());
}
template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N,int K,int LD>
matrix_const_reference_t<T,tb_mat<UPLO,NoT,DIAG,N,K,LD> >
trans(const matrix_const_reference_t<T,tb_mat<UPLO,Transpose,DIAG,N,K,LD> >& _a) {
  return matrix_const_reference_t<T,tb_mat<UPLO,NoT,DIAG,N,K,LD> >
    (tb_mat<UPLO,NoT,DIAG,N,K,LD>(_a.matrix().n(),_a.matrix().k(),
                                  _a.matrix().ld()),
     _a.data());
}


// TP

template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N>
matrix_const_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >
trans(const matrix_const_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >& _a) {
  return matrix_const_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >
    (tp_mat<UPLO,Transpose,DIAG,N>(_a.matrix().n()),
     _a.data());
}
template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N>
matrix_const_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >
trans(const matrix_const_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >& _a) {
  return matrix_const_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >
    (tp_mat<UPLO,NoT,DIAG,N>(_a.matrix().n()),
     _a.data());
}


template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N>
matrix_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >
trans(const matrix_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >& _a) {
  return matrix_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >
    (tp_mat<UPLO,Transpose,DIAG,N>(_a.matrix().n()),
     _a.data());
}
template <typename T,UpperLowerFlag UPLO,DiagonalFlag DIAG,int N>
matrix_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >
trans(const matrix_reference_t<T,tp_mat<UPLO,Transpose,DIAG,N> >& _a) {
  return matrix_reference_t<T,tp_mat<UPLO,NoT,DIAG,N> >
    (tp_mat<UPLO,NoT,DIAG,N>(_a.matrix().n()),
     _a.data());
}


# endif // DOXYGEN_SKIP

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC