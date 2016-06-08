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

    \internal

    \arg Defines expressions which evaluate to compile time error if
    certain types don't match. The idea is the have expression names
    describe the error.
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

#include <type_traits>

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** \todo error handling: provide high level predicates
    (`is_matrix_type<>`, `is_triangular<>`) for static_assrt();
    provide mixed static/ dynamic predicates for checking dimensions
    (`static_assert` for fixed dimensions, otherwise run-time
    `assert`)
 */


#  ifndef DOXYGEN_SKIP

//-----------------------------------------------------------------------------

// run time assertions

template <typename X> struct __assert_no_unit_diag {};

template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int LD>
struct __assert_no_unit_diag<tr_mat<UPLO,TRANS,DIAG,N,LD> > {
  __assert_no_unit_diag() { assert(DIAG!=UnitDiag && !"cannot modify diagonal"); }
};

template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int K,int LD>
struct __assert_no_unit_diag<tb_mat<UPLO,TRANS,DIAG,N,K,LD> > {
  __assert_no_unit_diag() { assert(DIAG!=UnitDiag && !"cannot modify diagonal"); }
};

template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N>
struct __assert_no_unit_diag<tp_mat<UPLO,TRANS,DIAG,N> > {
  __assert_no_unit_diag() { assert(DIAG!=UnitDiag && !"cannot modify diagonal"); }
};

//-----------------------------------------------------------------------------

#  endif

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
