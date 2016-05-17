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

    \arg Closures for vector/matrix computations (e.g., map operators to BLAS
    calls).
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** \defgroup vc_blas_expr Expression templates for evaluating matrix operations.
    \ingroup vc_blas_details

    Expressions and related operators on closures
    \sa vector_reference_t, vector_const_reference_t,
    matrix_reference_t, matrix_const_reference_t
 */

//-----------------------------------------------------------------------------

#if defined(_MSC_VER)
# pragma warning(push)
# pragma warning(disable:4100 4512)
# undef _BB // msc decides this must be a "constant"!?
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifndef DOXYGEN_SKIP
# include "blas_matrix_expr_generated_inc.hh"
#endif

#pragma GCC diagnostic pop

#if defined(_MSC_VER)
# pragma warning(pop)
#endif

# ifdef DOXYGEN_SKIP
/** \def _VC_BLAS_TRACE_EXPR
    \ingroup vc_blas
    Debugging aid: monitor evaluation of expression templates (compile time).
    \code
    # define _VC_BLAS_TRACE_EXPR 1 // debug: show what's going on
    \endcode
 */
# define _VC_BLAS_TRACE_EXPR
# endif

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

enum Expr_print { EXPR_print=200 };
enum Expr_prettyprint { EXPR_prettyprint=201 };
enum Expr_pmatlab { EXPR_pmatlab=202 };
enum Expr_pdata { EXPR_ptype=203 };

  //
  // UNIFY: print, etc. same for matrix and vector !!
  //

/// ostream << print() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_matrix,Expr_print,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << pp() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_matrix,Expr_prettyprint,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << matlab() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_matrix,Expr_pmatlab,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << matlab() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_matrix,Expr_pdata,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << print() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_vector,Expr_print,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << pp() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_vector,Expr_prettyprint,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};

/// ostream << matlab() \ingroup vc_blas_expr
template <typename A>
struct Expression<Type_vector,Expr_pmatlab,A> {
  Expression(const A& _a) : a(_a) {}
  A a;
};


//-----------------------------------------------------------------------------

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

/** Load zero matrix (ld_zero()) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    A=zeros;
    \endcode
*/
extern Expression<Type_initializer,Expr_zeros> zeros;

/** Load matrix of ones (ld_all()) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    A=ones;
    \endcode
*/
extern Expression<Type_initializer,Expr_ones> ones;

/** Load identity matrix (ld_eye()) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    A=eye;
    \endcode
*/
extern Expression<Type_initializer,Expr_eye> eye; // note: fails for vector assignment

/** Load random values (ld_rand()) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    A=random_values;
    \endcode
*/
extern Expression<Type_initializer,Expr_random_values> random_values;


//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
