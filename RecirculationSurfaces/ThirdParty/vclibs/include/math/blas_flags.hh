//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================


#ifndef __VC_MATH_BLAS_FLAGS_HH
#define __VC_MATH_BLAS_FLAGS_HH

#include "../base/debug.hh"

# ifndef _VC_DBG_BLAS_SCOPE
#  define _VC_DBG_BLAS_SCOPE()
//#  define _VC_DBG_BLAS_SCOPE() VC_DBG_SCOPE("blas") // debug: show what's going on
# endif

# ifdef DOXYGEN_SKIP
/** \def _VC_DBG_BLAS_SCOPE
    \ingroup vc_blas
    Debugging aid: define to track calls to BLAS routines (wrappers).
    \code
    # define _VC_DBG_BLAS_SCOPE VC_DBG_SCOPE("blas") // debug: show what's going on
    \endcode
 */
# endif

namespace VC {
namespace math {
namespace blas {

//=============================================================================

  //
  // C++ wrappers to low-level blas functions
  //

  /** Multiply matrix from which side?
      \ingroup vc_blas
      See symm(), trmm(), trsm().
  */
  enum SideFlag {
    Left='L',   //!< A*B, A is left
    Right='R'   //!< B*A, A is right
  };

  /** Upper or lower triangular?
      \ingroup vc_blas
   */
  enum UpperLowerFlag {
    Upper='U', //!< upper
    Lower='L'  //!< lower
  };

  /** Transpose matrix?
      \ingroup vc_blas
   */
  enum TransposeFlag {
    NoTranspose='N', //!< no
    NoT=NoTranspose, //!< same as NoTranspose
    Transpose='T'    //!< transpose
  };

  /** Have unit diagonal?
      \ingroup vc_blas
   */
  enum DiagonalFlag {
    NoUnitDiag='N', //!< no
    NoU=NoUnitDiag, //!< same as NoUnitDiag
    UnitDiag='U'    //!< assume unit diagonal
  };

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

#endif // __VC_MATH_BLAS_FLAGS_HH