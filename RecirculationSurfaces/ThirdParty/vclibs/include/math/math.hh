//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision: 6 $
//
//=============================================================================

//
// for documentation only
//

#ifndef VC_MATH_HH
#define VC_MATH_HH

#include <cmath>
#include <cassert>
#include <limits>

namespace VC {
/// math related [\ref vc_math] \ingroup vc_math
namespace math {

  /** \defgroup vc_math Mathematics related functionality
      \sa VC::math
  */

/// acos with clamping argument towards +-1 \ingroup vc_math
template <typename T>
inline T acos_c(T _x) {
  assert(fabs(_x)<T(1)+std::numeric_limits<T>::epsilon()*128.0 && "invalid argument");
  if      (_x<T(-1)) _x=T(-1);
  else if (_x>T(+1)) _x=T(+1);
  return ::acos(_x);
}
/// asin with clamping argument towards +-1 \ingroup vc_math
template <typename T>
inline T asin_c(T _x) {
  assert(fabs(_x)<T(1)+std::numeric_limits<T>::epsilon()*128.0 && "invalid argument");
  if      (_x<T(-1)) _x=T(-1);
  else if (_x>T(+1)) _x=T(+1);
  return ::asin(_x);
}


} // namespace math
} // namespace VC

#endif // VC_MATH_HH