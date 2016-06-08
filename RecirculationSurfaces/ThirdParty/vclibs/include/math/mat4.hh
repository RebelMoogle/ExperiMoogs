//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_MATH_MAT4_HH
#define VC_MATH_MAT4_HH

#include "../base/mat4.hh"

//=============================================================================


namespace VC {
namespace math {

//-----------------------------------------------------------------------------

# ifdef DOXYGEN_SKIP
/** Read/write MATLAB level 4 MAT file format
    \ingroup vc_math

    The `namespace VC::math::mat4` is an alias for VC::base::mat4io.
    \see vc_mat4
 */
namespace mat4 {}
# else

namespace mat4 = VC::base::mat4io;

# endif

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace math
} // namespace VC

# endif // VC_MATH_MAT4_HH
