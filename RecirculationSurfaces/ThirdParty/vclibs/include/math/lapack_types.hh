#ifndef VC_MATH_LAPACK_TYPES_HH
#define VC_MATH_LAPACK_TYPES_HH

/** \file lapack_types.hh Math: define types used by LAPACK
 */

namespace VC {
namespace math {
namespace lapack {

typedef int INTEGER;   //!< integer type used by LAPACK
typedef int LOGICAL;   //!< bool type used by LAPACK
typedef void (*EXTERNAL_PROCEDURE)(); //!< function pointer used by LAPACK (dummy type)

} // namespace lapack
} // namespace math
} // namespace VC


#endif // VC_MATH_LAPACK_TYPES_HH
