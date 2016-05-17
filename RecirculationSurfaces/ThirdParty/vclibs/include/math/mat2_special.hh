//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#include "Mat.hh"

#ifndef __VC_MATH_MAT2_HH__
#define __VC_MATH_MAT2_HH__

namespace VC {
namespace math {

//
// specializations for 2x2 matrices
//

# ifndef VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/// specialization \ingroup vc_math_lam_special
template <typename T> struct inverse_t<T,2,2> {

  static bool get(Mat<T,2,2>& _a) {
    T d=det_t<T,2,2>::get(_a);
    if (d!=T(0)) {
      T* raw=_a.data();
      T a0=raw[0];
      raw[0]= raw[3]/d;
      raw[1]=-raw[1]/d;
      raw[2]=-raw[2]/d;
      raw[3]= a0/d;
      return true;
    }
    return false;
  }
};

/// specialization \ingroup vc_math_lam_special
template <typename T>
struct det_t<T,2,2> {
  static T get(const Mat<T,2,2>& _a) {
    const T* raw=_a.data();
    return raw[0]*raw[3]-raw[1]*raw[2];
  }
};

/// specialization \ingroup vc_math_lam_special
template <typename T>
struct sy_eig_t<T,2,2> {
  static typename Mat<T,2,2>::col_t get(const Mat<T,2,2>& _a) {
    typename Mat<T,2,2>::col_t w;
    eig2x2(_a.data(),w.data());
    return w;
  }
  static typename Mat<T,2,2>::col_t get(const Mat<T,2,2>& _a,Mat<T,2,2>& _v) {
    typename Mat<T,2,2>::col_t w;
    eig2x2(_a.data(),w.data(),_v.data());
    return w;
  }
};

//-----------------------------------------------------------------------------
# endif // VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/** \defgroup vc_math_lam_2x2 Functions of 2x2 matrices
    \ingroup vc_math_lam
 */

/** Get rotation matrix.
    \ingroup vc_math_lam_2x2
    \param _rad angle (radiants)
    \return rotation matrix `R` such that `R*X` is a counter-clockwise rotation
 */
template <typename T>
Mat<T,2,2> rot2(T _rad) {
  T c=cos(_rad), s=sin(_rad);
  Mat<T,2,2> r;
  r.elts()=VecN<T,4>(c,s,-s,c);
  return r;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace math
} // namespace VC

#endif // __VC_MATH_MAT2_HH__
