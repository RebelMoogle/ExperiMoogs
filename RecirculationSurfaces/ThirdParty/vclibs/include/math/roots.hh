//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================


#ifndef __VC_MATH_ROOTS_HH
#define __VC_MATH_ROOTS_HH

# include <limits>
# include <cmath>

#include "../base/debug.hh"

namespace VC {
namespace math {

//=============================================================================

/** Compute real roots of quadratic equation _a*x^2+_b*x+_c == 0.
    \ingroup vc_math
    \tparam T argument and return values
    \tparam S used for internal computations (see, e.g., real_roots_2d())
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
*/
template <typename T,typename S>
int real_roots_2(T _a,T _b,T _c,T* _x) {
  S d=_b*_b-S(4)*_a*_c;
  if (d<S(0))
    return 0;

  d=sqrt(d);
  S q=S(-0.5)*(_b+(_b>=T(0) ? S(1) : S(-1))*d);

  if (d<std::numeric_limits<T>::epsilon()) {
    _x[0]=_x[1]=T(q/_a);
    return 1;
  }

  S x1=q/_a, x2=_c/q;

  if (x1<x2) {
    _x[0]=T(x1);
    _x[1]=T(x2);
  }
  else {
    _x[0]=T(x2);
    _x[1]=T(x1);
  }
  return 2;
}

//-----------------------------------------------------------------------------

/** Compute real roots of quadratic equation _a[0]*x^2+_a[1]*x+_a[2] == 0.
    \ingroup vc_math
    \tparam T argument and return values
    \tparam S used for internal computations (see, e.g., real_roots_2d())
    \param[in] _a 3 polynomial coefficients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
*/
template <typename T,typename S>
int real_roots_2(const T* _a,T* _x) {
  return real_roots_2<T,S>(_a[0],_a[1],_a[2],_x);
}

/** Compute real roots of quadratic equation _a[0]*x^2+_a[1]*x+_a[2] == 0.
    \ingroup vc_math
    \param[in] _a 3 polynomial coefficients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
*/
template <typename T>
int real_roots_2(const T* _a,T* _x) {
  return real_roots_2<T,T>(_a[0],_a[1],_a[2],_x);
}

/** Compute real roots of quadratic equation _a*x^2+_b*x+_c == 0.
    \ingroup vc_math
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
*/
template <typename T>
int real_roots_2(T _a,T _b,T _c,T* _x) {
  return real_roots_2<T,T>(_a,_b,_c,_x);
}

/** Compute real roots of quadratic equation always in double precision.
    \ingroup vc_math
    \param[in] _a 3 polynomial coefficients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
    \sa real_roots_2()
 */
template <typename T>
int real_roots_2d(const T* _a,T* _x) {
  return real_roots_2<T,double>(_a[0],_a[1],_a[2],_x);
}

/** Compute real roots of quadratic equation always in double precision.
    \ingroup vc_math
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (0,1,2)
    \sa real_roots_2()
 */
template <typename T>
int real_roots_2d(T _a,T _b,T _c,T* _x) {
  return real_roots_2<T,double>(_a,_b,_c,_x);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Compute real roots of cubic equation x^3+_a*x^2+_b*_x+_c == 0.
    \ingroup vc_math
    \tparam T argument and return values
    \tparam S used for internal computations (see, e.g., real_roots_3d())
    Note: assumes that coefficient ot x^3 equals 1.
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
 */
template <typename T,typename S>
int real_roots_3(T _a,T _b,T _c,T* _x) {
  // from Numerical Receipes
  S Q=(_a*_a-S(3)*_b)/S(9);
  S R=(S(2)*_a*_a*_a-S(9)*_a*_b+S(27)*_c)/S(54);

  if (R*R<Q*Q*Q) {
    S theta=R/sqrt(Q*Q*Q);
    if (theta<S(-1)) theta=-1; // clamp
    if (theta>S(+1)) theta=+1;
    theta=acos(theta);
    Q=sqrt(Q)*S(-2);
    S x0=Q*cos(theta/S(3))-_a/S(3);
    S x1=Q*cos((theta+S(2)*S(M_PI))/S(3))-_a/S(3);
    S x2=Q*cos((theta-S(2)*S(M_PI))/S(3))-_a/S(3);

    if (x0>x1) std::swap(x0,x1);
    if (x0>x2) std::swap(x0,x2);
    if (x1>x2) std::swap(x1,x2);

    _x[0]=T(x0); _x[1]=T(x1); _x[2]=T(x2);

    return 3;
  }

  S A=-((R>S(0)) ? S(+1) : S(-1))*(fabs(R)+sqrt(R*R-Q*Q*Q));
  A=cbrt(A); // pow(A,S(1)/S(3));

  S B=(fabs(A)>std::numeric_limits<S>::epsilon()) ? Q/A : S(0);
  S x=(A+B)-_a/S(3);

  _x[0]=_x[1]=_x[2]=T(x);

  return 1;
}

//-----------------------------------------------------------------------------

/** Compute real roots of cubic equation x^3+_a*x^2+_b*_x+_c == 0.
    \ingroup vc_math
    Note: assumes that coefficient ot x^3 equals 1.
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
 */
template <typename T>
int real_roots_3(T _a,T _b,T _c,T* _x) {
  return real_roots_3<T,T>(_a,_b,_c,_x);
}

/** Compute real roots of cubic equation always in double precision.
    \ingroup vc_math
    Note: assumes that coefficient ot x^3 equals 1.
    \param _a coeffcients to x^2
    \param _b coeffcients to x
    \param _c coeffcients to constant 1
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
    \sa real_roots_3()
 */
template <typename T>
int real_roots_3d(T _a,T _b,T _c,T* _x) {
  return real_roots_3<T,double>(_a,_b,_c,_x);
}

/** Compute real roots of cubic equation _a[0]*x^3+_a[1]*x^2+_a[2]*_x+_a[3] == 0.
    \ingroup vc_math
    \tparam T argument and return values
    \tparam S used for internal computations (see, e.g., real_roots_3d())
    \param[in] _a coeffcients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
 */
template <typename T,typename S>
int real_roots_3(const T* _a,T* _x) {
  return real_roots_3<T,S>(_a[1]/_a[0],_a[2]/_a[0],_a[3]/_a[0],_x);
}

/** Compute real roots of cubic equation _a[0]*x^3+_a[1]*x^2+_a[2]*_x+_a[3] == 0.
    \ingroup vc_math
    \param[in] _a coeffcients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
 */
template <typename T>
int real_roots_3(const T* _a,T* _x) {
  return real_roots_3<T,T>(_a[1]/_a[0],_a[2]/_a[0],_a[3]/_a[0],_x);
}

/** Compute real roots of cubic equation always in double precision.
    \ingroup vc_math
    \param[in] _a coeffcients
    \param[out] _x stores (0,1,2) solutions
    \return number of real roots (1 or 3)
    \sa real_roots_3()
 */
template <typename T>
int real_roots_3d(const T* _a,T* _x) {
  return real_roots_3<T,double>(_a[1]/_a[0],_a[2]/_a[0],_a[3]/_a[0],_x);
}

//-----------------------------------------------------------------------------

//=============================================================================



//=============================================================================
} // namespace math
} // namespace VC

#endif // __VC_MATH_ROOTS_HH
