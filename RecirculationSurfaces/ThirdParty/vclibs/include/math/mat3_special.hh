//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#include "Mat.hh"

#ifndef __VC_MATH_MAT3_HH__
#define __VC_MATH_MAT3_HH__

namespace VC {
namespace math {

//
// specializations for 3x3 matrices
//

# ifndef VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/// specialization \ingroup vc_math_lam_special
template <typename T>
struct inverse_t<T,3,3> {
  static bool get(Mat<T,3,3>& _a) {
    typename Mat<T,3,3>::row_t a=_a.row(0), b=_a.row(1), c=_a.row(2);
    typename Mat<T,3,3>::row_t ab=a%b, bc=b%c, ca=c%a;

    T d=(a|bc);
    if (d==T(0))
      return false;

    _a.col(0)=bc/d;
    _a.col(1)=ca/d;
    _a.col(2)=ab/d;

    return true;
  }
};

//-----------------------------------------------------------------------------

/// specialization \ingroup vc_math_lam_special
template <typename T>
struct det_t<T,3,3> {
  static T get(const Mat<T,3,3>& _a) {
    return det(_a.col(0),_a.col(1),_a.col(2));
  }
};

//-----------------------------------------------------------------------------

/** specialization
    \ingroup vc_math_lam_special

    \bug Check this in presence of aggressive math optimization of new
    gcc versions (>=4.8.1): behavior of optimized code differs to that
    of non-optimized code!
 */
template <typename T>
struct sy_eig_t<T,3,3> {

  /// compute eigenvalues as roots of characteristic polynomial
  static typename Mat<T,3,3>::col_t get(const Mat<T,3,3>& _a) {

    const T* pa=_a.data();

    // diagonal matrix
    if (fabs(pa[3])+fabs(pa[6])+fabs(pa[7])==T(0)) {
      return typename Mat<T,3,3>::col_t(pa[0],pa[4],pa[8]).sort();
    }

    double m[6]={pa[0], pa[3],pa[4], pa[6],pa[7],pa[8]};

    double a= m[0]+m[2]+m[5];
    double b= m[1]*m[1]  +m[3]*m[3] + m[4]*m[4] -
              m[0]*m[2] - m[0]*m[5] - m[2]*m[5];
    double c= -m[5]*m[1]*m[1] + 2.0*m[1]*m[3]*m[4] -
               m[2]*m[3]*m[3] -     m[0]*m[4]*m[4] + m[0]*m[2]*m[5];

    double x[3];
    int n=real_roots_3(-a,-b,-c, x);
    assert(n==3);
    use_nowarn(n);

    // see also http://en.wikipedia.org/wiki/Eigenvalue_algorithm

    return typename Mat<T,3,3>::col_t(x); // roots are sorted
  }

  /** Compute eigenvalues as roots of characteristic polynomial.
      Compute eigenvectors from Cayley-Hamilton theorem. Expect a
      slightly less accurate result compared to lapack::syev().
      In addition, signs and orientation may differ!
   */
  static typename Mat<T,3,3>::col_t get(const Mat<T,3,3>& _a,Mat<T,3,3>& _v) {

    const T* pa=_a.data();

    // diagonal matrix
    if (fabs(pa[3])+fabs(pa[6])+fabs(pa[7])==T(0)) {
      _v.load_eye();
      return typename Mat<T,3,3>::col_t(pa[0],pa[4],pa[8]).sort();
    }

    double m[6]={pa[0], pa[3],pa[4], pa[6],pa[7],pa[8]};

    double a= m[0]+m[2]+m[5];
    double b= m[1]*m[1]  +m[3]*m[3] + m[4]*m[4] -
              m[0]*m[2] - m[0]*m[5] - m[2]*m[5];
    double c= -m[5]*m[1]*m[1] + 2.0*m[1]*m[3]*m[4] -
               m[2]*m[3]*m[3] -     m[0]*m[4]*m[4] + m[0]*m[2]*m[5];

    double x[3];
    int n=real_roots_3(-a,-b,-c, x);
    assert(n==3);
    use_nowarn(n);

    // see also http://en.wikipedia.org/wiki/Eigenvalue_algorithm

    // Cayley-Hamilton

    Mat<double,3,3> a0=_a; a0-=diag(Mat<double,3,3>::col_t(x[0],x[0],x[0]));
    Mat<double,3,3> a1=_a; a1-=diag(Mat<double,3,3>::col_t(x[1],x[1],x[1]));
    Mat<double,3,3> a2=_a; a2-=diag(Mat<double,3,3>::col_t(x[2],x[2],x[2]));

    _v.col(0)=(a1*a2).col(0); // should work for multiple eigenvalues
    _v.col(1)=(a0*a2).col(1); //  (just too much work)
    _v.col(2)=(a0*a1).col(2);

    for (int j=0;j<3;++j)
      _v.col(j).normalize();

    return typename Mat<T,3,3>::col_t(x); // roots are sorted
  }
};

//-----------------------------------------------------------------------------
# endif // VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/** \defgroup vc_math_lam_3x3 Functions of 3x3 matrices
    \ingroup vc_math_lam
 */

/** Get rotation matrix that rotates about x-axis.
    \ingroup vc_math_lam_3x3
    \param _rad angle (radiants)
    \return rotation matrix `R` such that `R*X` is a counter-clockwise rotation
    about the x-axis
 */
template <typename T>
Mat<T,3,3> rot3x(T _rad) {
  T c=cos(_rad), s=sin(_rad);
  Mat<T,3,3> r;
  r.elts()=VecN<T,9>(1,0,0, 0,c,s, 0,-s,c);
  return r;
}
/** Get rotation matrix that rotates about y-axis.
    \ingroup vc_math_lam_3x3
    \param _rad angle (radiants)
    \return rotation matrix `R` such that `R*X` is a counter-clockwise rotation
    about the y-axis
 */
template <typename T>
Mat<T,3,3> rot3y(T _rad) {
  T c=cos(_rad), s=sin(_rad);
  Mat<T,3,3> r;
  r.elts()=VecN<T,9>(c,0,-s, 0,1,0, s,0,c);
  return r;
}
/** Get rotation matrix that rotates about z-axis.
    \ingroup vc_math_lam_3x3
    \param _rad angle (radiants)
    \return rotation matrix `R` such that `R*X` is a counter-clockwise rotation
    about the z-axis
 */
template <typename T>
Mat<T,3,3> rot3z(T _rad) {
  T c=cos(_rad), s=sin(_rad);
  Mat<T,3,3> r;
  r.elts()=VecN<T,9>(c,s,0, -s,c,0, 0,0,1);
  return r;
}

/** Get rotation matrix that rotates around `_x`.
    \ingroup vc_math_lam_3x3
    \param _rad angle (radiants)
    \param _x
    \return rotation matrix `R` such that `R*X` is a counter-clockwise rotation
    around the vector `_x`.
 */
template <typename T>
Mat<T,3,3> rot3(T _rad,const VecN<T,3>& _x) {
  T len=norm(_x), c=cos(_rad), s=sin(_rad);
  Mat<T,3,3> r;
  if (len!=T(0)) {
    T x=_x[0]/len, y=_x[1]/len, z=_x[2]/len, c1=T(1)-c;;
    r.elts()=VecN<T,9>(x*x*c1+c,   y*x*c1+z*s, x*z*c1-y*s,
                       x*y*c1-z*s, y*y*c1+c,   y*z*c1+x*s,
                       x*z*c1+y*s, y*z*c1-x*s, z*z*c1*c    );
  }
  else
    r.load_eye();

  return r;
}

/** Get rotation that rotates direction `_x` into `_y`.
    \param _x direction
    \param _y direction
    \return rotation matrix
 */
template <typename T>
Mat<T,3,3>& rot3(const VecN<T,3>& _x,const VecN<T,3>& _y) {
  VecN<T,3> x=_x/_x.norm(), y=_y/_y.norm();
  VecN<T,3> z=x%y;

  T angle=angle_u(x,y);
  return rot3(angle,z);
}

/** Get Euler angles from rotation matrix `Rz*Ry*Rx`.
    \ingroup vc_math_lam_3x3
    \param _r rotation matrix (no check!), assume order `_r=Rz*Ry*Rx`.
    \return Euler angles
 */
template <typename T>
VecN<T,3> euler_angles(const Mat<T,3,3>& _r) {
  const T* r=_r.data();
  T x,y,z;
  if (fabs(fabs(r[2])-1)>std::numeric_limits<T>::epsilon()) {
    T y1=-asin(r[2]);
    //T y2=M_PI-phi1;
    T x1=atan2(r[5]/cos(y1),r[8]/cos(y1));
    //T x2=atan2(r[5]/cos(y2),r[8]/cos(y2));
    T z1=atan2(r[1]/cos(y1),r[0]/cos(y1));
    //T z2=atan2(r[1]/cos(y2),r[0]/cos(y2));

    x=x1; // always take "first" solution
    y=y1;
    z=z1;
  }
  else {
    z=T(0);
    if (r[2]<T(0)) { // ==-1
      y=T(M_PI_2);
      z=z+atan2(r[3],r[6]);
    }
    else {
      y=-T(M_PI_2);
      z=-z+atan2(-r[3],-r[6]);
    }
  }
  return VecN<T,3>(x,y,z);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace math
} // namespace VC

#endif // __VC_MATH_MAT3_HH__
