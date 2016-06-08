//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MATH_QUATERNION
#define VC_MATH_QUATERNION

//== INCLUDES =================================================================

#include <complex>
#include <cassert>
#include <limits>
#include <iostream>

#include <vclibs/math/VecN.hh>
#include <vclibs/math/Mat.hh>


//== CLASS DEFINITION =========================================================

// TODO: remove this ASAP !!!
# if defined(_MSC_VER)
#  define constexpr
# endif

/** \defgroup vc_math_quaternion Quaternions
    \ingroup vc_math
    \sa [\ref vc_math_lam], [\ref vc_math_lav]
 */

namespace VC {
namespace math {

/** Quaternion
    \ingroup vc_math_quaternion

    This class is based on VC::math::VecN<T,4> elts() with real component
    real() referring to elts()[real_index()] and imaginary components imag()
    referring to `elts()[image_index():imag_index+2]`.

    Note that norm() is defined as `norm(q)==sqrt(q|q)==sqrt(sqr(q))` (some
    definitions don't apply the square root).

    \sa VC::math::VecN
 */
template <typename T>
class Quaternion {
public:
  typedef T value_type;
  typedef Quaternion<T> self_t;
  typedef T real_t;
  typedef VecN<T,3> imag_t;
  typedef VecN<T,4> elts_t;
  typedef VC::math::Mat<real_t,4,4> mat4x4_t;
  typedef VC::math::Mat<real_t,3,3> mat3x3_t;

  /// quaternion with zero real() and imag() components
  Quaternion() {}
  /// real quaternion with zero imag() component
  Quaternion(real_t _re) : Quaternion() { m_elts[real_index()]=_re; }
  /// quaternion from real() `_re` and imag() `_im`
  Quaternion(real_t _re,imag_t _im) { real()=_re; imag()=_im; }
  /// quaternion from 4-vector `_elts` (mind mapping of components)
  explicit Quaternion(const elts_t& _v) : m_elts(_v) {}
  /// quaternion from real() `_re` and imag() `(_x,_y,_z)`
  Quaternion(real_t _re,real_t _x,real_t _y,real_t _z)
    : Quaternion(_re,imag_t(_x,_y,_z)) {}
  /// pure imaginary quaternion
  Quaternion(const imag_t& _im) : Quaternion(0,_im) {}
  /// quaternion from complex number with `imag()==(_c.imag(),0,0)`
  Quaternion(const std::complex<T>& _c) : Quaternion(_c.real(),_c.imag(),0,0) {}
  /** Unit quaternion from vectors `_v1,_v2`.
      This represents a rotation around `_v1%_v2` by an angle of
      `2*acos(_v1|_v2)`. Note that `_v1` and `_v2` should be unit vectors;
      they are not normalized by the constructor!
  */
  Quaternion(const imag_t& _v1,const imag_t& _v2)
    : Quaternion(-(_v1|_v2),_v1%_v2) {}

  Quaternion(const self_t&) = default;
  self_t& operator=(const self_t&) = default;

  /** Return quaternion representing a rotation around `_axis`.
      \param _alpha rotation angle (rad)
      \param _axis direction with arbitrary but non-zero length
      \return `Quaternion(cos(_alpha/2),sin(alpha/2)*_axis)`
   */
  static self_t rotation(real_t _alpha,const imag_t& _axis) {
    real_t omega=_alpha/real_t(2);
    return self_t(cos(omega),_alpha*(sin(omega)/_axis.norm()));
  }

  /// cast to 4-vector (mind mapping of components)
  operator elts_t() const { return m_elts; }

  /// assign real Quaternion
  self_t& operator=(real_t _re) { return *this=Quaternion(_re); }
  /// assign purely imaginary Quaternion
  self_t& operator=(imag_t _im) { return *this=Quaternion(_im); }

  /// get elements as 4-vector (mind mapping of components)
  const elts_t elts() const { return m_elts; }
  /// get elements as 4-vector (mind mapping of components)
  elts_t elts() { return m_elts; }

  /// index of real component in elts()
  static constexpr int real_index() { return 3; }
  /// index of first imaginary component in elts()
  static constexpr int imag_index() { return 0; }

  real_t real() const { return m_elts[real_index()]; } //!< real part
  real_t& real() { return m_elts[real_index()]; } //!< real part

  /// imaginary part
  const imag_t& imag() const {
    return imag_t::to_cv(m_elts.begin()+imag_index());
  }
  /// imaginary part
  imag_t& imag() {
    return imag_t::to_v(m_elts.begin()+imag_index());
  }

  /// add quaternion
  self_t& operator+=(const self_t& _q){ m_elts+=_q.m_elts; return *this; }
  /// subtract quaternion
  self_t& operator-=(const self_t& _q){ m_elts-=_q.m_elts; return *this; }

  /// addition
  self_t operator+(const self_t& _q) const { return self_t(m_elts+_q.m_elts); }
  /// subtraction
  self_t operator-(const self_t& _q) const { return self_t(m_elts-_q.m_elts); }
  /// unary `+`
  const self_t& operator+() const { return *this; }
  /// unary `-`
  self_t operator-() const { return self_t(-m_elts); }

  /// add scalar
  self_t& operator+=(const real_t& _s) {
    m_elts[real_index()]+=_s; return *this;
  }
  /// subtract scalar
  self_t& operator-=(const real_t& _s) {
    m_elts[real_index()]-=_s; return *this;
  }
  /// addition
  self_t operator+(const real_t& _s) { return self_t(real()+_s,imag()); }
  /// subtraction
  self_t operator-(const real_t& _s) { return self_t(real()-_s,imag()); }

  /// multiply with scalar
  self_t& operator*=(const real_t& _s) { m_elts*=_s; return *this; }
  /// divide by scalar
  self_t& operator/=(const real_t& _s) { m_elts/=_s; return *this; }
  ///  multiplication with scalar
  self_t operator*(const real_t& _s) { return self_t(m_elts*_s); }
  //// division by scalar
  self_t operator/(const real_t& _s) { return self_t(m_elts/_s); }

  /// get conjugate
  self_t conj() const { return self_t(real(),-imag()); }
  /// get inverse
  self_t inverse() const { return conj()/sqr(); }
  /// get norm as `sqrt(q|q)`
  real_t norm() const { return m_elts.norm(); }
  /// get squared norm as `(q|q)`
  real_t sqr() const { return m_elts.sqr(); }
  /// make `*this` unit quaternion
  self_t& normalize() { m_elts.normalize(); return *this; }
  /// get unit quaternion from `*this`
  self_t unit_quaternion() const { return *this/norm(); }

  /// dot product
  real_t operator|(const self_t& _q) const { return (m_elts|_q.m_elts); }
  /// Hamilton product
  self_t operator*(const self_t& _q) {
    return self_t(real()*_q.real()-(imag()|_q.imag()),
                  _q.imag()*real()+imag()*_q.real()+(imag()%_q.imag())); // TODO: expand
  }
  /// Hamilton product
  self_t& operator*=(const self_t& _q){ return *this=*this*_q; }
  /// cross product
  self_t operator%(const self_t& _q) const {
    return self_t(imag()%_q.imag());
  }
  /// cross product
  self_t& operator%=(const self_t& _q) { return *this=*this%_q; }

  /// test for unit quaternion (uses some epsilon)
  bool is_unit_quaternion() const {
    return fabs(sqr()-real_t(1))<std::numeric_limits<real_t>::epsilon()*4;
  }
  /// test for purely imaginary quaternion (exact test `real()==0`)
  bool is_pure() const { return real()==real_t(0); } // exact test!

  /// spherical linear interpolation
  self_t slerp(const self_t& _q0,const self_t& _q1,real_t _t) {
    real_t m0=_q0.norm();
    real_t m1=_q1.norm();
    real_t m=(real_t(1)-_t)*m0+_t*m1;

    self_t p0=_q0/m0;
    self_t p1=_q1/m1;
    real_t w=(p0.conj()*p1).real();
    assert(real_t(-1)<=w && w<=real_t(1));
    real_t theta=acos(w);
    self_t p=p0*sin((real_t(1)-_t)*theta);
    p+=p1*sin(_t*theta);
    p*=m/sin(theta);

    return p;
  }

  /// get 4x4 matrix representation
  void to_real_matrix(mat4x4_t& _m) const {
    auto   m=_m.begin(); // column major
    real_t a=m_elts[real_index()];
    real_t b=m_elts[imag_index()];
    real_t c=m_elts[imag_index()+1];
    real_t d=m_elts[imag_index()+2];
    m[ 0]= a; m[ 1]=-b; m[ 2]=-c; m[ 3]=-d;
    m[ 4]= b; m[ 5]= a; m[ 6]= d; m[ 7]=-c;
    m[ 8]= c; m[ 9]=-d; m[10]= a; m[11]= b;
    m[12]= d; m[13]= c; m[14]=-b; m[15]= a;
  }

  /// get orthogonal matrix from **unit quaternion**
  void to_orthogonal_matrix(mat3x3_t& _m) const {
    to_orthogonal_matrix(_m.begin(),3);
  }
  /** Get 4x4 orthogonal matrix from **unit quaternion**.
      orthogonal 3x3 block `_m(1:3,1:3)` and identity map for the
      (homogeneous) 4-th coordinate (`_m(4,:)=_m(:,4)'=[0 0 0 1]`).
      \param[out] _m homogeneous transformation matrix
   */
  void to_orthogonal_matrix(mat4x4_t& _m) const {
    to_orthogonal_matrix(_m.begin(),4);
    real_t* m=_m.begin();
    m[ 3]=m[ 7]=m[11]=            // last row
    m[12]=m[13]=m[14]=real_t(0);  // last column
    m[15]            =real_t(1);
  }
  /** Get 3x3 orthogonal matrix from **unit quaternion**.
      \param[out] matrix entries (column major)
      \param _ld leading dimension of `_m`
   */
  void to_orthogonal_matrix(real_t* _m,int _ld) const {
    assert(_ld>=3);
    assert(is_unit_quaternion());
    real_t w=m_elts[real_index()],
           x=m_elts[imag_index()],
           y=m_elts[imag_index()+1],
           z=m_elts[imag_index()+2];

    _m[0]=T(1)-T(2)*y*y-T(2)*z*z;
    _m[1]=     T(2)*x*y+T(2)*w*z;
    _m[2]=     T(2)*x*z-T(2)*w*y;
    _m+=_ld;

    _m[0]=     T(2)*x*y-T(2)*w*z;
    _m[1]=T(1)-T(2)*x*x-T(2)*z*z;
    _m[2]=     T(2)*y*z+T(2)*w*x;
    _m+=_ld;

    _m[0]=     T(2)*x*z+T(2)*w*y;
    _m[1]=     T(2)*y*z-T(2)*w*x;
    _m[2]=T(1)-T(2)*x*x-T(2)*y*y;
  }

  // TODO: general quaternion

private:
  elts_t m_elts = elts_t(0,0,0,0);
};

/// Quaternion*Scalar \ingroup vc_math
template <typename T>
Quaternion<T> operator*(T _s,const Quaternion<T>& _q) {
  return _q*_s;
}
/// Quaternion/Scalar \ingroup vc_math
template <typename T>
Quaternion<T> operator/(T _s,const Quaternion<T>& _q) {
  return _q/_s;
}

/// same as Quaternion::norm() \ingroup vc_math
template <typename T>
T norm(const Quaternion<T>& _q) { return _q.norm(); }

/// same as Quaternion::sqr() \ingroup vc_math
template <typename T>
T sqr(const Quaternion<T>& _q) { return _q.sqrnorm(); }

/// same as Quaternion::unit_quaternion() \ingroup vc_math
template <typename T>
Quaternion<T> unit_quaternion(const Quaternion<T>& _q) {
  return _q.unit_quaternion();
}

/// output Quaternion as 4-vector \ingroup vc_math
template <typename T>
std::ostream& operator<<(std::ostream& _out,const Quaternion<T>& _q) {
  return _out << _q.elts();
}

//-----------------------------------------------------------------------------

# if defined(_MSC_VER)
#  undef constexpr // TODO: remove this ASAP !!!
# endif

//=============================================================================
} //namespace math
} //namespace VC
//=============================================================================
#endif // VC_MATH_QUATERNION defined
