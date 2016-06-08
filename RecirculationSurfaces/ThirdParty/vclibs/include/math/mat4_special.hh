//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#include "Mat.hh"

#ifndef __VC_MATH_MAT4_HH__
#define __VC_MATH_MAT4_HH__

namespace VC {
namespace math {

//
// specializations for 4x4 matrices
//

# ifndef VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/// specialization \ingroup vc_math_lam_special
template <typename T>
struct inverse_t<T,4,4> {
  static bool get(Mat<T,4,4>& _a) {

    // source: Intel

    T* dst=_a.data();
    T tmp[12]; /* temp array for pairs */
    T src[16]; /* array of transpose source matrix */
    T det; /* determinant */
    /* transpose matrix */
    for (int i = 0; i < 4; i++) {
      src[i] = dst[i*4];
      src[i + 4] = dst[i*4 + 1];
      src[i + 8] = dst[i*4 + 2];
      src[i + 12] = dst[i*4 + 3];
    }
    /* calculate pairs for first 8 elements (cofactors) */
    tmp[0] = src[10] * src[15];
    tmp[1] = src[11] * src[14];
    tmp[2] = src[9] * src[15];
    tmp[3] = src[11] * src[13];
    tmp[4] = src[9] * src[14];
    tmp[5] = src[10] * src[13];
    tmp[6] = src[8] * src[15];
    tmp[7] = src[11] * src[12];
    tmp[8] = src[8] * src[14];
    tmp[9] = src[10] * src[12];
    tmp[10] = src[8] * src[13];
    tmp[11] = src[9] * src[12];
    /* calculate first 8 elements (cofactors) */
    dst[0] = tmp[0]*src[5] + tmp[3]*src[6] + tmp[4]*src[7];
    dst[0] -= tmp[1]*src[5] + tmp[2]*src[6] + tmp[5]*src[7];
    dst[1] = tmp[1]*src[4] + tmp[6]*src[6] + tmp[9]*src[7];
    dst[1] -= tmp[0]*src[4] + tmp[7]*src[6] + tmp[8]*src[7];
    dst[2] = tmp[2]*src[4] + tmp[7]*src[5] + tmp[10]*src[7];
    dst[2] -= tmp[3]*src[4] + tmp[6]*src[5] + tmp[11]*src[7];
    dst[3] = tmp[5]*src[4] + tmp[8]*src[5] + tmp[11]*src[6];
    dst[3] -= tmp[4]*src[4] + tmp[9]*src[5] + tmp[10]*src[6];
    dst[4] = tmp[1]*src[1] + tmp[2]*src[2] + tmp[5]*src[3];
    dst[4] -= tmp[0]*src[1] + tmp[3]*src[2] + tmp[4]*src[3];
    dst[5] = tmp[0]*src[0] + tmp[7]*src[2] + tmp[8]*src[3];
    dst[5] -= tmp[1]*src[0] + tmp[6]*src[2] + tmp[9]*src[3];
    dst[6] = tmp[3]*src[0] + tmp[6]*src[1] + tmp[11]*src[3];
    dst[6] -= tmp[2]*src[0] + tmp[7]*src[1] + tmp[10]*src[3];
    dst[7] = tmp[4]*src[0] + tmp[9]*src[1] + tmp[10]*src[2];
    dst[7] -= tmp[5]*src[0] + tmp[8]*src[1] + tmp[11]*src[2];
    /* calculate pairs for second 8 elements (cofactors) */
    tmp[0] = src[2]*src[7];
    tmp[1] = src[3]*src[6];
    tmp[2] = src[1]*src[7];
    tmp[3] = src[3]*src[5];
    tmp[4] = src[1]*src[6];
    tmp[5] = src[2]*src[5];

    tmp[6] = src[0]*src[7];
    tmp[7] = src[3]*src[4];
    tmp[8] = src[0]*src[6];
    tmp[9] = src[2]*src[4];
    tmp[10] = src[0]*src[5];
    tmp[11] = src[1]*src[4];
    /* calculate second 8 elements (cofactors) */
    dst[8] = tmp[0]*src[13] + tmp[3]*src[14] + tmp[4]*src[15];
    dst[8] -= tmp[1]*src[13] + tmp[2]*src[14] + tmp[5]*src[15];
    dst[9] = tmp[1]*src[12] + tmp[6]*src[14] + tmp[9]*src[15];
    dst[9] -= tmp[0]*src[12] + tmp[7]*src[14] + tmp[8]*src[15];
    dst[10] = tmp[2]*src[12] + tmp[7]*src[13] + tmp[10]*src[15];
    dst[10]-= tmp[3]*src[12] + tmp[6]*src[13] + tmp[11]*src[15];
    dst[11] = tmp[5]*src[12] + tmp[8]*src[13] + tmp[11]*src[14];
    dst[11]-= tmp[4]*src[12] + tmp[9]*src[13] + tmp[10]*src[14];
    dst[12] = tmp[2]*src[10] + tmp[5]*src[11] + tmp[1]*src[9];
    dst[12]-= tmp[4]*src[11] + tmp[0]*src[9] + tmp[3]*src[10];
    dst[13] = tmp[8]*src[11] + tmp[0]*src[8] + tmp[7]*src[10];
    dst[13]-= tmp[6]*src[10] + tmp[9]*src[11] + tmp[1]*src[8];
    dst[14] = tmp[6]*src[9] + tmp[11]*src[11] + tmp[3]*src[8];
    dst[14]-= tmp[10]*src[11] + tmp[2]*src[8] + tmp[7]*src[9];
    dst[15] = tmp[10]*src[10] + tmp[4]*src[8] + tmp[9]*src[9];
    dst[15]-= tmp[8]*src[9] + tmp[11]*src[10] + tmp[5]*src[8];
    /* calculate determinant */
    det=src[0]*dst[0]+src[1]*dst[1]+src[2]*dst[2]+src[3]*dst[3];

    if (det==T(0))
      return false;

    /* calculate matrix inverse */
    det = 1/det;
    for (int j = 0; j < 16; j++)
      dst[j] *= det;

    return true;
  }
};

//-----------------------------------------------------------------------------
# endif // VC_MAT_NO_SPECIALIZATIONS
//-----------------------------------------------------------------------------

/** \defgroup vc_math_lam_4x4 Functions of 4x4 matrices
    \ingroup vc_math_lam
 */


/** Interpret 4x4 matrix as transformation of homogeneous coordinates.
    \ingroup vc_math_lam_4x4
    \todo check transformations
 */
template <typename T>
struct transformation_t<T,4> {
  typedef transformation_t<T,4> self_t; //!< transformation_t<T,4>
  typedef Mat<T,4,4> mat4_t; //!< matrix type
  typedef Mat<T,3,3> mat3_t; //!< 3x3 matrix (linear transformation part)
  typedef VecN<T,4>  vec4_t; //!< 3-vector or point in 3d
  typedef VecN<T,3>  vec3_t; //!< point/direction in homogeneous coordinates

  mat4_t& m;                  //!< the "wrapped" matrix
  transformation_t(mat4_t& _m) : m(_m) {}

  /// access first three entries of column `_j` (`m(1:3,_j)`)
  const vec3_t& col3(int _j) const { return vec3_t::to_cv(m.col(_j).data()); }
  /// access first three entries of column `_j` (`m(1:3,_j)`)
  vec3_t& col3(int _j) { return vec3_t::to_v(m.col(_j).data()); }
  // get first three entries of row `_j` (`m(_i,1:3)`)
  vec3_t row3(int _i) {
    return vec3_t(m.elts()[_i],m.elts()[_i+4],m.elts()[_i+8]);
  }
  // set first three entries of row `_j` (`m(_i,1:3)`)
  void set_row3(int _i,const vec3_t _x) {
    m.elts()[_i  ]=_x[0];
    m.elts()[_i+4]=_x[1];
    m.elts()[_i+8]=_x[2];
  }

  /// get matrix
  const mat4_t& get() const { return m; }

  /** Get linear operator from affine transformation matrix
      \return upper left block `_a(1:3,1:3)`
  */
  mat3_t get_linear() const { return mat3_t(col3(0),col3(1),col3(2)); }
  /** Set linear operator from affine transformation matrix.
      \param _a3 linear operator that replaces block `m(1:3,1:3)`
      \return `*this`
  */
  self_t& set_linear(const mat3_t& _a3) {
    col3(0)=_a3.col(0);
    col3(1)=_a3.col(1);
    col3(2)=_a3.col(2);
    return *this;
  }

  /** Access translation vector from affine transformation matrix.
      \return upper right block `_a(1:3,4)`
  */
  const vec3_t& translation() const { return col3(3); }
  /** Access translation vector from affine transformation matrix.
      \return upper right block `_a(1:3,4)`
  */
  vec3_t& translation() { return col3(3); }

  /** Apply homogeneous transformation to point.
      \param _x point
      \return transformed point
  */
  vec3_t map_point(const vec3_t& _x) const {
    vec4_t y(m.col(0)*_x[0]);
    y+=m.col(1)*_x[1];
    y+=m.col(2)*_x[2];
    y+=m.col(3);
    return vec3_t(y[0],y[1],y[2])/y[3];
  }

  /** Apply homogeneous transformation to vector.
      \param _x vector/direction
      \return transformed vector/direction
  */
  vec3_t map_vector(const vec3_t& _x) const {
    vec3_t y(col3(0)*_x[0]);
    y+=col3(1)*_x[1];
    y+=col3(2)*_x[2];
    return y;
  }

  /** Apply translation `_x` (same as `glTranslate`).
      See [glTranslate](http://www.opengl.org/sdk/docs/man2/xhtml/glTranslate.xml).
      \param _x translation vector
      \return `*this`
  */
  self_t& translate(const vec3_t& _x) {
    m.col(3)[0]+=(_x|row3(0));
    m.col(3)[1]+=(_x|row3(1));
    m.col(3)[2]+=(_x|row3(2));
    m.col(3)[3]+=(_x|row3(3));
    return *this;
  }

  /** Apply scaling `_x` (same as `glScale`).
      See [glScale](http://www.opengl.org/sdk/docs/man2/xhtml/glScale.xml).
      \param _x scaling in x,y,z direction
      \return `*this`
  */
  self_t& scale(const VecN<T,3>& _x) {
    m.col(0)*=_x[0];
    m.col(1)*=_x[1];
    m.col(2)*=_x[2];
    return *this;
  }

  /** Get rotation around axis `_x`.
      \param _angle in *degrees*
      \param _x rotation axis
      \return rotation matrix
  */
  static mat4_t rotation(T _angle,const vec3_t& _x) {
    mat4_t r;
    r.set_linear(rot3(_angle*T(M_PI/180.0),_x));
    r.col(3)=vec4_t(0,0,0,1);
    r.set_row3(3,vec3_t(0,0,0));
    return r;
  }

  /** Get rotation that rotates direction `_x` into `_y`.
      \param _x direction
      \param _y direction
      \return rotation matrix
  */
  static mat4_t rotation(const vec3_t& _x,const vec3_t& _y) {
    mat4_t r;
    r.set_linear(rot3(_x,_y));
    r.col(3)=vec4_t(0,0,0,1);
    r.set_row3(3,vec3_t(0,0,0));
    return r;
  }

  /** Apply rotation around axis `_x` (same as `glRotate`).
      See [glRotate](http://www.opengl.org/sdk/docs/man2/xhtml/glRotate.xml).
      \param _angle in *degrees*
      \param _x rotation axis
      \return `*this`
  */
  self_t& rotate(T _angle,const vec3_t& _x) {
    m=m*rotation(_angle,_x);
    return *this;
  }

  /** Apply rotation that rotates direction `_x` into `_y`.
      \param _x direction
      \param _y direction
      \return `*this`
  */
  self_t& rotate(const vec3_t& _x,const vec3_t& _y) {
    m=m*rotation(_x,_y);
    return *this;
  }

  /** Apply viewing transformation (same as `gluLookAt`).
      See [gluLookAt](http://www.opengl.org/sdk/docs/man2/xhtml/gluLookAt.xml).
      \param _eye eye point
      \param _center focus point
      \param _up up vector
      \return `*this`
  */
  self_t& lookAt(const vec3_t& _eye,const vec3_t& _center,const vec3_t& _up) {
    vec3_t z(_eye);
    z-=_center;

    z.normalize();

    vec3_t x(_up%z);
    vec3_t y(z%x);

    x.normalize();
    y.normalize();

    mat4_t a;
    a.col(0)=vec4_t(x[0],y[0],z[0],T(0));
    a.col(1)=vec4_t(x[1],y[1],z[1],T(0));
    a.col(2)=vec4_t(x[2],y[2],z[2],T(0));
    a.col(3)=vec4_t(0,0,0,1);

    m=m*a;
    return translate(-_eye);
  }

  /** Apply viewing transformation (same as `glFrustum`)
      See [glFrustum](http://www.opengl.org/sdk/docs/man2/xhtml/glFrustum.xml).
      \param _left
      \param _right
      \param _bottom
      \param _top
      \param _near
      \param _zfar
      \return `*this`
  */
  self_t& frustum(T _left,T _right,T _bottom, T _top,T _znear, T _zfar) {
    assert(_znear>T(0) && _zfar>_znear);

    T x, y, a, b, c, d;

    x = (T(2)*_znear) / (_right-_left);
    y = (T(2)*_znear) / (_top-_bottom);
    a = (_right+_left) / (_right-_left);
    b = (_top+_bottom) / (_top-_bottom);
    c = -(_zfar+_znear) / ( _zfar-_znear);
    d = -(T(2)*_zfar*_znear) / (_zfar-_znear);

    mat4_t X;
    X.col(0)=vec4_t(x,   T(0),T(0),T(0));
    X.col(1)=vec4_t(T(0),y,   T(0),T(0));
    X.col(2)=vec4_t(a,   b,   c,   T(-1));
    X.col(3)=vec4_t(T(0),T(0),d,   T(0));

    m=m*X;

    return *this;
  }

  /** Apply inverse viewing transformation (see frustum())
      \param _left
      \param _right
      \param _bottom
      \param _top
      \param _near
      \param _zfar
      \return `*this`
  */
  self_t& inverse_frustum(T _left,T _right,T _bottom, T _top,T _znear, T _zfar) {
    assert(_znear>T(0) && _zfar>_znear);

    T x, y, a, b, c, d;

    x = (_right-_left) / (T(2)*_znear);
    y = (_top-_bottom) / (T(2)*_znear);
    a = (_right+_left) / (T(2)*_znear);
    b = (_top+_bottom) / (T(2)*_znear);
    c = (_zfar+_znear) / (T(2)*_zfar*_znear);
    d = (_znear-_zfar) / (T(2)*_zfar*_znear);

    mat4_t X;
    X.col(0)=vec4_t(x,   T(0),T(0), T(0));
    X.col(1)=vec4_t(T(0),y,   T(0), T(0));
    X.col(2)=vec4_t(T(0),T(0),T(0), d);
    X.col(3)=vec4_t(a,   b,   T(-1),c);

    m=m*X;

    return *this;
  }

  /** Apply orthographic viewing transformation (same as `glOrtho`)
      See [glOrtho](http://www.opengl.org/sdk/docs/man2/xhtml/glOrtho.xml).
      \param _left
      \param _right
      \param _bottom
      \param _top
      \param _near
      \param _zfar
      \return `*this`
   */
  self_t& ortho(T _left,T _right,T _bottom, T _top,T _znear, T _zfar) {
    assert(_znear>T(0) && _zfar>_znear);

    T x,y,z,tx,ty,tz;

    x =  T(2)/(_right-_left);
    y =  T(2)/(_top-_bottom);
    z = -T(2)/(_zfar-_znear);
    tx = -(_right+_left)/(_right-_left);
    ty = -(_top+_bottom)/(_top-_bottom);
    tz = -(_zfar+_znear)/(_zfar-_znear);

    mat4_t X;
    X.col(0)=vec4_t(x,   T(0),T(0), T(0));
    X.col(1)=vec4_t(T(0),y,   T(0), T(0));
    X.col(2)=vec4_t(T(0),T(0),z,    T(0));
    X.col(3)=vec4_t(tx,  ty,  tz,   T(1));

    m=m*X;

    return *this;
  }

  /** Apply inverse orthographic projection (see ortho())
      \param _left
      \param _right
      \param _bottom
      \param _top
      \param _near
      \param _zfar
      \return `*this`
   */
  self_t& inverse_ortho(T _left,T _right,T _bottom, T _top,T _znear, T _zfar) {
    assert(_znear>T(0) && _zfar>_znear);

    T x, y, a, b, c, d;

    x = (_right-_left) / (T(2));
    y = (_top-_bottom) / (T(2));
    a = (_right+_left) / (T(2));
    b = (_top+_bottom) / (T(2));
    c = -(_zfar+_znear) / (T(2));
    d = (_znear-_zfar) / (T(2));

    mat4_t X;
    X.col(0)=vec4_t(x,   T(0),T(0), T(0));
    X.col(1)=vec4_t(T(0),y,   T(0), T(0));
    X.col(2)=vec4_t(T(0),T(0),d,    T(0));
    X.col(3)=vec4_t(a,   b,   c,    T(1));

    m=m*X;

    return *this;
  }


  /** Apply perspective projection (same as `gluPerspective`)
      See [gluPerspective](http://www.opengl.org/sdk/docs/man2/xhtml/gluPerspective.xml).
      \param _fovy field of view angle
      \param _aspect aspect ratio
      \param _near
      \param _zfar
      \return `*this`
  */
  self_t& perspective(T _fovy,T _aspect,T _znear,T _zfar) {
    T xmin, xmax, ymin, ymax;

    ymax=_znear*tan(_fovy*T(M_PI/360.0));
    ymin=-ymax;

    xmin=ymin*_aspect;
    xmax=ymax*_aspect;

    return frustum(xmin,xmax,ymin,ymax,_znear,_zfar);
  }

  /** Apply perspective projection (same as `gluPerspective`)
      \param _fovy field of view angle
      \param _aspect aspect ratio
      \param _near
      \param _zfar
      \return `*this`
  */
  self_t& inverse_perspective(T _fovy,T _aspect,T _znear,T _zfar) {
    T xmin, xmax, ymin, ymax;

    ymax=_znear*tan(_fovy*T(M_PI/360.0));
    ymin=-ymax;

    xmin=ymin*_aspect;
    xmax=ymax*_aspect;

    return inverse_frustum(xmin,xmax,ymin,ymax,_znear,_zfar);
  }

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} // namespace math
} // namespace VC


#endif // __VC_MATH_MAT4_HH__
