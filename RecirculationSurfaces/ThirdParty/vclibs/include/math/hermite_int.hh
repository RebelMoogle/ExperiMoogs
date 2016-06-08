
//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MATH_HERMITEINT_HH
#define VC_MATH_HERMITEINT_HH

#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>

//=============================================================================

namespace VC {
namespace math {

/** \defgroup vc_hermite_int Hermite interpolation
    \ingroup vc_math
    Provides Hermite interpolation of univariate data in arbitrary
    dimensions.

    - eval() evaluates a cubic piece

    - eval_at(), eval_linspace() evaluate a C^1 cubic spline

    - linspace() provides a uniform sampling of an interval

    - length(), accum_length(), length_bounds(), accum_arclength()
      estimate arc-length.

    - subdivide() subdivides the curve recursively to meet a spatial
      error tolerance. (This function is intended for *"small"* vectors,
      see below.)

    Knot vectors must be *strictly monotone* sequences. The preferred
    order is *ascending*. Some operations may require ascending order.

    The implementation is designed to require a minimum of supported
    vector operations and to minimize use of temporary vector
    variables. This means pre-allocation by caller and no use of de
    Casteljau algorithm: some functions take a `work` parameter for
    storing vector values temporary results. For convenience we
    provide also functions which construct temporary vectors and avoid
    the extra `work` parameter. This is indicated in the
    documentation and recommended only for small, fixed size vector
    objects like, e.g, VC::math::VecN.

    The vector class must support the following operations

    | operation  | specification               |
    |------------|-----------------------------|
    | assignment | `vec=vec`                   |
    | scaling    | `vec*=scalar (vec*scalar)`  |
    | axpy       | `vec+=vec*scalar, vec-=vec` |
    | 2-norm     | `norm2(vec)`                |

    In addition it must define `value_type` which is the scalar type.

    @todo find better stopping criterion for subdivide() (more uniform sampling)
    @todo approximate length up to tolerance tol<upper-lower by subdivision
    @todo reparametrization (e.g., arc length)
    @todo cardinal splines (with option to estimate and store dx)
 */

/// Hermite interpolation [\ref vc_hermite_int] \ingroup vc_hermite_int
namespace hermite_int {

//-----------------------------------------------------------------------------

/** Evaluate Hermite interpolant on [_ta,_tb] (position and tangent).
    \ingroup vc_hermite_int
    For the <em>special case</em> _ta==_tb eval() returns _x0 and _dx0!
    \param[out] _x result (position)
    \param[out] _dx result (tangent), no output if _dx==0
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _t evaluate at this parameter value
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
 */
template <typename vec>
void eval(vec& _x,vec* _dx,
          const vec& _x0,const vec& _dx0,
          const vec& _x1,const vec& _dx1,
          const typename vec::value_type& _t,
          typename vec::value_type _ta=typename vec::value_type(0),
          typename vec::value_type _tb=typename vec::value_type(1)) {

  typedef typename vec::value_type real;

  if (_ta==_tb)   { _x=_x0; if (_dx!=0) *_dx=_dx0; return; }

  real h=_tb-_ta;
  real t=(_t-_ta)/h;
#ifndef __CUDACC__
  assert(std::isfinite(t));
#else
  assert(isfinite(t));
#endif

  if (t==real(0)) { _x=_x0; if (_dx!=0) *_dx=_dx0; return; }
  if (t==real(1)) { _x=_x1; if (_dx!=0) *_dx=_dx1; return; }

  real t2=t*t;
  real t3=t2*t;

  real cx0 =real(2)*t3-real(3)*t2+real(1);   // (1+2*t)*)1-t)^2
  real cdx0=t3-real(2)*t2+t;                 // t*(1-t)^2
  real cx1 =real(3)*t2-real(2)*t3;           // t2*(3-2*t)
  real cdx1=t3-t2;                           // t2*(t-1)

  //_x =_x0*cx0 + _x1*cx1 + _dx0*(cdx0*h) + _dx1*(cdx1*h);
  _x =_x0;
  _x*=cx0;
  _x+= _x1*cx1;
  _x+=_dx0*(cdx0*h);
  _x+=_dx1*(cdx1*h);

  if (_dx!=0) {
    real ex0 =(t2-t)*real(6);               // 6*(t2-t)
    real edx0=real(3)*t2-real(4)*t+real(1); // (t-1)*(3*t-1)
    real ex1 =-ex0;                         // -6*(t2-t)
    real edx1=real(3)*t2-real(2)*t;         // (3*t-2)*t

    //*_dx=_x0*(ex0/h) + _x1*(ex1/h) + _dx0*edx0 + _dx1*edx1;
    *_dx=_x0;
    *_dx*=(ex0/h);
    *_dx+=_x1*(ex1/h);
    *_dx+=_dx0*edx0;
    *_dx+=_dx1*edx1;
  }
}

/** Evaluate Hermite interpolant on [_ta,_tb] (position and tangent).
    \ingroup vc_hermite_int
    \param[out] _x result (position)
    \param[out] _dx result (tangent)
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _t evaluate at this parameter value
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
 */
template <typename vec>
void eval(vec& _x,vec& _dx,
          const vec& _x0,const vec& _dx0,
          const vec& _x1,const vec& _dx1,
          typename vec::value_type _t,
          typename vec::value_type _ta=typename vec::value_type(0),
          typename vec::value_type _tb=typename vec::value_type(1)) {
  eval(_x,&_dx,_x0,_dx0,_x1,_dx1,_t,_ta,_tb);
}

/** Find parameter interval of _t in knot vector using linear forward search.
    \ingroup vc_hermite_int
    \param _knots base of knot vector
    \param _n number of knots
    \param _t parameter value
    \return index 0<=i<_n-1 of first knot of interval or -1 or _n-1 if out of bounds
 */
template <typename real>
int find_interval_forward(const real* _knots,size_t _n,real _t) {
  size_t i;
  assert(_n>1);
  if (_t<_knots[0])
    return -1;
  for (i=0;i<_n-1;++i)
    if (_knots[i]<=_t && _t<=_knots[i+1]) // stick to <= for both
      return int(i);
  return int(i);
}

/** Find parameter interval of _t in knot vector using linear backward search.
    \ingroup vc_hermite_int
    \param _knots base of knot vector
    \param _n number of knots
    \param _t parameter value
    \return index 0<=i<_n-1 of first knot of interval or -1 or _n-1 if out of bounds
 */
template <typename real>
int find_interval_backward(const real* _knots,size_t _n,real _t) {
  size_t i;
  assert(_n>1);
  if (_t>_knots[_n-1])
    return _n;
  for (i=_n-2;i>0;--i)
    if (_knots[i]<=_t && _t<=_knots[i+1]) // stick to <= for both
      return int(i);
  return int(i);
}

/// functor that wraps find_interval_forward() \ingroup vc_hermite_int
template <typename real>
struct find_interval_forward_t {
  int operator()(const real* _knots,size_t _n,real _t) {
    return find_interval_forward(_knots,_n,_t);
  }
};

/// functor that wraps find_interval_backward() \ingroup vc_hermite_int
template <typename real>
struct find_interval_backward_t {
  int operator()(const real* _knots,size_t _n,real _t) {
    return find_interval_forward(_knots,_n,_t);
  }
};

/** Sample Hermite interpolating cubic C^1 spline at parameter values _t.
    \ingroup vc_hermite_int
    Extrapolation for any parameters  out of range [knots[0],_knots[_nk-1]].
    \tparam FIND functor to find knot intervals; the default setting
    LinearSearch uses a specialization with linear search using `FIND=`
    find_interval_forward() or find_interval_backward() depending on the
    bounds of the first knot interval
    \param[out] _xout write _nt position samples
    \param[out] _dxout _nt write tangent samples (no output if _dxout==0)
    \param _x _nk positions interpolated at _knots
    \param _dx _nk tangents interpolated at _knots
    \param _knots knot vector with _kn entries
    \param _nk >=2  number of _knots (and _x,_dx entries)
    \param _find functor to find knot intervals
    \param _t array of parameter values for sampling of size _nt (should be sorted)
    \param _nt number of sample points (size of _t and _xout, _dxout)
 */
template <typename vec,typename FIND>
void eval_at(vec* _xout,vec* _dxout,
             const vec* _x,const vec* _dx,
             const typename vec::value_type* _knots,size_t _nk,
             FIND _find,
             const typename vec::value_type* _t,size_t _nt) {
  assert(_nk>=2);

  int    idx=0;  // restart search for knot interval from current index

  for (size_t i=0;i<_nt;++i) {
    if (_t[i]<_knots[idx])
      idx=0;  // restart search from beginning
    int idxnew=_find(_knots+idx,_nk-idx,_t[i]);
    if (idxnew<0 || idx+idxnew+1>=int(_nk)) {
      // extrapolate
      int j=idxnew<0 ? 0 : int(_nk-2);
      eval(_xout[i],_dxout!=0 ? _dxout+i : 0,
           _x[j],_dx[j],_x[j+1],_dx[j+1], _t[i], _knots[j],_knots[j+1]);
      continue;
    }

    idx+=idxnew;

    eval(_xout[i],_dxout!=0 ? _dxout+i : 0,
         _x[idx],_dx[idx],_x[idx+1],_dx[idx+1],
         _t[i], _knots[idx],_knots[idx+1]);
  }
}

/** Sample Hermite interpolating cubic C^1 spline at parameter values _t.
    \ingroup vc_hermite_int
    Extrapolation for any parameters  out of range [knots[0],_knots[_nk-1]].
    \param[out] _xout write _nt position samples
    \param[out] _dxout _nt write tangent samples (no output if _dxout==0)
    \param _x _nk positions interpolated at _knots
    \param _dx _nk tangents interpolated at _knots
    \param _knots knot vector with _kn entries
    \param _nk >=2  number of _knots (and _x,_dx entries)
    \param _t array of parameter values for sampling of size _nt (should be sorted)
    \param _nt number of sample points (size of _t and _xout, _dxout)
    Uses find_interval_forward() or find_interval_backward() to find knot
    intervals.
 */
template <typename vec>
void eval_at(vec* _xout,vec* _dxout,
             const vec* _x,const vec* _dx,
             const typename vec::value_type* _knots,size_t _nk,
             const typename vec::value_type* _t,size_t _nt) {
  assert(_nk>=2);

  if (_knots[0]<_knots[1])
    eval_at(_xout,_dxout,_x,_dx,_knots,_nk,
            find_interval_forward_t<typename vec::value_type>(),_t,_nt);
  else {
    assert(_knots[0]!=_knots[1]);
    eval_at(_xout,_dxout,_x,_dx,_knots,_nk,
            find_interval_backward_t<typename vec::value_type>(),_t,_nt);
  }
}

/** Sample interval [_ta,_tb] uniformly.
    \ingroup vc_hermite_int
    Linear interpolation between _ta and _tb (including bounds) with uniformly
    distributed samples.
    \param[out] _t write _n values
    \param _ta lower interval bound
    \param _tb upper interval bound, _ta<=_tb
    \param _n >=2 number of samples
 */
template <typename real>
void linspace(real* _t,real _ta,real _tb,size_t _n) {
  assert(_n>=2);

  real rn=real(_n-1);

  for (real ri=real(0);ri<=rn;ri+=real(1)) {
    real s=ri/rn;
    *_t++=_ta*(real(1)-s)+_tb*s;
  }
}

/** Sample Hermite interpolating cubic C^1 spline in linspace() spanned by [_ta,_tb]
    \ingroup vc_hermite_int
    Extrapolation for any parameters  out of range [knots[0],_knots[_nk-1]].
    \tparam FIND functor to find knot intervals; the default setting
    LinearSearch uses a specialization with linear search using `FIND=`
    find_interval_forward() or find_interval_backward() depending on the
    bounds of the first knot interval
    \param[out] _xout write _nt position samples
    \param[out] _dxout _nt write tangent samples (no output if _dxout==0)
    \param _x _nk positions interpolated at _knots
    \param _dx _nk tangents interpolated at _knots
    \param _knots knot vector with _kn entries
    \param _nk number of _knots (and _x,_dx entries)
    \param _ta <= knots[0] lower bound of sample interval
    \param _tb <= _knots[_nk-1] upper bound of sample interval, _tb>=_ta
    \param _nt >=2 number of samples (including _ta and _tb)
 */
template <typename vec,typename FIND>
void
eval_linspace(vec* _xout,vec* _dxout,
              const vec* _x,const vec* _dx,
              const typename vec::value_type* _knots,size_t _nk,
              FIND _find,
              typename vec::value_type _ta,typename vec::value_type _tb,size_t _nt) {

  typedef typename vec::value_type real;

  assert(_nk>=2);
  assert(_nt>=2);

  size_t k=0;
  int    idx=0;

  real rn=real(_nt-1);

  for (real ri=real(0);ri<=rn;ri+=real(1),++k) {

    real s=ri/rn;
    real t=_ta*(real(1)-s)+_tb*s;

    int idxnew=_find(_knots+idx,_nk-idx,t);

    if (idxnew<0 || idx+idxnew+1>=int(_nk)) {
      // extrapolate
      int j=idxnew<0 ? 0 : int(_nk-2);
      eval(_xout[k],_dxout!=0 ? _dxout+k : 0,
           _x[j],_dx[j],_x[j+1],_dx[j+1], t, _knots[j],_knots[j+1]);
      continue;
    }

    idx+=idxnew;

    eval(_xout[k],_dxout!=0 ? _dxout+k : 0,
         _x[idx],_dx[idx],_x[idx+1],_dx[idx+1], t,_knots[idx],_knots[idx+1]);

  }
}

/** Sample Hermite interpolating cubic C^1 spline in linspace() spanned by [_ta,_tb]
    \ingroup vc_hermite_int
    Extrapolation for any parameters  out of range [knots[0],_knots[_nk-1]].
    \param[out] _xout write _nt position samples
    \param[out] _dxout _nt write tangent samples (no output if _dxout==0)
    \param _x _nk positions interpolated at _knots
    \param _dx _nk tangents interpolated at _knots
    \param _knots knot vector with _kn entries
    \param _nk number of _knots (and _x,_dx entries)
    \param _ta <= knots[0] lower bound of sample interval
    \param _tb <= _knots[_nk-1] upper bound of sample interval, _tb>=_ta
    \param _nt >=2 number of samples (including _ta and _tb)
    Uses find_interval_forward() or find_interval_backward() to find
    knot intervals.
 */
template <typename vec>
void
eval_linspace(vec* _xout,vec* _dxout,
              const vec* _x,const vec* _dx,
              const typename vec::value_type* _knots,size_t _nk,
              typename vec::value_type _ta,typename vec::value_type _tb,
              size_t _nt) {
  assert(_nk>=2);

  if (_knots[0]<_knots[1])
    eval_linspace(_xout,_dxout,_x,_dx,_knots,_nk,
                  find_interval_forward_t<typename vec::value_type>(),
                  _ta,_tb,_nt);
  else {
    assert(_knots[0]!=_knots[1]);
    eval_linspace(_xout,_dxout,_x,_dx,_knots,_nk,
                  find_interval_backward_t<typename vec::value_type>(),
                  _ta,_tb,_nt);
  }
}

/** Approximate arc lengths.
    \ingroup vc_hermite_int
    Approximates arc length of parameter intervals as distance between
    control points and accumulates lengths in _s.
    \param[out] _s write _n accumulated arc length values sum(||x[i]-x[i-1]||)+_s0
    (ignored if _s==0)
    \param _work temporary storage for a single vector
    \param _x _n>1 control points (positions interpolated at _s[i])
    \param _n number of control points and parameter values
    \param _sa lower bound/offset to _s
    \return arc length
 */
template <typename vec>
typename vec::value_type
accum_length(typename vec::value_type* _s,vec& _work,const vec* _x,size_t _n,
             typename vec::value_type _sa=typename vec::value_type(0)) {

  typedef typename vec::value_type real;

  assert(_n>1);
  real   s=_sa,sum(0);

  if (_s!=0)
    _s[0]=s;

  for (size_t i=1;i<_n;++i) {
    // _work=_x[i]-_x[i-1];
    _work=_x[i];
    _work-=_x[i-1];
    sum+=norm2(_work);
    if (_s!=0)
      _s[i]=sum;
  }
  return sum-_sa;
}

/** Approximate arc lengths.
    \ingroup vc_hermite_int
    Approximates arc length of parameter intervals as distance between
    control points and accumulates lengths in _s.
    <b>Constructs \c vec object for temporary storage!</b>
    \param[out] _s write _n accumulated arc length values sum(||x[i]-x[i-1]||)+_s0
    (ignored if _s==0)
    \param _x _n>1 control points (positions interpolated at _s[i])
    \param _n number of control points and parameter values
    \param _sa lower bound/offset to _s
    \return arc length
 */
template <typename vec>
typename vec::value_type
accum_length(typename vec::value_type* _s,const vec* _x,size_t _n,
             typename vec::value_type _sa=typename vec::value_type(0)) {
  vec work;
  return accum_length(_s,work,_x,_n,_sa);
}

/** Get lower and upper bound on arc length of cubic Hermite interpolant over [_ta,_tb].
    \param _lower bound difference _x1-_x2
    \param _upper bound as length of control polygon of cubic Bezier curve
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param[out] _work temporary storage for two vectors
    \sa length()
 */
template <typename vec>
void
length_bounds(typename vec::value_type& _lower,typename vec::value_type& _upper,
              const vec& _x0,const vec& _dx0,
              const vec& _x1,const vec& _dx1,
              typename vec::value_type _ta,typename vec::value_type _tb,vec _work[2]) {

  typedef typename vec::value_type real;

  _work[0]=_x1;
  _work[0]-=_x0;
  _lower=norm2(_work[0]);

  real h=_tb-_ta;

  _work[0]=_dx0;
  _work[0]*=(h/real(3));  // b1-b0

  real s=norm2(_work[0]);

  _work[1]=_dx1;
  _work[1]*=(-h/real(3)); // b2-b3
  s+=norm2(_work[1]);

  _work[1]+=_x1;
  _work[1]-=_x0;
  _work[1]-=_work[0];
  s+=norm2(_work[1]);

  _upper=s;
}

/** Get lower and upper bound on arc length of cubic Hermite interpolant over [_ta,_tb].
    <b>Constructs \c vec objects for temporary storage!</b>
    \param _lower bound difference _x1-_x2
    \param _upper bound as length of control polygon of cubic Bezier curve
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \sa length()
 */
template <typename vec>
void length_bounds(typename vec::value_type& _lower,typename vec::value_type& _upper,
                   const vec& _x0,const vec& _dx0,
                   const vec& _x1,const vec& _dx1,
                   typename vec::value_type _ta,typename vec::value_type _tb) {
  vec work[2];
  length_bounds(_lower,_upper,_x0,_dx0,_x1,_dx1,_ta,_tb,work);
}

/** Approximate arc length h([_ta,_t]) of cubic Hermite interpolant over [_ta,_tb].
    \ingroup vc_hermite_int
    Compute length from linear interpolation using _k+2 samples.
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _t parameter value
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param[out] _work temporary storage for two vectors, on output _work[0] equals value x(_t)
    \param _k number of samples in interval
    \return length
    \sa accum_arclength()
*/
template <typename vec>
typename vec::value_type
length(const vec& _x0,const vec& _dx0,
       const vec& _x1,const vec& _dx1,
       typename vec::value_type _t,
       typename vec::value_type _ta,
       typename vec::value_type _tb,vec _work[2],size_t _k) {

  typedef typename vec::value_type real;

  if (_k==0) {
    _work[0]=_x1;
    _work[0]-=_x0;
    return norm2(_work[0]);
  }
  real s(0);
  int  iw=0;

  real rn=real(_k+1);
  real h=_tb-_ta, delta=_t-_ta;


  _work[iw]=_x0;

  for (real ri=real(1);ri<=rn;ri+=real(1)) {
    real t=ri/rn*delta;
    int  jw=(iw+1)%2;

    eval(_work[jw],(vec*) 0,_x0,_dx0,_x1,_dx1,t,real(0),h);
    _work[iw]-=_work[jw];
    s+=norm2(_work[iw]);

    iw=jw;
  }

  if (iw!=0)
    _work[0]=_work[1];

  return s;
}

/** Approximate arc length h([_ta,_t]) of cubic Hermite interpolant over [_ta,_tb].
    \ingroup vc_hermite_int
    Compute length from linear interpolation using _k+2 samples.
    <b>Constructs \c vec objects for temporary storage!</b>
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _t parameter value
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param _k number of samples in interval
    \return length
    \sa accum_arclength()
*/
template <typename vec>
typename vec::value_type
length(const vec& _x0,const vec& _dx0,
       const vec& _x1,const vec& _dx1,
       typename vec::value_type _t,
       typename vec::value_type _ta,typename vec::value_type _tb,size_t _k) {
  vec work[2];
  return length(_x0,_dx0,_x1,_dx1,_t,_ta,_tb,work,_k);
}

/** Approximate arc length of cubic Hermite interpolant.
    \ingroup vc_hermite_int
    Compute length from linear interpolating using _k+2 samples.
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param _work temporary storage for two vectors
    \param _k number of samples in interval
    \return length
    \sa accum_arclength()
*/
template <typename vec>
typename vec::value_type
length(const vec& _x0,const vec& _dx0,
       const vec& _x1,const vec& _dx1,
       typename vec::value_type _ta,typename vec::value_type _tb,
       vec _work[2],size_t _k) {
  return length(_x0,_dx0,_x1,_dx1,_tb,_ta,_tb,_work,_k);
}

/** Approximate arc length of cubic Hermite interpolant.
    \ingroup vc_hermite_int
    Compute length from linear interpolating using _k+2 samples.
    <b>Constructs \c vec objects for temporary storage!</b>
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param _k number of samples in interval
    \return length
    \sa accum_arclength()
*/
template <typename vec>
typename vec::value_type
length(const vec& _x0,const vec& _dx0,
       const vec& _x1,const vec& _dx1,
       typename vec::value_type _ta,typename vec::value_type _tb,size_t _k) {
  vec work[2];
  return length(_x0,_dx0,_x1,_dx1,_ta,_tb,work,_k);
}

/** Accumulate arc length of polynomial pieces.
    \ingroup vc_hermite_int
    Approximates arc length of polynomial pieces from interpolation at
    _k samples (using length()) and accumulates length in _s.
    \param[out] _s write _n accumulated values sum(length(PIECE[i]))+_sa
    \param _work temporary storage for a single vector
    \param _x _n interpolation points (at _s[i] for output)
    \param _dx _n tangents at _s[i]
    \param _knots original knot vector (_n knots)
    \param _n >1 number of control points and parameter values
    \param _k >0 number of samples in each knot interval (parameter to length())
    \param _sa lower bound/offset to _s
    \return arc length
    \sa accum_length()
 */
template <typename vec>
typename vec::value_type
accum_arclength(typename vec::value_type* _s,
                vec _work[2],
                const vec* _x,const vec* _dx,
                const typename vec::value_type* _knots,size_t _n,size_t _k=3,
                typename vec::value_type _sa=typename vec::value_type(0)) {

  typedef typename vec::value_type real;

  if (_k==0)
    return accum_length(_s,_work[0],_x,_n,_sa);

  assert(_n>1);
  assert(_k>0);

  real   s=_sa;

  _s[0]=s;
  for (size_t i=1;i<_n;++i)
    _s[i]=_s[i-1]+
      length(_x[i-1],_dx[i-1],_x[i],_dx[i],_knots[i-1],_knots[i],_work,_k);

  return _s[_n-1]-_sa;
}

/** Accumulate arc length of polynomial pieces.
    \ingroup vc_hermite_int
    Approximates arc length of polynomial pieces from interpolation at
    _k samples (using length()) and accumulates length in _s.
    \param[out] _s write _n accumulated values sum(length(PIECE[i]))+_sa
    \param _x _n interpolation points (at _s[i] for output)
    \param _dx _n tangents at _s[i]
    \param _knots original knot vector (_n knots)
    \param _n >1 number of control points and parameter values
    \param _k >0 number of samples in each knot interval (parameter to length())
    \param _sa lower bound/offset to _s
    \return arc length
    \sa accum_length()
 */
template <typename vec>
typename vec::value_type
accum_arclength(typename vec::value_type* _s,
                const vec* _x,const vec* _dx,
                const typename vec::value_type* _knots,size_t _n,size_t _k=3,
                typename vec::value_type _sa=typename vec::value_type(0)) {

  vec work[2];
  return accum_arclength(_s,work,_x,_dx,_knots,_n,_k,_sa);
}

//-----------------------------------------------------------------------------

/** Helper to subdivide().
    \ingroup vc_hermite_int
    Uses Bernstein-Bezier form and de Casteljau's algorithm for recursive
    (mid-interval) subdivision.
    Subdivides until
    \code
    lower<_tol && (upper-lower)<_tol*typename vec::value_type(0.5)
    \endcode
    where \c upper and \c lower as from length_bounds().
    <b>Constructs \c vec objects for temporary storage!</b>
 */
template <typename vec,typename output_iter>
output_iter _subdivide(typename vec::value_type _tol,
                       output_iter _xout,const vec _b[4],
                       size_t _maxdepth) {

  typedef typename vec::value_type real;

  real lower=norm2(_b[3]-_b[1]);
  real upper=norm2(_b[1]-_b[0])+norm2(_b[2]-_b[1])+norm2(_b[3]-_b[2]);

  if (_maxdepth==0) {
#ifndef __CUDACC__
    assert(!"reached maximum recursion level");
#endif

    *_xout=_b[0];
    ++_xout;
    return _xout;
  }

  if (lower<_tol && (upper-lower)<_tol*real(0.5)) {
    *_xout=_b[0];
    ++_xout;
    return _xout;
  }

  vec b10=(_b[0]+_b[1])*real(0.5);
  vec b11=(_b[1]+_b[2])*real(0.5);
  vec b12=(_b[2]+_b[3])*real(0.5);

  vec b20=(b10+b11)*real(0.5);
  vec b21=(b11+b12)*real(0.5);

  vec b30=(b20+b21)*real(0.5);

  vec c[4];
  c[0]=_b[0]; c[1]=b10; c[2]=b20; c[3]=b30;

  _xout=_subdivide(_tol,_xout,c,_maxdepth-1);

  c[0]=b30;   c[1]=b21; c[2]=b12; c[3]=_b[3];

  _xout=_subdivide(_tol,_xout,c,_maxdepth-1);

  return _xout;
}

/** Subdivide for sampling positions.
    \ingroup vc_hermite_int
    Subdivide polynomial until tolerance _tol is reached (see _subdivide()).
    The output is written to _xout. The output does \a not include the final
    position _x1!
    <b>Uses de Casteljau's algorithm and constructs potentially many(!)
    temporary \c vec objects!</b>
    \param _tol tolerance on maximum segment length and maximum length error,
    see _subdivide()
    \param _xout store samples \a not including the final point _x1
    \param _x0 value at t=_ta
    \param _dx0 derivative d/dt at t=_ta
    \param _x1 value at t=_tb
    \param _dx1 derivative d/dt at t=_tb
    \param _ta lower bound of parameter interval
    \param _tb upper bound of parameter interval
    \param _maxdepth restrict recursion depth
 */
template <typename vec,typename output_iter>
output_iter subdivide(typename vec::value_type _tol,output_iter _xout,
                      const vec& _x0,const vec& _dx0,const vec& _x1,const vec& _dx1,
                      typename vec::value_type _ta,typename vec::value_type _tb,
                      size_t _maxdepth=32) {

  typedef typename vec::value_type real;

  assert(_tol>0);
  assert(_ta<_tb);
  real h=_tb-_ta;
#ifndef __CUDACC__
  assert(std::isfinite(h));
#else
  assert(isfinite(h));
#endif

  vec  b[4];
  b[0]=_x0;
  b[1]=b[0]+_dx0*(h/real(3));
  b[3]=_x1;
  b[2]=b[3]-_dx1*(h/real(3));

  return _subdivide(_tol,_xout,b,_maxdepth);
}

/** Subdivide for sampling positions.
    \ingroup vc_hermite_int
    Subdivide spline until tolerance _tol is reached (see _subdivide()).
    The output is written to _xout. The output does \a not include the final
    position _x1!
    <b>Uses de Casteljau's algorithm and constructs potentially many(!)
    temporary \c vec objects!</b>
    \param _tol tolerance on maximum segment length and maximum length error,
    see _subdivide()
    \param _xout store samples \a not including the final point _x1
    \param _x _n interpolation points (at _s[i] for output)
    \param _dx _n tangents at _s[i]
    \param _knots original knot vector (_n knots)
    \param _n >1 number of control points and parameter values
    \param _maxdepth restrict recursion depth
 */
template <typename vec,typename output_iter>
output_iter subdivide(typename vec::value_type _tol,output_iter _xout,
                      const vec* _x,const vec* _dx,
                      const typename vec::value_type* _knots,size_t _n,
                      size_t _maxdepth=32) {
  assert(_tol>0);
  assert(_n>1);
  for (size_t i=1;i<_n;++i)
    _xout=subdivide(_tol,_xout,
                    _x[i-1],_dx[i-1],_x[i],_dx[i],_knots[i-1],_knots[i],
                    _maxdepth);

  return _xout;
}

//-----------------------------------------------------------------------------

} // namespace hermite_int
} // namespace math
} // namespace VC

//-----------------------------------------------------------------------------


//=============================================================================

#endif // VC_MATH_HERMITEINT_HH
