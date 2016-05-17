//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MATH_ODE_HH
#define VC_MATH_ODE_HH

#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

#include "hermite_int.hh"

//=============================================================================

namespace VC {
namespace math {

/** \defgroup vc_ode Ordinary differential equations (ODE)
    Solve initial value problems for ordinary differential equations (ODE)
    in the form \f[ \frac{\partial}{\partial t} y = f(t,y)\f].
    \ingroup vc_math
 */

/** \example vclibs/math/example_rk43.cc
    \ingroup vc_ode
    Shows how to configure and use VC::math::RK43.
 */


/// Solve ordinary differential equations  [\ref vc_ode] \ingroup vc_ode
namespace ode {

/** State of evaluation in integrator (e.g., RK43) or user defined evaluator.
    Describes integration state (RK43::eval_state()) but is also used as
    exceptions that are thrown on particular events (e.g., OutOfDomain by
    evaluator).

    For exceptions, the evaluator instance (RK43::evaluator_t) throws
    ArithmeticError, OutOfDomain, or ForceStop.

    - Generally, the integrator object (RK43) catches any exception
      `state`, stops integration, sets the integration state to
      `state` and may rethrows `state` (e.g., RK43::initialize() and
      RK43::step(), in contrast the high level routine
      RK43::integrate() always returns `state` and does not throw).

    - OutOfDomain is not rethrown. Instead the integrator decreases
      the step size until the boundary is hit, i.e., current step size
      <= RK43::hmin, and the final state HitBoundary is thrown. (The
      evaluator object may throw HitBoundary as well in
      base_evaluator::output(), then the integrator assumes the
      boundary was hit with sufficient accuracy.)

    StateStr [EvalState] provides a textual description of a state.

    \sa RK43, StateStr
 */
enum EvalState {
  Undefined=-1,
  Success=0,       //!< integration (step) successful
  InvalidOption,   //!< invalid option/parameter (e.g., missing maximum step size)
  ArithmeticError, //!< evaluation fails (arithmetic yields undefined value)
  OutOfDomain,     //!< evaluation out of domain bounds fails
  HitBoundary,     //!< integration (step) stopped at domain boundary
  CriticalPoint,   //!< integration (step) stopped in critical point
  ForceStop ,      //!< integration was forced to stop
  StepUpdate       //!< internal state: repeat last step with updated h
};

/// textual description of EvalState (e.g., StateStr[OutOfDomain])
extern const char** StateStr;

//-----------------------------------------------------------------------------

/** Evaluate dy and communicate with ODE solver.
    \ingroup vc_ode
    \tparam real_t scalar type (time t)
    \tparam vec_t vector type (value y and derivative dy)
    \sa RK43
 */
template <typename real,typename vec>
struct base_evaluator {
  typedef real real_t; //!< scalar
  typedef vec  vec_t;  //!< vector

  /** Initialize vector _x to have same dimensions as _template.
      Default implementation executes _x=_template.
   */
  void initialize_vector(vec_t& _x,const vec_t& _template) {
    _x=_template;
  }

  /** Initialize solving: called by ODE_SOLVER::initialize().
      \sa RK43, RK43::initialize()
      This method may throw
      - ArithmeticError
      - ForceStop
      - OutOfDomain
      Default implementation is empty.
  */
  template <typename ODE_SOLVER>
  void initialize(real_t _t,const vec_t& _y,ODE_SOLVER*)
    throw(VC::math::ode::EvalState);

  /** Evaluate `d/dy f(_t,_y) and store result in `_dy.`
      This method may throw
      - ArithmeticError
      - OutOfDomain
      Default implementation is undefined.
  */
  void dy(real_t _t,const vec_t& _y,vec_t& _dy) throw(VC::math::ode::EvalState);

  /** Output after successful step.

      This method is
      - first called after initialize(), in this case there was
        no step yet and `ODE_SOVER::t()==ODE_SOVER::tprev()`.
      - called before storing (Solution::push()) results of the step,
        results are stored afterward regardless of a throw.

      The solver's state may be Success, CriticalPoint or HitBoundary.
      This method may throw
      - ForceStop
      - CriticalPoint
      - HitBoundary
   */
  void output(real_t _t,const vec_t& _y,const vec_t& _dy) {
    use_nowarn(_t); use_nowarn(_y); use_nowarn(_dy);
  }
};

//-----------------------------------------------------------------------------

/** Solution to ODE.
    \ingroup vc_ode
    \tparam real_t scalar type (time t)
    \tparam vec_t vector type (value y and derivative dy)
    Note: Use reserve() for optimal memory performance.
    \sa RK43
 */
template <typename real,typename vec>
struct Solution {
  typedef real real_t; //!< scalar
  typedef vec  vec_t;  //!< vector
  typedef Solution<real,vec> self_t; //!< self

  std::vector<real_t> t;  //!< time
  std::vector<vec_t>  y;  //!< value/position
  std::vector<vec_t>  dy; //!< derivative/direction

  /// push triple `(_t,_y,_dy)` from integration step
  void push(real_t _t,const vec_t& _y,const vec_t& _dy) {
    t.push_back(_t); y.push_back(_y); dy.push_back(_dy);
  }
  /// pop last triple
  void pop() {
    t.pop_back(); y.pop_back(); dy.pop_back();
  }
  /// clear solution
  void clear() { t.clear(); y.clear(); dy.clear();  }

  /// get size (same as t.size() or y.size(); dy.size() may equal 0)
  size_t size() const {
    assert(t.size()==y.size());
    assert(t.size()==dy.size() || dy.size()==0);
    return t.size();
  }

  /// Is solution empty?
  bool empty() const { return size()==0; }

  /// reserve storage for _n steps
  void reserve(unsigned _n) {
    t.reserve(_n); y.reserve(_n); dy.reserve(_n);
  }

  /// swap solutions
  void swap(Solution<real_t,vec_t>& _other) {
    t.swap(_other.t); y.swap(_other.y); dy.swap(_other.dy);
  }

  /// reverse solution (e.g., after backwards integration)
  void reverse() {
    std::reverse(t.begin(),t.end());
    std::reverse(y.begin(),y.end());
    std::reverse(dy.begin(),dy.end());
  }

  /** Append an _other solution.
      \param _other will be appended
      \param _toffset time offset will be subtracted from `_other.t`
      \param _ibegin start index
  */
  void append(const Solution& _other,
              real_t _toffset=real_t(0),int _ibegin=0) {
    int n=_other.t.size();
    if (n<=_ibegin) return;
    int m=this->t.size();
    this->t.reserve(m+n);
    this->y.reserve(m+n);
    this->dy.reserve(m+n);

    for (int i=_ibegin;i<n;++i)
      push(_other.t[i]-_toffset,_other.y[i],_other.dy[i]);
  }
  /** Reverse _other and append() it.
      Reverse orientation of `_other.dy`
      \param _other will be appended
      \param _toffset time offset will be subtracted from `_other.t`
  */
  void append_reverse(const Solution& _other,real_t _toffset=real_t(0)) {
    int n=_other.t.size();
    if (n==0) return;
    int m=this->t.size();
    this->t.reserve(m+n);
    this->y.reserve(m+n);
    this->dy.reserve(m+n);

    for (int i=n-1;i>=0;--i) {
      push(_other.t[i]-_toffset,_other.y[i],_other.dy[i]);
      dy.back()*=real_t(-1);
    }
  }

  /** @name Evaluation using C^1 cubic Hermite interpolation

      The given methods evaluate the current solution to obtain a
      resampled or refined version using function in namespace
      VC::math::hermite_int.

      - Note that cubic Hermite interpolation is used *regardless* of
        the ODE solver which produced the solution and which might
        show higher (or lower) approximation order!

      - Cubic interpolation is appropriate for the RK43 solver.

      Some methods evaluate the derivatives dy optionally (or not at
      all).

      \sa [\ref vc_hermite_int]

      @{
   */


  /** Evaluate at given _t.
      \param[in] _sol solution from RK43
      \param[in,out] _rsol destination for refined solution, input `_rsol.t`
      specifies interpolation points positions are output to `_rsol.y` and
      `_rsol.dy` (see below)
      \param _dy no output to `_rsol.dy` if `false` (then
      `_rsol.dy.size()==0` on exit)
   */
  void eval_at(self_t& _rsol,bool _dy=true) {
    if (size()==0) {
      _rsol.clear();
      return;
    }

    size_t n=_rsol.t.size();
    _rsol.y.resize(n);
    _rsol.dy.resize(_dy ? n : size_t(0));

    hermite_int::eval_at
      (&_rsol.y.front(),_rsol.dy.empty() ? (vec_t*) 0 : &_rsol.dy.front(),
       &y.front(),&dy.front(),&t.front(),size(),
       &_rsol.t.front(),n);
  }

  /** Refine solution: evaluate in hermite_int::linspace(_t0,_t1,_n).
      \param[in] _sol solution from RK43

      \param[out] _rsol destination for refined solution computes
      `_rsol.y` and `_rsol.`t and `_rsol.dy` (see below)
      \param _t no output to `_rsol.t` if `false` (then
      `_rsol.t.size()==0` on exit)
      \param _dy no output to `_rsol.dy` if `false` (then
      `_rsol.dy.size()==0` on exit)
   */
  void eval_in(real_t _t0,real_t _t1,size_t _n,self_t& _rsol,
               bool _t=true,bool _dy=true) {
    if (size()==0 || _n==0 || _t0==_t1) {
      _rsol.clear();
      return;
    }

    _rsol.y.resize(_n);
    _rsol.t.resize(_t ? _n : size_t(0));
    _rsol.dy.resize(_dy ? _n : size_t(0));

    if (_rsol.t.empty()) {
      hermite_int::eval_linspace
        (&_rsol.y.front(),_rsol.dy.empty() ? (vec_t*) 0 : &_rsol.dy.front(),
         &y.front(),&dy.front(),&t.front(),size(),
         _t0,_t1,_n);
    }
    else {
      hermite_int::linspace(&_rsol.t.front(),_t0,_t1,_n);
      eval_at(_rsol,_dy);
    }
  }

  /** Refine solution: evaluate in hermite_int::linspace(tmin,tmax,_n).
      Here, tmin, tmax are t.front() and t.back(), respectively.
      \param[in] _sol solution from RK43
      \param[out] _rsol destination for refined solution computes `_rsol.y` and
      `_rsol.t` and `_rsol.dy` (see below)
      \param _t no output to `_rsol.t` if `false` (then
      `_rsol.t.size()==0` on exit)
      \param _dy no output to `_rsol.dy` if `false` (then
      `_rsol.dy.size()==0` on exit)
  */
  void eval(size_t _n,self_t& _rsol,bool _t=true,bool _dy=true) {
    if (size()==0 || _n==0)
      _rsol.clear();
    else
      eval_in(t.front(),t.back(),_n,_rsol,_t,_dy);
  }

  /** Refine solution by recursive subdivision.
      Calls VC::math::hermite_int::subdivide().
      \param _tol arc length tolerance (see VC::math::hermite_int::subdivide())
      \param _rsol output iterator for positions y *including* final point
      \return output iterator
   */
  template <typename output_iter>
  output_iter
  subdivide(real_t _tol,output_iter _rsol) {
    if (size()==0)
      return _rsol;

    output_iter ii=
      hermite_int::subdivide(_tol,_rsol,&y.front(),&dy.front(),&t.front(),size());
    *ii=y.back();
    ++ii;

    return ii;
  }

  /// @}

};

/// Test functions \ingroup vc_ode
namespace TestFunctions {
  /** For A1,A2,A4,F4 see
      <em>W.H.Enright et al. Interpolant for Runge-Kutta Formulas.
      ACM Trans. Math. Soft., vol 12, no. 3, 1986, 193--218</em>
   */
  enum Id {
    A1, A2, A4, F4, P123, P234, P345, T2, T3, LINEAR, SINK, NumFunctions
  };

  /// get name (or 0)
  const char* get_name(Id _id);
  /// query problem dimension
  int get_dim(Id _id);
  /// query initial value
  void get_y0(Id _id,double* _y);
  /// get derivative
  void get_dy(Id _id,double _t,const double* _y,double* _dy);
  /// get ground truth
  void get_y(Id _id,double _t,double* _y);

} // namespace TestFunctions

//-----------------------------------------------------------------------------

} // namespace ode
} // namespace math
} // namespace VC

//-----------------------------------------------------------------------------


//=============================================================================

#endif // VC_MATH_ODE_HH
