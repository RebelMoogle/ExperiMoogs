//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MATH_ODE_RK43_HH
#define VC_MATH_ODE_RK43_HH

#include <cmath>

#include "../system/platform.hh" // alloca
#include "../base/debug.hh"

#include "ode.hh"
#include "hermite_int.hh"

#include <cassert>
#include <vector>
#include <limits>
#include <iostream>

//=============================================================================

namespace VC {
namespace math {

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4290)
# endif

namespace ode {

/** Options to RK43.
    \ingroup vc_ode
*/
template <typename real>
struct rk43_options_t {
  typedef real value_type;
  typedef real real_t;

  /// default constructor (values see below)
  rk43_options_t() : hmax(0.0), h0(0.0), rsmin(0.0),
                     atol(1e-6), rtol(1e-3), rho(0.8) {}

  /// maximum step size >0 (must be set: the default value is 0, which is \a invalid)
  real hmax;
  /// initial step size (determined automatically in <=0.0, default is 0.0)
  real h0;
  /** Relative minimum step size/minimum norm of dy.
      rsmin>=0 determines state changes to HitBoundary (decreased step size)
      or CriticalPoint (norm(dy)/max(norm(y),atol/rtol)*h<rsmin).
      rsmin must be set: the default value is 0, which is \a invalid.
  */
  real rsmin;
  /// absolute tolerance >0 (controls adaptive step size, default value is 1e-6)
  real atol;
  /// relative tolerance >0 (controls adaptive step size, default value is 1e-3)
  real rtol;
  /// damping factor >0 for adaptive size (default value is 0.8)
  real rho;
};

/** Runge-Kutta integrator.
    \ingroup vc_ode

    4th order Runge-Kutta integrator with 3rd order error term for
    step size control.

    Implements solver described in
    \code
    Detlev Stalling and Hans-Christian Hege.
    Fast and resolution independent line integral convolution.
    SIGGRAPH 1995, pp. 249--256
    \endcode
    RK43 uses the same weights, but there are a couple of modifications
    included.

    The evaluator defines methods as in ode::base_evaluator for
    \arg initialization,
    \arg computing derivatives dy=f(t,y),
    \arg reporting results.
    These methods may throw exceptions of type EvalState to signal
    events.

    The evaluator state eval_state() is explained in detail with enum
    EvalState. States/events are communicated as C++ exceptions. This is
    hidden by integrate(). However, initialize() and step() may throw
    exceptions!

    RK43 can integrate \a backwards in time, for both integrate()
    (_t0>_t1) and _step() (_t1<0). In this case, orientation of
    derivatives from EVAL::dy() is reversed for integration (RK43
    takes care of that!). Note that RK43 passes the "original"
    derivatives (as from EVAL::dy()) to EVAL::output() and solution().

    \todo unit test for RK43 (compare to matlab ode45)
    \todo swap_move (can be implemented as swap(a,b) or a=b) -- specialization
    \todo vector initialization as global function with specialization?!
 */
template <typename EVAL>
class RK43 {
public:

  typedef EVAL                         evaluator_t; //!< evaluator
  typedef typename evaluator_t::real_t real_t;      //!< scalar type
  typedef typename evaluator_t::vec_t  vec_t;       //!< vector type

  typedef Solution<real_t,vec_t>       solution_t;  //!< solution

  typedef rk43_options_t<real_t>       options_t;   //!< options

  //// constructor sets default options (see options_t)
  RK43(const options_t& _options=options_t())
# ifndef DOXYGEN_SKIP
    : options(_options),
      m_eval(0),m_sol(0),
      m_t(0), m_t0(0),  m_h(0), m_hnew(0), m_signh(0),
      m_normy(-1), m_normdy(-1),
      m_eval_state(Undefined),
      m_nsteps(0), m_nevals(0), m_near_boundary(false)
# endif
  {}

  /// options
  options_t options;

  /// get current solution
  const solution_t& solution() const {
    assert(m_sol!=0);
    return *m_sol;
  }
  /// get current evaluator
  evaluator_t& evaluator() const {
    assert(m_eval!=0);
    return *m_eval;
  }
  /// get current state
  EvalState eval_state() const { return m_eval_state; }

  /** Caught OutOfDomain during last step().
      This flag indicates that we got "near" the domain boundary
      (where "near" is relative to the last step's magnitude).
      step() decreased the step size h() such that eval_state() is
      Success (or HitBoundary). Use this method if an early stop,
      potentially followed by extrapolation/bisection, is preferred
      to iteratively decreasing the step size. Note that the flag
      should be used only in combination with a (relative) step
      size threshold!
   */
  bool near_boundary() const { return m_near_boundary; }

  /// user sets updated steps size _hnew>0 before throwing StepUpdate
  void step_update(real_t _hnew) {
    m_hnew=_hnew; // read back after StepUpdate
    assert(_hnew>real_t(0));
  }

  /** Initialize integrator.
      Initialization calls _evaluator->dy() and may hence fail
      and \a throw an exception.
      \param _eval set evaluator object
      \param _y0 set initial value (i.e, y()==_y0 after call)
      \param _t0 set initail time (i.e., t()==_t0 after call)
      \param _solution set solution object, steps are added by
      Solution::push(); if _solution==0 then steps are \a not
      recorded
      \param _append don't call Solution::clear() if _true
      (no effect if _solution==0)
   */
  void initialize(evaluator_t* _eval,const vec_t& _y0,real_t _t0,
                  solution_t* _solution=0,bool _append=false)
    throw(EvalState);

  /// get current time
  real_t t() const { return m_t; }
  // get time of previous step
  real_t tprev() const { return m_t0; }
  /// get current value/position
  const vec_t& y() const { return m_v[RY]; }
  /// get current derivative/direction
  const vec_t& dy() const { return m_v[RDY]; }
  /// get current (absolute) step size
  real_t h() const { return m_h; }
  /// get current orientation (-1|+1,0=undefined state)
  real_t sign_h() const { return m_signh; }

  /** Interpolate value at _t from current state.
      This method does not call <tt>evaluator().dy()</tt>.
      It can be called safely from <tt>Evaluator::output()</tt>.
      \param _t interpolation point must be in interval [tprev(),t()]
      \param[out] _y store interpolated position
      \param[out] _dy store interpolated derivative (ignored if ==0)
      \return \c false on failure (extrapolation due to invalid _t or
      invalid eval_state())
   */
  bool interpolate(real_t _t,vec_t& _y,vec_t* _dy=0);

  /// get number of steps since initialize() was called
  int nsteps() const { return m_nsteps; }
  /** Get number of calls to evaluator().eval().
      Count is rest by initialize() (as for nsteps()). We count
      only calls to EVAL::eval() which did \a not thrown an exception!
  */
  int nevals() const { return m_nevals; }

  /** Take one step.
      Assumes initialize() was called. An exception is \a thrown on failure or
      certain events.
      \param _t1 maximum time, integrate \a backwards in time if _t1<t()
      \return true if _t1 was reached (decreases last step size to reach _t1)
      and false else
   */
  bool step(real_t _t1=std::numeric_limits<real_t>::infinity()) throw(EvalState);

  /** Integrate from _(_t0,_y0) to _t1.
      Calls initialize() and step() until _t1 is reached (within a tolerance)
      or boundary or critical point was reached.

      \arg Integrate \a backwards in time if _t0>_t1
      \arg Guesses \a valid options from parameters (if _t1<inf), see
      guess_options(), The old options are restored on return.
      \arg This method does \a not throw  exceptions!

      \param _eval set evaluator object
      \param _y0 set initial value (i.e, y()==_y0 after call)
      \param _t0 set initail time (i.e., t()==_t0 after call)
      \param _t1 maximum time
      \param _solution set solution object
      \param _append same as for initialize()
      \param _maxsteps ignored if ==0, otherwise return ForceStop
      if nsteps()>_maxsteps
      \return Success if _t1 was reached, the current state else.
      Note:
      \arg This method does not throw exeptions!
      \arg The use of the _evaluator, _y0, _t0, _solution parameters is as for
      initialize(). Time _t1 is used to control the loop calling step(_t1).
   */
  EvalState integrate(evaluator_t* _eval,const vec_t& _y0,real_t _t0,real_t _t1,
                      solution_t* _solution=0,bool _append=false,int _maxsteps=0);

  /** Guess default options and update _opts.
      This function is called by integrate() to provide valid options.
      Update such that
      \arg _opts.rsmin >= 1e3*eps(_t1-_t0)
      \arg _opts.hmax = (_t1-_t0)/10 if hmax<=0
      \arg _opts.h0=_opts.hmax/2 if h0<=0
  */
  static void guess_options(options_t& _opts,real_t _t0,real_t _t1);

  /** @name Refinement of solution.

      The solution of RK43 can be refined using a C^1 Hermite interpolation
      between pairs (y[i],dy[i]) and (y[i+1],dy[i+1]). We provide functions
      in namespace VC::math::hermite_int which are called by respective
      methods of RK43.

      The methods of VC::math::ode::Solution are appropriate for this
      purpose. They are wrapped by methods of RK43 below for
      convenience.

      \sa [\ref vc_hermite_int].
      @{
  */

  /// refine solution: evaluate at given _t \sa Solution::eval_at()
  RK43& eval_at(solution_t& _rsol,bool _dy=true) {
    assert(m_sol!=0);
    m_sol->eval_at(_rsol,_dy);
    return *this;
  }

  /// refine solution: evaluate in interval \sa Solution::eval_in()
  RK43& eval_in(real_t _t0,real_t _t1,size_t _n,
                solution_t& _rsol,bool _t=true,bool _dy=true) {
    assert(m_sol!=0);
    m_sol->eval_in(_t0,_t1,_n,_rsol,_t,_dy);
    return *this;
  }

  /// refine solution: evaluate in interval \sa Solution::eval()
  RK43& eval(size_t _n,solution_t& _rsol,bool _t=true,bool _dy=true) {
    assert(m_sol!=0);
    m_sol->eval(_n,_rsol,_t,_dy);
    return *this;
  }

  /// subdivide solution \sa  Solution::subdivide()
  template <typename output_iter>
  output_iter
  subdivide(real_t _tol,output_iter _rsol) {
    assert(m_sol!=0);
    return m_sol->subdivide(_tol,_rsol);
  }

  /// @}

private:
  RK43(const RK43&);

  /// determine initial step size from _y0, _dy0, and options (does not throw)
  real_t _initial_stepsize(const vec_t& _y0,const vec_t& _dy0) const;
  /** Core of step().
      May be iterated (called by _step1()). Does not throw or set m_eval_state.
      Keeps previous y,dx in RK3,RK3.
  */
  EvalState _step0(bool _update_h);
  /// core of step() (calls _step0(), does not throw)
  void _step1(bool _update_h);

  /// indices into m_v (define vector values state variables/temporaries)
  enum VReg {
    RY,RDY,
    RTMP,RTMP1,
    RK1,RK2,RK3,RK4,
    VRegEnd
  };

  evaluator_t*        m_eval;
  solution_t*         m_sol;

  vec_t               m_v[VRegEnd]; //!< vector valued state variables (see VReg)
  real_t              m_t,m_t0;     //!< time and previous steps' time
  real_t              m_h,m_hnew;   //!< absolute step size, step_update() sets m_hnew
  real_t              m_signh;      //!< sign of h +|-1 (direction)
  real_t              m_normy;      //!< norm of current y
  real_t              m_normdy;     //!< norm of current dy
  EvalState           m_eval_state; //!< current state
  int                 m_nsteps;     //!< number of steps since initialize())
  int                 m_nevals;     //!< number of dy evaluations sice initialize()
  bool                m_near_boundary; //!< set if OutOfDomain detected during step()
};

//-----------------------------------------------------------------------------

template <typename EVAL>
typename RK43<EVAL>::real_t
RK43<EVAL>::_initial_stepsize(const vec_t& /*_y0*/,const vec_t& /*_dy0*/) const {
  // "inspired" by Matlab ode45 (line 280 + odeargument)
  assert(m_normy>=real_t(0));
  assert(m_normdy>=real_t(0));
  real_t h=options.hmax;
  real_t threshold=options.atol/options.rtol;
  real_t rh=(m_normdy/std::max(m_normy,threshold))/
    (options.rho*sqrt(sqrt((options.rtol))));
  if (h*rh>real_t(1))
    h=real_t(1)/rh;

  if (options.h0>real_t(0))
    h=std::min(h,options.h0);

  if (h<options.rsmin) h=options.rsmin;
  if (h>options.hmax) h=options.hmax;

  return h;
}

//-----------------------------------------------------------------------------

template <typename EVAL>
void RK43<EVAL>::initialize(evaluator_t* _eval,const vec_t& _y0,real_t _t0,
                            solution_t* _solution,bool _append)
  throw(EvalState) {
  if (options.rsmin<=real_t(0) ||
      options.hmax    < real_t(0) || // will be set
      options.h0      < real_t(0) || // will be set
      options.atol <=real_t(0) ||
      options.rtol <=real_t(0) ||
      options.rho  <=real_t(0)) {
    VC_DBG_P(options);
    throw InvalidOption;
  }

  m_eval=_eval;
  m_sol=_solution;

  for (int i=0;i<VRegEnd;++i)
    m_eval->initialize_vector(m_v[i],_y0);

  m_eval->initialize(m_t=m_t0=_t0,_y0,this);

  if (m_sol!=0 && !_append)
    m_sol->clear(); // we don't m->sol.reserve() storage

  m_eval_state=Success;

  m_nsteps=m_nevals=0;

  m_v[RY]=_y0;
  m_normy=norm2(_y0);

  try {
    m_eval->dy(_t0,_y0,m_v[RDY]);
    ++m_nevals;
  }
  catch (EvalState state) {
    throw (m_eval_state=state);
  }
  m_normdy=norm2(m_v[RDY]);

  try {
    m_eval->output(_t0,_y0,m_v[RDY]);
  } catch (EvalState state) {
    assert(state==Success || state==ForceStop ||
           state==HitBoundary || state==CriticalPoint);
    m_eval_state=state;
  }

  if (m_sol!=0)
    m_sol->push(_t0,_y0,m_v[RDY]);

  m_h=_initial_stepsize(_y0,m_v[RDY]);
  m_signh=real_t(0); // undefined

  if (m_eval_state!=Success)
    throw m_eval_state;

  m_near_boundary=(m_eval_state==OutOfDomain);
}

//-----------------------------------------------------------------------------

template <typename EVAL>
bool RK43<EVAL>::interpolate(real_t _t,vec_t& _y,vec_t* _dy) {
  if (m_eval_state!=Success && m_eval_state!=HitBoundary &&
      m_eval_state!=CriticalPoint)
    return false;

  if (m_t0==m_t)
    return false;

  if (m_t0<m_t) {
    if (!(m_t0<=_t && _t<=m_t))
      return false;
  }
  else {
    if (!(m_t<=_t && _t<=m_t0))
      return false;
  }
  assert(fabs(m_t0-m_t)/(m_t0-m_t)*m_signh==real_t(1));

  hermite_int::eval(_y,_dy,
                    m_v[RK3],m_v[RK4],m_v[RY],m_v[RDY],_t,m_t0,m_t);

  return true;
}


//-----------------------------------------------------------------------------

template <typename EVAL>
EvalState RK43<EVAL>::_step0(bool _update_h) {
  assert(m_eval!=0);
  assert(m_eval_state==Success ||
         m_eval_state==OutOfDomain ||
         m_eval_state==StepUpdate);
  assert(real_t(0)<m_h && m_h<=options.hmax);

  // assume we have y(t) and dy(t,y) evaluated in m_v[RY] and m_v[RDY]

  assert(real_t(fabs(m_signh))==real_t(1));
  real_t h=m_h*m_signh;

  try {
    m_v[RTMP]=m_v[RY];        // Save for retries:
    m_v[RTMP1]=m_v[RDY];      // RY and RDY are modified only on Success.

    //
    // RK-step
    //

    m_v[RK1]=m_v[RTMP1]*h;    // k1=dy(t,y)*h (m_v[RDY] had "original" orientation
                              // assume we already have last evaluation point

    m_v[RTMP1]=m_v[RY];
    m_v[RTMP1]+=m_v[RK1]*real_t(0.5);
    m_eval->dy(m_t+h*real_t(0.5),m_v[RTMP1], m_v[RK2]);
    m_v[RK2]*=h;                                         // k2=dy(t+h/2,y+k1/2)*h
    ++m_nevals;

    m_v[RTMP1]=m_v[RY];
    m_v[RTMP1]+=m_v[RK2]*real_t(0.5);
    m_eval->dy(m_t+h*real_t(0.5),m_v[RTMP1], m_v[RK3]);
    m_v[RK3]*=h;                                         // k3=dy(t+h/2,y+k2/2)*h
    ++m_nevals;

    m_v[RTMP1]=m_v[RY];
    m_v[RTMP1]+=m_v[RK3];
    m_eval->dy(m_t+h,m_v[RTMP1], m_v[RK4]);
    m_v[RK4]*=h;                                         // k4=dy(t+h,y+k3)*h
    ++m_nevals;

    //m_v[RTMP]=m_v[RY]; // did not touch
    m_v[RTMP]+=m_v[RK1]/real_t(6);
    m_v[RTMP]+=m_v[RK2]/real_t(3);
    m_v[RTMP]+=m_v[RK3]/real_t(3);
    m_v[RTMP]+=m_v[RK4]/real_t(6);

  } catch (EvalState state) {
    assert(state==OutOfDomain || state==ArithmeticError);

    return state; // ===>
  }

  // RTMP contains new position y

  EvalState state=Success;

  try {
    m_eval->dy(m_t+h,m_v[RTMP], m_v[RTMP1]); // RTMP1 contains new derivative dy
    ++m_nevals;
  } catch (EvalState state) {
    assert(state==OutOfDomain || state==ArithmeticError);
    if (state==OutOfDomain) {
      return state; // retry with smaller h
    }

    // ==> LEAVE
  }

  real_t normynew=norm2(m_v[RTMP]);

  m_hnew=m_h;
  real_t threshold=options.atol/options.rtol;

  if (_update_h) {
    //
    // error control: decide whether we have to repeat the step
    //

    m_v[RK4]-=m_v[RTMP1]*h; // sign
    real_t delta=norm2(m_v[RK4])/real_t(6);

    // experimental: should have absolute and relative error
    delta/=std::max(threshold,std::max(m_normy,normynew));

#ifndef __CUDACC__
    if (!std::isfinite(delta)) {
#else
    if (!isfinite(delta)) {
#endif
      VC_DBG_TRACE("unflagged ArithmeticError");
      return ArithmeticError; // ===>
    }

    real_t hnew=
        std::min(options.hmax,
                 m_h*real_t(sqrt(sqrt(options.rho*options.rtol/delta))) );
    // explicit cast to real_t?!

    if (delta>options.rtol) {

      if (hnew!=m_h) {
        m_hnew=hnew;

        return StepUpdate; // ===>
      }
    }
    else
      m_hnew=std::min(hnew,options.hmax);

  } // adaptive step control

  //
  // We are done.
  //

  // --- rotating register indices?! ---

  m_v[RK3]=m_v[RY];             // previous result
  m_v[RK4]=m_v[RDY];

  m_v[RY]=m_v[RTMP];            // copy result from temporary
  m_v[RDY]=m_v[RTMP1];          // keep "original" orientation

                                //  (don't multiply by m_signh)
  real_t normy0=std::max(std::max(m_normy,normynew),threshold);
                                // potential problem is CP in 0

  m_normy=normynew;
  m_normdy=norm2(m_v[RDY]);

  assert(state==Success);

  if (normy0==real_t(0) || // avoid division by zero
      m_hnew==real_t(0) ||
      (m_hnew*m_normdy/normy0<=options.rsmin)) {
    //VC_DBG_P(normy0);
    state=CriticalPoint;
  }

  return state;
}

//-----------------------------------------------------------------------------

template <typename EVAL>
void RK43<EVAL>::_step1(bool _update_h) {
  assert(m_eval!=0);
  assert(m_eval_state==Success);
  assert(real_t(0)<m_h && m_h<=options.hmax);

  const int MAXTRIES=64;
  int       ntries;
  bool      inside=true;
  bool      updateh=_update_h;
  int       nevals0=m_nevals;
  EvalState state;

  m_near_boundary=false;

  for (ntries=MAXTRIES;ntries>0;--ntries) {

    state=_step0(updateh);

    switch (state) {

    case Success:          // proceed
    case HitBoundary:
    case ForceStop:
    case ArithmeticError:
      goto _FINISH_STEP;

    case CriticalPoint:    // proceed w/ HitBoundary
      if (!inside)
        state=HitBoundary; // stopped near boundary
      goto _FINISH_STEP;

    case StepUpdate:       // redo with new update step size
      m_eval_state=StepUpdate;
      m_h=m_hnew;
      updateh=false;
      //VC_DBG_P(m_h);
      continue;

    case OutOfDomain:      // redo with new update step size
      inside=false;
      updateh=false;
      m_near_boundary=true;
      m_h*=real_t(0.25);
      if (m_h*m_normdy/std::max(m_normy,options.atol/options.rtol)<options.rsmin) {
        state=HitBoundary;
        goto _FINISH_STEP;
      }
      else
        m_eval_state=OutOfDomain;
      continue;

    default:
      VC_DBG_P(state);
      assert("invalid state"!=0);
      break;
    }
  }

  _FINISH_STEP:

  if (ntries==0) {
    VC_DBG_TRACE("exceeded MAXTRIES (StepUpdate)");
    VC_DBG_P(nevals()-nevals0);
    use_nowarn(nevals0);

    assert(m_eval_state==StepUpdate || m_eval_state==OutOfDomain);
    state=Success; // we hope the best
  }

  assert(state!=OutOfDomain);
  assert(state!=StepUpdate);

  m_eval_state=state;
}

//-----------------------------------------------------------------------------

template <typename EVAL>
bool RK43<EVAL>::step(real_t _t1) throw(EvalState) {

  m_signh=(_t1>=m_t) ? real_t(+1) : real_t(-1);

  bool updateh=true;

 _REPEAT_STEP:

  assert(m_eval!=0);
  assert(m_signh>=real_t(0) ? (m_t<=_t1) : (m_t>=_t1));

  if (m_eval_state!=Success)
    throw m_eval_state;

  if (m_t==_t1)
    return true;

  _step1(updateh);

  real_t h=m_h*m_signh;

  switch (m_eval_state) {

  case Success:
    m_t0=m_t;
    if ((m_signh>real_t(0)) ? (m_t+h>=_t1) : (m_t+h<=_t1)) {
      // got too far: interpolate
      m_v[RK1]=m_v[RY];
      m_v[RK2]=m_v[RDY]; // RK3,RK4: y(t), dy(t,y); RK1,RK2: y(t+h), dy(t+h,y)

      hermite_int::eval
        (m_v[RY],&m_v[RDY],
         m_v[RK3],m_v[RK4], m_v[RK1],m_v[RK2], _t1, m_t,m_t+h);

      m_t=_t1; // and we are done!
    }
    else
      m_t+=h;

    m_h=m_hnew;
    break;

  case ForceStop:
  case ArithmeticError:
    return false;

  default:
    assert(m_eval_state==CriticalPoint ||
           m_eval_state==HitBoundary);
    m_t+=h;
  }

  ++m_nsteps;

  try {
    m_eval->output(m_t,m_v[RY],m_v[RDY]);
  } catch (EvalState state) {

    if (state==StepUpdate) {

      m_eval_state=Success;

      real_t hnew=m_hnew; // step_update() stores step here

      hnew=std::min(hnew,options.hmax);

      if (hnew!=m_h) {

        m_hnew=hnew;
        m_t-=h;
        m_v[RY]=m_v[RK3];
        m_v[RDY]=m_v[RK4];
        m_normy=norm2(m_v[RY]);
        m_normdy=norm2(m_v[RDY]);
        updateh=false; // forbid update

        goto _REPEAT_STEP;
      }
    }
    assert((state==Success ||
            state==ForceStop ||
            state==CriticalPoint ||
            state==HitBoundary) && "output() may not throw another state" != 0);

    if (m_eval_state==Success) // output() state overrides only Success
      m_eval_state=state;
  }

  if (m_sol!=0)
    m_sol->push(m_t,m_v[RY],m_v[RDY]);

  if (m_eval_state!=Success)
    throw m_eval_state;

  return m_signh>real_t(0) ? (m_t<_t1) : (m_t>_t1);
}

//-----------------------------------------------------------------------------

template <typename EVAL>
EvalState
RK43<EVAL>::integrate(evaluator_t* _eval,const vec_t& _y0,real_t _t0,real_t _t1,
                      solution_t* _solution,bool _append,int _maxsteps) {

  options_t opt=options;
  guess_options(options,_t0,_t1);

  try {
    this->initialize(_eval,_y0,_t0,_solution,_append);

    if (fabs(_t0-_t1)*
        m_normdy/std::max(m_normdy,options.atol/options.rtol)<options.rsmin)
      return Success;

    while (this->step(_t1)) {
      if (_maxsteps>0 && nsteps()==_maxsteps)
        throw ForceStop;
    }

  } catch (EvalState state) {
    m_eval_state=state;
  }
  options=opt;

  return m_eval_state;
}

//-----------------------------------------------------------------------------

template <typename EVAL>
void RK43<EVAL>::guess_options(options_t& _opts,real_t _t0,real_t _t1) {
  real_t span=fabs(_t1-_t0);

#ifndef __CUDACC__
  assert(std::isfinite(span));
#else
  assert(isfinite(span));
#endif

  if (_opts.rsmin<=real_t(0))  // fixme: taking absolute value
    _opts.rsmin=std::numeric_limits<real_t>::epsilon()*span*real_t(1e3);
  if (_opts.hmax<=real_t(0))
    _opts.hmax=span/real_t(10.0);
  if (_opts.h0<=real_t(0))
    _opts.h0=_opts.hmax*0.5;
  if (_opts.atol<=real_t(0))
    _opts.atol=1e-6;
  if (_opts.rtol<=real_t(0))
    _opts.rtol=1e-3;
  if (_opts.rho<=real_t(0))
    _opts.rho=0.8;
}

//-----------------------------------------------------------------------------

/** Output _options to _s.
    \ingroup vc_ode
 */
template <typename real_t>
std::ostream& operator<<(std::ostream& _s,
                         const rk43_options_t<real_t>& _options) {
  return _s << "hmax=" << _options.hmax
            << ", rsmin=" << _options.rsmin
            << ", h0=" << _options.h0
            << ", atol=" << _options.atol
            << ", rtol=" << _options.rtol
            << ", rho=" << _options.rho;
}

/** Output Solution _sol (values y only) to _s.
    \ingroup vc_ode
 */
template <typename real_t,typename vec_t>
std::ostream& operator<<(std::ostream& _s,
                         const Solution<real_t,vec_t>& _sol) {
  for (typename std::vector<vec_t>::const_iterator ii=_sol.y.begin();
       ii!=_sol.y.end();++ii)
    _s << *ii << std::endl;
  return _s;
}

//-----------------------------------------------------------------------------

# if defined(_MSC_VER)
#  pragma warning (pop)
# endif

} // namespace ode
} // namespace math
} // namespace VC

//-----------------------------------------------------------------------------


//=============================================================================

#endif // VC_MATH_ODE_RK43_HH
