//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_PROTOTHREADS_HH
#define VC_BASE_PROTOTHREADS_HH


//== INCLUDES =================================================================

#include <type_traits>

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

/** \defgroup vc_base_pt Protothreads
    \ingroup vc_base

   [Protothreads](http://dunkels.com/adam/pt/) are stack-less
   [coroutines](http://en.wikipedia.org/wiki/Coroutine).

   This is a slightly modified version of the original work referenced
   above.

   Note that all `pt_()` are **macros**!

   ## Example

   ~~~{.cpp}
   struct : pt_state_t { int result=2; } pt;
   auto primes=[pt]() mutable {
     pt_begin(pt);
     pt_yield(pt); // 2 is even
     ++pt.result;  // 3 and all other primes are odd
     while (true) {
       pt_yield(pt);
       int i,n=pt.result;
       do {
         n+=2;
         for (i=3;i*i<=n && n%i!=0;i+=2) {}
       } while (i*i<=n);
       pt.result=n;
     }
     pt_end(pt);
   };

   for (int i=0;i<100;++i)
     cout << primes()->result << ' ';
   cout << endl;
   ~~~

   ## Limitations

   - Local variables are **not** preserved between calls. Use
     `mutable` captures in lambdas instead.
   - The portable [protothreads
     implementation](http://dunkels.com/adam/pt/expansion.html) relies
     on (a particular abuse of) the `switch` statement. You **cannot**
     use `switch` in a protothread.
*/

# ifndef DOXYGEN_SKIP

#  define LC_RESUME(c) switch (c) { case 0:
#  define LC_SET(c) c=__LINE__; case __LINE__:
#  define LC_END(c) }


#  define PT_RETURN(pt,s)                       \
  do {                                          \
    pt.ps_=VC::base::pt::STATE::s;              \
    return VC::base::pt::return_state(pt)(pt);  \
  } while (false)                               \

# define PT_WAIT_UNTIL(pt,c,rv) \
  do {                                           \
  LC_SET(pt.lc_);                                \
  if (!(c))                                      \
    PT_RETURN(pt,rv);                            \
  } while (true);

# endif

/** Reset `pt` to initial state.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_reset(pt)                                           \
  do {                                                          \
    typedef std::remove_reference<decltype((pt))>::type X;      \
    X* p=&(pt);                                                 \
    p->~X();                                                    \
    new(p) X;                                                   \
  } while (false)

/** Declare start of a protothread.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_begin(pt)                           \
  { bool _yield_flag=true;                      \
    LC_RESUME(pt.lc_)

/** Yield from current protothread.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_yield(pt)                           \
  do {                                          \
    _yield_flag=false;                          \
    LC_SET(pt.lc_);                             \
    if (!_yield_flag)                           \
      PT_RETURN(pt,YIELDED);                    \
  } while (false)

/** Yield until `c` becomes true.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_yield_until(pt,c)                  \
  do {                                         \
    _yield_flag=false;                         \
    LC_SET(pt.lc_);                            \
    if (!_yield_flag || !(c))                  \
      PT_RETURN(pt,YIELDED);                   \
  } while (false)

/** Wait until `c` becomes true.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_wait_until(pt,c)                   \
  LC_SET(pt.lc_);                              \
  if (!(c))                                    \
    PT_RETURN(pt,WAITING);

/** Wait while `c` is true.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_wait_while(pt,c)                    \
  pt_wait_until(pt,!(c))

/** Exit `pt` and reset it state.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_exit(pt)                            \
  do {                                          \
    pt_reset(pt);                               \
    PT_RETURN(pt,EXITED);                       \
  } while (false)

/** Reset `pt` it state and return.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_restart(pt)                         \
  do {                                          \
    pt_reset(pt);                               \
    PT_RETURN(pt,WAITING);                      \
  } while (false)

/** Declare end of a protothread.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_end(pt)                             \
  LC_END(pt.lc_);                               \
  _yield_flag=false;                            \
  pt_reset(pt);                                 \
  PT_RETURN(pt,ENDED);                          \
  }

/** Make `pt` wait for `child_call` to finish.
    \ingroup vc_base_pt
    \hideinitializer
 */
# define pt_wait_thread(pt,child_call)          \
  pt_wait_while(pt,VC::base::pt::is_active((child_call)))

//-----------------------------------------------------------------------------
/// Protothreads data structures \ingroup vc_base_pt
namespace pt {
//-----------------------------------------------------------------------------

/// state of a protothread \ingroup vc_base_pt
enum class STATE {
  WAITING=0,
  YIELDED=1,
  EXITED=2,
  ENDED=3
};

/** State of a protothread
    \ingroup vc_base_pt

   In contrast to the [original
   implementation](http://dunkels.com/adam/pt/), I store not only the
   local contuination information but the thread_state(): a
   protothread returns a pointer to `state_t`. Feel free to use
   sub-classes, e.g., to store results. The **default CTOR** is called
   for resetting the thread.

   **Caveal:** The members `lc_` and `ps_` are declared `public` for
   the sake of a simple implementation. **Don't** use them directly!
*/
struct state_t {
  int    lc_=0;              //!< local continuation state
  STATE  ps_=STATE::WAITING; //!< protothread state

  /// get thread state
  STATE thread_state() const { return ps_; }
  /// Is thread waiting from `pt_wait_until()`?
  bool is_waiting() const { return ps_==STATE::WAITING; }
  /// Has thread yielded by `pt_yield()`?
  bool has_yielded() const { return ps_==STATE::YIELDED; }
  /// Has thread reached `pt_end()`?
  bool has_ended() const { return ps_==STATE::ENDED; }
  /// Has thread stopped due to `pt_exit()`?
  bool has_exited() const { return ps_==STATE::EXITED; }
  /// is_waiting() or has_yielded()
  bool is_active() const { return int(ps_)<int(STATE::ENDED); }
};

/// Is thread waiting from `pt_wait_until()`? \ingroup vc_base_pt
inline bool is_waiting(STATE _ps) { return _ps==STATE::WAITING; }
/// Has thread yielded by `pt_yield()`? \ingroup vc_base_pt
inline bool has_yielded(STATE _ps) { return _ps==STATE::YIELDED; }
/// Has thread reached `pt_end()`? \ingroup vc_base_pt
inline bool has_ended(STATE _ps) { return _ps==STATE::ENDED; }
/// Has thread stopped due to `pt_exit()`? \ingroup vc_base_pt
inline bool has_exited(STATE _ps) { return _ps==STATE::EXITED; }
/// is_waiting() or has_yielded() \ingroup vc_base_pt
inline bool is_active(STATE _ps) { return int(_ps)<int(STATE::ENDED); }
/// is_waiting() or has_yielded() \ingroup vc_base_pt
inline bool is_active(const state_t& _ps) { return _ps.is_active(); }

/** Return copy of persistent state.
    \ingroup vc_base_pt

    Utility used by PT_RETURN (and thus any form of `pt_yield` or
    `pt_wait_*()`): defines what part of `persistent_state_t`
    is copied on yield/return.

    Define **specialization** to change copy behavior.

    \tparam STATE persistent state used by protothread
 */
template <typename STATE>
struct return_state_t {
  typedef STATE persistent_state_t;
  persistent_state_t operator()(const persistent_state_t& _state) const {
    static_assert(std::is_base_of
                  <VC::base::pt::state_t,persistent_state_t>::value,
                  "return_state_t<> requires subclass of state_t");
    return _state;
  }
};

/// see return_state_t \ingroup vc_base_pt
template <typename STATE>
return_state_t<STATE> return_state(const STATE&) {
  return return_state_t<STATE>();
}

//-----------------------------------------------------------------------------
} // namespace pt
//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_PROTOTHREADS_HH defined
