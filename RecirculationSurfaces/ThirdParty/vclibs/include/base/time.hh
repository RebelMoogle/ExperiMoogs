//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_TIME_HH
#define VC_BASE_TIME_HH


//== INCLUDES =================================================================

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <string>

#include "../system/platform.hh" // THREAD_LOCAL

//== CLASS DEFINITION =========================================================

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4068)
# endif

# define __VC_BASE_POSIX_RT_CLOCK 1
# define __VC_BASE_POSIX_TIMES 2
# define __VC_BASE_POSIX_CLOCK 0
# define __VC_BASE_WIN_HPC 3
# define __VC_BASE_OSX_MACH 4
# define __VC_BASE_STD_CHRONO 5

# if (defined(_POSIX_TIMERS) && _POSIX_TIMERS>0) && !defined(__VC_DONT_USE_RT_CLOCK)
#   define __VC_BASE_CPU_TIME_CLOCK  __VC_BASE_POSIX_RT_CLOCK
# else
#  if !defined(_WIN32) && !defined(_WIN64) && 0 // AVOID THIS! resolution 1/100s
#    include <sys/times.h>
#    define __VC_BASE_CPU_TIME_CLOCK  __VC_BASE_POSIX_TIMES
#  endif
# endif

# if (defined(_WIN32) || defined(_WIN64)) && !defined(__VC_DONT_USE_RT_CLOCK)
#   if !(defined(VC_PLATFORM_WINDOWS_MINGW))
#     include <windows.h>
#     define __VC_BASE_CPU_TIME_CLOCK __VC_BASE_WIN_HPC
#   endif
# endif

#if defined(__APPLE__) && !defined(__VC_DONT_USE_RT_CLOCK)
#  include <mach/mach_time.h>
#  define __VC_BASE_CPU_TIME_CLOCK __VC_BASE_OSX_MACH
#endif

# ifndef __VC_BASE_CPU_TIME_CLOCK
#  if  (__cplusplus > 199711L)
#    define __VC_BASE_CPU_TIME_CLOCK __VC_BASE_STD_CHRONO
#    include <chrono>
#  else
#    define __VC_BASE_CPU_TIME_CLOCK __VC_BASE_POSIX_CLOCK
#  endif
# endif

//-----------------------------------------------------------------------------

namespace VC {
namespace base {

/** \defgroup vc_base_time Measuring time
    \ingroup vc_base
*/

class cpu_time_diff_t;

/** \class cpu_time_t time.hh
    \brief Yields per-process CPU time.
    \ingroup vc_base_time


    cpu_time_t measures timestamps (CPU time, not wall clock), it is
    designed for timing short intervals (\sa cpu_time_diff_t).

    The static function cpu_time_t::seconds() returns the time in
    seconds since an unspecified starting time.

    - High resolution timers are used if possible. For POSIX systems
      supporting `_POSIX_TIMERS`,
      `clock_gettime(CLOCK_PROCESS_CPUTIME_ID)` is called. (Mind the
      notes on SMP systems.)

    - On all other systems either the (low resolution) POSIX `clock()`
      or `std::chrono::high_resolution_clock` (C++11) is used.  Note
      that `clock()` is prone to overflow/wrap around even for descent
      time intervals!

    Use cpu_time_diff_t for measuring time *differences*. Note that
    differences of seconds() (`double`) may suffer from lack of
    numerical precision, wheras cpu_time_diff_t() does not.

    ## Note on `std::chrono`

    This code avoids the C++11
    [chrono](http://en.cppreference.com/w/cpp/chrono) library if
    possible. The main reason is that the `chono` implementation of
    [high_resolution_clock](http://en.cppreference.com/w/cpp/chrono/high_resolution_clock)
    is *not necessarily* high resolution.

    The "fallback" implementation that uses `std::chrono` converts
    `chrono`'s internal representation of `std::time_point<>` to
    nanoseconds. This is, because `cpu_time_diff_t` is a subclass of
    `cpu_time_t`, and they are supposed share the *same state* which
    is interpreted as "time point" or "duration". Note that suing
    `chrono` comes with an overhead in particular for non-optimized
    code.

    ## Credits

    - Windows HPC timers added by Michael Schiefer.
    - OS X timers added by Martin Kirst.
 */
class cpu_time_t {
public:

  enum { CLOCK = __VC_BASE_CPU_TIME_CLOCK };

  /// create timestamp
  cpu_time_t();

  /// resultion() is meaningful only if high resolution timers are used
  static bool is_highres() {
    return (CLOCK==__VC_BASE_POSIX_RT_CLOCK) || (CLOCK==__VC_BASE_WIN_HPC) || (CLOCK==__VC_BASE_OSX_MACH);
  }
  /// get resolution in fraction of a second (`0.0` is *unknown*)
  static double resolution() { return sm_resolution; }

  /// get time
  double seconds() const;

  /// same as seconds
  operator double() const { return seconds(); }

  /** Get a human readable string representation.
      Selects the best unit
      - nanoseconds `"ns"`
      - milliseconds `"ms"`
      - seconds `"s"` (up to 300 before switching to minutes:seconds)
      - minutes `"min"` (up to 120 before switching to hours)
      - hours `"h"` (up to 48 before switching to days)
      - days `"d"`

      and constructs a string `"##d:#h:#m:#s"` or `"#{ns|us|ms}"`
      (`#` stands for numbers), zero entries will not be suppressed.
      \param _fpformat sprintf-like format string for floating point
      numbers (`double)` and their unit (`string`) shown for s,ns,ms.
      Default is "%.4lg%s".
      \param _unit select best unit on base of `max(_unit,(double) *this)`,
      use to obtain same units/string formats, the default value is `0.0`
      (no influence)
      \return formatted string
   */
  std::string to_s(const char* _fpformat=0,double _unit=0.0) const;

  /** @name global stop watch

      tic() and toc() implement a simple global stop watch as in
      Matlab. (If supported by the compiler, each thread has it's own
      stop watch.)

      Calls are delegated to StopWatch::global().

      \sa StopWatch

      @{
   */

  /// start global stop watch (StopWatch::tic())
  static void tic();

  /// stop global stop watch (StopWatch::toc())
  static cpu_time_diff_t toc();
  /// same, additonally base::oprintf() `_msg` to `_os` with `%s` as to_s()
  static cpu_time_diff_t toc(std::ostream& _os,
                             const char* _msg="elapsed time %s\n");

  /// @}

private:
  friend class cpu_time_diff_t;

# if __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_POSIX_RT_CLOCK
  struct timespec m_tp;    //!< timestamp
# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_POSIX_TIMES
  clock_t         m_clock;
  struct tms      m_tms;      //!< timestamp
  static long     sm_clk_tck; //!< clock ticks
# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_POSIX_CLOCK
  clock_t         m_clock;    //!< real time as returned by clock
# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_WIN_HPC
  static LARGE_INTEGER sm_frequency;
  LARGE_INTEGER        m_count;
# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_OSX_MACH
  uint64_t        m_time;
#  elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_STD_CHRONO
  typedef std::chrono::high_resolution_clock clock_t;
  size_t          m_time;     //!< nanoseconds
# else
#  error "invalid __VC_BASE_CPU_TIME_CLOCK"
# endif

  static double     sm_resolution; //!< clock resolution
};

//-----------------------------------------------------------------------------

/** \class cpu_time_diff_t time.hh
    \brief Difference of cpu_time_t time-stamps
    \ingroup vc_base_time
 */
class cpu_time_diff_t : public cpu_time_t {
public:
  /// get difference to _reference
  cpu_time_diff_t(const cpu_time_t& _reference);
  /// zero time difference
  cpu_time_diff_t();

  /// add _dt
  cpu_time_diff_t& add(const cpu_time_diff_t& _dt);
  /// add _dt
  cpu_time_diff_t& operator+=(const cpu_time_diff_t& _dt) { return add(_dt);  }
  /// add _dt
  cpu_time_diff_t operator+(const cpu_time_diff_t& _dt) const {
    return cpu_time_diff_t(*this).add(_dt);
  }
};

//-----------------------------------------------------------------------------

# if __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_POSIX_RT_CLOCK

  //
  // POSIX real-time clock
  //

inline cpu_time_t::cpu_time_t() {
  //int rv=clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&m_tp);
  int rv=clock_gettime(CLOCK_THREAD_CPUTIME_ID,&m_tp);
  //int rv=clock_gettime(CLOCK_REALTIME,&m_tp);
  assert(rv==0 && "EINVAL?! clk_id specified is not supported on this system?!"!=0);
  use_nowarn(rv);
}
inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_tp.tv_sec-=_reference.m_tp.tv_sec;
  m_tp.tv_nsec-=_reference.m_tp.tv_nsec;
}
inline cpu_time_diff_t::cpu_time_diff_t() {
  m_tp.tv_sec=time_t(0);
  m_tp.tv_nsec=long(0);
}
inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_tp.tv_sec+=_dt.m_tp.tv_sec;
  m_tp.tv_nsec+=_dt.m_tp.tv_nsec;
  return *this;
}
inline double cpu_time_t::seconds() const {
  return double(m_tp.tv_sec)+double(m_tp.tv_nsec)*1e-9;
}

# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_POSIX_TIMES

  //
  // POSIX times()
  //

inline cpu_time_t::cpu_time_t() {
  m_clock=times(&m_tms);
}
inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_clock-=_reference.m_clock;
  m_tms.tms_utime -=_reference.m_tms.tms_cutime;
  m_tms.tms_stime -=_reference.m_tms.tms_stime;
  m_tms.tms_cutime-=_reference.m_tms.tms_cutime;
  m_tms.tms_cstime-=_reference.m_tms.tms_cstime;
}
inline cpu_time_diff_t::cpu_time_diff_t() {
  m_clock=
    m_tms.tms_utime=m_tms.tms_stime=
    m_tms.tms_cutime=m_tms.tms_cstime= clock_t(0);
}
inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_clock+=_dt.m_clock;
  m_tms.tms_utime +=_dt.m_tms.tms_cutime;
  m_tms.tms_stime +=_dt.m_tms.tms_stime;
  m_tms.tms_cutime+=_dt.m_tms.tms_cutime;
  m_tms.tms_cstime+=_dt.m_tms.tms_cstime;
  return *this;
}
inline double cpu_time_t::seconds() const {
  return double(m_tms.tms_utime+m_tms.tms_stime)+double(sm_clk_tck);
}
# elif __VC_BASE_CPU_TIME_CLOCK == __VC_BASE_POSIX_CLOCK

  //
  // POSIX clock()
  //

inline cpu_time_t::cpu_time_t() {
  m_clock=std::clock();
}
inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_clock-=_reference.m_clock;
}
inline cpu_time_diff_t::cpu_time_diff_t() {
  m_clock=clock_t(0);
}
inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_clock+=_dt.m_clock;
  return *this;
}
inline double cpu_time_t::seconds() const {
  return double(m_clock)/double(CLOCKS_PER_SEC);
}

# elif __VC_BASE_CPU_TIME_CLOCK==__VC_BASE_WIN_HPC

  //
  // Windows performance counters
  //

inline cpu_time_t::cpu_time_t()  {
  int rv=QueryPerformanceCounter(&m_count);
  assert(rv!=0 && "failed to query performance counter");
}
inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_count.QuadPart -= _reference.m_count.QuadPart;
}
inline cpu_time_diff_t::cpu_time_diff_t() {
  m_count.QuadPart=0;
}
inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_count.QuadPart += _dt.m_count.QuadPart;
  return *this;
}
inline double cpu_time_t::seconds() const {
  return double(m_count.QuadPart)/double(sm_frequency.QuadPart);
}

# elif __VC_BASE_CPU_TIME_CLOCK == __VC_BASE_OSX_MACH

  //
  // OS X: Mach Absolute Time Units
  //

inline cpu_time_t::cpu_time_t()  { m_time = mach_absolute_time(); }

inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_time -= _reference.m_time;
}

inline cpu_time_diff_t::cpu_time_diff_t()  { m_time=0; }

inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_time+=_dt.m_time;
  return *this;
}

inline double cpu_time_t::seconds() const {
  return double(m_time)*sm_resolution;
}

# elif __VC_BASE_CPU_TIME_CLOCK == __VC_BASE_STD_CHRONO

  //
  // std::chrono
  //

inline cpu_time_t::cpu_time_t()  {
  m_time=std::chrono::time_point<clock_t,std::chrono::nanoseconds>
    (clock_t::now()).time_since_epoch().count();
}

inline cpu_time_diff_t::cpu_time_diff_t(const cpu_time_t& _reference) {
  m_time-=_reference.m_time;
}

inline cpu_time_diff_t::cpu_time_diff_t()  { m_time=0; }

inline cpu_time_diff_t& cpu_time_diff_t::add(const cpu_time_diff_t& _dt) {
  m_time+=_dt.m_time;
  return *this;
}

inline double cpu_time_t::seconds() const { return double(m_time)*1e-9; }

# else
#  error "invalid __VC_BASE_CPU_TIME_CLOCK"
# endif

//-----------------------------------------------------------------------------

/** \class StopWatch time.hh
    \brief Simple stop watch based on cpu_time_t and cpu_time_diff_t.
    \ingroup vc_base_time
    \sa ScopedStopWatch
 */
class StopWatch {
public:
  /// start stop watch
  void tic() { m_time=cpu_time_t(); }
  /// stop watch
  cpu_time_diff_t toc() const {
    return cpu_time_diff_t(m_time);
  }
  /// same, and base::oprintf() _msg to `_os` with `%s` as cpu_time_t::to_s()
  cpu_time_diff_t toc(std::ostream& _os,
                      const char* _msg="elapsed time %s\n");

  /// get global stop watch used by cpu_time_t::tic() (per thread instance)
  static StopWatch& global();

private:
  /// called by toc() (not inline)
  static cpu_time_diff_t _toc(std::ostream& _os,
                              const cpu_time_diff_t& _dt,const char* _msg);
  cpu_time_t m_time;
};

//-----------------------------------------------------------------------------

/** \class ScopedStopWatch time.hh
    \brief StopWatch that measures and outputs object life time.
    \ingroup vc_base_time
 */
class ScopedStopWatch : public StopWatch {
public:
  /// same as for StopWatch::toc(), output by object destructor
  ScopedStopWatch(std::ostream& _os=std::cerr,
                  const char* _msg="elapsed time %s\n")
    : out(_os), msg(_msg) { this->tic(); }
  ~ScopedStopWatch() { this->toc(out,msg.c_str()); }

  std::ostream&     out;
  const std::string msg;
private:
  /// no restart
  void tic() { StopWatch::tic(); }
};

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP
  //
  // force POD type
  //
  extern THREAD_LOCAL char __global_StopWatch_global[sizeof(StopWatch)];
# endif


# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wstrict-aliasing"

inline StopWatch& StopWatch::global() {
  return *((StopWatch*) __global_StopWatch_global);
}

# pragma GCC diagnostic pop

//-----------------------------------------------------------------------------

inline void cpu_time_t::tic() { StopWatch::global().tic(); }
inline cpu_time_diff_t cpu_time_t::toc() { return StopWatch::global().toc(); }
inline cpu_time_diff_t cpu_time_t::toc(std::ostream& _os,const char* _msg) {
  return StopWatch::global().toc(_os,_msg);
}
inline cpu_time_diff_t StopWatch::toc(std::ostream& _os,const char* _msg) {
  return _toc(_os,toc(),_msg);
}

//-----------------------------------------------------------------------------

/** Stream output of VC::base::cpu_time_t
    \ingroup vc_base_time
*/
  template <typename T,unsigned K>
  std::ostream& operator<<(std::ostream& _s,const ::VC::base::cpu_time_t& _t) {
    return _s << _t.to_s();
}

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
# if defined(_MSC_VER)
#  pragma warning (pop)
# endif
//=============================================================================
#endif // VC_BASE_TIME_HH defined
