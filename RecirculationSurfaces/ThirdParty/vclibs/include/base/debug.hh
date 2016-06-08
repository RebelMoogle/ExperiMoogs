//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_DEBUG_HH
#define VC_BASE_DEBUG_HH


//== INCLUDES =================================================================

#include <string>
#include <cassert>
#include <iostream>

#include "../system/platform.hh"

#ifdef __VC_HAS_SIGTRAP
# include <csignal>
#endif


#include "backtrace.hh"
#include "time.hh"

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

/** \defgroup vc_base_debug Logging and debugging tools
    \ingroup vc_base

    VC::base::Logger provides a logging at different levels that are
    set and evaluated at *runtime*. Log events are *always*
    evaluated, regardless of the code version (also for `NDEBUG`
    defined).

    The debugging macros in debug.hh are evaluated at *compile* time
    in that *no* code is produced if `NDEBUG` is defined! Generally,
    debug macros provide more information like source file and line
    number (one essential reason why macros are used.).

    All macro definitions within debug.hh depend on `NDEBUG`
    - Macros have *no effect* if `NDEBUG` is defined
    - The evalutaion of macros (and their arguments!) applies only
      to debug code.<br>
      This implies argument (evaluation) of these macros **must not**
      have any side effects!
    - All macros which "execute" statements, i.e., produce output, are
      defined as a code *block* within braces `{}`.

    Note that all `const char*` arguments to `dbg_*`
    functions must refer to strings which remain *persistent* at the
    specified location. (*Hint*: use only constant strings, e.g.,
    `"string"`.)

    **Important note:**

    Defining `NDEBUG` changes the definition of the *internally* used
    type `_DbgErrorLocation`. This, libraries may no more be binary
    compatible or "link-able" when mixing debug code with optimized
    (`NDEBUG`) code.

    Therefore, any method that takes a `_DbgErrorLocation` should be
    defined *alternatively* in `NDEBUG` code using `DbgLocation`.
 */

# ifndef DOXYGEN_SKIP

  //
  // Keep these here to have the following functions inlined.
  // Don't mess around with these globals! (Although one cannot
  // harm anything.)
  //

extern std::ostream*              _dbg_stream;
extern THREAD_LOCAL std::ostream* _dbg_stream_thread;
extern bool                       _dbg_colors;
extern THREAD_LOCAL bool          _dbg_colors_thread;
extern const std::ostream*        _dbg_null;

extern const char*                _dbg_line_prefix;
extern THREAD_LOCAL const char*   _dbg_line_prefix_thread;

extern const char*                _dbg_here_short_format;

extern const char*                _dbg_base_path;

extern bool                       _dbg_being_debugged;

// check if `gdb` is attached (`==1`) or not (`==0`) or don't know (`==-1`)
int _dbg_is_gdb_present();

// "removes" leading ".*#{_basepath}/" from _path, return pointer into _path
const char* _dbg_strip_path(const char* _path,const char* _basepath);
// prints short location information
std::ostream& _dbg_here_short(const char* _file,int _line,const char* _function);
// print dbg_line_prefix()
std::ostream& _dbg_write_line_prefix();
// get FIXME prefix
const char* _dbg_fixme();
# endif

/** Set output stream for debugging messages.
    \ingroup vc_base_debug

    The default is `std::cerr`. Consider calling
    VC::base::Logger::set_properties() for simple initialization,
    e.g., from environment variable.

    \param _os output stream
    \param _colors Use ANSI escape sequences to render colors?
    `0`=no, `1`=yes/forced, any other value implies isa_tty().

    \sa dbg_set_stream_thread(), dbg_stream(), dbg_null(),
    VC::base::Logger::set_properties(), [\ref vc_base_ansiesc]
 */
void dbg_set_stream(std::ostream& _os,int _colors=-1);

/** Set use of ANSI escape sequences to render colors.
    \ingroup vc_base_debug
    \sa dbg_set_stream(), [\ref vc_base_ansiesc]
 */
inline void dbg_set_ansi_colors(bool _colors) {
  _dbg_colors=_colors;
}

/** Overrrides dbg_set_stream() for current thread.
    \ingroup vc_base_debug

    Enables use of per-thread streams to avoid force of
    synchronization for debug output.

    Changes back to global setting for `_os==dbg_null()`.

    Arguments as for dbg_set_stream()
 */
void dbg_set_stream_thread(std::ostream& _os,int _colors=-1);

/** Set use of ANSI escape sequences to render colors.
    \ingroup vc_base_debug
    \sa dbg_set_stream_thread(), [\ref vc_base_ansiesc]
 */
inline void dbg_set_ansi_colors_thread(bool _colors) {
  _dbg_colors_thread=_colors;
}

/** Suppress any debug output to this stream.
    \ingroup vc_base_debug
    - Use with dbg_set_stream() to suppress debug output.
    - **Don't ever use this stream!**
    \sa dbg_set_stream()
 */
inline std::ostream& dbg_null() { return (std::ostream&) *_dbg_null; }

/** Debugging message are output to this stream.
    \ingroup vc_base_debug
    Note: dbg_stream() can be set by VC::base::Logger::set_properties().
    \sa dbg_set_stream(), VC::base::Logger::set_properties()
*/
inline std::ostream& dbg_stream() {
  if (_dbg_stream_thread!=0)
    return *_dbg_stream_thread;
  return *_dbg_stream;
}

/** Any line output from debug macros will be prefixed (default is "").
    \ingroup vc_base_debug
    String may contain format codes, see [\ref vc_base_ansiesc].
    \sa dbg_line_prefix(), dbg_set_line_prefix_thread()
 */
inline void dbg_set_line_prefix(const char* _prefix) {
  _dbg_line_prefix=_prefix;
}
/** Overrides dbg_line_prefix() for current thread.
    \ingroup vc_base_debug
    Changes back to prefix set by dbg_set_prefix() for `_prefix==0`.
    \ingroup vc_base_debug
    \sa dbg_set_line_prefix()
 */
inline void dbg_set_line_prefix_thread(const char* _prefix) {
  _dbg_line_prefix_thread=_prefix;
}
/** Any line output from debug macros will be prefixed (default is `""`).
    \ingroup vc_base_debug
    \sa dbg_set_line_prefix()
 */
inline const char* dbg_line_prefix() {
  if (_dbg_line_prefix_thread!=0)
    return _dbg_line_prefix_thread;
  return _dbg_line_prefix!=0 ? _dbg_line_prefix : "";
}

/** Set format string for short output of location.
    \ingroup vc_base_debug
    printf-like format string with substitutions for

    | format | substitute                                                  |
    |--------|-------------------------------------------------------------|
    | `"%F"` | source file name                                            |
    | `"%L"` | line number                                                 |
    | `"%f"` | function name                                               |
    | `"%t"` | `pthread_self()` if available; same as `"%T"` else          |
    | `"%T"` | `std::this_thread::get_id()`, should equal `pthread_self()` |
    | `"%%"` | evaluating to `"%"`.                                        |

    String may contain format codes, see [\ref vc_base_ansiesc].

    The default is `"^!r%F [%L]:!^."`.

    \sa dbg_short_format()
 */
inline void dbg_set_short_format(const char* _format) {
  _dbg_here_short_format=_format;
}
/** Get format string for short output of location.
    \ingroup vc_base_debug
    \sa dbg_set_short_format()
 */
inline const char* dbg_short_format() {
  return _dbg_here_short_format;
}

/** Set base path which will be stripped from `__FILE__`.
    \ingroup vc_base_debug
    Provides shorter output if `__FILE__` refers to an absolute, potentially
    long path name. The default base path is taken from `ENV[VC_DBG_BASEPATH]`.
 */
inline void dbg_strip_basepath(const char* _path) {
  _dbg_base_path=_path;
}

/** Check if we are being debugged.
    - Current version works on Linux only (see `_dbg_is_gdb_present()`).
    - For perfomance reasons, check only at process
      initialization. Will not work when attaching `gdb` to a running
      process!
 */
inline bool dbg_being_debugged() {
  return _dbg_being_debugged;
}

//-----------------------------------------------------------------------------

/** \class DbgLocation debug.hh
    \ingroup vc_base_debug
    \brief Describes location in source and stack trace.

    Instances of DbgLocation are used by a couple of debug macros to
    indicate where a a particular event has occurred. An example is
    \ref VC_DBG_LOCATION.
 */
struct DbgLocation {
  enum {
    MaxStackDepth=6 //!< maxmium depth for stack
  };

  DbgLocation(const char* _file=0,int _line=0,const char* _function_name=0)
    : file(_dbg_strip_path(_file,_dbg_base_path)),
      line(_line), function_name(_function_name),
      stack(backtrace_t(MaxStackDepth)) {
  }

  /// make this human readable (multi-line) string
  std::string to_s(const char* _delimiter="\n",bool _stack=true) const;

  const char*                      file;           //!< source `__FILE__`
  int                              line;           //!< source `__LINE__`
  const char*                      function_name;  //!< e.g., `__FUNCTION__`
  brief_backtrace_t<MaxStackDepth> stack;          //!< from backtrace_t
};

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4512)
# endif

/** \class DbgScopeMsg debug.hh
    \ingroup vc_base_debug
    \brief Print message on construction and desctruction.
 */
class DbgScopeMsg {
public:
  DbgScopeMsg(std::ostream& _os,const std::string& _msg,const DbgLocation& _loc);
  ~DbgScopeMsg();
private:
  std::ostream& m_os;
  std::string   m_msg;
  DbgLocation   m_location;
  cpu_time_t    m_time;
};

# if defined(_MSC_VER)
#  pragma warning (pop)
# endif

//-----------------------------------------------------------------------------

# define VC_DBG_STREAM ::VC::base::dbg_stream()

# ifndef NDEBUG

/** \def VC_DBG_LOCATION
    \brief Construct DbgLocation object.
    \ingroup vc_base_debug

    Uses preprocessor definition `__FILE__` to set DbgLocation::file, etc.
*/
#  if defined(__GNUC__)
#   define VC_DBG_LOCATION()                                            \
  ::VC::base::DbgLocation(__FILE__,__LINE__,__PRETTY_FUNCTION__)
#  elif defined(_MSC_VER)
#   define VC_DBG_LOCATION()                                    \
  ::VC::base::DbgLocation(__FILE__,__LINE__,__FUNCTION__)
#  else
#   define VC_DBG_LOCATION()                    \
  ::VC::base::DbgLocation(__FILE__,__LINE__,0)
#  endif

typedef DbgLocation _DbgErrorLocation;

#  define _VC_DBG_IS_NULL (&(VC_DBG_STREAM)==&::VC::base::dbg_null())

#  define _VC_DBG_PREFIX() ::VC::base::_dbg_write_line_prefix()

#  define _VC_DBG_HERE(msg)                                             \
  ((VC_DBG_STREAM << msg << "\n"), ::VC::base::_dbg_write_line_prefix() \
   << ' '                                                               \
   << VC_DBG_LOCATION().to_s                                            \
   ((std::string("\n")+::VC::base::dbg_line_prefix()+' ').c_str()))

#  if defined(__GNUC__)
#   define _VC_DBG_HERE_SHORT(msg)                                      \
  (::VC::base::_dbg_here_short(__FILE__,__LINE__,__PRETTY_FUNCTION__) << msg)
#  elif defined(_MSC_VER)
#   define _VC_DBG_HERE_SHORT(msg)                                      \
  (::VC::base::_dbg_here_short(__FILE__,__LINE__,__FUNCTION__) << msg)
#  else
#   define _VC_DBG_HERE_SHORT(msg)                              \
  (::VC::base::_dbg_here_short(__FILE__,__LINE__,0) << msg)
#  endif

//
// Refers to "if (_VC_DBG_IS_NULL) { ... debug output ... }"
//
// We do it this way because there should be no attempt to write to
// dbg_null(), even if it does nothing. This means there should be
// as few as possible performance penalty when switching to
// dbg_null().
// In the same sense, std::flush should be used only once.
//

/** \def VC_DBG_HERE
    \brief dump message _msg and DbgLocation
    \ingroup vc_base_debug
*/
#  define VC_DBG_HERE(msg)                                              \
  do { if (!_VC_DBG_IS_NULL) {                                          \
      _VC_DBG_PREFIX(); _VC_DBG_HERE(msg) << std::endl << std::flush;   \
    } } while (false)

/** \def VC_DBG_HERE_SHORT
    \brief dump __FILE__ [ __LINE__]: _msg
    \ingroup vc_base_debug
*/
#  define VC_DBG_HERE_SHORT(msg)                                \
  do { if (!_VC_DBG_IS_NULL) { _VC_DBG_PREFIX();                \
      _VC_DBG_HERE_SHORT(msg) << std::endl << std::flush; }     \
  } while (false)

/** \def VC_DBG_TRACE(msg)
    \brief dump msg
    \ingroup vc_base_debug
*/
#  define VC_DBG_TRACE(msg)                                     \
  do { if (!_VC_DBG_IS_NULL) {                                  \
      _VC_DBG_HERE_SHORT(msg) << std::endl << std::flush;       \
    } } while (false)

/** \def VC_DBG_P
    \brief dump variable _var using ostream.operator>>
    \ingroup vc_base_debug
*/
#  define VC_DBG_P(var)                                         \
  do { VC_DBG_TRACE(#var " = " << (var)); } while (false)

/** \def VC_DBG_PS
    \brief dump "string" variable _var using ostream.operator>>
    Same as VC_DBG_P but enclose value in '"'.
    \ingroup vc_base_debug
*/
#  define VC_DBG_PS(var)                                                \
  do { VC_DBG_TRACE(#var " = \"" << (var) << '"'); } while (false)

/** \def VC_DBG_PQS
    \brief dump `QString` variable _var using ostream.operator>>
    Same as VC_DBG_P but enclose value in '"'.
    \ingroup vc_base_debug
*/
#  define VC_DBG_PQS(var)                                               \
  do { VC_DBG_TRACE(#var " = \"" <<                                     \
                    (var).toLocal8Bit().data() << '"'); } while (false)


/** \def VC_DBG_P_IF
    \brief conditional VC_DBG_P
    \ingroup vc_base_debug
*/
#  define VC_DBG_P_IF(condition,var)                    \
  do { if (condition) VC_DBG_P(var); } while (false)

/** \def VC_DBG_COUNT
    \brief dump message and counter
    \ingroup vc_base_debug
*/
#  define VC_DBG_COUNT(msg)                             \
  do { static unsigned long count=0ul;                  \
    VC_DBG_TRACE("#(" << count++ << ") " << msg);      \
  } while (false)

/** \def VC_DBG_TRACE_ONCE
    \brief dump message once on first time of execution
    \ingroup vc_base_debug
*/
#  define VC_DBG_TRACE_ONCE(msg)                                \
  do { static bool __was_here(false);                           \
    if (!__was_here) { __was_here=true; VC_DBG_TRACE(msg); }    \
  } while (false)

/** \def VC_DBG_FIXME
    \brief dump FIXME message once (on first time of execution)
    \ingroup vc_base_debug
*/
#  define VC_DBG_FIXME(msg)                                     \
  do { VC_DBG_TRACE_ONCE(::VC::base::_dbg_fixme() << msg); } while(false)

/** \def VC_DBG_BACKTRACE
    \brief dump backtrace_t (deeper stacktrace then VC_DBG_HERE)
    \ingroup vc_base_debug
*/
#  define VC_DBG_BACKTRACE(msg)                                 \
  do { if (!_VC_DBG_IS_NULL) {                                  \
      _VC_DBG_PREFIX();                                         \
      _VC_DBG_HERE_SHORT(msg) << std::endl;                     \
      ::VC::base::_dbg_write_line_prefix();                     \
      VC_DBG_STREAM << (::VC::base::backtrace_t().symbolic      \
                        (::VC::base::dbg_line_prefix()))        \
                    << std::endl << std::flush;                 \
    } } while (false)

#  ifdef __VC_HAS_SIGTRAP
/** \def VC_DBG_BREAK
    \brief trigger breakpoint (if platform supports SIGTRAP)
    \ingroup vc_base_debug
*/
#   define VC_DBG_BREAK(msg)                                            \
  do {   if (!_VC_DBG_IS_NULL) {                                        \
      _VC_DBG_HERE(msg);                                                \
      VC_DBG_STREAM << "\n --- breakpoint ---\n" << std::flush;         \
    }                                                                   \
    else std::cerr << msg << "\n --- breakpoint ---\n" << std::flush;   \
    std::raise(SIGTRAP);                                                \
  } while (false)
#  else
/** \def VC_DBG_BREAK
    \brief trigger breakpoint (if platform supports SIGTRAP)
    \ingroup vc_base_debug
*/
#   define VC_DBG_BREAK(msg)                                            \
  do { if (!_VC_DBG_IS_NULL) {                                          \
      _VC_DBG_HERE(msg); VC_DBG_STREAM << std::endl << std::flush;      \
    }                                                                   \
    assert(!"SIGTRAP is not supported by platform. Raising SIGABRT..."); \
  } while (false)
#  endif

/** \def VC_DBG_BREAK_IF
    \brief trigger conditional breakpoint
    \ingroup vc_base_debug
*/
#  define VC_DBG_BREAK_IF(cond)                                 \
    do { if (cond) {                                            \
        VC_DBG_TRACE(# cond);                                   \
        VC_DBG_BREAK(">>> Hit conditional breakpoint <<<");     \
      } } while (false)

/** \def VC_DBG_DBREAK
    \brief same as VC_DBG_BREAK but break only if dbg_being_debugged()
    \ingroup vc_base_debug
*/
#  define VC_DBG_DBREAK(msg)                                          \
  do { if (::VC::base::dbg_being_debugged()) { VC_DBG_BREAK(msg); }   \
         else { VC_DBG_TRACE_ONCE("skipping break point"); }          \
     } while (false)

/** \def VC_DBG_DBREAK_IF
    \brief same as VC_DBG_BREAK_IF but break only if dbg_being_debugged()
    \ingroup vc_base_debug
*/
#    define VC_DBG_DBREAK_IF(cond)                                      \
  do { if (::VC::base::dbg_being_debugged()) { VC_DBG_BREAK_IF(cond); } \
       else { VC_DBG_TRACE_ONCE("skipping conditional break point"); }  \
     } while (false)

/** \def VC_DBG_SCOPE
    \brief print msg on entry and exit of block
    \ingroup vc_base_debug
    \sa VC::base::DbgScopeMsg
*/
#  define VC_DBG_SCOPE(msg)                     \
  ::VC::base::DbgScopeMsg __scope ## __LINE__   \
  (VC_DBG_STREAM,(msg),VC_DBG_LOCATION());

extern DbgLocation __dbg_tic_location;
extern cpu_time_t  __dbg_tic_time;

/** \def VC_DBG_TIC
    \brief start simple stop watch (a la Matlab)
    \ingroup vc_base_debug
    \sa VC_DBG_TIC, cpu_time_t
*/
#  define VC_DBG_TIC()                                          \
  do { ::VC::base::__dbg_tic_time=::VC::base::cpu_time_t();     \
    ::VC::base::__dbg_tic_location=VC_DBG_LOCATION();           \
  } while (false)

/** \def VC_DBG_TOC
    \brief stop stop-watch and print elapsed time since last VC_DBG_TIC()
    \ingroup vc_base_debug
    \sa VC_DBG_TIC, cpu_time_t
*/
#  define VC_DBG_TOC()                                                  \
  do { ::VC::base::__dbg_tic_location=VC_DBG_LOCATION();                \
    if (!_VC_DBG_IS_NULL) {                                             \
      ::VC::base::_dbg_write_line_prefix()                              \
        << VC_DBG_LOCATION().to_s(" : ",false) << " - elapsed time: "   \
        << ::VC::base::cpu_time_diff_t                                  \
        (::VC::base::__dbg_tic_time).to_s()<< "\n";                     \
      ::VC::base::_dbg_write_line_prefix()                              \
          << "\tstarted @ " <<                                          \
          ::VC::base::__dbg_tic_location.to_s(" ",false)                \
          << std::endl << std::flush;                                   \
    }                                                                   \
  } while (false)

# else // NDEBUG

struct _DbgErrorLocation {
  _DbgErrorLocation() {}
};

#  define VC_DBG_LOCATION() ::VC::base::_DbgErrorLocation()
#  define VC_DBG_HERE(msg)
#  define VC_DBG_HERE_SHORT(msg)
#  define VC_DBG_TRACE(msg)
#  define VC_DBG_P(var)
#  define VC_DBG_PS(var)
#  define VC_DBG_PQS(var)
#  define VC_DBG_P_IF(condition,var)
#  define VC_DBG_COUNT(msg)
#  define VC_DBG_TRACE_ONCE(msg)
#  define VC_DBG_FIXME(msg)
#  define VC_DBG_BACKTRACE(msg)
#  define VC_DBG_BREAK(msg)
#  define VC_DBG_BREAK_IF(cond)
#  define VC_DBG_DBREAK(msg)
#  define VC_DBG_DBREAK_IF(cond)
#  define VC_DBG_SCOPE(msg)
#  define VC_DBG_TIC()
#  define VC_DBG_TOC()

# endif // NDEBUG

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_DEBUG_HH defined
