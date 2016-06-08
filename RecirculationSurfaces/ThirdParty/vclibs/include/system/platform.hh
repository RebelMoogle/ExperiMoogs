//==============================================================================
// $TEMPLATE_HEADLINE$
//-----------------------------------------------------------------------------
// $Id: debug.hh,v 1.14 2001/11/19 09:01:57 roessl Exp $
// $Date: 2001/11/19 09:01:57 $
// $Revision: 1.14 $
//=============================================================================

#ifndef VC_SYSTEM_PLATFORM_HH
#define VC_SYSTEM_PLATFORM_HH

//=============================================================================

/** \file platform.hh
    \ingroup vc_system
    Platform related definitions.
 */

//== INCLUDE FILES ============================================================

//== DEFINITIONS ==============================================================

/** \defgroup vc_platform Platform dependent definitions.

    Includes definitions, maps system calls (underscore form of POSIX
    calls in Windows), implements missing functions.
 */

#ifdef DOXYGEN_SKIP

/** \def __RESTRICT
    \brief  Compiler's version of \c restrict keyword (if available).
    \ingroup vc_platform
 */
# define __RESTRICT

/** \def THREAD_LOCAL
    \ingroup vc_platform
    Compiler's version of \c __thread keyword (if available).
 */

/** \def __VC_HAS_SIGTRAP
    \ingroup vc_platform
    Defined if POSIX signals and SIGTRAP are available and SIGTRAP can be
    used as breakpoint.
*/
# define __VC_HAS_SIGTRAP

/** \def VC_ALLOCA_LIMIT
    \ingroup vc_platform
    Preferred maximum size (in bytes) for memory allocated by `alloca()`.

    Define this value externally (e.g., `-DVC_ALLOCA_LIMIT=0x10000`)
    to control the limit. A value of `0` should _void_ usage of
    `alloca` in algorithms. The default setting allows use of
    `alloca()`.
 */
# define VC_ALLOCA_LIMIT

#endif

//------------------------------------------------------------------------------

#if !defined (VC_PLATFORM)
# if defined(_WIN64)
#  define VC_PLATFORM "Win64"
#  define VC_PLATFORM_WINDOWS 64
# elif defined(_WIN32)
#  define VC_PLATFORM "Win32"
#  define VC_PLATFORM_WINDOWS 32
# else
#  define VC_PLATFORM "Unix"
#  define VC_PLATFORM_UNIX 1
# endif
#endif

#if defined(VC_PLATFORM_UNIX)

//-UNIX{-----------------------------
# include <unistd.h>
# include <fcntl.h>

# define _sys(call) call
# define __VC_HAS_SIGTRAP // debug.hh
//-}UNIX-----------------------------

#elif defined(VC_PLATFORM_WINDOWS)

//-WINDOWS{--------------------------

// POSIX (system) calls
# include <fcntl.h>
# include <io.h>
# include <stdlib.h>
# include <process.h>
# include <malloc.h>

# define _sys(call) _##call

# define fileno _sys(fileno)
# define tempnam _sys(tempnam)
# define getpid _sys(getpid)
# define isatty _sys(isatty)

# if !(defined(__MINGW32__) || defined(__MINGW64__))
#  define popen _sys(popen)
#  define pclose _sys(pclose)
#  define snprintf _sys(snprintf)

template <typename T> inline T log2(T _x) { return log(_x)/log(T(2)); }
template <typename T> inline T cbrt(T _x) { return pow(_x,T(1)/T(3)); }

# else
#   define VC_PLATFORM_WINDOWS_MINGW
#   include <stdio.h>

/********* port __GNUC_PREREQ macro to mingw *********/
# if !defined __GNUC_PREREQ

# if !defined __MINGW_H
#  include <_mingw.h>
#  define __GNUC_PREREQ(major, minor)       __MINGW_GNUC_PREREQ(major, minor)
# else
#  if defined (__GNUC_MINOR__)
#   define __GNUC_PREREQ(major, minor)      __GNUC__ > (major) || (__GNUC__ == (major) && __GNUC_MINOR__ >= (minor)))
#  else
#   define __GNUC_PREREQ(major, minor)      0
#  endif
# endif

# endif /* __GNUC_PREREQ */

# endif

// math
#  if !defined(_CRT_RAND_S)
#   define drand48() \
    (double(::rand())/double(RAND_MAX)) // WEAK (potential flaw)
#  else
inline double
drand48() {
 unsigned int v; do rand_s(&v); while(v==UINT_MAX); return double(v)/double(UINT_MAX);
}
#  endif

# include <float.h>

# if !defined(_USE_MATH_DEFINES)
#  define M_E           2.7182818284590452354   /* e */
#  define M_LOG2E       1.4426950408889634074   /* log_2 e */
#  define M_LOG10E      0.43429448190325182765  /* log_10 e */
#  define M_LN2         0.69314718055994530942  /* log_e 2 */
#  define M_LN10        2.30258509299404568402  /* log_e 10 */
#  define M_PI          3.14159265358979323846  /* pi */
#  define M_PI_2        1.57079632679489661923  /* pi/2 */
#  define M_PI_4        0.78539816339744830962  /* pi/4 */
#  define M_1_PI        0.31830988618379067154  /* 1/pi */
#  define M_2_PI        0.63661977236758134308  /* 2/pi */
#  define M_2_SQRTPI    1.12837916709551257390  /* 2/sqrt(pi) */
#  define M_SQRT2       1.41421356237309504880  /* sqrt(2) */
#  define M_SQRT1_2     0.70710678118654752440  /* 1/sqrt(2) */
# endif

namespace std {
#  if __cplusplus>199711L
  template <typename T> inline bool isfinite(T _x) { return _finite(_x) != 0; }
  template <typename T> inline bool isnan(T _x) { return _isnan(_x) != 0; }
  template <typename T> inline bool isinf(T _x) { return _finite(_x)==0 || _isnan(_x)!=0; }

  //template <typename T> inline T log2(T _x) { return log(_x)/log(T(2)); }
  //template <typename T> inline T cbrt(T _x) { return pow(_x,T(1)/T(3)); }
#  endif
}
//-}WINDOWS--------------------------
#else
# error "unknown platform"
#endif


#if defined(__GNUC__)

# if !defined(VC_PLATFORM_WINDOWS_MINGW)
#   include <alloca.h>
# endif

# define __RESTRICT __restrict
# define THREAD_LOCAL __thread

#elif defined(__clang__)

// check this!
# define __RESTRICT __restrict
# define THREAD_LOCAL __thread

#elif defined(_MSC_VER)

# define __RESTRICT __restrict
# define THREAD_LOCAL __declspec(thread)

#else

# define __RESTRICT
# define THREAD_LOCAL

#endif

#ifndef VC_ALLOCA_LIMIT
# define VC_ALLOCA_LIMIT 0x10000
#endif

//== CLASS DEFINITIONS ========================================================


//== INLINE IMPLEMENTATION ====================================================

/** Avoid warning 'unused variable'.
    \ingroup vc_platform
 */
template <typename T> void use_nowarn(const T&) {}

//=============================================================================
#endif // VC_SYSTEM_PLATFORM_HH
