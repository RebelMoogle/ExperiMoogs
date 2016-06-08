//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_BACKTRACE_HH
#define VC_BASE_BACKTRACE_HH


//== INCLUDES =================================================================

# if defined(ARCH_LINUX)
#  include <execinfo.h>
#  define __VC_BASE_BACKTRACE_SUPPORTED
# endif

#include <cassert>
#include <vector>
#include <string>

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

template <unsigned MaxDepth>
class brief_backtrace_t;

/** \class backtrace_t backtrace.hh
    \ingroup vc_base_debug
    \brief Backtrace: stack of callers for debugging purposes.

    The implementation is platform dependent. backtrace_t might or might
    not be fully implemented.

    If is_supported()==false then any instance will yield is_defined()=false.

    MaxDepth specifies the maxmium depth. A "dynamic" _max_depth<=MaxDepth
    may be requested on construction.

    backtrace_t stores backtrace information only. It does not store
    any information on queried symbols (insert_symbols()).

    brief_backtrace_t may be used as a "lighter-weight" object to store
    a backtrace with prescribed maxmimum depth <= MaxDepth
    (typically << MaxDepth). Use provided copy constructors and note that
    a backtrace may be "truncated" this way.

    \sa brief_backtrace_t
 */
# if defined(__VC_BASE_BACKTRACE_SUPPORTED) || defined(DOXYGEN_SKIP)
class backtrace_t {
public:
  enum { MaxDepth = 64 }; //!< maximum stack depth

  /// constructor queries backtrace
  backtrace_t(unsigned _max_depth=MaxDepth) : m_depth(0) {
    for (unsigned i=0;i<MaxDepth;++i)
      m_addr[i]=(void*) 0;

    if (_max_depth>MaxDepth)
      _max_depth=MaxDepth;

    m_depth=backtrace((void**) m_addr,_max_depth);
  }
  /// copy from brief_backtrace_t _brief
  template <unsigned N>
  backtrace_t(const brief_backtrace_t<N>& _brief);
  /// Is backtrace supported on this platform?
  static bool is_supported() { return true; }
  /// get depth of call stack
  unsigned depth() const { return m_depth; }

  /// returns depth()==0 (essentially if is_supported()==false)
  bool is_defined() const { return depth()>0; }
  /// get maximum depth
  static unsigned max_depth() { return MaxDepth; }

  typedef const void* caller_t; //!< stored pointers
  caller_t* begin() const { return (caller_t*) m_addr; } //!< callers
  caller_t* end() const { return   (caller_t*) (m_addr+m_depth); } //!< callers

  /** Get symbols as strings from caller_t, demangle if possible.
      \param _symbols symbols are inserter here
      \param _demangle demangle C++ symbols (if available)
      \param _suppress_paths removed paths from shared objects

      \a Notes:
      \arg symbols() internally calls malloc() (a Linux flaw,
      unfortunately). Hence if the heap of the process is corrupted,
      then a call to symbols() may cause a segmentation fault! Don't
      use this, when you really want to use a debugger ;-) or use fork().
      \arg Symbol information might not be available for external
      libraries, functions declared "static", and symbols in the same
      object files as "main()".
  */
  void insert_symbols(std::vector<std::string>& _symbols,
                      bool _demangle=true,bool _suppress_paths=true) const;
  /** Provided for convenience: calls get_symbols() and concatenates _symbols.
      \param _delimiter is appended between symbol strings
      \return string containing all symbols (or empty string if !is_defined())

      See also notes on insert_symbols()!
   */
  std::string
  symbolic(const char* _delimiter="\n") const;

private:
  template <unsigned N> friend class brief_backtrace_t;

  caller_t    m_addr[MaxDepth];  //!< stack frame
  unsigned    m_depth;           //!< actual depth (elements in \c d_addr)
};

# else

class backtrace_t {
public:
  enum { MaxDepth = 0 };

  backtrace_t() {}
  backtrace_t(unsigned /*_max_depth=MaxDepth*/) {}
  template <unsigned N>
  backtrace_t(const brief_backtrace_t<N>& _brief);
  static bool is_supported() { return false; }
  unsigned depth() const { return 0; }
  bool is_defined() const { return false; }
  static unsigned max_depth() { return MaxDepth; }

  typedef const void* caller_t; //!< stored pointers
  caller_t* begin() const { return (caller_t*) 0; } //!< callers
  caller_t* end() const { return   (caller_t*) 0; } //!< callers

  void insert_symbols(std::vector<std::string>& _symbols,
                      bool _demangle=true,bool _suppress_paths=true) const;
  std::string
  symbolic(const char* _delimiter="\n") const;
};
# endif // __VC_BASE_BACKTRACE_SUPPORTED


//-----------------------------------------------------------------------------

/** \class backtrace_t backtrace.hh
    \ingroup vc_base_debug
    \brief Brief backtrace.

    This class is provided only to store (partial) backtrace_t information
    if a fixed maxmimum depth. It is intended for passing backtrace
    information with exceptions.

    Note that brief_backtrace_t is intendet to be \a lightweight, i.e.,
    _MaxDepth should be small!

    \sa backtrace_t
 */
# if defined(__VC_BASE_BACKTRACE_SUPPORTED) || defined(DOXYGEN_SKIP)
template <unsigned _MaxDepth=8>
class brief_backtrace_t {
public:
  enum { MaxDepth = _MaxDepth };  //!< maximum stack depth

  /// copy from _backtrace
  brief_backtrace_t(const backtrace_t& _backtrace) {
    unsigned i;
    for (i=0;i<_MaxDepth && _backtrace.depth();++i)
      m_addr[i]=_backtrace.m_addr[i];
    for (;i<MaxDepth;++i)
      m_addr[i]=(backtrace_t::caller_t) 0;
  }

private:
  friend class backtrace_t;

  backtrace_t::caller_t m_addr[MaxDepth]; //!< call stack
};

template <unsigned N>
backtrace_t::backtrace_t(const brief_backtrace_t<N>& _brief) {
  unsigned maxdepth=N, i;
  if (maxdepth>MaxDepth)
    maxdepth=MaxDepth;

  for (i=0;i<maxdepth && _brief.m_addr[i]!=0;++i)
    m_addr[i]=_brief.m_addr[i];

  m_depth=i;

  for (;i<MaxDepth;++i)
    m_addr[i]=0;
}

# else

template <unsigned _MaxDepth=8>
class brief_backtrace_t {
public:
  enum { MaxDepth = _MaxDepth };  //!< maximum stack depth
  brief_backtrace_t(const backtrace_t& /*_backtrace*/) {}
};

template <unsigned N>
backtrace_t::backtrace_t(const brief_backtrace_t<N>& /*_brief*/) {}


# endif // __VC_BASE_BACKTRACE_SUPPORTED

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_BACKTRACE_HH defined
