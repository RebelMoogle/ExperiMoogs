//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_PRINTF_HH
#define VC_BASE_PRINTF_HH


//== INCLUDES =================================================================

#include <string>
#include <iostream>

#include "exception.hh"

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4290)
# endif

# if defined(__GNUC__)
#  define _IS_PRINTF_LIKE(m,n) __attribute__((format(printf,m,n)))
# else
#  define _IS_PRINTF_LIKE(m,n)
# endif

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** \defgroup vc_base_printf sprintf-like text formatting
    \ingroup vc_base

    This implementation is based on `std::sprintf`, i.e., the same rules
    apply, and there is *no* typechecking , etc.

    An alternative may be <a href="http://www.boost.org">boost</a>'s
    format library.
*/

/** \brief \c std::sprintf into std::string
    \ingroup vc_base_printf
    \param _str output (_str.clear() is called on failure)
    \param _format format string and values ... (see `std::sprintf`)
    \return number of characters printed or -1 on failure, se `std::sprintf`
 */
int ssprintf(std::string& _str,const char* _format,...)
  _IS_PRINTF_LIKE(2,3);

//-----------------------------------------------------------------------------

/** \brief std::ssprintf() and return result
    \ingroup vc_base_printf
    \param _format format string and values ... (see `std::sprintf`)
    \return string
    \throw runtime_error if sprintf returned -1
 */
std::string formatf(const char* _format,...) throw(vc_runtime_error)
  _IS_PRINTF_LIKE(1,2);

//-----------------------------------------------------------------------------

/** \brief \c std::fprintf to \c std::ostream.
    \ingroup vc_base_printf
    \param  _os input/output
    \param _format format string and values ... (see `std::sprintf`)
    \return _os
    \throw runtime_error if sprintf returned -1
 */
std::ostream& oprintf(std::ostream& _os,
                      const char* _format,...) throw(vc_runtime_error)
  _IS_PRINTF_LIKE(2,3);

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

namespace detail {
  template <typename T,typename S> class iosformat {
  public:
    iosformat(const T& _data,const char* _format)
      : data(S(_data)), format(_format) {}
    S           data;
    const char* format;
  };
  template <typename T,typename S>
  std::ostream& operator<<(std::ostream& _os,const iosformat<T,S>& _data) {
    return oprintf(_os,_data.format,_data.data);
  }
} // detail

# endif // DOXYGEN_SKIP

/** \brief Output formatted floating point number.
    \ingroup vc_base_printf
    <code>
    cout << fmt_fp(1.23456789,"%1.2g");
    </code>
    \tparam T `float` or `double` (or any type that can be converted to `double`)
    \param _data number to output
    \param _format format string
    \sa oprintf
 */
template <typename T>
detail::iosformat<T,double> fmt_fp(const T& _data,const char* _format="%g") {
  return detail::iosformat<T,double>(_data,_format);
}

/** \brief Output formatted integer number.
    \ingroup vc_base_printf
    <code>
    cout << fmt_int(1234,"0x%08d");
    </code>
    \tparam T any type that can be converted to `int`
    \param _data number to output
    \param _format format string
    \sa oprintf
 */
template <typename T>
detail::iosformat<T,int> fmt_int(const T& _data,const char* _format="%d") {
  return detail::iosformat<T,int>(_data,_format);
}

/** \brief Output formatted unsigned integer number.
    \ingroup vc_base_printf
    <code>
    cout << fmt_uint(1234u,"0x%08u");
    </code>
    \tparam T any type that can be converted to `int`
    \param _data number to output
    \param _format format string
    \sa oprintf
 */
template <typename T>
detail::iosformat<T,unsigned> fmt_uint(const T& _data,const char* _format="%u") {
  return detail::iosformat<T,unsigned int>(_data,_format);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

# undef _IS_PRINTF_LIKE

# if defined(_MSC_VER)
#  pragma warning (pop)
# endif

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_PRINTF_HH defined
