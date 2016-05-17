//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_EXCEPTION_HH
#define VC_BASE_EXCEPTION_HH


//== INCLUDES =================================================================

#include <exception>
#include <string>

#include "debug.hh"

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

/** \def VC_RUNTIME_ERROR
    \brief Construct vc_runtime_error with additonal debug information.
    \ingroup vc_base
    Construct a vc_runtime_error with stack backtrace information on souce
    file, line number and function.
    \sa vc_runtime_error, VC_HERE, [\ref vc_base_debug]
 */
# define VC_RUNTIME_ERROR(msg) \
  ::VC::base::vc_runtime_error((msg),VC_DBG_LOCATION())

/** \class vc_runtime_error exception.hh
    \brief Runtime error exception.


    May store additonal debug information compared to `std::runtime_error`.
    \ingroup vc_base
*/
class vc_runtime_error : public std::exception {
public:
  /// construct from error message _msg
  vc_runtime_error(const std::string& _msg) throw()
    : m_msg(_msg) {}
  /** Construct from error _msg and DbgLocation.
      DbgLocation is meaningless for non-debug (`defined(NDEBUG)`) code.
      The macro VC_RUNTIME_ERROR() will construct a vc_runtime_error
      with information on souce file, line number and function.
   */
  vc_runtime_error(const std::string& _msg,const _DbgErrorLocation& _loc) throw()
    : m_msg(_msg), m_location(_loc) {}

# ifdef NDEBUG
  vc_runtime_error(const std::string& _msg,const DbgLocation&) throw()
    : m_msg(_msg) {}
# endif

  virtual ~vc_runtime_error() throw() {}

  /** Similar to message() but may be modified.
      \arg return message().c_str() non-debug code
      \arg sets up a new text including debug information for debug code
   */
  virtual const char* what() const throw();

  /// get message as passed to constructor
  const std::string& message() const { return m_msg; }

private:
    std::string        m_msg;      //!< error message
    _DbgErrorLocation  m_location; //!< DbgLocation if `!defined(NDEBUG)`
#ifndef NDEBUG
  mutable std::string  m_buf;      //!< buffer => backtrace_t::symbolic() once
#endif
};

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_EXCEPTION_HH defined
