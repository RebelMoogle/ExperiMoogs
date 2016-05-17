//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_DEMANGLE_HH
#define VC_BASE_DEMANGLE_HH


//== INCLUDES =================================================================

#include <string>

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

/** Demangle C++ symbols (if supported by platform).
    \ingroup vc_base_debug
 */
std::string demangle(char const * _mangled_name);

/** Demangle C++ symbols (if supported by platform).
    \ingroup vc_base_debug
 */
inline
std::string demangle(const std::string& _mangled_name) {
  return demangle(_mangled_name.c_str());
}

//-----------------------------------------------------------------------------

/** Demangle result of `typeid(T).name()`.
    \ingroup vc_base_debug
    Depending on the compiler, _name should already be demangled. However,
    e.g., GNU g++ returns a mangled name.
 */
inline
std::string demangle_typeid_name(const char* _name) {
#ifdef __GNUC__
    return demangle(_name);
#else
    return std::string(_name);
#endif
}

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_DEMANGLE_HH defined
