//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_TTY_HH
#define VC_BASE_TTY_HH


//== INCLUDES =================================================================

#include <string>
#include <iostream>

/** \defgroup vc_base_tty Query information about TTY.
    \ingroup vc_base

    The global flag force_not_a_tty forces isa_tty() to always return
    `false`. The flag will be set automatically if the environment
    variable `VC_FORCE_NOT_A_TTY` is defined.
*/

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {

//-----------------------------------------------------------------------------

/// global flag to isa_tty() (default is false) \ingroup vc_base_tty
extern bool force_not_a_tty;

//-----------------------------------------------------------------------------

/** Try to figure out if input _is comes from a TTY.
    \ingroup vc_base_tty

    Note: there is no portable way to do this in C++ (require the file
    descriptor fd for `_is`). In doubt, we call ::isatty(fd) only for
    `fd=0` (`cin`).

    \returns true if _os is a tty and false if not, or if we cannot tell,
    or if global force_not_a_tty==true
 */
bool isa_tty(const std::istream& _is);


/** Try to figure out if _os outputs to a TTY.
    \ingroup vc_base_tty

    Note: there is no portable way to do this in C++ (require the file
    descriptor fd for `_os`). In doubt, we call `::isatty(fd)` only for
    `fd=1,2` (`cout`,`cerr`,`clog`).

    \returns true if _os is a tty and false if not, or if we cannot tell,
    or if global force_not_a_tty==true
 */
bool isa_tty(const std::ostream& _os);

//-----------------------------------------------------------------------------

/** Query current size of the terminal window.
    \ingroup vc_base_tty

    \arg *Unix*: Queries terminal size for either `stdin`,
    `stdout,` or `stderr`, whichever is a TTY, or use the fallback.

    \arg Fallback for *any platform*: Query environment variables
    `"COLUMNS"` and `"LINES"` (as set and updated, e.g., by
    `bash`). If any of these variables is not set or cannot be
    converted to a positive integer `get_terminal_size()` provides
    standard values and returns `false`.


    \param[out] _width width (number of characters)
    \param[out] _height height (number of characters)
    \return `true` if the information could be queried from the host system
    (including environment variables), otherwise `false`.
    In the second case \c _width and \c _height were set to "reasonable"
    standard values.
 */
bool get_terminal_size(int& _width,int& _height);

//-----------------------------------------------------------------------------


//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_TTY_HH defined
