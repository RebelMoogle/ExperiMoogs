//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_ANSI_ESCAPE_HH
#define VC_BASE_ANSI_ESCAPE_HH


//== INCLUDES =================================================================

#include <string>
#include <iostream>

#include <functional>

/** \defgroup vc_base_ansiesc Simple text rendering with ANSI escape sequences.
    \ingroup vc_base

    The use of ANSI escape sequences may (terminal) or may not (not isatty())
    be supported. We provide few functions for string substitution and
    stream output which

    - substitute format codes by ANSI escape sequences, or
    - suppress formatting (and filter out any ANSI escape sequences).

    The following *format codes* are supported

    | code    | fornmat                                                         |
    |---------|-----------------------------------------------------------------|
    | `"^!*"` | bold                                                            |
    | `"^!_"` | underline                                                       |
    | `"^!:"` | faint                                                           |
    | `"^!-"` | negative/inverse                                                |
    | `"^!+"` | positive                                                        |
    | `"^!."` | normal / normal style                                           |
    | `"^!x"` | with foreground color `x=r,g,b,c,m,y,k,w `(red,...,black,white) |
    | `"^!X"` | with bright foreground color `X=R,G,B,C,M,Y,K,W`                |
    | `"^!^"` | emphasis style                                                  |
    | `"^!!"` | error style                                                     |
    | `"^!?"` | warning style                                                   |

    Any sequence `"\x1b[*@"` with `'*'=sequence of characters <64, and `'@'` is
    any character in `64..126` is assumed to be an ANSI escape sequence.
*/

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {

/** @name String constants for ANSI escape sequences
    \ingroup vc_base_ansiesc
    @{
 */

extern const char* ESC_BOLD;
extern const char* ESC_UNDERLINE;
extern const char* ESC_FAINT;
extern const char* ESC_NEGATIVE;
extern const char* ESC_POSITIVE;
extern const char* ESC_NORMAL;

extern const char* ESC_RED;
extern const char* ESC_GREEN;
extern const char* ESC_BLUE;
extern const char* ESC_CYAN;
extern const char* ESC_MAGENTA;
extern const char* ESC_YELLOW;
extern const char* ESC_BLACK;
extern const char* ESC_WHITE;

extern const char* ESC_HI_RED;
extern const char* ESC_HI_GREEN;
extern const char* ESC_HI_BLUE;
extern const char* ESC_HI_CYAN;
extern const char* ESC_HI_MAGENTA;
extern const char* ESC_HI_YELLOW;
extern const char* ESC_HI_BLACK;
extern const char* ESC_HI_WHITE;

extern const char* ESC_EMPH;
extern const char* ESC_ERROR;

/// @}

//-----------------------------------------------------------------------------

/** Get escape sequence to _code.
    \param _code character `"X"` in `"^!X"`, see [\ref vc_base_ansiesc].
    \return sequence or empty string (always returns a valid C string)
 */
const char* ansi_escape(int _code);

/** Suppress/filter any format codes or escape sequences from _s.
    \param _s input string
    \return std::string with format codes and ANSI escape sequences removed
 */
std::string suppress_ansi_escape(const char* _s);

/// Same for std::string \ingroup vc_base_ansiesc
inline std::string suppress_ansi_escape(const std::string& _s) {
  return suppress_ansi_escape(_s.c_str());
}

/** Substitute format codes in _s by ASCII escape sequences.
    \param _s input string
    \param _render render or suppress format codes/ANSI sequences
    \return std::string with format codes replaced by ANSI escape sequences
    or suppressed (suppress_ansi_escape() if `_render==false`)

 */
std::string render_ansi_escape(const char* _s,bool _render=true);

/// Same for std::string \ingroup vc_base_ansiesc
inline std::string
render_ansi_escape(const std::string& _s,bool _render=true) {
  return render_ansi_escape(_s.c_str(),_render);
}

/** Render or suppress ANSI escape sequences and write to _os.
    \ingroup vc_base_ansiesc

    \param _os output stream
    \param _s input string
    \param _render true: write render_ansi_escape(),
    false: write suppress_ansi_escape()
    return _os
 */
std::ostream&
write_ansi_escape(std::ostream& _os,const char* _s,bool _render=true);

/// Same for std::string \ingroup vc_base_ansiesc
inline std::ostream&
write_ansi_escape(std::ostream& _os,const std::string& _s,bool _render=true) {
  return write_ansi_escape(_os,_s.c_str(),_render);
}

//-----------------------------------------------------------------------------

/** Get string containing escape sequences from a textual description.
    \ingroup vc_base_ansiesc
    The specification is a *comma separated* list of tokens. Either the
    corresponding escape sequence is added to the result, or the `_unknown`
    callback is invoked for unknown (=invalid) tokens.

    The parser is based on regular expressions and rather simple. There is
    no error checking other than the callback. Double specification,
    e.g., of colors lead to redundant escape sequences!

    The following tokens are recognized

    - Colors `r,red, g,green, b,blue, c,cyan, m,magenta, y,yellow,
    k,black, w,white`
    - Modifier `hi` if specified the highlight variant (of *any* color in
    the specification is used)
    - `bf, *, bold`
    - `u, _, underline`
    - `f, faint`
    - `0, normal`
    - `p, pos, positive, +`
    - `n, neg, negative, -`
    - `em, emph` (predefined ESC_EMPH)
    - `err, error` (predefined ESC_ERROR)

    The recognition is *not* case sensitive.

    \param _begin start of specification
    \param _end end of specification
    \param _unknown callback function for unknown primitives, use, e.g., for
    extending the parser
 */
std::string
get_ansi_escape_sequence(const std::string::const_iterator& _begin,
                         const std::string::const_iterator& _end,
                         std::function<void(std::string)> _unknown=
                         std::function<void(std::string)>());

//-----------------------------------------------------------------------------

/** Same as above but description is given as string \c _spec
    \ingroup vc_base_ansiesc
 */
inline std::string
get_ansi_escape_sequence(const std::string& _spec,
                         std::function<void(std::string)> _unknown=
                         std::function<void(std::string)>()) {
  return get_ansi_escape_sequence(_spec.begin(),_spec.end(),_unknown);
}

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_ANSI_ESCAPE_HH defined
