//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_MARKUP_HH
#define VC_BASE_MARKUP_HH


//== INCLUDES =================================================================

#include <string>
#include <memory>

#include <boost/regex.hpp>

/** \defgroup vc_base_markup Simple markup for text rendering.
    \ingroup vc_base

    VC::base::markup::Scanner and VC::base::markup::Formatter
    implement a minimal markup engine for simple text rendering.

    \sa [\ref VC::base::markup], [\ref vc_base_ansiesc]

    \bug get_terminal_size() seems to return wrong values (sometimes)
    \bug *FIXED!?! (wrong input)* check for missing spaces in text
    \bug *FIXED!?!* `"*text*,"` fails (fix also string in 
    VC::appdb::args::markup_summary())
*/

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {

/// minimal markup engine \ingroup vc_base_markup
namespace markup {

class Formatter;

/// text recognized by Scanner \ingroup vc_base_markup
typedef std::pair<std::string::const_iterator,
                  std::string::const_iterator> text_t;

/// outline state of Scanner \ingroup vc_base_markup
enum OutlineState { Paragraph,  List };

/// text style state of Scanner \ingroup vc_base_markup
enum TextStyle { Normal, Emphasis, Bold, Typewriter };

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Implement simple text markup.

    The scanner breaks a text into tokens, which are output by an instance of
    Formatter. The scanner does not provide any sophisticated error
    checking.

    Currently, the following markup tags are supported (`N` and `M` are
    decimal digits>=1, and <tt>@bol</tt> refers to the beginning of a line).

    Text structure

    \arg <tt>.hN</tt> (<tt>@bol</tt>) defines a __heading__ of level `N` (until
    end of line)
    \arg <tt>*</tt> (<tt>@bol</tt>) defines an item in a __bullet list__, the 
    level is given by the number of `*`.
    \arg <tt>.p</tt> or <tt>.pN</tt> or <tt>.pN:M</tt> start a new 
    __paragraph__ with indentation level `N` and `M` for all following lines.

    Text styles

    \arg <tt>*bold*</tt>
    \arg <tt>_emphasis_</tt>
    \arg `@typewriter@`
    \arg `.style`  defines a style sheet (until end of line), the format 
    depends on the Formatter

    \sa [\ref Formatter]
 */
class Scanner {
public:
  Scanner(std::shared_ptr<Formatter> _formatter=
                std::shared_ptr<Formatter>((Formatter*) 0))
    : m_formatter(_formatter),
      m_os(Paragraph), m_ts(Normal), m_level(0) {
    m_start=m_stop;
  }

  /// set formatter
  void setFormatter(std::shared_ptr<Formatter> _formatter) {
    m_formatter=_formatter;
  }

  /// scan _text and report to current Formatter
  void scan(const std::string& _text);

  /// current text style
  TextStyle style() const { return m_ts; }

  /// current outline state
  OutlineState outline() const { return m_os; }

  /// current outline level
  int level() const { return m_level; }

private:

  bool scan_heading();
  bool scan_bullet();
  bool scan_newline();
  bool scan_line();
  bool scan_hrule();
  bool scan_paragraph();
  bool scan_style();

  bool scan_words(std::string::const_iterator _start,
                  std::string::const_iterator _stop);
  bool switch_style(int _);

  std::string::const_iterator m_start, m_stop;
  boost::smatch               m_match;
  std::shared_ptr<Formatter> m_formatter;
  OutlineState                m_os;
  TextStyle                   m_ts;
  int                         m_level;
};

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Output tokens recognized by Scanner

    This class provides an interface for the actual text rendering.

    \sa [\ref Scanner]
 */
class Formatter {
public:
  Formatter() : m_scanner(0) {}

  /** @name callbacks
      @{
   */
  virtual void startScan(const Scanner* _scanner);
  virtual void heading(int _level,text_t _text) = 0;
  virtual void bullet(int _level) = 0;
  virtual void paragraph(int _level,bool _force) = 0;
  virtual void word(text_t __text) = 0;
  virtual void style(TextStyle _prev,TextStyle _cur) = 0;
  virtual void hrule() = 0;
  virtual void styleSheet(const std::string& _style);
  virtual void stopScan();

  /// @}

protected:
  const Scanner* m_scanner;
};

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Output debug messages, strictly not a formatter
    \sa [\ref vc_base_debug]
 */
class FormatDump : public Formatter {
public:
  FormatDump() {}
  virtual ~FormatDump();

  virtual void heading(int _level,text_t _text);
  virtual void bullet(int _level);
  virtual void paragraph(int _level,bool _force);
  virtual void word(text_t _text);
  virtual void style(TextStyle _prev,TextStyle _cur);
  virtual void styleSheet(const std::string& _style);
  virtual void hrule();
};

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Output reformatted text with markup codes.
 */
class FormatSelf : public Formatter {
public:
  /// write output to \a _out
  FormatSelf(std::ostream& _out) : m_out(_out), m_needspace(false) {}
  virtual ~FormatSelf();

  virtual void startScan(const Scanner* _scanner);
  virtual void heading(int _level,text_t _text);
  virtual void bullet(int _level);
  virtual void paragraph(int _level,bool _force);
  virtual void word(text_t _text);
  virtual void style(TextStyle _prev,TextStyle _cur);
  virtual void hrule();
  virtual void styleSheet(const std::string& _style);
  virtual void stopScan();

protected:
  std::ostream& m_out;
  bool m_needspace;
};

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Render ASCII text to terminal window using ANSI escape sequences.

    Tries to query the width of the terminal window (`AutoWidth`, use fixed
    width otherwise) and format text using white space indentation and
    ANSI escape sequences (if `AlwaysColors` or `AutoColors` and
    VC::base::isa_tty()).

    Various properties of the output can be set in a
    MarkupformatText::StyleSheet instance, which can be provided via
    setStyleSheet() or from within the text using the <tt>.style</tt>
    primitive.

    Examples for defining styles 
    (see also VC::base::get_ansi_escape_sequence() for definitions in `{}`)

    \arg define style
    \code
    .style h1:{u,*,yellow},bf:{white,*},em:{green,hi},tt:{m,hi}
    \endcode

    \arg reset/load default style
    \code
    .style reset
    \endcode

    The __environment variable__ `VC_MARKUP_STYLE` can be used to define 
    a default style, e.g.
    \code
    export VC_MARKUP_STYLE="h1:{u,*,yellow},bf:{white,*},em:{green,hi},tt:{m,hi}"
    \endcode

    \sa [\ref vc_base_ansiesc]
 */
class FormatText : public Formatter {
public:
  /// argument to setColors() [\ref vc_base_ansiesc]
  enum UseColors { NoColors=0, AlwaysColors=1, AutoColors=-1 };
  /// argument to fixWidth(): determine width from terminal settings
  enum Width { AutoWidth=-1 };

  /** \ingroup vc_base_markup
      \brief Style sheet for text rendering.

      The syntax for specification to apply() (or the constructor) is
      as follows

      \arg Tokens are given in a comma separated list.
      \arg Each token has the form <tt>TOKEN:{DATA}</tt>.
      \arg `TOKEN` is one of <tt>em,emph,_, b,bold,bf,*, tt,@, hD, ulD</tt>, where
      `D` is a decimal digit.
      \arg `DATA` is a comma separated list as for
      VC::base::get_ansi_escape_sequence() with the following extensions:
      \arg `'C'` defines the bullet character in `ulD` as a single character
      `C`.
      \arg `upcase` in `hD` renders headings in upper case letters.
      \arg <em>Unknown tokens are ignored.</em>
      \sa [\ref FormatText]
   */
  class StyleSheet {
  public:
    /// use default styles
    StyleSheet() { reset(); }
    /// load style from `_style`
    StyleSheet(const std::string& _style);
    /// load default style
    void reset();
    /// apply `_style` to current settings
    void apply(const std::string& _style);

    enum { MaxHeadings=3, MaxBullets=3 };

    std::string h[MaxHeadings];   //!< headings style
    std::string ul[MaxBullets];   //!< bullets style
    std::string emph, bold, tt;
    bool        hup[MaxHeadings]; //!< headings upper case?
    char        ulch[MaxBullets]; //!< bullet character
  private:
    int  m_level;      //!< temporary state during string parsing
    void seth(const std::string& _item);  //!< callback setting hup
    void setul(const std::string& _item); //!< callback setting ulch
  };

  /** Constructor
      \param _out output stream
      \param _colors use ANSI escape sequences? `AutoColors` queries
      VC::base::isa_tty()
      \param _width output width, `AutoWidth` queries the terminal
      setting if possible
   */
  FormatText(std::ostream& _out,
# if defined VC_PLATFORM_WINDOWS
             UseColors _colors=NoColors,
# else
             UseColors _colors=AutoColors,
# endif

             int _width=AutoWidth);
  virtual ~FormatText();

  /// enable or disable use of ANSI escape sequences
  void setColors(UseColors _colors);
  /// set output width (`_w` may be set to `AutoWidth`)
  void fixWidth(int _w);
  /// set right margin for all text
  void setRightMargin(int _n);
  /// set left margin for all text except <tt>.h1</tt> headings
  void setLeftMargin(int _n);
  /// access the style sheet
  StyleSheet getStyleSheet() { return m_style; }

  virtual void startScan(const Scanner* _scanner);
  virtual void heading(int _level,text_t _text);
  virtual void bullet(int _level);
  virtual void paragraph(int _level,bool _force);
  virtual void word(text_t _text);
  virtual void style(TextStyle _prev,TextStyle _cur);
  virtual void hrule();
  virtual void styleSheet(const std::string& _style);
  virtual void stopScan();

protected:

  std::ostream& color(const std::string& _color);
  std::ostream& basic_indent();

  std::ostream&  m_out;
  StyleSheet     m_style;

  int  m_width;     //!< specified width (may be AutoWidth)
  int  m_curWidth;  //!< "real" with for current scan
  int  m_col;       //!< current column
  int  m_right;     //!< right margin
  int  m_left;      //!< left margin (except .h1 headings)
  bool m_tty;
  bool m_colors;
  bool m_neednewline;
  bool m_needspace;
};

//-----------------------------------------------------------------------------

/** \ingroup vc_base_markup
    \brief Output (possibly incorrect) HTML code.

    Very simplistic HTML output without checks. Intended to create small
    snippets of, e.g., help text.
 */
class FormatHTML : public Formatter {
public:
  /// write output to \a _out
  FormatHTML(std::ostream& _out) : m_out(_out), m_needspace(false) {}
  virtual ~FormatHTML();

  virtual void startScan(const Scanner* _scanner);
  virtual void heading(int _level,text_t _text);
  virtual void bullet(int _level);
  virtual void paragraph(int _level,bool _force);
  virtual void word(text_t _text);
  virtual void style(TextStyle _prev,TextStyle _cur);
  virtual void hrule();
  //virtual void styleSheet(const std::string& _style);
  virtual void stopScan();

protected:
  void openUL(int _level);
  void closeUL(int _level);

  std::ostream& m_out;
  int  m_lastLevel;
  int  m_ul;
  bool m_needspace;
};

//=============================================================================
} // namespace markup
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_MARKUP_HH defined
