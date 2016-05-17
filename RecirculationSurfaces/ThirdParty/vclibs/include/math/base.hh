//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision: 6 $
//
//=============================================================================

#ifndef VC_BASE_HH
#define VC_BASE_HH

#include <string>
#include <memory>
#include <vector>
#include <istream>
#include <functional>

/** \def _VC_NO_BOOST_REGEX
    \brief Don't compile features that require linking Boost regex library.
    \ingroup vc_base
 */

/** \def _VC_NO_BOOST_STREAMS
    \brief Don't compile features that require linking Boost stream i/o libraries.
    \ingroup vc_base
 */

namespace VC {
/// basic tools [\ref vc_base] \ingroup vc_base
namespace base {

/** \defgroup vc_base Basic tools
    \sa VC::base

    ## Environment variables

    | variable             | effect                                                        |
    |----------------------|---------------------------------------------------------------|
    | `LOG`                | sets levels,... for \ref vc_base_debug "logging/debug output" |
    | `VC_FORCE_NOT_A_TTY` | forces no \ref vc_base_tty "tty"                              |
    | `VC_MARKUP_STYLE`    | define style for \ref vc_base_markup "markup"                 |
    | `EDITOR`             | to set config::editor::tty::editor                            |
    | `VC_EDITOR`          | to set config::editor::gui::editor                            |
    | `TERM`               | to set config::terminal::term                                 |
    |                      |                                                               |

    Specify "no tty", e.g., to control output of
    \ref vc_base_ansiesc "escape codes".

    \todo helpers like not_copyable

    \example vclibs/base/test_log.cc
    Illustrates basic usage of VC::base::Logger.

    \example vclibs/base/test_backtrace.cc
    Applies debugging tools, VC::base::cpu_time, and
    VC::base::backtrace_t. (Note that debug symbols are only
    shown for backtraces from calls to libraries, hence the
    \c extern functions.

    \example vclibs/base/test_memory.cc
    Compares use of small_vector and bounded_small_vector to
    std::vector.

    \example vclibs/base/test_printf.cc
    [\ref vc_base_printf]

    \example vclibs/base/test_compressed_streams.cc
    [\ref vc_base_utility]

    \example vclibs/base/test_ansi_escape.cc
    [\ref vc_base_ansiesc]

    \example vclibs/base/test_markup.cc
    [\ref vc_base_markup]

    \example vclibs/base/test_tty.cc
    [\ref vc_base_tty]

    \example vclibs/base/test_primtype.cc
    [\ref vc_base_primtype]
*/

/** \defgroup vc_base_config Basic configuration.
    \ingroup vc_base
*/

/// Basic configuration [\ref vc_base_config] \ingroup vc_base_config
namespace config {
  /// preferred editor command \ingroup vc_base_config
  namespace editor {
    namespace tty {

      /** Get editor command. \ingroup vc_base_config
          Use value of `editor`, try to set `editor` if unset.
          \return  Returns command _without_ resoling the editor's path,
          `%s` the placeholder is inserted as file name.
       */
      std::string command();

      extern const char* editor; //!< default is `getenv(EDITOR)`
      extern const char NANO[];
      extern const char VI[];
      extern const char VIM[];
      extern const char EMACS[];
      extern const char EMACSCLIENT[];
    } // namespace tty
   namespace gui {
      extern const char* editor; //!< default is `getenv(VC_EDITOR)`
      /** If available same as `editor`.
          Don't detach process automatically but wait until editor is closed.
          This feature may not be available (i.e., `editor_wait=0`.).
      */
      extern const char* editor_wait;
      extern const char GEDIT[];
      extern const char GVIM[];
      extern const char GVIM_WAIT[];
      extern const char EMACS[];
      extern const char EMACSCLIENT[];
      extern const char EMACSCLIENT_WAIT[];
      extern const char NOTEPAD[];
      extern const char NOTEPADPP[];

      /** Get editor command. \ingroup vc_base_config
          Use value of `editor`, try to set `editor` if unset.
          \return  Returns command _without_ resoling the editor's path,
          `%s` the placeholder is inserted as file name.
          \param _wait use `editor_wait` if available
       */
      std::string command(bool _wait=false);

   } // namespace gui
  } // namespace editor

  /// preferred terminal \ingroup vc_base_config
  namespace terminal {
     extern const char* term; //!< terminal command (default: `getenv(TERM)`)
     extern const char XTERM[];
     extern const char GNOME_TERMINAL[];
     extern const char KONSOLE[];

     /** Get terminal command. \ingroup vc_base_config
         Use value of `term`, try to set `term` if unset.
         \return  Returns command _without_ resoling the terminal's path,
         `%s` the placeholder is inserted as command to be executed by the
         terminal.
       */
      std::string command();
  } // namespace terminal

} // namespace config

/** \defgroup vc_base_utility Collection of small utilities.
    \ingroup vc_base

    - ensure_binary_stream() deal with Windows: force standard streams
      opened in `ios::binary` mode

    - get_istream(), get_ostream() provide (open) streams from file
      names with extra functionality like handling of standard
      streams, compression,...

    - compressed(), decompressed() provide compressed streams
*/

/// Utilities [\ref vc_base_utility] \ingroup vc_base_utility
namespace utility {

/// standard integer file descriptors \ingroup vc_base_utility
enum standard_fileno_t {
  STDIN=0,
  STDOUT=1,
  STDERR=2
};

/** Ensure that `_fd` is opened in `ios::binary` mode.
    \ingroup vc_base_utility

    Applies only to Windows, which uses different conventions for text
    mode and binary mode, see, e.g., <a
    href="stdcxx.apache.org/doc/stdlibug/30-4.html">here</a>.

    \throw VC::base::vc_runtime_error on failure
 */
void ensure_binary_stream(standard_fileno_t _fd);

/** Construct new input stream from `_filename`.
    \ingroup vc_base_utility

    - Files are opened as `std::ifstream(_filename,ios::binary)`.
    - `_filename=="-"` refers to `stdin`
    - `_filename=="string:DATA"` creates an
      `std::istringstream("DATA")`, where some standard escape
      sequences in `DATA` will be substituted.
    - `_filename=="cmd:COMMAND"` will read from `popen("COMMAND")`.
    - The returned `shared_ptr` takes care of closing the stream,
      except for `std::cin` which will not be closed.
    - Detects and handles compressed streams by decompressed() (where
      compression_from_stream() decides on `compression_t`).

    \throw    VC::base::vc_runtime_error on failure
*/
std::shared_ptr<std::istream> get_istream(const std::string& _filename);

/** Construct new output stream from `_filename`.
    \ingroup vc_base_utility

    - Files are opened as `std::ofstream(_filename,ios::binary)`.
    - `_filename=="-"` refers to `stdout`
    - The returned `shared_ptr` takes care of closing the stream,
      except for `std::cout` which will not be closed.
    - Handles compression of output based on file name suffix
      (compression_from_suffix() decides on compression, handled by
      compressed()).

    \throw VC::base::vc_runtime_error on failure
*/
std::shared_ptr<std::ostream> get_ostream(const std::string& _filename);

//-----------------------------------------------------------------------------

/// generator for `istream` factory
typedef std::function<std::shared_ptr<std::istream>(const std::string&)>
    istream_generator_t;
/// register generator with factory
bool register_istream_generator(const istream_generator_t& _g);

//-----------------------------------------------------------------------------

/** Cache `istream` contents.

    Read data once and store as string. The **prefix** `cache:` is
    recognized by

    - get_istream() which returns an `istringstream`
    - read_file() which returns the stored string.

    All `_filename` arguments to methods of `istream_cache_t` are without
    prefix `cache:`.
 */
class istream_cache_t {
  struct map_t;
  istream_cache_t();
  std::unique_ptr<map_t> m_map;
public:
  istream_cache_t(const istream_cache_t&) = delete;
  istream_cache_t& operator=(const istream_cache_t&) = delete;

  /// get singleton instance
  static istream_cache_t& instance();

  /// clear cache
  void clear();
  /// remove `_filename` and contents from cache; no effect if `!is_cached()`
  void erase(const std::string& _filename);
  /// test if `_filename` is cached
  bool is_cached(const std::string& _filename) const {
    return get(_filename)!=nullptr;
  }
  /// get cached data or `nullptr`
  const std::string* get(const std::string& _filename) const;
  /// insert `_filename`
  const std::string& insert(const std::string& _filename);
  /// get `istream`, calls `insert`
  std::shared_ptr<std::istream> get_istream(const std::string& _filename);
};

//-----------------------------------------------------------------------------

/** Read all remaining data from `_in` as `string`.
    \ingroup vc_base_utility

    - _in input, will `unsetf()` `std::ios::skipws`
    \return text (no error handling, must check `_in`)
 */
inline std::string read_all(std::istream& _in) {
  _in.unsetf(std::ios::skipws);      // No white space skipping!
  return std::string(std::istreambuf_iterator<char>(_in.rdbuf()),
                     std::istreambuf_iterator<char>());
  // no need to reset flag
}

/** Call read_all() for get_istream().
    \ingroup vc_base_utility

    Note: this function may
    \param _filename
    \return data
    \throw VC::base::vc_runtime_error on failure (tests `iso::bad()`)
 */
std::string read_file(const std::string& _filename);

//-----------------------------------------------------------------------------

/// stream compression algorithms \ingroup vc_base_utility
enum compression_t {
  NONE=0,
  GZIP,
  BZIP2,
  ZLIB
};

/** Get `compression_t` from `_filename` suffix.
    \ingroup vc_base_utility

    Suffixes are treated *case sensitive*.

    | suffix   | `compression_t` |
    |----------|-----------------|
    | `".gz"`  | GZIP            |
    | `".bz2"` | BZIP2           |
    | `".z"`   | ZLIB            |

    \param _filename file name (eventually with suffix/extension)
    \param _stripped file name with compression suffix removed (e.g.,
    `"a.txt.gz" => "a.txt"`) or original `_filename` if `NONE` was
    returned

    \return compression type or `NONE`
 */
compression_t compression_from_suffix(const std::string& _filename,
                                      std::string* _stripped=nullptr);

/** Remove compression suffix from `_filename`.
    \ingroup vc_base_utility

    \param _filename
    \return `_stripped` from call to compression_from_suffix()
*/
std::string remove_compression_suffix(const std::string& _filename);

/** Get `compression_t` from stream `_in`.
    \ingroup vc_base_utility

    \param _in stream will be examined by reading magic pattern and
    followed by `putback()`.

    Note that a call to compression_from_suffix() will _block_ if
    there is no data available and not VC::base::isa_tty()!

    \return compression type or `NONE`
 */
compression_t compression_from_stream(std::istream& _in);

/** Construct new input stream from `_in` with decompression.
    \ingroup vc_base_utility
    \param _in original input stream, compressed data will be read from here
    \param _c compression type
    \return new stream or `_in` if `_c==NONE`
 */
std::shared_ptr<std::istream> decompressed(std::shared_ptr<std::istream> _in,
                                           compression_t _c);

/** Construct new input stream from `_in` with decompression.
    \ingroup vc_base_utility

    Calls compression_from_stream() to detect `compression_t`.
    \param _in original input stream, compressed data will be read from here
    \return new stream or `_in` if `_c==NONE`
 */
std::shared_ptr<std::istream> decompressed(std::shared_ptr<std::istream> _in);

/** Construct new output stream from `_out` with compression.
    \ingroup vc_base_utility
    \param _out original output stream, compressed data will be written here
    \param _c compression type
    \return new stream or `_out` if `_c==NONE`
 */
std::shared_ptr<std::ostream> compressed(std::shared_ptr<std::ostream> _out,
                                         compression_t _c);

//-----------------------------------------------------------------------------

/** Read CVS data from `_in` and pass each record to `_process_row`.

    - There is no strict syntax matching. Errors are silently ignored, and
      we try to do something reasonable.
    - Parsing fields and related error checking is delegated to `_process_row`.
    - Whitespace is ignored other than in quoted strings.
    - There is no character escaping.

    `parse_csv()` splits each input line, or data *row*, into text
    fields that are passed to `_process_row`.

    - The first argument to `_process_row` is the current *line number*
      (starting at 1 and counting empty and comment lines).
    - The second argument contains the fields of the current row. Note
      that the number of fields is not checked and may vary!
    - Returning `false` stops parsing.

    \param _in input stream
    \param _process_row consumes data one row at a time
    \param _delimiters list of *characters* that act as delimiters
    \param _quotes list of *characters* that act as quotes
    \param _comments list of *characters* that start a comment line: must be
    first non-whitespace character in line
 */
void
parse_csv(std::istream& _in,
          std::function<bool(size_t,
                             const std::vector<std::string>&)> _process_row,
          const std::string& _delimiters=",",
          const std::string& _quotes="\"",
          const std::string& _comments="#");

//-----------------------------------------------------------------------------

/** Test if `_text` starts with `_pattern`.
    \ingroup vc_base_utility
*/
bool starts_with(const std::string& _text,const std::string& _pattern);

/** Test if `_text` ends with `_pattern` and `_text.size()>_pattern.size()`.
    \ingroup vc_base_utility
*/
bool ends_with(const std::string& _text,const std::string& _pattern);

//-----------------------------------------------------------------------------

} // namespace utility
} // namespace base
} // namespace VC

#endif // VC_BASE_HH
