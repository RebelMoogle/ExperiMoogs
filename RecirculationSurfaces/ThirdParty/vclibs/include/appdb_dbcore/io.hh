//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_IO_HH
#define VC_APPDB_IO_HH

// included by appdb.hh

//== INCLUDES =================================================================

#include <iostream>

//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

/** \defgroup vc_appdb_io Read and write databases (AppDB)
    \ingroup vc_appdb

    Read and write data base contents from and to files.

    **Supported file formats**

    \arg \ref label_appdb_io_conf "Unix Conf-like file format"

    \arg Note that new file formats can be supported easily by using
    bindings to scripting languages, e.g., [\ref vc_appdb_rb].

    \todo convert to/from boost::options
    \todo convert to `boost::property_tree` (provides XML)
 */

/** \page vc_appdb_io_conf AppDB Conf file format
    \anchor label_appdb_io_conf
    \ingroup vc_appdb_io

    ### Lines

    - One _record_ per line, exception: line ending with `'\'`.

    - White space acts as _separator_, extra white spaces is ignored.

    - `#` starts a comment, the remainder of the line is ignored,
      use `\#` to escape.

    - A line `__END__` denotes the end of the file, subsequent lines
      are ignored.

    - A line `include{FILE}` will include `FILE` (currently no
      checking for circular includes, other than a simple counter), see
      also io::NoInclude, io::IgnoreInclude, io::SloppyInclude

    - Lines that contain `R"cxx"(` or `)cxx""` (exclusively,
      and each starting with beginning of line) are _ignored_. Thus
      the file could be included in a C++ file and compiled to a C++11
      raw string.

      (Note that the _prefix_ is fixed, and we don't do any further
      syntax checking. The extra quote `"` in the prefix help for
      proper syntax highlighting in conf files.)

    ### Records

    A _record_ is a pair `KEY = VALUE`, where `KEY` is a key starting
     with a slash `/` (see [\ref vc_appdb]). Additionally, a record
     may include a type specifier: `KEY = {{TYPE}} VALUE`. This is
     used

    - _only_ for reading to create a Value instance of a given type _and_
    - _only_ if the node does not yet have a value attached _and_
    - _only_ the io::CreateValues is specified.

    Here, `TYPE` refers to a type as defined by make::register_type() or
    make::register_alias().

    ### Values

    The remainder of the line is interpreted as `VALUE` with _spaces_ acting
    as _separators_.

    - The number of separate items determines the dimension of the value.
      It should match the given dimension in the Database.
      See also io::SloppyDimensions.

    - Use `{...}` or `"..."` or <tt>'...'</tt> for grouping and `\"`
      or `\\}` or `\'` (and `\\`) for escaping.

    - The standard C escape sequences (e.g., `\n`) are supported.

    ### Groups

    - A line `[PREFIX]` defines a path used as a prefix to all
      subsequent keys that do not begin with a slash `/`.

    - `PREFIX` defines a "group". It should be a valid key. (For
      convenience, a missing leading slash `/` is added automatically,
      and a trailing slash `/` is ignored.

    - After a `PREFIX` was set, all keys that do not start with `/`
      are treated _relative_ to `PREFIX`, which keys that start with
      `/` are just treated as before (as _absolute_ paths).

    - Defining a prefix _without_ any relative key following will
      _not_ automatically create nodes for the CreateNodes option.


    ### More

    - A string value of the form `=>KEY` (no space) can be used to
      **reference** another node `KEY`.

    - Use resolve_references() to create references to other nodes
      _after_ successfully reading input.

    \todo Multiline string with indentation ignored << ... >>
 */


/// Read and write databases \ingroup vc_appdb_io
namespace io {

//-----------------------------------------------------------------------------

/// specify Conf file format
enum ConfFormat {
  Conf  //!< argument to read
} ;

/// options to read()
enum Options {
  None=0,
  CreateNodes=1,         //!< create non-existent nodes on reading
  CreateValues=2,        //!< create non-existent values on reading (Value::String)
  SloppyDimensions=4,    //!< no error for _missing_ dimensions
  AtomicRead=8,          //!< read all or nothing on failure
  NoInclude=0x10,        //!< including a file results in an error
  IgnoreInclude=0x20,    //!< ignore any included files (no error)
  SloppyInclude=0x40,    //!< ignore errors from included files
  WriteDotNodes=0x80,    //!< write nodes with `key_tail()=~/\.*/` (like Unix dot files)
  WriteDotChilds=0x100,  //!< write children of "hidden dot nodes"
  WriteReferences=0x200, //!< write nodes that reference another node
  WriteGroups=0x400,     //!< write_if() and others use groups/prefixes
  Append=0x800,          //!< write_file() opens file for appending
  ReadSetsAlways=0x1000  //!< reading modifies even of old equals new value
};

/** compose an Option value from `_str`.
    \ingroup vc_appdb_io

    Every character in `_str` (case-sensitive) sets one
    option. Multiple settings are ignored. Invalid characters produce
    a warning on VC::base::Logger `appdb`.

    | character    | Constant             |
    |--------------|----------------------|
    | `'n'`        | io::CreateNodes      |
    | `'v'`        | io::CreateValues     |
    | `'D'`        | io::SloppyDimensions |
    | `'a'`        | io::AtomicRead       |
    | `'X'`        | io::NoInclude        |
    | `'i'`        | io::IgnoreInclude    |
    | `'I'`        | io::SloppyInclude    |
    | <tt>'.'</tt> | io::WriteDotNodes    |
    | <tt>':'</tt> | io::WriteDotChilds   |
    | `'r'`        | io::WriteReferences  |
    | ``+'`        | io::Append           |
    | `'g'`        | io::WriteGroups      |
    | `'s'`        | io::ReadSetsAlways   |

    - E.g., `"nv."` refers to `(CreateNodes|CreateValues|WriteDotNodes)`.
    - An empty string refers to io::None.

    \param _str describes options
    \return value
    \sa options_string()
 */
Options compose_options(const char* _str);

/// same as above \ingroup vc_appdb_io
inline Options compose_options(const std::string& _str) {
  return compose_options(_str.c_str());
}

/** Get string from `_options`.
    \param _options output warning on VC::base::Logger `appdb` for
    invalid value
    \return string as for compose_options()
    \sa compose_options()
 */
std::string options_string(Options _options);

/** \brief Read Conf file.
    \ingroup vc_appdb_io

    - No nodes are created unless io::CreateNodes applies, this option
      creates all nodes required to represent the given path.

    - No values are created unless io::CreateValues. This option is
      potentially **dangerous**: it uses make::sv() with Value::String
      to create values with the "right" dimension. This may not be
      what you intend!

    - Dimensions must match exactly unless io::SloppyInclude applies.

    - `read()` will set data read from each line in the conf file
      immediately unless io::AtomicRead applies. For an atomic read,
      the database remains untouched on any error. (This significantly
      slower.)

    - `read()` will output error and warning messages to
      VC::base::Logger `"appdb"`.

    \sa [\ref vc_appdb_io_conf]
    \param _in input stream
    \param _db database
    \param _options options
    \return success
 */
bool read(ConfFormat,std::istream& _in,
         AbstractDatabase* _db,Options _options=None);

/** Same as read() with key `_patterns` as predicates
    \ingroup vc_appdb_io
    \param _in input stream
    \param _db database
    \param _pattern defines a regular expressions for keys
    separated. Only nodes with matching keys are read (calls `regex_match()`!).
    There will be no syntax checks on the value for mismatched keys.
    \param _options options
    \return success
    \sa [\ref vc_appdb_io_conf]
 */
bool read_matching(ConfFormat,std::istream& _in,AbstractDatabase* _db,
                   const std::string& _pattern,
                   Options _options=None);

/** Same as read() but open file `_filename` for reading.
    \ingroup vc_appdb_io
    \param _filename file is opened by  base::utility::get_istream().
    \param _db database
    \param _options options
    \return success
*/
bool read_file(ConfFormat,const std::string& _filename,
               AbstractDatabase* _db,Options _options=None);

/** Same as read_matching() but open file `_filename` for reading.
    \ingroup vc_appdb_io
    \param _filename file is opened by  base::utility::get_istream().
    \param _db database
    \param _pattern
    \param _options options
    \return success
*/
bool read_matching_from_file(ConfFormat,const std::string& _filename,
                             AbstractDatabase* _db,
                             const std::string& _pattern,Options _options=None);

/** Read a single `KEY = VALUE` line \ingroup vc_appdb_io.
    `_line` **must not** include comments, `__END__` or `include{FILE}`
    directives.
    \param _line the data
    \param _db database
    \param _options options
    \return success
    \sa [\ref vc_appdb_io_conf]
 */
bool read_line(ConfFormat,const std::string& _line,
               AbstractDatabase* _db,Options _options=None);

/** Parse `VALUE` string `_value` and set `_node` value.
    \param _node node
    \param _value `VALUE` part of `KEY = VALUE` record
    \param _options options
    \return success
 */
bool read_value(ConfFormat,Node* _node,
		const std::string& _value,Options _options=None);

/** Write Conf file. \ingroup vc_appdb_io
    database `_db` as Conf file to `_out`.
    \param _out output stream
    \param _db database
    \param _options options (only some are recognized, e.g.,
    WriteDotNodes, WriteDotChilds)
    \return `_out`
    \sa [\ref vc_appdb_io_conf]
 */
std::ostream&
write(ConfFormat,std::ostream& _out,
      AbstractDatabase* _db,Options _options=None);

/** Same as write() with additional `_predicate`
    \ingroup vc_appdb_io
    \param _out output stream
    \param _db database
    \param _predicate write only values of nodes for which `_predicate`
    evaluates to `true`
    \param _options options (only some are recognized, e.g.,
    WriteDotNodes, WriteDotChilds)
    \return `_out`
    \sa [\ref vc_appdb_io_conf]
    \todo option WriteGroups does not sort entries efficiently ("preorder")
 */
std::ostream&
write_if(ConfFormat,std::ostream& _out,
         AbstractDatabase* _db,std::function<bool(const Node*)> _predicate,
         Options _options=None);

/** Same as write_if() with key `_patterns` as predicates
    \ingroup vc_appdb_io
    \param _out output stream
    \param _db database
    \param _patterns define a list of regular expressions for keys
    separated by `':'`; only nodes with keys matching _any_ expression are
    written.
    \param _options options (only some are recognized, e.g.,
    WriteDotNodes, WriteDotChilds)
    \return `_out`
    \sa [\ref vc_appdb_io_conf]
 */
std::ostream&
write_matching(ConfFormat,std::ostream& _out,
               AbstractDatabase* _db,const std::string& _patterns,
               Options _options=None);

/** Combines write_if() and write_matching() for convenience.
    \ingroup vc_appdb_io
    Form predicate by logical *and* `_predicate` (first test) and test for
    matching `_patterns`.
 */
std::ostream&
write_matching_if(ConfFormat,std::ostream& _out,
                  AbstractDatabase* _db,const std::string& _patterns,
                  std::function<bool(const Node*)> _predicate,
                  Options _options=None);

/** Same as write() but open `_filename` for writing.
    \ingroup vc_appdb_io
    \param _filename file
    \param _db database
    \param _options options (WriteDotNodes, Append)
    \return success
*/
bool write_file(ConfFormat,const std::string& _filename,
                AbstractDatabase* _db,Options _options=None);

/** Same as write_matching() but open `_filename` for writing.
    \ingroup vc_appdb_io
    \param _filename file
    \param _db database
    \param _patterns
    \param _options options (WriteDotNodes, Append)
    \return success
*/
bool write_matching_to_file(ConfFormat,const std::string& _filename,
                            AbstractDatabase* _db,const std::string& _patterns,
                            Options _options=None);

/** Get string representation of `_value`.
    \ingroup vc_appdb_io
    Returns a string representation for `_value` as output by write_if() and
    readable by read_value():
    \arg combines dimensions, and
    \arg insert separators and/or escape characters if required.
    \param _value may have multiple dimensions
    \return string representation
  */
std::string format(ConfFormat,Value* _value);

/** Reinterpret values to create references. \ingroup vc_appdb_io
    \arg Searches `_db` for nodes with `string` values (of size 1) of the form
    `"=>KEY"`. _Note that there is no whitespace!_
    (The node may not reference another node!)
    \arg If `KEY` refers to a node create a reference to this node, which will
    "shadow" the value. (The value itself remains untouched.)
    \arg Otherwise report an error.
    \return `true` if there were no errors
 */
bool resolve_references(AbstractDatabase* _db);

//=============================================================================
} // namespace io
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_IO_HH defined
