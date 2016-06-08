//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_ARGS_HH
#define VC_APPDB_ARGS_HH

// included by appdb.hh

//== INCLUDES =================================================================



//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

//-----------------------------------------------------------------------------

/** \defgroup vc_appdb_args Command line parsing (AppDB)
    \ingroup vc_appdb

    \arg Use parse() for parsing a command line.
    \arg A usage summary is extracted by markup_summary(). The display of this
    summary is handled by parse() (calling print_usage_and_exit(), which calls
    markup_summary(), which uses get_help()).
    \arg The utility function complain_and_exit() can be used to display a
    ParseArgsException.
    \arg Diagnostic messages are may be written to VC::base::Logger `"appdb"`.

    \example vclibs/appdb/examples/example_args.cc
 */

/// Command line parsing (AppDB) \ingroup vc_appdb_args
namespace args {

/// default argument to get_help(), markup_summary() \ingroup vc_appdb_args
extern std::string DEFAULT_DESCRIPTION;

/** Get description for `_node`.
    \ingroup vc_appdb_args
    The following substitutions are applied (if values exist):
    \arg `"%T"` is replaced by the value's type
    \arg `"%S"` is replaced by `"[N]"`, where `N` is the size _or_
    by `""` if `N<2`
    \arg `"%V"` is replaced by the current value

    \param _node look for help in node `_node->key()+"/.help"`
    \param _brief get brief description if `true`
    \param _default_desc used if no description given
    \return description with substitutions applied, may be empty
    \sa see VC::appdb::make::help()
 */ 
std::string get_help(const Node& _node,bool _brief=false,
                     const std::string& _default_desc=DEFAULT_DESCRIPTION);

/** Get  summary for children of `_args`. 
    \ingroup vc_appdb_args
    \param _args find descriptions of children
    \param _brief use brief descriptions (VC::appdb::make::help())
    \param _default_desc used if no description given
    \arg Uses get_help() to obtain  descriptions
    \arg Prints headings and message if there exists help for `_args`
    \return text to be processed by VC::base::markup::Scanner
    \sa VC::appdb::import_args(), VC::appdb::make::help()
 */
std::string markup_summary(Node* _args,bool _brief=false,
                           const std::string& _default_desc=DEFAULT_DESCRIPTION);

/// Same as above but specify node with key `_args` \ingroup vc_appdb_args
std::string markup_summary(AbstractDatabase* _db,const std::string& _args,
                           bool _brief=false,
                           const std::string& _default_desc=DEFAULT_DESCRIPTION);


/** Parse command line arguments from `_argc` and `_argv`.
    \ingroup vc_appdb_args
   
    \arg Arguments refer to nodes that are children of
    `_parent`. Their values are set, and their information is used,
    e.g., to generate help.

    \arg Creates child node `"help"` (appdb::Event) which triggers 
    print_usage_and_exit(). The utility function make::help() can be
    used to create help.

    \arg The option `--help` calls print_usage_and_exit() to display
    help.

    \arg Creates child node `".argv"` which stores a copy of all arguments 
    as passed.

    \arg Creates child node `".leading"` which takes all of `_argc+1`
    before the first argument (leading \c '-') \a or the \c @ directive. A
    single dash \c '-' is _not_ considered an argument.

    \arg Creates child node `".trailing"` which takes all of `argc`
    after \c '--' (double dash).  These are not interpreted by parse().

    \arg A _complete_ database can be read from a file by io::read() using
    \code
    @ filename
    \endcode
    or
    \code 
    @filename
    \endcode
    on the command line.

    \arg The special option `--set-value KEY=VALUE` sets _any_ value in the
    database (provide full key). This feature may be a potential a potential
    annoyance!

    \arg The special option `--dump-database=FILE` writes the current
    contents of the database to `FILE`. This is useful for creating a 
    "template" for a response file but may be a potential annoyance!

    \arg The default _long_ options use the syntax `--option=VALUE`.
    Only for _boolean_ values and _events_ (appdb::Event), the
    `=VALUE` part may be omitted; without specification the values is
    set to `true`. The format of `VALUE` must conform to
    io::read_value(); this applies especially to vector valued data.

    \arg _Short_ options are realized as child nodes name `".X"` where
    `X` is a single character. These nodes are supposed to _reference_
    the respective "long option nodes"; short_opt() may be used to
    create short options conveniently. Short options may be combined,
    e.g. as `-xzf filename` for `-x` and `-z` referring to a boolean
    value (effect: set `true`) or an event (effect: _trigger_ event).

    \arg On failure diagnostic messages are output to VC::base::Logger
    `"appdb"` and a ParseArgsException is thrown. The utility function
    complain_and_exit() can be used to process the exception.

    \arg Note that io:read_value() is used for reading and
    converting arguments. This means that whitespace _splits_ strings, 
    even though _shell quotes_ are used, e.g.,
    \code
    --data="1 2 3"
    \endcode
    is treated as a value of Value::size() 3. If a instead single 
    string is expected as input (`size()==1`), use additional quotes or 
    braces for `read_value()`, e.g.,
    \code
    --data="\"1 2 3\""
    \endcode
    or
    \code   
    --data="{1 2 3}"
    \endcode
    _Please keep in mind that the outer level quotation is substituted by 
    the shell._    

    Example for calling parse()
    \code
    try {
      args::parse(argc,(const char**) argv,&db,"/local/args",true);
    }
    catch (ParseArgsException& e) {
      args::complain_and_exit(e);
    }
    \endcode

    Example for a command line
    \code
    example_args a b c --atol=111 --method c --eval -em l \
    --set-value /local/param/list=1 -- d e f
    \endcode
    \arg sets `.leading={a,b,c}` (size 3)
    \arg sets `atol`
    \arg sets `method`
    \arg triggers `eval`
    \arg triggers `eval` again (short option `-e`)
    \arg sets method again (short option `-m`)
    \arg sets `/local//param/list=1`
    \arg sets `.trailing={d,e,f}` (size 3)

    \param _argc # of arguments (as passed to `main()`)
    \param _argv # arguments (as passed to `main()`)
    \param _parent nodes storing/referencing arguments reside below
    \param _noadd don't add nodes; use for simple argument checking

    \sa markup_summary()
 */
void parse(int _argc,const char** _argv,Node* _parent,bool _noadd=true)
  throw (ParseArgsException);

/** Same as above, `_parent` node's path is created if necessary. 
    \ingroup vc_appdb_args
 */
void parse(int _argc,const char** _argv,
           AbstractDatabase* _db,const std::string& _parent,bool _noadd=true)
  throw (ParseArgsException);

/** Print usage information. \ingroup vc_appdb_args

    Print information for `_parent` to `cerr` and
    exits with exit code `-1`. Calls markup_summary().
 */
void print_usage_and_exit(Node* _parent);

/** Prints `_e` and exits. \ingroup vc_appdb_args
    Prints to `cerr` and exists with exit code `-1`.
    Calls VC::base::write_ansi_escape().
 */
void complain_and_exit(const ParseArgsException& _e);

/** Create short option `_opt` for `_node`. \ingroup vc_appdb_args
    \sa parse()
 */
void short_opt(Node& _node,char _opt);

//=============================================================================
} // namespace args
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_ARGS_HH defined
