//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_LOG_HH
#define VC_BASE_LOG_HH


//== INCLUDES =================================================================

#include <string>
#include <iostream>
#include <sstream>
#include <cassert>

# if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4244)
# endif

#include <mutex>

# if defined(_MSC_VER)
#  pragma warning (pop)
# endif

#include "backtrace.hh"

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {
//-----------------------------------------------------------------------------

/** \class Logger
    \brief Log events depending on level().
    \ingroup vc_base_debug

    \arg Minimizes run-time overhead in case of logging is switched off (level()
    does not match).

    \arg Note: data is copied in expression templates (otherwise we can't
    output temporary objects like "a+b"). Wrap data by Logger::ref() to force
    use of const references instead.

    \arg Synchonizes output to same Logger between threads.

    \arg Inspired by \c <a
    href="http://log4r.rubyforge.org">Log4r</a>, the manual section
    part on <a href="http://log4r.rubyforge.org/manual.html#art"> the
    art of logging</a> is worth reading.

    \arg Setup of loggers from external parameters can eb done easily using
    Logger::set_properties().

    Usage
    \code
    Logger mylog("mylog");
    Logger* log=Logger::get("mylog");      // query registry by name
    assert(log==&mylog);

    mylog.set_output(Logger::null());      // no output
    mylog.set_output(cerr);                // default

    mylog.set_level(Logger::OFF);          // no logging
    mylog.set_level(Logger::ALL);          // maximum logging
    mylog.set_level(Logger::WARN);         // default

    //
    // Log events are output using "debug()=" ("warn()=",...) and
    // concatenated by ",". Trailing newline is inserted automatically,
    // and the output is flushed.
    //

    mylog.debug()="No output for level==WARN";

    mylog.warn()="This is a warning message ",234,' ',"composed of several",
                 "output elements";

    mylog.info()="This is a",Logger::endl,"multiline message";

    mylog.set_trace(Logger::DEBUG);        // provides stack trace up to level DEBUG


    // You may want to force use of const references:

    NonCopyable   obj1; // we can't copy this
    MyBulkyObject obj2; // we don't want to copy this

    mylog.debug()=s,Logger::ref(obj1),Logger::ref(obj2);

    \endcode

    Note:
 */
class Logger {
public:
  /// log level
  enum Level {
    OFF=0, //!< no logging
    FATAL, //!< fatal()
    ERROR, //!< error()
    WARN,  //!< warn()
    INFO,  //!< info()
    DEBUG, //!< debug()
    ALL=DEBUG //!< maximum logging
  };

  /// construct new logger (calls set_output())
  Logger(const std::string& _name,Level _level=WARN,
         std::ostream* _out=&std::cerr,int _ansi_colors=ms_color_mode);

  /// construct new logger (calls set_output())
  Logger(const std::string& _name,std::ostream& _out,Level _level=WARN,
         int _ansi_colors=ms_color_mode);

  ~Logger() { unregister_me(); }

  /** @name attibutes and setting
      @{
  */

  /// get name
  const std::string& name() const { return m_name; }

  /// get level
  Level level() const { return m_level; }
  /// set level
  Logger& set_level(Level _level) {
    m_level=_level;
    return *this;
  }

  /// null stream, no output
  static std::ostream& null() { return *ms_null_stream; }

  /** Set output stream.
      \param _os Output stream: no output in case of _os==null()
      \param _colors Use ANSI escape sequences to render colors?
      0=no, 1=yes/forced, any other value implies \c isatty() for
      \c cerr, \c clog, \c cout or no colors.
  */
  Logger& set_output(std::ostream& _os,int _colors=ms_color_mode);

  /// get output stream (use only as argument to other.set_output())
  std::ostream& output() const {
    return m_out!=0 ? *m_out : *ms_null_stream;
  }
  /// Render colors with ANSI escape sequences?
  bool colors() const { return m_colors; }

  /// Render colors with ANSI escape sequences?
  Logger& set_colors(bool _colors) { m_colors=_colors; return *this; }

  /// set default value for _color argument in set_output() (default is -1)
  static void set_default_color_mode(int _mode) {
    ms_color_mode=_mode;
  }

  /// enable stack trace for every output up to _level
  Logger& set_trace(Level _level) { m_trace_level=_level; return *this; }

  /// get logger by name
  static Logger* get(const std::string& _name);

  /** Get log level from _string representation _str.
      Abbreviations are allowed, matching is not case sensitive, e.g.,
      `"o"` matches OFF, `"warn"` matches WARN (`"WARNING"` does not match).
      \param _str string describing level
      \param _name if `!=0` return name as from strlevel() above (use as
      error indicator)
      \return Level or OFF in case of no match (then `*_name==0`)
   */
  static Level strlevel(const std::string& _str,const char** _name);

  /// get string representation of _level
  static const char* strlevel(Level _level);

  /** Set output and level for loggers as described in _properties.

      The _properties string has the syntax
      \code
      key1:[name1 name2] key2:[...] ...
      \endcode
      where `key` describes either a Level or an output or another
      *keyword*, `name1`, etc. are the names of loggers. Names are
      separated by white space. The following special names are
      recognized:

      - `"*"` refers to \a all logger instances
      - `"$"` refers to debug output VC::base::dbg_stream() and is
        valid only if key is an output descriptor (calls
        VC::base::dbg_set_stream()).

      Valid output descriptors are
      \code
      stdout cout 1
      stderr cerr 2
      null           // Logger::null() (or VC::base::dbg_null() for "$")
      \endcode
      Otherwise output is treated as a file name. The following
      *suffixes* are substituted file name:

      - `'+'` open file for appending, remove suffix (a seperator
        line `"=== PID ===\n"` is inserted to the file where PID is the
        process id)
      - '$' is replaced by ".PID" to file name (where PID is the process
        id from `getpid()`)

      Valid keywords are

      - `COLOR_MODE` with values `"off"`,`"on"`,`"auto"`: call
        set_default_color_mode() with argument `0,1,-1`
      - `COLORS` force color mode for loggers
      - `NO_COLORS` disable color mode for loggers

      Example
      \code
      Logger::set_properties(
       "OFF:[*] stderr:[*] "               // default for all
       "WARN:[mylog yourlog] "             // levels for specific
       "cout:[mylog] /tmp/log+:[yourlog] " // output for specific, append
       "stderr:[$] "                       // call dbg_set_stream(cerr)
       "COLOR:[mylog] "                    // force color mode
       "NO_COLOR:[yourlog] "               // disable color mode
       "COLOR_MODE:[auto]"                 // default color mode for
                                           //  loggers *to be* constructed
      );
      \endcode

      ## Notes

      - This function is intended as (default) initialization, e.g,
        in your \c main() function. It acts only on existent logger
        instances (i.e. all of them if they are declared \c static).

      - Files opened by set_properties() are closed on exit of
        the process. Don't reopen these files from anywhere else.

      - Hint: You use a command line argument or an environment
        variable (from getenv() see below) as _properties string.

      \param _properties string
      \return false on failure (syntax error)

      \sa dbg_stream()
   */
  static bool set_properties(const std::string& _properties);

  /// same as above, ignores NULL (e.g., as from getenv())
  static bool set_properties(const char* _properties) {
    if (_properties!=0)
      return set_properties(std::string(_properties));
    return false;
  }

  /// @}

  /** (used internally: store pointer/output object => define ostream<<Reference)
      \ingroup vc_base_debug
      \internal
  */
  template <typename T>
  struct Reference {
    Reference(const T& _data) : data(&_data) {}
    const T* data; //!< pointer to data
  };

  /** Force use of reference.
      \code
      std::string      s; // we don't care if it's copied (reference counting)
      NonCopyable   obj1; // we can't copy this
      MyBulkyObject obj2; // we don't want to copy this

      mylog.debug()=s,Logger::ref(obj1),Logger::ref(obj2);
      \endcode
   */
  template <typename T>
  inline static Reference<T> ref(const T& _data) { return _data; }


  template <typename T> struct EventLeafNode;
  template <typename T,typename NEXT> struct EventNode;

# ifndef DOXYGEN_SKIP

  /** (used internally: returned by Logger::operator()(), etc.)
      \ingroup vc_base_debug
      \internal
  */
  struct Event {
    Event(const Logger* _logger,Level _level) : logger(_logger), level(_level) {}
    ~Event() {
      if (logger!=0)
        logger->processEvent(level);
    }
    std::ostringstream& buffer() const {
      assert(logger!=0);
      return logger->m_buffer;
    }

    template <typename T>
    EventLeafNode<T> operator=(T _t) {
      return EventLeafNode<T>(this,_t);
    }

    const Logger*  logger;     //!< logger
    Level          level;      //!< log level (set by Logger::operator()())
  };


  /** (used internally: closure to Logger::operator=())
      \ingroup vc_base_debug
      \internal
  */
  template <typename T,typename NEXT>
  struct EventNode {
    Event*  event;
    T       data;
    NEXT    next;

    EventNode(Event* _event,T _data,const NEXT& _next)
      : event(_event), data(_data),next(_next) {}

    ~EventNode() {
      if (event!=0 && event->logger!=0)
        event->buffer() << data << next;
    }

    template <typename T1>
    EventNode<T,EventNode<NEXT,EventLeafNode<T1> > >
    operator,(T1 _data) {
      Event* e=event;
      assert(e!=0);

      event=0;
      return EventNode<T,EventNode<NEXT,EventLeafNode<T1> > >
        (e,this->data,EventNode<NEXT,EventLeafNode<T1> >
         (0,next,EventLeafNode<T1>(0,_data)));
      // Note: NEXT is of type EventLeafNode<...> -- we don't care
    }
  };

  /** (used internally: closure to Event::operator=())
    \ingroup vc_base_debug
    \internal
  */
  template <typename T>
  struct EventLeafNode {
    Event* event;
    T      data;

    EventLeafNode(Event* _event,T _data)
      : event(_event), data(_data) {}

    ~EventLeafNode() {
      if (event!=0 && event->logger!=0)
        event->buffer() << data;
    }

    template <typename T1>
    EventNode<T,EventLeafNode<T1> >
    operator,(T1 _data) {
      Event* e=event;
      assert(e!=0);
      event=0;

      return EventNode<T,EventLeafNode<T1> >
        (e,this->data,EventLeafNode<T1>(0,_data));
    }

  };

# endif // DOXYGEN_SKIP


  /** @name methods for logging
      Usage
      \code
      Logger mylog("mylog",Logger::INFO);

      // Messages are output by Event::operator=() and concatenated by
      // the comma operator.

      mylog.debug()="Debug message will not be printed for level INFO.";
      mylog.info()="This is a DEBUG message.", "And here is more text";

      // Note that newline and flush will be appended to message
      // automatically.

      \endcode
      @{
  */

  static const char* endl; //!< end of line

  Event operator()(Level _level) const {
    if (_level<=m_level && m_out!=ms_null_stream) {
      m_access.lock();
      m_buffer.clear();
      m_buffer.str("");
      return Event(this,_level);
    }
    return Event(0,OFF);
  }
  /// log DEBUG event: usage \code log.debug() << message << ...; \endcode
  Event debug() const { return (*this)(DEBUG); }
  /// log INFO event
  Event info() const { return (*this)(INFO); }
  /// log WARN event
  Event warn() const { return (*this)(WARN); }
  /// log ERROR event
  Event error() const { return (*this)(ERROR); }
  /// log FATAL event
  Event fatal() const { return (*this)(FATAL); }

  /// @}

protected:
  /// output buffer (may apply additonal formatting)
  virtual void processEvent(Level _level) const;

private:
  friend struct Event;

  Logger(const Logger&); // better something like "inherit not_copyable"?!

  void register_me();          //!< register logger
  void unregister_me();        //!< unregister

  /// define or invoke synchronization handler
  static void sync_stream(std::ostream* _stream,
                          void (*_handler)(std::ostream*)=nullptr);

  static const Logger::Level   ms_level[]; //!< levels
  static const char*           ms_level_str[]; //!< levels as strings
  static const char*           ms_level_str_esc[]; //!< same w/ escape sequences
  static const char*           ms_level_color_esc[]; //!< escape sequences
  static std::ostream*         ms_null_stream; //!< null stream
  static int                   ms_color_mode; //!< color mode 0=off,1=on,-1=auto

  std::string                  m_name;   //!< name of log
  Level                        m_level;  //!< current log level
  std::ostream*                m_out;    //!< output stream
  Level                        m_trace_level; //!< stack trace depth
  bool                         m_colors; //!< ANSI colors
  mutable std::ostringstream   m_buffer; //!< buffer to output Event

  mutable std::mutex           m_access; //!< synchronize access
};

//-----------------------------------------------------------------------------

template <typename T,typename NEXT>
std::ostream& operator<<(std::ostream& _out,
                         const Logger::EventNode<T,NEXT>& _node) {
  return _out << _node.data << _node.next;
}
template <typename T>
std::ostream& operator<<(std::ostream& _out,
                         const Logger::EventLeafNode<T>& _node) {
  return _out << _node.data;
}
template <typename T>
std::ostream& operator<<(std::ostream& _out,
                         const Logger::Reference<T>& _ref) {
  return _out << *(_ref.data);
}

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_LOG_HH defined
