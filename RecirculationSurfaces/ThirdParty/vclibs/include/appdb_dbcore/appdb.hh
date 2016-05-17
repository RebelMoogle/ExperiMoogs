//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_HH
#define VC_APPDB_HH

#include <climits>
#include <cstdlib>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <functional>

#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/cast.hpp>
#include <boost/signals2/signal.hpp> // boost::signals2 (thread safe)
namespace boost_signals = boost::signals2;

#include "exceptions.hh"
#include "conv.hh"
#include "buffer.hh"

//== INCLUDES =================================================================



//== CLASS DEFINITION =========================================================

namespace VC {

namespace base {
  class Logger;
}

/// An application database \ingroup vc-appdb
namespace appdb {

/** Logger for \ref vc_appdb "AppDB" related output */
extern VC::base::Logger dblog;

//-----------------------------------------------------------------------------

class AbstractDatabase;

class Value;
class Node;
class Database;

# ifndef DOXYGEN_SKIP

struct NodePtrHash {
  size_t operator()(const Node* _node) const;
};

struct NodePtrEqual {
  bool operator()(const Node* _a,const Node* _b) const;
};

struct NodePtrLess {
  bool operator()(const Node* _a,const Node* _b) const;
};

# endif

//-----------------------------------------------------------------------------

/** \defgroup vc_appdb AppDB: an application database

    The application database stores key-value pairs, where the keys
    define an hierarchical name space (similar to the file system).

    \arg Keys are strings of the form `"/[.a-zA-Z0-9_-]+/..." without
    immediate repetition of '/' or '.' and without '-' as first
    character, e.g., `"/a/b/c"`. Each key identifies a Node in the
    Database. The keys define parent-child relations for nodes.

    \arg Each Node may have a Value. Values are data of arbitrary
    type, which might have a single (fixed) dimension/size. Values are
    accessed by a getter/setter interface, which is defined for
    certain data types (currently `int`, `double`, and `std::string`).
    More on \ref vc_appdb_values "values..."

    \arg Nodes are stored in a Database. Every Database implements the
    AbstractDatabase interface.

    The general idea is to store a relatively _small_ amount of values,
    which are in turn of relatively _small_ size. We aim at input
    parameters for applications, including such parameters that are
    modified _interactively_.

    The database provides means to _access_ data (via the Value
    interface) and to _observe_ changes. Notifications are implemented
    using the <a
    href="http://www.boost.org/doc/html/signals2.html">boost/signals2</a>
    library.

    A Node refers to a key in a database. Node instances are _never_
    deleted, whereas Value instances can be attached to or detached
    from nodes (using `std::shared_ptr<Value>`).

    **More information**

    \arg \ref vc_appdb_values "Predefined implementations of" Value,
    their \ref vc_appdb_make "simple usage", and defining
    \ref vc_appdb_sym "Symbolic values" for `enum` data
    \arg VC::appdb makes use of \ref vc_appdb_exc "exceptions"
    \arg Database \ref vc_appdb_io "input/output"
    \arg Dealing with \ref vc_appdb_args "command line arguments"
    \arg \ref vc_appdb_qt "Qt integration" including
    <a href="http://www.freedesktop.org/wiki/Software/dbus">DBUS</a>
    interface
    \arg <a href="http://www.ruby-lang.org/en/">Ruby</a>
    \ref vc_appdb_rb "bindings" including an interactive shell (the
    latter requires Qt)
    \arg \ref vc_appdb_conv "Data conversion" and more \ref vc_appdb_detail
    "details"

    \todo inefficient Value::get_binary(), Value::set_binary():
    supports only Value::get_string(), Value::set_string(), thus
    requires always specializations from `genericvalues.hh` -- could
    setup buffer dynamically according to Value::preferred_accessor()

    \todo Prefer `std::function` (and `std::bind`, `std::regex`) vs
    `boost` counterpart for C++11 (`__cplusplus>201100L`): `namespace
    cxx=std` or `namespace cxx=boost`.

    \example vclibs/appdb/examples/example.cc
 */

/** \defgroup vc_appdb_detail Internally used classes (AppDB)
    \ingroup vc_appdb
 */

//-----------------------------------------------------------------------------


/** \defgroup vc_appdb_values Predefined implementations of Value interface (AppDB)
    \ingroup vc_appdb

    **Special values**

    \arg NilValue no value, access throws InvalidIndexException,
    global constant VC::appdb::NIL
    \arg Event no state, only emits signal

    **Generic values**

    Generic values provide solutions to accessing externally defined
    values via pointers or getter/setter functions. The type is a
    template parameter, there is a variant for `enum` types with a
    sym::SymbolTable (see \ref vc_appdb_sym "here"), and there are
    values that manage their own storage.

    \arg PValue, PSymbolicValue, PEnumValue access via pointer to user data
    \arg PSValue, SPSymbolicValue, PSEnumValue access via `std::shared_ptr` to user data
    \arg SValue, SSymbolicValue, SEnumValue same as PValue but manage own storage
    \arg FValue, FSymbolicValue, FEnumValue access via user provided getter/setter
    functions

    **Utilities**

    Note that values are generally referenced through
    `std::shared_ptr<Value>`. One way to construct a new values is
    <a href="http://en.cppreference.com/w/cpp/memory/shared_ptr/make_shared">`make_shared<>`</a>.

    There are also couple of \ref vc_appdb_make "utility functions" to
    ease the construction of Value instances.
 */

/** \class Value
    \brief Interface to accessing a value in an AbstractDatabase
    \ingroup vc_appdb
    \sa [\ref vc_appdb_values]
 */
class Value {
public:
  /// hint on best access
  enum TypeHint {
    Int='i', Double='d', String='s', Binary='y'
  };

  Value(const std::type_info& _type,unsigned _size,
        TypeHint _preferred_accessor=String);
  virtual ~Value();

  /** @name properties
      @{
  */

  /// get type
  const std::type_info& type() const { return m_type; };
  /// get type hint for preferred access (e.g., `double` for `type()=float`)
  TypeHint preferred_accessor() const { return m_preferred_accessor; }
  /// synonym for size()==0
  bool is_empty() const { return m_size==0; }
  /// get size/dimension
  unsigned size() const { return m_size; }
  /// get node, may be 0
  Node* node() const { return m_node; }
  /// get database, may be zero
  AbstractDatabase* db() const;
  /// logical clock, advanced on notify_changed()
  size_t time() const { return m_time; }

  /// translate TypeHint to Buffere::Type
  static TypeHint bufferType2ValueTypeHint(Buffer::Type _type);

  /// @}

  /** @name basic interface
      @{
  */

  virtual void set_string(const std::string& _data,unsigned _i) = 0;
  virtual void set_int(int _data,unsigned _i) {
    std::string sdata;
    ConvThrow<int,std::string>()(_data,sdata);
    set_string(sdata,_i);
  }
  virtual void set_double(double _data,unsigned _i) {
    std::string sdata;
    ConvThrow<double,std::string>()(_data,sdata);
    set_string(sdata,_i);
  }
  virtual void set_binary(const Buffer& _buf);

  virtual std::string get_string(unsigned _i) const = 0;
  virtual int get_int(unsigned _i) const {
    std::string sdata=get_string(_i);
    int data;
    ConvThrow<std::string,int>()(sdata,data);
    return data;
  }
  virtual double get_double(unsigned _i) const {
    std::string sdata=get_string(_i);
    double data;
    ConvThrow<std::string,double>()(sdata,data);
    return data;
  }
  virtual void get_binary(Buffer& _buf) const;

  // @}

  /** @name signals
      @{
  */

  /// suppress notification (increase counter)
  void silence() { ++m_notify; }
  /// re-enable notification after silence() (decrease counter)
  void unsilence() { --m_notify; }
  /// Are notifications suppressed?
  bool is_silent() const { assert(m_notify>=0); return m_notify>0; }

  /// Automate calls to Value::silence(), Value::unsilence() \ingroup vc_appdb
  class DelayNotification {
  public:
    Value& value;
    bool   suppress;
    size_t time;
    /** Delay or suppress notification for _value during lifetime
        \param _value
        \param _suppress if false, notification is
        delayed until destructor is called; else, notification is suppressed
     */
    DelayNotification(Value& _value,bool _suppress=false)
      : value(_value), suppress(_suppress), time(_value.time()) {
      _value.silence();
    }
    ~DelayNotification() {
      value.unsilence();
      if (!value.is_silent() && !suppress && time!=value.time()) {
        value.notify_changed();
      }
    }
  };

  /// trigger notification
  virtual void notify_changed();

  typedef boost_signals::signal<void(Value*)> signal_changed_t; //!< signal type

  signal_changed_t changed;  //!< signal emitted on change

  /// @}


  /** @name access
      @{
   */

  void set_int(int _val) { set_int(_val,0); }
  void set_double(double _val) { set_double(_val,0); }
  void set_string(std::string _val) { set_string(_val,0); }

  void set(int _val,unsigned _i=0) { this->set_int(_val,_i); }
  void set(double _val,unsigned _i=0) { this->set_double(_val,_i); }
  void set(const std::string& _val,unsigned _i=0) { this->set_string(_val,_i); }
  void set(const char* _val,unsigned _i=0) { this->set_string(_val,_i); }

  template <typename T>
  void set_array(const T* _val) {
    DelayNotification s(*this);
    for (unsigned i=0;i<this->size();++i)
      this->set(_val[i],i);
  }
  template <typename T>
  void set2(const T& _a,const T& _b) {
    assert(this->size()==2);
    T data[2]; data[0]=_a,data[1]=_b;
    this->set_array(data);
  }
  template <typename T>
  void set3(const T& _a,const T& _b,const T& _c) {
    assert(this->size()==3);
    T data[3]; data[0]=_a,data[1]=_b; data[2]=_c;
    this->set_array(data);
  }
  template <typename T>
  void set4(const T& _a,const T& _b,const T& _c,const T& _d) {
    assert(this->size()==3);
    T data[3]; data[0]=_a,data[1]=_b; data[2]=_c; data[3]=_d;
    this->set_array(data);
  }

  int get_int() const { return this->get_int(0); }
  double get_double() const { return this->get_double(0); }
  std::string get_string() const { return this->get_string(0); }

  void _get(int& _val,unsigned _i=0) const { _val=get_int(_i); }
  void _get(double& _val,unsigned _i=0) const { _val=get_double(_i); }
  void _get(std::string& _val,unsigned _i=0) const { _val=get_string(_i); }

  template <typename T>
  void get(T& _val,unsigned _i=0) const {
    typename BestType<T>::type val;
    _get(val,_i);
    ConvThrow<typename BestType<T>::type,T>()(val,_val);
  }
  template <typename T> T get(unsigned _i=0) const {
    T val; get(val,_i); return val;
  }

  template <typename T>
  void get_array(T* _data) const {
    for (unsigned i=0;i<size();++i)
      get(_data[i],i);
  }
  template <typename T>
  void get2(T& _a,T& _b) const { get(_a,0); get(_b,1); }
  template <typename T>
  void get3(T& _a,T& _b,T& _c) const { get(_a,0); get(_b,1); get(_c,2); }
  template <typename T>
  void get4(T& _a,T& _b,T& _c,T& _d) const {
    get(_a,0); get(_b,1); get(_c,2); get(_d,3);
  }

  /// compose string of size() values
  std::string to_string(const std::string& _delimiter=" ");

  /// @}

  /** @name Assignment and comparison

      @{
  */

  /// load value of `_other` into `this`
  virtual Value& assign(const Value& _other)
    throw (InvalidIndexException,BadConversionException);

  /// load value of `_other[_i]` into `this[_j]`
  virtual Value& assign(const Value& _other,unsigned _i,unsigned _j)
    throw (InvalidIndexException,BadConversionException);


  /** Compare values.

      \arg Values are converted automatically for comparison, i.e.,
      values may be equal even though their type()s are not.

      \arg Values are compared lexicographically: a missing dimensions
      refers to a smaller values, i.e., size() mismatch results in
      a valid comparison.

      \param _other compare `this` <=> `_other`
      \return `-1` if `this<_other`, `0` if values are equal, and `+1`
      if `this>_other`
   */
  virtual int compare(const Value& _other) throw (BadConversionException);

  /// @}

  /** @name Observing values.
      @{
  **/

  /** Create "buddy" that obverses `this` value.

      The buddy is a new Value instance that behaves like `this` one
      by *forwarding* all get() and set() requests to the buddy and
      observing `this` instance's `changed` signal.

      - create_buddy() is useful in situations where
        Node::create_reference() cannot be used.

      - The *user is responsible* to ensure that `this` remains valid
        for the *lifetime* of the returned buddy!

      \return buddy
   */
  std::shared_ptr<Value> create_buddy();

  /// @}

protected:
  friend class Node;

  void check_index(unsigned _i) const {
    if (_i>=m_size) //(!(0<=_i && _i<m_size))
      throw InvalidIndexException((Value*) this,_i,VC_DBG_LOCATION());
  }

  Node*                 m_node;
  const std::type_info& m_type;
  size_t                m_time;
  unsigned              m_size;
  int                   m_notify;
  TypeHint              m_preferred_accessor;
};

/// write Value \ingroup vc_appdb
std::ostream& operator<<(std::ostream& _s,const Value& _value);

//-----------------------------------------------------------------------------

/** Undefined Value
    \ingroup vc_appdb_values

    Use NilValue if your prefer assigning a Value rather than leaving
    a Node without value. Accessing a NilValue always raises an
    InvalidIndexException.
 */
class NilValue : public Value {
public:
  NilValue() : Value(typeid(void),0) {}
  virtual ~NilValue() {}
  virtual void set_string(const std::string& /*_data*/,unsigned _i) {
    throw InvalidIndexException(this,_i,VC_DBG_LOCATION());
  }
  virtual std::string get_string(unsigned _i) const {
    throw InvalidIndexException((NilValue*) this,_i,VC_DBG_LOCATION());
    return std::string();
  }
  virtual void notify_changed() {
    throw InvalidIndexException((NilValue*) this,0,VC_DBG_LOCATION());
  }
};

extern const std::shared_ptr<NilValue> NIL; //!< global NilValue \ingroup vc_appdb

//-----------------------------------------------------------------------------

/** Special value to trigger events.
    \ingroup vc_appdb_values

    Event does not store any data. It is uses as a dummy value to
    trigger events on notify_changed().

    The only argument passed is the index of the data.
 */
class Event : public Value {
public:
  Event(unsigned _size=1) : Value(typeid(void),_size), m_what(0) {}
  virtual ~Event() {}

  virtual void set_int(int,unsigned _i) {
    check_index(_i);
    fire(_i);
  }
  virtual void set_double(double,unsigned _i) {
    check_index(_i);
    fire(_i);
  }
  virtual void set_string(const std::string&,unsigned _i) {
    check_index(_i);
    fire(_i);
  }
  virtual std::string get_string(unsigned _i) const {
    check_index(_i);
    return std::string("<event>");
  }

  virtual void notify_changed() { Value::notify_changed(); }
  /// notify a change on signal triggered, `_i` is obtained by what()
  virtual void fire(unsigned _i) { m_what=_i; notify_changed(); }
  /// short for `fire(0)`
  void fire() { fire(0); }

  signal_changed_t& triggered() { return changed; }
  const signal_changed_t& triggered() const { return changed; }

  /// get index that was passed to fire()
  unsigned what() const { return m_what; }

protected:
  unsigned m_what; //!< index provided to set_XXX
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Interface used for definition of a NodeIterator
    \ingroup vc_appdb_detail

    You won't require this class unless you implement your own
    database and iterators. (In this case you may want to have a look
    at AbstractNodeIterator_STL.)
 */
class AbstractNodeIterator {
public:
  virtual AbstractNodeIterator* clone() const = 0;
  virtual ~AbstractNodeIterator() {}
  virtual Node* deref() const = 0;
  virtual void incr() = 0;
  virtual void decr() = 0;
  virtual bool is_equal(const AbstractNodeIterator* _other) const = 0;
};


/** Generic node iterator (similar to STL iterators).
    \ingroup vc_appdb.

    NodeIterator is, e.g., returned by AbstractDatabase::begin(). Its
    implementation hidden in an AbstractNodeIterator.
 */
class NodeIterator {
public:
  typedef Node* value_type;
  NodeIterator(AbstractNodeIterator* _i) : m_i(_i) {}
  NodeIterator(const NodeIterator& _other) : m_i(_other.m_i->clone()) {}
  virtual ~NodeIterator();
  NodeIterator& operator=(const NodeIterator& _other) {
    delete m_i; m_i=_other.m_i->clone(); return *this;
  }
  AbstractNodeIterator* abstract_iterator() const { return m_i; }

  Node* operator*() const { return m_i->deref(); }
  NodeIterator& operator++() { m_i->incr(); return *this; }
  NodeIterator& operator--() { m_i->decr(); return *this; }
  bool operator==(const NodeIterator& _other) const { return m_i->is_equal(_other.m_i); }
  bool operator!=(const NodeIterator& _other) const { return !m_i->is_equal(_other.m_i); }
protected:
  AbstractNodeIterator* m_i;
};

/** Implements an AbstractNodeIterator for an STL container.
    \ingroup vc_appdb_detail

    You won't require this class unless you implement your own
    database and iterators.
*/
template <typename Iterator>
class AbstractNodeIterator_STL : public AbstractNodeIterator {
public:
  Iterator impl;

  AbstractNodeIterator_STL(const Iterator& _i) : impl(_i) {}
  virtual AbstractNodeIterator* clone() const {
    return new AbstractNodeIterator_STL(impl);
  }
  virtual Node* deref() const { return *impl; }
  virtual void incr() { ++impl; }
  virtual void decr() { --impl; }
  virtual bool is_equal(const AbstractNodeIterator* _other) const {
    return impl==dcast((AbstractNodeIterator*) _other)->impl;
  }

  static AbstractNodeIterator_STL<Iterator>*
  dcast(AbstractNodeIterator* _other) {
    AbstractNodeIterator_STL<Iterator>* other=
      dynamic_cast<AbstractNodeIterator_STL*>(_other);
    assert(other!=0);
    return other;
  }
};

//-----------------------------------------------------------------------------

/// pair of NodeIterator objects defining a range \ingroup vc_appdb
typedef std::pair<NodeIterator,NodeIterator> NodeRange;

/** \class AbstractDatabase
    \brief Interface to an application database.
    \ingroup vc_appdb
 */
class AbstractDatabase {
public:

  AbstractDatabase(const std::string& _name);
  virtual ~AbstractDatabase();

  const std::string& name() const { return m_name; }
  virtual size_t size() const = 0;

  /// add a single node
  virtual Node* add_node(const std::string& _key)
    throw(InvalidKeyException,NodeExistsException) = 0;

  /** Calls add_node() for any partial path of _path (if node does not
      exist).
      \sa import_args()
   */
  void add_nodes(const std::string& _path) throw(InvalidKeyException);

  /// find Node, may return 0
  virtual Node* find_node(const std::string& _key) const {
    NodeIterator ii=find_node2(_key);
    return *ii;
  }
  /// find Node, may return end()
  virtual NodeIterator find_node2(const std::string& _key) const = 0;
  /// find Node, raise NoSuchNodeException on failure
  Node& operator[](const std::string& _key) const throw(NoSuchNodeException);
  /// find Node or add_node() if it does not exist yet
  Node* find_or_add_node(const std::string& _key) throw(InvalidKeyException);

  /** @name ranges
      @{
  */
  /// NodeRange(begin(),end()) represents all nodes
  virtual NodeIterator begin() const = 0;
  /// NodeRange(begin(),end()) represents all nodes
  virtual NodeIterator end() const = 0;
  /// advance _iterator and skip subtrees, different level denotes end
  virtual void advance_on_same_level(NodeIterator& _iterator) const = 0;
  /// get subtree of _node, may be an emty range
  virtual NodeRange subtree(const Node* _node) const = 0;

  ///@}

  /** @name hierarchy
      @{
  */

  /// get parent of _node, may be 0
  virtual Node* parent(const Node* _node) const = 0;
  /// get _children of _node
  virtual std::vector<Node*> children(const Node* _node) const = 0;
  /// Determine of _node has any children.
  virtual bool has_children(const Node* _node) const;

  /// @}


  /** @name Events
      @{
  */

  /// event identifiers for signal node_changed and handle_node_event()
  enum Event {
    CreateNode,       //!< data=0
    CreateReference,  //!< data=pointer to previously referenced Node (or 0)
    RemoveReference,  //!< data=pointer to referenced Node (or 0)
    SetValue,         //!< data=pointer to previously attached Value (or 0)
    RemoveValue,      //!< data=pointer to attached Value (or 0)
    ValueChanged,     //!< data=0
    DestroyDatabase   //!< node=0,data=0, called by destructor
  };
  /// signaltype
  typedef boost_signals::signal<void(Event,Node*,void*)> signal_node_changed_t;

  signal_node_changed_t node_changed; //!< signal emitted on changed

  /// trigger notification
  virtual void handle_node_event(Event _event,Node* _node,void* _data);

  /// @}

  /// debugging: dump all node keys to VC::base::Logger "appdb"
  virtual void log_list_nodes() const;
  /// debugging: dump all nodes keys and node values to VC::base::Logger "appdb"
  virtual void log_list_nodes_values() const;

  /// Is _key a valid key? (s.a. InvalidKeyException)
  static bool is_valid_key(const std::string& _key);
  /// split key /a/b/c" => "/a", _remainder="/b/c"
  static std::string key_first(const std::string& _key,std::string* _remainder=0);
  /// split key /a/b/c" => "/a/b", _tail="c"
  static std::string key_head(const std::string& _key,std::string* _tail=0);
  /// split key /a/b/c" => "c", _remainder="/a/b"
  static std::string key_tail(const std::string& _key,std::string* _head=0);

protected:

  virtual void new_node_event(Node* _node) { handle_node_event(CreateNode,_node,0); }

  std::string    m_name;       //!< db name
};

//-----------------------------------------------------------------------------

/** \class Database
    \brief Implements an application database.
    \ingroup vc_appdb
 */
class Database : public AbstractDatabase {
protected:
  typedef std::unordered_set<Node*,NodePtrHash,NodePtrEqual> hash_t;
  typedef std::set<Node*,NodePtrLess> set_t;
  typedef std::vector<Node*> array_t;

  typedef set_t::const_iterator range_iterator_t;
  typedef std::pair<range_iterator_t,range_iterator_t> range_t;

  typedef AbstractNodeIterator_STL<range_iterator_t> abstract_iterator_t;

public:
  Database(const std::string& _name);
  virtual ~Database();

  virtual size_t size() const { return m_hash.size(); }

  virtual Node* add_node(const std::string& _key)
    throw(InvalidKeyException,NodeExistsException);

  virtual Node* find_node(const std::string& _key) const;
  virtual NodeIterator find_node2(const std::string& _key) const
    throw(NoSuchNodeException);

  virtual NodeIterator begin() const {
    return NodeIterator(new abstract_iterator_t(m_set.begin()));
  }
  virtual NodeIterator end() const {
    return NodeIterator(new abstract_iterator_t(m_set.end()));
  }
  virtual NodeRange subtree(const Node* _node) const {
    range_t r=_subtree(_node);
    return NodeRange(NodeIterator(new abstract_iterator_t(r.first)),
                     NodeIterator(new abstract_iterator_t(r.second)));
  }

  virtual void advance_on_same_level(NodeIterator& _iterator) const {
    abstract_iterator_t* ii=abstract_iterator_t::dcast(_iterator.abstract_iterator());
    advance_on_same_level(ii->impl);
  }

  virtual Node* parent(const Node* _node) const;
  virtual std::vector<Node*> children(const Node* _noden) const;
  virtual bool has_children(const Node* _node) const;

  //virtual void log_list_nodes() const;

protected:
  array_t        m_nodes;      //!< stores all nodes
  hash_t         m_hash;       //!< key => node
  set_t          m_set;        //!< ordered map key => node

  range_t _subtree(const Node* _node) const;
  void advance_on_same_level(range_iterator_t& ii) const;
  range_iterator_t find_next_on_same_level(const Node* _node) const;
};

//-----------------------------------------------------------------------------

/** \class Node
    \brief Node in an AbstractDatabase.
    \ingroup vc_appdb
 */
class Node {
public:
  virtual ~Node(); // prefer to hide this

  /// get Database
  AbstractDatabase* db() const { return m_db; }
  /// get key
  const std::string& key() const { return m_key; }
  /// get hierarchy level
  unsigned level() const { return m_level; }

  /// get parent Node, may be 0 (short for AbstractDatabase::parent(this))
  Node* parent() const;
  /// get _children (short for AbstractDatabase::children(this,_children))
  std::vector<Node*> children() const {
    return m_db->children(this);
  }
  /// Does `this` Node have children? (short for AbstractDatabase::has_children(this))
  bool has_children() const { return m_db->has_children(this); }

  /// find child key()+_key_suffix, may return 0
  Node* find_child(const std::string& _key_suffix) const;
  /// same as find_child(), but throws NoSuchNodeException on failure
  Node& operator[](const std::string& _key_suffix) const
    throw(NoSuchNodeException);
  /// get iterator to `this` Node (short for Database::find_node2(key()))
  NodeIterator iterator() const { return m_db->find_node2(m_key); }

  /// add child by calling AbstractDatabase::add_child() for db()
  Node* add_child(const std::string& _key_suffix)
    throw(InvalidKeyException,NodeExistsException);
  /// return add_child() if find_child() returns `0`
  Node* find_or_add_child(const std::string& _key_suffix)
    throw(InvalidKeyException);

  /** @name References

      - Nodes may act as references to other nodes. Access of the
        node's value() is redirected to the referenced node.

      - References only apply to values, they do not change child
        relations!

      - There is no check for dependency cycles!

      - Why references? -- You can use references to define "virtual"
        groups of nodes/values. For instance, you define a tree of
        nodes that mirrors program logic/parameters of algoritms. Then
        you define another tree that refers to a dialog, which is
        presented to the user. Nodes in this dialog-tree can be
        realized as references.

      - *Limitation:* A Value attached to a node `B` that references
        `A` will always report `A` as Value::node(). There is no way
        to retrieve `B`. If this is required, you may want to use
        Value::create_buddy().

      @{
  */

  /// create a reference to _node
  void create_reference(Node* _node);
  /// remove reference
  void remove_reference();
  /// Is `this` node referencing another node?
  bool is_reference() const { return m_ref!=0; }
  /// get references node, may be 0
  Node* reference() const { return m_ref; }

  /// @}

  /** @name Values

      Access of values respects redirection from references.

      @{
  */

  /// attach _value, return `this`
  Node& set_value(const std::shared_ptr<Value>& _value);
  /// remove _value
  std::shared_ptr<Value> remove_value();
  /// Is a value attached to ]c this node?
  bool has_value() const { return value()!=0; }
  /// get value, may be 0 if !has_value()
  std::shared_ptr<Value> value() const {
    return m_ref!=0 ? m_ref->value() : m_value;
  }
  /// get value, return NIL if there is no value attached
  Value& value2() const {
    Value* v=value().get();
    return (v!=0 ? *v : *std::static_pointer_cast<Value>(NIL).get()); }
  /// get value, throw if !has_value()
  Value& value3() const throw(NoValueException);

  /// @}

  /** @name Events
      @{
  */

  /// switch notification of parent nodes (signal child_value_changed) on or off
  void set_notify_parent(bool _state) { m_notify_parent=_state; }
  /// Is notification of parent node enabled?
  bool notify_parent() const { return m_notify_parent; }

  typedef boost_signals::signal<void(Node*)> signal_value_changed_t; //!< signal type
  typedef boost_signals::signal<void(Node*,Node*,unsigned _relative_level)>
  signal_child_value_changed_t; //!< signal type
  typedef boost_signals::signal<void(Node*,Value*)> signal_value_replaced_t; //!< signal type

  signal_value_changed_t       value_changed; //!< signal emitted on change
  signal_child_value_changed_t child_value_changed; //!< signal emitted on change of child
  signal_value_replaced_t      value_replaced;

  /// trigger changed signal
  virtual void handle_value_changed(Node*);

  /// @}


protected:
  friend class Database;

  Node(Database* _db,const std::string& _key) throw(InvalidKeyException);
  Node(const std::string& _key) : m_value(), m_key(_key) { /* comparable dummy */ }

  AbstractDatabase*         m_db;
  std::shared_ptr<Value>    m_value;
  Node*                     m_ref;
  std::string               m_key;
  boost_signals::connection m_connection; // required for disconnect???
  unsigned                  m_level;
  bool                      m_notify_parent;

private:
  Node(const Node&);
};

//-----------------------------------------------------------------------------
#ifndef DOXYGEN_SKIP
//-----------------------------------------------------------------------------

inline AbstractDatabase* Value::db() const {
  assert(m_node!=0);
  return m_node->db();
}

//-----------------------------------------------------------------------------

inline Node*
AbstractDatabase::find_or_add_node(const std::string& _key)
  throw(InvalidKeyException) {
  Node* node=find_node(_key);
  return node!=0 ? node : add_node(_key);
}

//-----------------------------------------------------------------------------

inline Node& AbstractDatabase::operator[](const std::string& _key) const
  throw(NoSuchNodeException) {
  Node* node=find_node(_key);
  if (node==0)
    throw NoSuchNodeException((Database*) this,_key,VC_DBG_LOCATION());

  return *node;
}

inline std::string
AbstractDatabase::key_first(const std::string& _key,std::string* _remainder) {
  assert(_key.size()>1);
  size_t i=_key.find('/',1);
  assert(i!=std::string::npos);
  if (_remainder!=0)
    *_remainder=_key.substr(i+1);
  return _key.substr(0,i);
}

inline std::string
AbstractDatabase::key_head(const std::string& _key,std::string* _tail) {
  assert(_key.size()>1);
  size_t i=_key.rfind('/',std::string::npos);
  assert(i!=std::string::npos);
  if (_tail!=0)
    *_tail=_key.substr(i+1);
  return _key.substr(0,i);
}

inline std::string
AbstractDatabase::key_tail(const std::string& _key,std::string* _remainder) {
  assert(_key.size()>1);
  size_t i=_key.rfind('/',std::string::npos);
  assert(i!=std::string::npos);
  if (_remainder!=0)
    *_remainder=_key.substr(0,i);
  return _key.substr(i+1);
}

//-----------------------------------------------------------------------------

inline Node* Database::parent(const Node* _node) const {
  return _node->parent();
}

inline void
Database::advance_on_same_level(range_iterator_t& ii) const {
  if (ii==m_set.end())
    return;

  range_iterator_t jj=ii;

  if (++jj==m_set.end()) {
    ii=m_set.end();
    return;
  }

  if ((*jj)->level()>(*ii)->level())       // did we descend?
    jj=find_next_on_same_level(*ii);       // => then ignore

  ii=jj;
}

inline Database::range_iterator_t
Database::find_next_on_same_level(const Node* _node) const {
  Node dummy(_node->key()+"/\xff"); // invalid key; "maximum"
  range_iterator_t ii=m_set.upper_bound(&dummy);
  return ii;
}

inline bool Database::has_children(const Node* _node) const {
  assert(_node->db()==this);
  set_t::const_iterator ii=m_set.find((Node*) _node);
  assert(ii!=m_set.end());
  if (++ii==m_set.end())
    return false;

  return (*ii)->key().find(_node->key())==0;
}

//-----------------------------------------------------------------------------

inline Node* Node::find_child(const std::string& _key_suffix) const {
  return m_db->find_node(m_key+'/'+_key_suffix);
}

inline Node&
Node::operator[](const std::string& _key_suffix) const
  throw(NoSuchNodeException) {
  std::string key=m_key+'/'+_key_suffix;
  Node* node=m_db->find_node(key);
  if (node==0)
    throw NoSuchNodeException(m_db,key,VC_DBG_LOCATION());
  return *node;
}

inline Node* Node::add_child(const std::string& _key_suffix)
  throw(InvalidKeyException,NodeExistsException) {
  return m_db->add_node(key()+'/'+_key_suffix);
}

inline Node* Node::find_or_add_child(const std::string& _key_suffix)
  throw(InvalidKeyException) {
  Node* ch=find_child(_key_suffix);
  if (ch==0)
    ch=add_child(_key_suffix);
  return ch;
}

//-----------------------------------------------------------------------------

inline size_t NodePtrHash::operator()(const Node* _node) const {
  std::hash<std::string> h;
  return h(_node->key());
}

inline bool NodePtrEqual::operator()(const Node* _a,const Node* _b) const {
  return _a->key()==_b->key();
}

inline bool NodePtrLess::operator()(const Node* _a,const Node* _b) const {
  return _a->key()<_b->key();
}

//-----------------------------------------------------------------------------

inline Value::TypeHint Value::bufferType2ValueTypeHint(Buffer::Type _type) {
  switch (_type) {
  case Buffer::BYTE:
  case Buffer::BOOLEAN:
  case Buffer::INT16: case Buffer::UINT16:
  case Buffer::INT32: case Buffer::UINT32:
  case Buffer::INT64: case Buffer::UINT64:
    return Value::Int;
  case Buffer::DOUBLE:
  case Buffer::FLOAT:
    return Value::Double;
  case Buffer::STRING:
    return Value::String;
  default: ;
  }
  return Value::Binary;
}

//-----------------------------------------------------------------------------
#endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------

//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
# include "genericvalues.hh"
# include "sym.hh"
# include "make.hh"
# include "io.hh"
# include "args.hh"
//=============================================================================
#endif // VC_APPDB_HH defined
