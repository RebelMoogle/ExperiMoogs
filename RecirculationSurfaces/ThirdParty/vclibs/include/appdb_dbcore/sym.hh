//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_SYM_HH
#define VC_APPDB_SYM_HH

//== INCLUDES =================================================================

# include <unordered_map>
# include <map>
# include <string>
# include <memory>
# include <cctype>

// included by appdb.hh

//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

/** \defgroup vc_appdb_sym Symbolic values (AppDB).
    \ingroup vc_appdb

    Symbolic values use a SymbolTable instance to provide a mapping
    from a symbol (`std::string`) to a value. This is useful for
    providing an external representation of values, e.g., for saving
    in a human readable file format.

    Examples are `bool` values with external representations `"true",
    "false", ...` or `enum` values with user defined mappings.

    Preferred symbol tables can be defined for a specific `enum` type
    by specializing VC::appdb::sym::symbols, i.e.,
    `VC::appdb::sym::symbols<bool>::symbol_table()` refers to a
    predefined symbol table.

    The base class sym::SymbolTableInfo is not type-dependent and provides
    meta-information on sym::SymbolTable, e.g., for providing help on
    symbols.

    Any symbolic Value should inherit and implement the
    sym::SymbolicValue interface, which gives access to
    sym::SymbolTableInfo. Use sym::info() to decide if _any_ Value is
    a symbolic value and to access meta-information in case.

    \sa PSymbolicValue, SPSymbolicValue, SSymbolicValue,
    FSymbolicValue, PEnumValue, SPEnumValue, SEnumValue, FEnumValue
 */

/// definitions for SymbolTable \ingroup vc_appdb_sym
namespace sym {
/// options to ctor \ingroup vc_appdb_sym
enum Options {
  None=0,
  IgnoreCase=1,   //!< compare always lower case
  NoDefault=2,    //!< don't use default value
  VerifyStrict=4  //!< throw exception if there is no mapping
};

/// base class for SymbolTable \ingroup vc_appdb_sym
class SymbolTableInfo {
public:
  /// mapping `symbol => description`
  typedef std::map<std::string,std::string> descriptions_type;

  SymbolTableInfo(const std::string& _name);
  ~SymbolTableInfo();

  /// get name
  const std::string& name() const { return m_name; }
  /// get mapping (e.g., to list descriptions)
  const descriptions_type& describe() const { return m_dmap; }

protected:
  std::string       m_name; //!< name of mapping
  descriptions_type m_dmap; //!< map descriptions
};

/** Interface for symbol substitution table.
    \ingroup vc_appdb_sym

    \arg SymbolTable maps `std::string` to `T` and can be used, e.g.,
    to implement `enum` types with string constants.

    \arg If the mapping is not injective, the reverse map considers the
    _first_ mapping that was add()ed.

    \tparam T value type
    \tparam TM type that is used internally in map (e.g., pefer `int` for any
    `enum`)

    \sa Enum

    \todo fix `TM` to `int` and have add(), lookup(0 and rlookup() in base class
    (no inlining; template just casts type)
*/
template <typename T,typename TM=int>
class SymbolTable : public SymbolTableInfo {
public:
  typedef T value_type; //!< value_type

  /** ctor
      \param _default default Value if no mapping is defined (may be 0)
      \param _name name of mapping
      \param _options OR-combination of Options
   */
  SymbolTable(const T& _default,const std::string& _name,
              unsigned _options=sym::None)
    : SymbolTableInfo(_name), m_default(_default), m_options(_options) {
    m_map.reserve(32);
    m_rmap.reserve(32);
  }
  /** ctor
      \param _name name of mapping
      \param _options OR-combination of Options, sym::NoDefault will be set
   */
  SymbolTable(const std::string& _name,unsigned _options=sym::None)
    : SymbolTableInfo(_name), m_options(_options|sym::NoDefault) {
  }

  ~SymbolTable() {}

  /// get default Value
  const T& defaultValue() const { return m_default; }

  /// add mapping
  SymbolTable<T,TM>* add(const std::string& _symbol,const T& _value,
                         const std::string& _description=std::string()) {
    m_map[make_symbol(_symbol)]=(TM) _value;
    if (m_rmap.find((TM) _value)==m_rmap.end()) {
      m_rmap[(TM) _value]=_symbol;    // store original case
      m_dmap[_symbol]=_description;   // store only one description ...
    }
    else if (!_description.empty())
      m_dmap[_symbol]=_description;   // ... unless requested

    return this;
  }

  /** Lookup `_symbol`.
      Throws if `_symbol` is undefined _and_ there is no default value
      _or_ option `sym::VerifyStrict` applies.
      \param _symbol
      \return value
   */
  const T lookup(const std::string& _symbol) const
    throw(UndefinedSymbolException) {
    auto vi=m_map.find(make_symbol(_symbol));
    if (vi==m_map.end()) {
      if (!(m_options&sym::VerifyStrict) && !(m_options&sym::NoDefault))
        return m_default;
      throw UndefinedSymbolException
        (m_name,_symbol,std::string(),VC_DBG_LOCATION());
    }
    return detail::appdb_cast<T,TM>::value(vi->second);
  }

  /** Reverse lookup.
      Throw if no symbols is found and `sym::VerifyStrict` applies.
      \param _value find symbol for `_value`
      \return symbol or empty string (if no throw)
   */
  const std::string& rlookup(const T& _value) const
    throw(UndefinedSymbolException) {
    static const std::string NIL;
    auto si=m_rmap.find((TM) _value);
    if (si==m_rmap.end()) {
      if (!(m_options&sym::VerifyStrict))
        return NIL;
      std::string vstr("<?>");
      Conv<TM,std::string> conv; // fallback to numeric conversion
      conv((TM) _value,vstr);
      throw UndefinedSymbolException
        (m_name,std::string(),vstr,VC_DBG_LOCATION());
    }
    return si->second;
  }

  /** @name helpers used by subclasses of Value
      @{
   */

  /// set value: assume `VALUE` defines `void set_native(const T&,unsigned)`
  template <typename VALUE>
  void set_value(VALUE& _value,unsigned _i,const std::string& _symbol)
    throw(UndefinedSymbolException,BadConversionException) {
    auto self=this;
    if (self==0)
      throw BadConversionException(" (missing SymbolTable)",
                                   typeid(std::string),typeid(T),VC_DBG_LOCATION());
    _value.set_native(lookup(_symbol),_i);
  }

  /// get value: assume `VALUE` defines `T get_native(unsigned)`
  template <typename VALUE>
  std::string get_string_from_value(const VALUE& _value,unsigned _i) const
    throw(UndefinedSymbolException,BadConversionException){
    auto self=this;
    if (self==0)
      throw BadConversionException(" (missing SymbolTable)",
                                   typeid(std::string),typeid(T),VC_DBG_LOCATION());
    return rlookup(_value.get_native(_i));
  }

  /// @}

private:
  /// mapping `symbol => value
  typedef std::unordered_map<std::string,TM> map_t;
  /// reverse mapping
  typedef std::unordered_map<TM,std::string> rmap_t;

  /// lower case if IgnoreCase
  std::string make_symbol(const std::string& _symbol) const {
    if (m_options&sym::IgnoreCase) {
      std::string s=_symbol;
      for (auto ii=s.begin();ii!=s.end();++ii)
        *ii=std::tolower(*ii);
      return s;
    }
    else
      return _symbol;
  }

  map_t       m_map;        //!< mapping `symbol => (value,description)
  rmap_t      m_rmap;       //!< reverse mapping
  T           m_default;    //!< default value
  unsigned    m_options;    //!< options
};

//-----------------------------------------------------------------------------

/// Symbol table for `enum` types
template <typename T>
class Enum : public SymbolTable<T,int> {
public:
  typedef SymbolTable<T,int> base_type;
  typedef T value_type; //!< value_type

  /** ctor
      \param _default default Value if no mapping is defined (may be 0)
      \param _name name of mapping
      \param _options OR-combination of Options
   */
  Enum(const T& _default,const std::string& _name,unsigned _options=sym::None)
    : base_type(_default,_name,_options) {
  }
  /** ctor
      \param _name name of mapping
      \param _options OR-combination of Options, sym::NoDefault will be set
   */
  Enum(const std::string& _name,unsigned _options=sym::None)
    : base_type(_name,_options) {
  }

  Enum<T>* add(const std::string& _symbol,const T& _value,
               const std::string& _description=std::string()) {
    return (Enum<T>*) base_type::add(_symbol,_value,_description);
  }
};

//-----------------------------------------------------------------------------

/// Try to get a predefined symbol table for `T`. \ingroup vc_appdb
template <typename T>
struct symbols {
  /// The default implementation returns `0` (no table).
  static std::shared_ptr<SymbolTable<T> > symbol_table() {
    return std::shared_ptr<SymbolTable<T> >();
  }
};

template <>
struct symbols<bool> {
  /** Symbol table for `bool`.
    \ingroup vc_appdb
    - `"true","t","yes","on" map to `true`
    - `"false","f","no","off","nil" map to `false`
    The mapping is _not_ case sensitive, and there is _no_ default mapping.
*/
  static std::shared_ptr<Enum<bool> > symbol_table();
};

//-----------------------------------------------------------------------------

/// base class to "tag" symbolic values \ingroup vc_appdb_sym
class SymbolicValue {
public:
  SymbolicValue();
  virtual ~SymbolicValue();
  /// get access to SymbolTableInfo::name() and SymbolTableInfo::describe()
  virtual const SymbolTableInfo* symbols() const = 0;
};

/// get SymbolicValue or `0` for any other Value \ingroup vc_appdb_sym
inline const SymbolicValue* info(const Value* _value) {
  return dynamic_cast<const SymbolicValue*>(_value);
}

/// get SymbolicValue or `0` for any other Value \ingroup vc_appdb_sym
inline const SymbolicValue* info(const std::shared_ptr<Value>& _value) {
  return info(_value.get());
}

//-----------------------------------------------------------------------------

/** Register human readable description `_alias` for `_type`.
    \ingroup vc_appdb_sym
    \param _type some type
    \param _alias a human readable and/or pretty printed alias
    \return `0` if there exists already an alias, `1` else
    \sa describe_type()
 */
int register_type_alias(const std::type_info& _type,const std::string& _alias);

/** Register human readable description `_alias` for `T`.
    \ingroup vc_appdb_sym
    \tparam T some type
    \param _alias a human readable and/or pretty printed alias
    \return `0` if there exists already an alias, `1` else
    \sa describe_type()
 */
template <typename T>
int register_type_alias(const std::string& _alias) {
  T* dummy=nullptr;
  return register_type_alias(typeid(*dummy),_alias);
}

/** Get human readable description of `_type`.
    \ingroup vc_appdb_sym
    Get alias or demangled type. Eventually, namespaces and template arguments
    are stripped.
    \param _type some type
    \sa register_type_alias(), VC::base::demangle()
*/
std::string describe_type(const std::type_info& _type);

/** Get human readable description of `_type`.
    \ingroup vc_appdb_sym
    Get alias or demangled type. Eventually, namespaces and template arguments
    are stripped.
    \tparam T some type
    \param _dummy is of type `T`, will not be accessed
    \sa register_type_alias(), VC::base::demangle()
*/
template <typename T>
std::string describe_type(const T& _dummy) {
  return describe_type(typeid(_dummy));
}

/** Get human readable description of `_type`.
    \ingroup vc_appdb_sym
    Get alias or demangled type. Eventually, namespaces and template arguments
    are stripped.
    \tparam T some type
    \sa register_type_alias(), VC::base::demangle()
*/
template <typename T>
std::string describe_type() {
  T* dummy=nullptr;
  return describe_type(typeid(*dummy));
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
} // namespace sym
} // namespace appdb
} // namespace VC

# include "genericvalues.hh"

namespace VC {
namespace appdb {
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#ifndef DOXYGEN_SKIP
//-----------------------------------------------------------------------------
namespace detail {
//-----------------------------------------------------------------------------

template <typename T,typename ACCESSORMIXIN,typename TM>
class TSymbolicValueBase : public sym::SymbolicValue, public Value,
                           public ACCESSORMIXIN {
public:
  typedef T value_type;
  typedef sym::SymbolTable<T,TM> map_type;

  std::shared_ptr<map_type> map;

  TSymbolicValueBase(unsigned _size,
                     const std::shared_ptr<map_type>& _map,
                     TypeHint _type=detail::preferred<T>::accessor()) :
    sym::SymbolicValue(), Value(typeid(T),_size,_type),
    map(_map) {}
  virtual ~TSymbolicValueBase() {}

  virtual void set_string(const std::string& _data,unsigned _i) {
    map->set_value(*this,_i,_data);
    notify_changed();
  }
  virtual std::string get_string(unsigned _i) const {
    return map->get_string_from_value(*this,_i);
  }

  void set_native(const T& _data,unsigned _i) { this->set_data(_data,_i); }
  T get_native(unsigned _i) const { return this->get_data(_i); }

  /// get information on symbol table
  virtual const sym::SymbolTableInfo* symbols() const { return map.get(); }

protected:
  void verify(const T& _data) {
    if (map.get()==0)
      throw BadConversionException(" (missing SymbolTable)",
                                   typeid(std::string),typeid(T),VC_DBG_LOCATION());
    map->rlookup(_data); // verify: throw on mismatch
  }

  template <typename TYPE>
  void _set_binary(const Buffer& _buf,const TYPE& _dummy) {
    if (_buf.type()!=Buffer::type_of(_dummy) ||
        _buf.type()==Buffer::INVALID)
      throw BadConversionException(" (set_binary)",typeid(Buffer),typeid(TYPE),
                                   VC_DBG_LOCATION());

    size_t index=Buffer::HeaderSize;
    for (unsigned i=0;i<this->size();++i) {
      TYPE data;
      index=_buf.unpack(index,data);
      verify(appdb_cast<T,TYPE>::value(data));
      this->set_data(appdb_cast<T,TYPE>::value(data),i);
    }
    assert(index==_buf.size());
    this->notify_changed();
  }
  template <typename TYPE>
  void _get_binary(Buffer& _buf,const TYPE& _dummy) const {
    _buf.clear(Buffer::type_of(_dummy));
    for (unsigned i=0;i<this->size();++i)
      _buf.pack(TYPE(this->get_data(i)));
  }
};

//-----------------------------------------------------------------------------

template <typename T,typename ACCESSORMIXIN,typename TM>
class TSymbolicValue : public TSymbolicValueBase<T,ACCESSORMIXIN,TM> {
public:
  typedef T value_type;
  typedef TSymbolicValueBase<T,ACCESSORMIXIN,TM> base_type;
  typedef typename base_type::map_type map_type;

  TSymbolicValue(unsigned _size,const std::shared_ptr<map_type>& _map,
                 Value::TypeHint _type=detail::preferred<T>::accessor()) :
    base_type(_size,_map,_type) {}
  virtual ~TSymbolicValue() {}
};

template <typename ACCESSORMIXIN>
class TSymbolicValue<int,ACCESSORMIXIN,int>
  : public TSymbolicValueBase<int,ACCESSORMIXIN,int> {
public:
  typedef int value_type;
  typedef TSymbolicValueBase<int,ACCESSORMIXIN,int> base_type;
  typedef typename base_type::map_type map_type;

  TSymbolicValue(unsigned _size,const std::shared_ptr<map_type>& _map) :
    base_type(_size,_map) {}
  virtual ~TSymbolicValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    this->set_data(_data,_i);
    this->notify_changed();
  }
  virtual void set_double(double /*_data*/,unsigned _i) {
    this->check_index(_i);
    throw BadConversionException(" (set_double())",
                                 typeid(double),typeid(int),VC_DBG_LOCATION());
  }
  virtual void set_binary(const Buffer& _buf) {
    this->_set_binary(_buf,int(0));
  }
  virtual int get_int(unsigned _i) const {
    this->check_index(_i);
    return this->get_data(_i);
  }
  virtual double get_double(unsigned _i) const {
    return double(get_int(_i));
  }
  virtual void get_binary(Buffer& _buf) const {
    this->_get_binary(_buf,int(0));
  }
};

template <typename ACCESSORMIXIN>
class TSymbolicValue<double,ACCESSORMIXIN,double>
  : public TSymbolicValueBase<double,ACCESSORMIXIN,double> {
public:
  typedef double value_type;
  typedef TSymbolicValueBase<double,ACCESSORMIXIN,double> base_type;
  typedef typename base_type::map_type map_type;

  TSymbolicValue(unsigned _size,const std::shared_ptr<map_type>& _map) :
    base_type(_size,_map) {}
  virtual ~TSymbolicValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    this->set_data(double(_data),_i);
    this->notify_changed();
  }
  virtual void set_double(double _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    this->set_data(_data,_i);
    this->notify_changed();
  }
  virtual void set_binary(const Buffer& _buf) {
    this->_set_binary(_buf,double(0));
  }
  virtual int get_int(unsigned _i) const {
    this->check_index(_i);
    throw BadConversionException(" (get_int())",
                                 typeid(double),typeid(int),VC_DBG_LOCATION());
  }
  virtual double get_double(unsigned _i) const {
    this->check_index(_i);
    return this->get_data(_i);
  }
  virtual void get_binary(Buffer& _buf) const {
    this->_get_binary(_buf,double(0));
  }
};

template <typename ACCESSORMIXIN>
class TSymbolicValue<float,ACCESSORMIXIN,float>
  : public TSymbolicValueBase<float,ACCESSORMIXIN,float> {
public:
  typedef float value_type;
  typedef TSymbolicValueBase<float,ACCESSORMIXIN,float> base_type;
  typedef typename base_type::map_type map_type;

  TSymbolicValue(unsigned _size,const std::shared_ptr<map_type>& _map) :
    base_type(_size,_map)
  {}
  virtual ~TSymbolicValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    this->set_data(float(_data),_i);
    this->notify_changed();
  }
  virtual void set_double(double _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    float data;
    ConvThrow<double,float>()(_data,data);
    this->set_data(data);
    this->notify_changed();
  }
  virtual void set_binary(const Buffer& _buf) {
    this->_set_binary(_buf,double(0));
  }
  virtual int get_int(unsigned _i) const {
    this->check_index(_i);
    throw BadConversionException(" (get_int())",
                                 typeid(double),typeid(int),VC_DBG_LOCATION());
  }
  virtual double get_double(unsigned _i) const {
    this->check_index(_i);
    return double(this->get_data(_i));
  }
  virtual void get_binary(Buffer& _buf) const {
    this->_get_binary(_buf,double(0));
  }
};

//-----------------------------------------------------------------------------

template <typename SYMBOLICVALUEMIXIN>
class EnumValueBase : public SYMBOLICVALUEMIXIN {
public:
  typedef typename SYMBOLICVALUEMIXIN::value_type value_type;
  typedef SYMBOLICVALUEMIXIN base_type;
  typedef typename base_type::map_type map_type;

  EnumValueBase(unsigned _size,const std::shared_ptr<map_type>& _map)
    : base_type(_size,_map) {}

  EnumValueBase(value_type* _p,unsigned _size,const std::shared_ptr<map_type>& _map)
    : base_type(_p,_size,_map) {}

  EnumValueBase(const std::shared_ptr<value_type>& _p,unsigned _size,
                const std::shared_ptr<map_type>& _map)
    : base_type(_p,_size,_map) {}

  EnumValueBase(const value_type* _p,unsigned _size,const std::shared_ptr<map_type>& _map)
    : base_type(_p,_size,_map) {}

  template <typename GETTER,typename SETTER>
  EnumValueBase(const GETTER& _getter,const SETTER& _setter,unsigned _size,
                const std::shared_ptr<map_type>& _map)
    : base_type(_getter,_setter,_size,_map) {}

  virtual ~EnumValueBase() {}

  virtual value_type get_enum(unsigned _i) const {
    this->check_index(_i);
    return this->get_data(_i);
  }
  virtual void set_enum(const value_type& _data,unsigned _i) {
    this->check_index(_i);
    this->verify(_data);
    this->set_data(_data,_i);
    this->notify_changed();
  }

  void set_native(const value_type& _data,unsigned _i) { this->set_data(_data,_i); }
  value_type get_native(unsigned _i) { return this->get_data(_i); }

  void set(const value_type& _val,unsigned _i=0) { this->set_enum(_val,_i); }
  void _get(value_type& _val,unsigned _i=0) { _val=this->get_enum(_i); }

  virtual void set_int(int _data,unsigned _i) {
    set_enum(appdb_cast<value_type,int>::value(_data),_i);
  }
  virtual void set_double(double /*_data*/,unsigned _i) {
    this->check_index(_i);
    throw BadConversionException(" (set_double())",
                                 typeid(double),typeid(int),VC_DBG_LOCATION());
  }
  virtual void set_binary(const Buffer& _buf) {
    this->_set_binary(_buf,int(0)); // convert to int
  }
  virtual int get_int(unsigned _i) const {
    return int(get_enum(_i));
  }
  virtual double get_double(unsigned _i) const {
    return double(get_enum(_i));
  }
  virtual void get_binary(Buffer& _buf) const {
    this->_get_binary(_buf,int(0)); // convert to int
  }
};

//-----------------------------------------------------------------------------
} // namespace detail
# endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// Symbolic Value accessed via <b>P</b>ointer \ingroup vc_appdb_values
template <typename T,typename TM=T>
class PSymbolicValue
  : public detail::TSymbolicValue<T,detail::PointerAccessor<T>,TM> {
public:
  typedef T value_type;
  typedef detail::TSymbolicValue<T,detail::PointerAccessor<T>,TM> base_type;
  typedef typename base_type::map_type map_type;

  PSymbolicValue(T* _data,unsigned _size=1,
                 const std::shared_ptr<map_type>& _map=sym::symbols<T>::symbol_table()) :
    base_type(_size,_map) {
    this->initialize(_data);
  }
  virtual ~PSymbolicValue() {}
};

/// Symbolic Value accessed via <b>S</b>shared <b>P</b>ointer (`std::shared_ptr`) \ingroup vc_appdb_values
template <typename T,typename TM=T>
class SPSymbolicValue
  : public detail::TSymbolicValue<T,detail::SharedPointerAccessor<T>,TM> {
public:
  typedef T value_type;
  typedef detail::TSymbolicValue<T,detail::SharedPointerAccessor<T>,TM> base_type;
  typedef typename base_type::map_type map_type;

  SPSymbolicValue(const std::shared_ptr<T>& _data,unsigned _size=1,
                  const std::shared_ptr<map_type>& _map=sym::symbols<T>::symbol_table()) :
    base_type(_size,_map) {
    this->initialize(_data);
  }
  virtual ~SPSymbolicValue() {}
};

/// Symbolic Value with its own private <b>S</b>torage \ingroup vc_appdb_values
template <typename T,typename TM=T>
class SSymbolicValue
  : public detail::TSymbolicValue<T,detail::PrivateStorageAccessor<T>,TM> {
public:
  typedef T value_type;
  typedef detail::TSymbolicValue<T,detail::PrivateStorageAccessor<T>,TM> base_type;
  typedef typename base_type::map_type map_type;

  SSymbolicValue(const T* _init,unsigned _size=1,
                 const std::shared_ptr<map_type>& _map=sym::symbols<T>::symbol_table()) :
    base_type(_size,_map) {
    this->initialize(_init,_size);
  }
  SSymbolicValue(const T& _init) :
    base_type(&_init,1) {
    this->initialize(&_init,1);
  }
  virtual ~SSymbolicValue() {}
};

/** Symbolic Value accessed by getter/setter <b>F</b>unctions
    \ingroup vc_appdb_values
    \arg Getter function: `std::function<T()>` or `[](unsigned){}`
    \arg Setter function: `std::function<void(const T&,unsigned)`
    \tparam T value_type
    \tparam SCALAR specialization `SCALAR==true` strips `unsigned` index
    argument from getter and setter functions
    \sa FValue
 */
template <typename T,bool SCALAR=false,typename TM=T>
class FSymbolicValue
  : public detail::TSymbolicValue<T,detail::FunctionCallAccessor<T>,TM> {
public:
  typedef T value_type;
  typedef detail::TSymbolicValue<T,detail::FunctionCallAccessor<T>,TM> base_type;
  typedef typename base_type::map_type map_type;
  typedef typename detail::FunctionCallAccessor<T>::getter_t getter_t;
  typedef typename detail::FunctionCallAccessor<T>::setter_t setter_t;

  FSymbolicValue(const getter_t& _getter,const setter_t& _setter,
         unsigned _size=1,
         const std::shared_ptr<map_type>& _map=sym::symbols<T>::symbol_table()) :
    base_type(_size,_map) {
    this->initialize(_getter,_setter);
  }
  virtual ~FSymbolicValue() {}
};

template <typename T,typename TM>
class FSymbolicValue<T,true,TM>
  : public detail::TSymbolicValue<T,detail::FunctionCallAccessor<T,true>,TM> {
public:
  typedef T value_type;
  typedef detail::TSymbolicValue<T,detail::FunctionCallAccessor<T,true>,TM> base_type;
  typedef typename base_type::map_type map_type;
  typedef typename detail::FunctionCallAccessor<T,true>::getter_t getter_t;
  typedef typename detail::FunctionCallAccessor<T,true>::setter_t setter_t;

  FSymbolicValue(const getter_t& _getter,const setter_t& _setter,unsigned _size,
         const std::shared_ptr<map_type>& _map=sym::symbols<T>::symbol_table()) :
    base_type(1,_map) {
    assert(_size==1);
    this->initialize(_getter,_setter);
  }
  virtual ~FSymbolicValue() {}
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// Enum Value accessed via <b>P</b>ointer \ingroup vc_appdb_values
template <typename T>
class PEnumValue : public detail::EnumValueBase<PSymbolicValue<T,int> > {
public:
  typedef T value_type;
  typedef detail::EnumValueBase<PSymbolicValue<T,int> > base_type;
  typedef typename base_type::map_type map_type;

  PEnumValue(T* _data,unsigned _size=1,
             const std::shared_ptr<map_type>& _map=
             sym::symbols<T>::symbol_table()) :
    base_type(_data,_size,_map) {
    this->initialize(_data);
  }
  virtual ~PEnumValue() {}
};

/** Enum Value accessed via <b>S</b>shared <b>P</b>ointer (`std::shared_ptr`)
    \ingroup vc_appdb_values
*/
template <typename T>
class SPEnumValue : public detail::EnumValueBase<SPSymbolicValue<T,int> > {
public:
  typedef T value_type;
  typedef detail::EnumValueBase<SPSymbolicValue<T,int> > base_type;
  typedef typename base_type::map_type map_type;

  SPEnumValue(const std::shared_ptr<T>& _data,unsigned _size=1,
              const std::shared_ptr<map_type>& _map=
              sym::symbols<T>::symbol_table()) :
    base_type(_data,_size,_map) {
  }
  virtual ~SPEnumValue() {}
};

/// Enum Value with its own private <b>S</b>torage \ingroup vc_appdb_values
template <typename T>
class SEnumValue : public detail::EnumValueBase<SSymbolicValue<T,int> > {
public:
  typedef T value_type;
  typedef detail::EnumValueBase<SSymbolicValue<T,int> > base_type;
  typedef typename base_type::map_type map_type;

  SEnumValue(const T* _data,unsigned _size=1,
             const std::shared_ptr<map_type>& _map=
             sym::symbols<T>::symbol_table()) :
    base_type(_data,_size,_map) {
  }
  SEnumValue(const T& _data,
             const std::shared_ptr<map_type>& _map=
             sym::symbols<T>::symbol_table()) :
    base_type(&_data,1,_map) {
  }
  virtual ~SEnumValue() {}
};

/** Value accessed by getter/setter <b>F</b>unctions
    \ingroup vc_appdb_values

    getter and setter functions provided with the constructor can be
    constructed using `std::function<>`.

    \sa FValue, FSymbolicValue
 */
template <typename T,bool SCALAR=false>
class FEnumValue
  : public detail::EnumValueBase<FSymbolicValue<T,SCALAR,int> > {
public:
  typedef T value_type;
  typedef detail::EnumValueBase<FSymbolicValue<T,SCALAR,int> > base_type;
  typedef typename base_type::map_type map_type;
  typedef typename FSymbolicValue<T,SCALAR,int>::getter_t getter_t;
  typedef typename FSymbolicValue<T,SCALAR,int>::setter_t setter_t;

  FEnumValue(const getter_t& _getter,const setter_t& _setter,
             unsigned _size,
             const std::shared_ptr<map_type>& _map=
             sym::symbols<T>::symbol_table()) :
    base_type(_getter,_setter,_size,_map) {
  }
  virtual ~FEnumValue() {}
};

template <typename T>
class FEnumValue<T,true>
  : public detail::EnumValueBase<FSymbolicValue<T,true,int> > {
public:
  typedef T value_type;
  typedef detail::EnumValueBase<FSymbolicValue<T,true,int> > base_type;
  typedef typename base_type::map_type map_type;
  typedef typename FSymbolicValue<T,true,int>::getter_t getter_t;
  typedef typename FSymbolicValue<T,true,int>::setter_t setter_t;

  FEnumValue(const getter_t& _getter,const setter_t& _setter,
             const std::shared_ptr<map_type>& _map=
             sym::symbols<T>::symbol_table()) :
    base_type(_getter,_setter,1,_map) {
  }
  virtual ~FEnumValue() {}
};

//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_SYM_HH defined
