//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_MAKE_HH
#define VC_APPDB_MAKE_HH

// included by appdb.hh

//== INCLUDES =================================================================


//== CLASS DEFINITION =========================================================

namespace std {
class type_info;
}

namespace VC {
namespace appdb {

/** \defgroup vc_appdb_make Utility functions to create Value instances (APPDB)
    \ingroup vc_appdb_values

    - Creation of Value instances.
    - Type registry (register_type()) allows creation of user defined
      types.
 */

/// utility functions to create Value instances \ingroup vc_appdb_make
namespace make {

/** Create an SValue `<T>` for a given `_type`
    \ingroup vc_appdb_make
    \param _size dimension
    \param _type type
    \return SValue `<T>` created by `std::make_shared`.

    \todo extra functions for creating SEnumValue objects
    (and same for FValue,PValue)
 */
std::shared_ptr<Value> sv(unsigned _size,Value::TypeHint _type);

/** Create a scalar SValue `<T>` for a given _type.
    \ingroup vc_appdb_make
    \param _type type
    \return SValue `<T>` created by `std::make_shared`.
 */
inline std::shared_ptr<Value> sv(Value::TypeHint _type) {
  return sv(1,_type);
}

/// create scalar SValue
template <typename T>
inline std::shared_ptr<Value> sv(const T& _init) {
  return std::static_pointer_cast<Value>
    (std::make_shared<SValue<T> >(&_init,1));
}
template <>
inline std::shared_ptr<Value> sv(const bool& _init) {
  return std::static_pointer_cast<Value>
    (std::make_shared<SEnumValue<bool> >(&_init,1));
}
inline std::shared_ptr<Value> sv(const char* _init) {
  std::string init(_init);
  return std::static_pointer_cast<Value>
    (std::make_shared<SValue<std::string> >(&init,1));
}
template <typename T,typename MAP>
inline std::shared_ptr<Value>
sv(const T& _init,const std::shared_ptr<MAP>& _map) {
  return std::static_pointer_cast<Value>
    (std::make_shared<SEnumValue<T> >(&_init,1,_map));
}

/// create an SValue of size 2
template <typename T>
inline std::shared_ptr<Value> sv(const T& _a,const T& _b) {
  T init[]={_a,_b};
  return std::static_pointer_cast<Value>
    (std::make_shared<SValue<T> >(init,2));
}
template <>
inline std::shared_ptr<Value> sv(const bool& _a,const bool& _b) {
  bool init[]={_a,_b};
  return  std::static_pointer_cast<Value>
    (std::make_shared<SEnumValue<bool> >(init,2));
}

/// create an SValue of size 3
template <typename T>
inline std::shared_ptr<Value> sv(const T& _a,const T& _b,const T& _c) {
  T init[]={_a,_b,_c};
  return  std::static_pointer_cast<Value>
    (std::make_shared<SValue<T> >(init,3));
}
template <>
inline std::shared_ptr<Value> sv(const bool& _a,const bool& _b,const bool& _c) {
  bool init[]={_a,_b,_c};
  return  std::static_pointer_cast<Value>
    (std::make_shared<SEnumValue<bool> >(init,3));
}

/// create an SValue
template <typename T>
inline std::shared_ptr<Value> sv(const T* _init,unsigned _size) {
  static_assert(!std::is_same<T,char>::value,
                "ambiguous call to make::sv(), use std::string()");
  return std::make_shared<SValue<T> >(_init,_size);
}
template <>
inline std::shared_ptr<Value> sv(const bool* _init,unsigned _size) {
  return std::make_shared<SEnumValue<bool> >(_init,_size);
}

//-----------------------------------------------------------------------------

/** Create a PValue `<T>` for some `_data`.
    \ingroup vc_appdb_make
    \param _v data
    \param _size dimension
    \return PValue `<T>` created by `std::make_shared`.
 */
template <typename T>
inline std::shared_ptr<Value> pv(T* _v,unsigned _size=1) {
  return std::static_pointer_cast<Value>
    (std::make_shared<PValue<T> >(_v,_size));
}
template <>
inline std::shared_ptr<Value> pv(bool* _v,unsigned _size) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<PEnumValue<bool> >(_v,_size));
}
inline std::shared_ptr<Value> pv(bool* _v) { return pv(_v,1); }

template <typename T,typename MAP>
inline std::shared_ptr<Value>
pv(T* _v,const std::shared_ptr<MAP>& _map) {
  return std::static_pointer_cast<Value>
    (std::make_shared<PEnumValue<T> >(_v,1,_map));
}

/** Create an SPValue `<T>` for some `_data`.
    \ingroup vc_appdb_make
    \param _v data
    \param _size dimension
    \return SPValue `<T>` created by `std::make_shared`.
 */
template <typename T>
inline std::shared_ptr<Value> spv(const std::shared_ptr<T>& _v,unsigned _size=1) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<SPValue>(_v,_size));
}
template <>
inline std::shared_ptr<Value> spv(const std::shared_ptr<bool>& _v,unsigned _size) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<SPEnumValue<bool> >(_v,_size));
}
inline std::shared_ptr<Value> spv(const std::shared_ptr<bool>& _v) {
  return spv(_v,1);
}
template <typename T,typename MAP>
inline std::shared_ptr<Value>
spv(const std::shared_ptr<T>& _v,const std::shared_ptr<MAP>& _map) {
  return std::static_pointer_cast<Value>
    (std::make_shared<SPEnumValue<T> >(_v,1,_map));
}

/** Create an FValue `<T>` for some data.
    \ingroup vc_appdb_make
    Getter and setter may be _lambda_ expressions.
    \param _getter `std::function<T(unsigned)>` or `[](unsigned){}`
    \param _setter `std::function<void(const T&,unsigned)` or
    `[](const T&,unsigned){}`
    \param _size dimension
    \return FValue `<T>` created by `std::make_shared`.
 */
template <typename T>
inline std::shared_ptr<Value>
fv(const std::function<T(unsigned)>& _getter,
   const std::function<void(const T&,unsigned)>& _setter,unsigned _size) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<FValue<T> >(_getter,_setter,_size));
}

template <>
inline std::shared_ptr<Value>
fv(const std::function<bool(unsigned)>& _getter,
   const std::function<void(const bool&,unsigned)>& _setter,unsigned _size) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<FEnumValue<bool> >(_getter,_setter,_size));
}

template <typename T,typename MAP>
inline std::shared_ptr<Value>
fv(const std::function<T(unsigned)>& _getter,
   const std::function<void(const T&,unsigned)>& _setter,unsigned _size,
   const std::shared_ptr<MAP>& _map) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<FEnumValue<T> >(_getter,_setter,_size,_map));
}

/** Create an FValue `<T,false>` for some _scalar_ data.
    \ingroup vc_appdb_make
    Getter and setter may be _lambda_ expressions.
    \param _getter `std::function<T()>` or `[](){}`
    \param _setter `std::function<void(const T&)` or `[](const T&){}`
    \return FValue `<T>` created by `std::make_shared`.
 */
template <typename T>
inline std::shared_ptr<Value>
fv(const std::function<T(void)>& _getter,
   const std::function<void(const T&)>& _setter) {
  return std::static_pointer_cast<Value>
    (std::make_shared<FValue<T,true> >(_getter,_setter));
}

inline std::shared_ptr<Value>
fv(const std::function<bool(void)>& _getter,
   const std::function<void(const bool&)>& _setter) {
  return std::static_pointer_cast<Value>
    (std::make_shared<FEnumValue<bool,true> >(_getter,_setter));
}

template <typename T,typename MAP>
inline std::shared_ptr<Value>
fv(const std::function<T(unsigned)>& _getter,
   const std::function<void(const T&,unsigned)>& _setter,
   const std::shared_ptr<MAP>& _map) {
  return  std::static_pointer_cast<Value>
    (std::make_shared<FEnumValue<T,true> >(_getter,_setter,_map));
}

//-----------------------------------------------------------------------------

/** Create new node containing help.

    - Creates new node `_node->key()+"/.help"` of size 2 with
      `_text`and `brief`.
    - Assume that `%T`,`%S`,`%V` in text will be replaced by _type_,
      _size_, and _value_, respectively, on output.
    \param _node create help for node
    \param _text description
    \param _brief short description
    \sa VC::appdb::args::get_help()
 */
inline void help(const Node& _node,
                 const std::string& _text,
                 const std::string& _brief=std::string()) {
  const std::string key=_node.key();
  _node.db()->find_or_add_node(key+"/.help")
    ->set_value(make::sv(_text,_brief));
}

//-----------------------------------------------------------------------------

/** Register `_type` with sv(unsigned,const std::string&).
    \ingroup vc_appdb_make

    The types `std::string`, `int`, and `double` are predefined.

    \param _type describes a type; there may be multiple descriptions for
    one type (see register_alias())
    \param _create function that creates a `T` value of a given dimension
    \return `true`
 */
bool register_type(const std::type_info& _type,
                   std::function<std::shared_ptr<Value>(unsigned)> _create);
/** Register alias name for `_type`
    \ingroup vc_appdb_make

    The following alias names are predefined

    | type                                             | name(s)           |
    |--------------------------------------------------|-------------------|
    | NilValue (set to global constant VC::appdb::NIL) | `0`, `nil`, `NIL` |
    | `std::string`                                    | `s`, `string`     |
    | `int`                                            | `i`, `int32`      |
    | `double`                                         | `d`, `float`, `f` |

    \param _type must have been registered using register_type()
    \param _alias alias name
    \return `true`
 */
bool register_alias(const std::type_info& _type,const std::string& _alias)
  throw(NoSuchTypeException);

/** Get constructor Value of `_type` as given to register_type().
    \ingroup vc_appdb_make
 */
std::function<std::shared_ptr<Value>(unsigned)>
get_constructor(const std::type_info& _type) throw(NoSuchTypeException);

/** Get constructor Value of `_type` as given to register_type().
    \ingroup vc_appdb_make
 */
std::function<std::shared_ptr<Value>(unsigned)>
get_constructor(const std::string& _type) throw(NoSuchTypeException);


//-----------------------------------------------------------------------------

/** Create an SValue `<T>` for a given `_type`.
    \ingroup vc_appdb_make
    \param _size dimension
    \param _type type as defined by register_type()
    \return SValue `<T>` created by `std::make_shared`.
 */
inline std::shared_ptr<Value> sv(unsigned _size,const std::type_info& _type)
  throw(NoSuchTypeException) {
  return get_constructor(_type)(_size);
}

/** Create an SValue `<T>` for a given `_type`.
    \ingroup vc_appdb_make
    \param _type type as defined by register_type()
    \return SValue `<T>` created by `std::make_shared`.
 */
inline std::shared_ptr<Value> sv(const std::type_info&_type)
  throw(NoSuchTypeException) {
  return sv(1,_type);
}

/** Create an SValue `<T>` for a given `_type`.
    \ingroup vc_appdb_make
    \param _size dimension
    \param _type type as defined by register_type()
    \return SValue `<T>` created by `std::make_shared`.
 */
inline std::shared_ptr<Value> sv(unsigned _size,const std::string& _type)
  throw(NoSuchTypeException) {
  return get_constructor(_type)(_size);
}

/** Create an SValue `<T>` for a given type `T`.
    \ingroup vc_appdb_make
    \param _size dimension
    \return SValue `<T>` created by `std::make_shared`.
 */
template <typename T>
std::shared_ptr<Value> sv(unsigned _size) throw(NoSuchTypeException) {
  return sv(_size,typeid(T));
}

/** Create an SValue `<T>` for a given type `T`.
    \ingroup vc_appdb_make
    \return SValue `<T>` created by `std::make_shared`.
 */
template <typename T>
std::shared_ptr<Value> sv() throw(NoSuchTypeException) {
  return sv(typeid(T));
}


//=============================================================================
} // namespace make
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_MAKE_HH defined
