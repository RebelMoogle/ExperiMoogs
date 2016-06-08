//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_GENERICVALUES_HH
#define VC_APPDB_GENERICVALUES_HH

// included by appdb.hh

//== INCLUDES =================================================================



//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

#ifndef DOXYGEN_SKIP
//-----------------------------------------------------------------------------
namespace detail {
//-----------------------------------------------------------------------------

// avoid Visual C++ warning about casting (int/enum to bool)
template <typename A,typename B>
struct appdb_cast {
  static A value(const B& _b) { return A(_b); }
};
template <typename B>
struct appdb_cast<bool,B> {
  static bool value(const B& _b) { return _b!=B(0); }
};

//-----------------------------------------------------------------------------

template <typename T>
struct preferred {
  static Value::TypeHint accessor() { return Value::String; }
};

template <>
struct preferred<int> {
  static Value::TypeHint accessor() { return Value::Int; }
};
template <>
struct preferred<unsigned> {
  static Value::TypeHint accessor() { return Value::Int; }
};
template <>
struct preferred<long> {
  static Value::TypeHint accessor() { return Value::Int; }
};
template <>
struct preferred<unsigned long> {
  static Value::TypeHint accessor() { return Value::Int; }
};
template <>
struct preferred<char> {
  static Value::TypeHint accessor() { return Value::Int; }
};
template <>
struct preferred<unsigned char> {
  static Value::TypeHint accessor() { return Value::Int; }
};
// It's an "enum" type and should use symbol table
// template <>
// struct preferred<bool> {
//   static Value::TypeHint accessor() { return Value::Int; }
// };
template <>
struct preferred<double> {
  static Value::TypeHint accessor() { return Value::Double; }
};
template <>
struct preferred<float> {
  static Value::TypeHint accessor() { return Value::Double; }
};

//-----------------------------------------------------------------------------

template <typename T>
class PointerAccessor {
public:
  typedef T value_type;

  PointerAccessor() : m_data(0) {}

  T* data() { return m_data; }
  const T* data() const { return m_data; }

protected:
  void check_pointer() const {
    if (m_data==0)
      throw InvalidIndexException((Value*) this,0,VC_DBG_LOCATION());
  }
  void initialize(T* _data) { m_data=_data; }

  T get_data(unsigned _i) const { check_pointer(); return m_data[_i]; }
  void set_data(const T& _data,unsigned _i) { check_pointer(); m_data[_i]=_data; }

  T* m_data;
};

template <typename T>
class SharedPointerAccessor {
public:
  typedef T value_type;

  SharedPointerAccessor() : m_data() {}

  const std::shared_ptr<T>& data() const { return m_data; }

protected:
  void check_pointer() const {
    if (m_data.get()==0)
      throw InvalidIndexException((Value*) this,0,VC_DBG_LOCATION());
  }
  void initialize(const std::shared_ptr<T>& _data) { m_data=_data; }

  T get_data(unsigned _i) const {
    check_pointer();
    return m_data.get()[_i];
  }
  void set_data(const T& _data,unsigned _i) {
    check_pointer();
    m_data.get()[_i]=_data;
  }

  std::shared_ptr<T> m_data;
};

template <typename T>
class PrivateStorageAccessor : public PointerAccessor<T> {
public:
  typedef T value_type;

  PrivateStorageAccessor() {}
  ~PrivateStorageAccessor() { delete[] this->m_data; }

protected:
  void initialize(const T* _data,unsigned _size) {
    assert(this->m_data==0);
    this->m_data=new T[_size];
    if (_data!=0)
      for (unsigned i=0;i<_size;++i)
        this->m_data[i]=_data[i];
  }
};

template <typename T,bool SCALAR=false>
class FunctionCallAccessor {
public:
  typedef T value_type;
  typedef std::function<T(unsigned)> getter_t;
  typedef std::function<void(const T&,unsigned)> setter_t;

  getter_t getter;
  setter_t setter;

protected:
  void initialize(const getter_t _getter,const setter_t& _setter) {
    getter=_getter;
    setter=_setter;
  }

  T get_data(unsigned _i) const { return getter(_i); }
  void set_data(const T& _data,unsigned _i) { setter(_data,_i); }
};

template <typename T>
class FunctionCallAccessor<T,true> {
public:
  typedef T value_type;
  typedef std::function<T()> getter_t;
  typedef std::function<void(const T&)> setter_t;

  getter_t getter;
  setter_t setter;

protected:
  void initialize(const getter_t _getter,const setter_t& _setter) {
    getter=_getter;
    setter=_setter;
  }

  T get_data(unsigned _i) const { assert(_i==0); return getter(); }
  void set_data(const T& _data,unsigned _i) { assert(_i==0); setter(_data); }
};

//-----------------------------------------------------------------------------

template <typename T,typename ACCESSORMIXIN>
class TValueBase : public Value, public ACCESSORMIXIN {
public:
  typedef T value_type;

  TValueBase(unsigned _size) :
    Value(typeid(T),_size,detail::preferred<T>::accessor()) {}
  virtual ~TValueBase() {}

  virtual void set_string(const std::string& _data,unsigned _i) {
    check_index(_i);
    T data;
    ConvThrow<std::string,T>()(_data,data);
    this->set_data(data,_i);
    notify_changed();
  }
  virtual std::string get_string(unsigned _i) const {
    check_index(_i);
    std::string sdata;
    ConvThrow<T,std::string>()(this->get_data(_i),sdata);
    return sdata;
  }

protected:
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
      this->set_data(data,i);
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

template <typename T,typename ACCESSORMIXIN>
class TValue : public TValueBase<T,ACCESSORMIXIN> {
public:
  typedef T value_type;

  TValue(unsigned _size) :
    TValueBase<T,ACCESSORMIXIN>(_size) {}
  virtual ~TValue() {}
};

template <typename ACCESSORMIXIN>
class TValue<int,ACCESSORMIXIN> : public TValueBase<int,ACCESSORMIXIN> {
public:
  typedef int value_type;

  TValue(unsigned _size) :
    TValueBase<int,ACCESSORMIXIN>(_size) {}
  virtual ~TValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
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
class TValue<double,ACCESSORMIXIN> : public TValueBase<double,ACCESSORMIXIN> {
public:
  typedef double value_type;

  TValue(unsigned _size) :
    TValueBase<double,ACCESSORMIXIN>(_size) {}
  virtual ~TValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
    this->set_data(double(_data),_i);
    this->notify_changed();
  }
  virtual void set_double(double _data,unsigned _i) {
    this->check_index(_i);
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
class TValue<float,ACCESSORMIXIN> : public TValueBase<float,ACCESSORMIXIN> {
public:
  typedef float value_type;

  TValue(unsigned _size) :
    TValueBase<float,ACCESSORMIXIN>(_size)
  {}
  virtual ~TValue() {}

  virtual void set_int(int _data,unsigned _i) {
    this->check_index(_i);
    this->set_data(float(_data),_i);
    this->notify_changed();
  }
  virtual void set_double(double _data,unsigned _i) {
    this->check_index(_i);
    float data;
    ConvThrow<double,float>()(_data,data);
    this->set_data(data,_i);
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
} // namespace detail
# endif // DOXYGEN_SKIP
//-----------------------------------------------------------------------------

/// Value accessed via <b>P</b>ointer \ingroup vc_appdb_values
template <typename T>
class PValue : public detail::TValue<T,detail::PointerAccessor<T> > {
public:
  typedef T value_type;

  PValue(T* _data,unsigned _size=1) :
    detail::TValue<T,detail::PointerAccessor<T> >(_size) {
    this->initialize(_data);
  }
  virtual ~PValue() {}
};

/// Value accessed via <b>S</b>shared <b>P</b>ointer (`std::shared_ptr`) \ingroup vc_appdb_values
template <typename T>
class SPValue : public detail::TValue<T,detail::SharedPointerAccessor<T> > {
public:
  typedef T value_type;

  SPValue(const std::shared_ptr<T>& _data,unsigned _size=1) :
    detail::TValue<T,detail::PointerAccessor<T> >(_size) {
    this->initialize(_data.get());
  }
  virtual ~SPValue() {}
};

/// Value with its own private <b>S</b>torage \ingroup vc_appdb_values
template <typename T>
class SValue : public detail::TValue<T,detail::PrivateStorageAccessor<T> > {
public:
  typedef T value_type;

  SValue(const T* _init,unsigned _size=1) :
    detail::TValue<T,detail::PrivateStorageAccessor<T> >(_size) {
    static_assert(!std::is_same<T,char>::value,
                  "ambiguous call to make::sv(const char*), use std::string()");
    this->initialize(_init,_size);
  }
  virtual ~SValue() {}
};

/** Value accessed by getter/setter <b>F</b>unctions
    \ingroup vc_appdb_values

    getter and setter functions provided with the constructor can be
    constructed using `std::function<>`.
 */
template <typename T,bool SCALAR=false>
class FValue : public detail::TValue<T,detail::FunctionCallAccessor<T> > {
public:
  typedef T value_type;
  typedef typename detail::FunctionCallAccessor<T>::getter_t getter_t;
  typedef typename detail::FunctionCallAccessor<T>::setter_t setter_t;

  FValue(const getter_t& _getter,const setter_t& _setter,
         unsigned _size) :
    detail::TValue<T,detail::FunctionCallAccessor<T> >(_size) {
    this->initialize(_getter,_setter);
  }
  virtual ~FValue() {}
};

/** Value accessed by getter/setter <b>F</b>unctions
    \ingroup vc_appdb_values

    getter and setter functions provided with the constructor can be
    constructed using `std::function<>`.
 */
template <typename T>
class FValue<T,true> : public detail::TValue<T,detail::FunctionCallAccessor<T,true> > {
public:
  typedef T value_type;
  typedef typename detail::FunctionCallAccessor<T,true>::getter_t getter_t;
  typedef typename detail::FunctionCallAccessor<T,true>::setter_t setter_t;

  FValue(const getter_t& _getter,const setter_t& _setter) :
    detail::TValue<T,detail::FunctionCallAccessor<T,true> >(1) {
    this->initialize(_getter,_setter);
  }
  virtual ~FValue() {}
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_GENERICVALUES_HH defined
