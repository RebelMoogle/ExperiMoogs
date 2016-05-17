//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_BUFFER_HH
#define VC_APPDB_BUFFER_HH

#include <climits>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>

#include <vclibs/base/endian.hh>

#include "exceptions.hh"

//== INCLUDES =================================================================



//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

//-----------------------------------------------------------------------------

/** Buffer for storing "binary" Value data.
    \ingroup vc_appdb_conv
 */
class Buffer {
public:
  typedef unsigned char byte_t;
  typedef char* cstr_t;

  enum Type {
    // use DBus constants
    INVALID=0,
    BYTE='y',
    BOOLEAN='b',
    INT16='n', UINT16='q',
    INT32='i', UINT32='u',
    INT64='x', UINT64='t',
    DOUBLE='d',
    STRING='s',
    STRUCT='r',
    FLOAT='f' // undefined in DBus protocol
  };

  static bool is_valid_type(int _type) {
    switch (_type) {
    case INVALID: case BYTE:   case BOOLEAN:
    case INT16:   case INT32:  case INT64:
    case UINT16:  case UINT32: case UINT64:
    case DOUBLE: case FLOAT:
    case STRUCT: case STRING:
      return true;
    default:
      return false;
    }
  }

  enum { HeaderSize=4, InitialIndex=HeaderSize };

  Buffer()  { clear(INVALID); }

  void resize(size_t _size) { m_buf.resize(_size+HeaderSize); }
  void clear(Type _type) {
    m_buf.resize(HeaderSize); *((int32_t*) begin())=_type;
  }
  void cp(const byte_t* _data,size_t _n);

  size_t size() const { return m_buf.size(); }
  Type type() const {
    assert(m_buf.size()>=HeaderSize); return Type(*((int32_t*) begin()));
  }

  byte_t* begin() { return &*m_buf.begin(); }
  byte_t* end() { return &*m_buf.end(); }

  const byte_t* begin() const { return &*m_buf.begin(); }
  const byte_t* end() const { return &*m_buf.end(); }

  Buffer& pack(const byte_t* _data,size_t _size) {
    m_buf.resize(m_buf.size()+_size);
    std::memcpy(end()-_size,_data,_size);
    return *this;
  }
  template <typename T>
  Buffer& pack(const T& /*_data*/) {
    throw BadConversionException(" (pack)",
                                 typeid(T),typeid(Buffer),VC_DBG_LOCATION());
  }

  size_t unpack(size_t _index,byte_t* _data,size_t _size) const {
    assert(/*0<=_index &&*/ _index+_size<=size());
    std::memcpy(_data,begin()+_index,_size);
    return _index+_size;
  }
  template <typename T>
  size_t unpack(size_t /*_index*/,T& /*_data*/) const {
    throw BadConversionException(" (unpack)",
                                 typeid(Buffer),typeid(T),VC_DBG_LOCATION());
  }

  template <typename T>
  static Type type_of(const T&) { return INVALID; }

private:
  std::vector<byte_t> m_buf;
};

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP

template<> inline Buffer& Buffer::pack<uint8_t>(const uint8_t& _data) {
  m_buf.push_back(_data); return *this;
}
template<> inline Buffer& Buffer::pack<int8_t>(const int8_t& _data) {
  m_buf.push_back((unsigned char) _data); return *this;
}

template<> inline Buffer& Buffer::pack<int16_t>(const int16_t& _data) {
  using namespace VC::base::endian;
  int16_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}
template<> inline Buffer& Buffer::pack<uint16_t>(const uint16_t& _data) {
  using namespace VC::base::endian;
  uint16_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}

template<> inline Buffer& Buffer::pack<int32_t>(const int32_t& _data) {
  using namespace VC::base::endian;
  int32_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}
template<> inline Buffer& Buffer::pack<uint32_t>(const uint32_t& _data) {
  using namespace VC::base::endian;
  uint32_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}

template<> inline Buffer& Buffer::pack<int64_t>(const int64_t& _data) {
  using namespace VC::base::endian;
  int64_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}
template<> inline Buffer& Buffer::pack<uint64_t>(const uint64_t& _data) {
  using namespace VC::base::endian;
  uint64_t data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}

template<> inline Buffer& Buffer::pack<float>(const float& _data) {
  using namespace VC::base::endian;
  float data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}
template<> inline Buffer& Buffer::pack<double>(const double& _data) {
  using namespace VC::base::endian;
  double data=convert<Little>::native(_data);
  return pack((const byte_t*) &data,sizeof(data));
}

template<> inline Buffer& Buffer::pack<Buffer::cstr_t>(const Buffer::cstr_t& _data) {
  using namespace VC::base::endian;
  uint32_t n=uint32_t(std::strlen(_data));
  pack(n);
  return pack((const byte_t*) _data,n+1);
}
template<> inline Buffer& Buffer::pack<std::string>(const std::string& _data) {
  using namespace VC::base::endian;
  uint32_t n=uint32_t(_data.size());
  pack(n);
  return pack((const byte_t*) _data.c_str(),n+1);
}

//-----------------------------------------------------------------------------

template<> inline size_t
Buffer::unpack<int8_t>(size_t _index,int8_t& _data) const {
  _data=*((int8_t*) begin()); return _index+1;
}
template<> inline size_t
Buffer::unpack<uint8_t>(size_t _index,uint8_t& _data) const {
  _data=*begin(); return _index+1;
}

template<> inline size_t
Buffer::unpack<int16_t>(size_t _index,int16_t& _data) const {
  using namespace VC::base::endian;
  int16_t data;
  unpack(_index,(byte_t*) &data,sizeof(int16_t));
  _data=convert<Native>::from<int16_t,Little>(data);
  return _index+sizeof(int16_t);
}
template<> inline size_t
Buffer::unpack<uint16_t>(size_t _index,uint16_t& _data) const {
  using namespace VC::base::endian;
  uint16_t data;
  unpack(_index,(byte_t*) &data,sizeof(uint16_t));
  _data=convert<Native>::from<uint16_t,Little>(data);
  return _index+sizeof(uint16_t);
}

template<> inline size_t
Buffer::unpack<int32_t>(size_t _index,int32_t& _data) const {
  using namespace VC::base::endian;
  int32_t data;
  unpack(_index,(byte_t*) &data,sizeof(int32_t));
  _data=convert<Native>::from<int32_t,Little>(data);
  return _index+sizeof(int32_t);
}
template<> inline size_t
Buffer::unpack<uint32_t>(size_t _index,uint32_t& _data) const {
  using namespace VC::base::endian;
  uint32_t data;
  unpack(_index,(byte_t*) &data,sizeof(uint32_t));
  _data=convert<Native>::from<uint32_t,Little>(data);
  return _index+sizeof(uint32_t);
}

template<> inline size_t
Buffer::unpack<int64_t>(size_t _index,int64_t& _data) const {
  using namespace VC::base::endian;
  int64_t data;
  unpack(_index,(byte_t*) &data,sizeof(int64_t));
  _data=convert<Native>::from<int64_t,Little>(data);
  return _index+sizeof(int64_t);
}
template<> inline size_t
Buffer::unpack<uint64_t>(size_t _index,uint64_t& _data) const {
  using namespace VC::base::endian;
  uint64_t data;
  unpack(_index,(byte_t*) &data,sizeof(uint64_t));
  _data=convert<Native>::from<uint64_t,Little>(data);
  return _index+sizeof(uint64_t);
}

template<> inline size_t
Buffer::unpack<float>(size_t _index,float& _data) const {
  using namespace VC::base::endian;
  float data;
  unpack(_index,(byte_t*) &data,sizeof(float));
  _data=convert<Native>::from<float,Little>(data);
  return _index+sizeof(float);
}
template<> inline size_t
Buffer::unpack<double>(size_t _index,double& _data) const {
  using namespace VC::base::endian;
  double data;
  unpack(_index,(byte_t*) &data,sizeof(double));
  _data=convert<Native>::from<double,Little>(data);
  return _index+sizeof(double);
}

template<> inline size_t
Buffer::unpack<std::string>(size_t _index,std::string& _data) const {
  uint32_t n;
  _index=unpack(_index,n);
  _data=std::string(((char*) begin())+_index,n);
  assert(m_buf[_index+n]==0);
  return _index+n+1;
}

//-----------------------------------------------------------------------------

template<> inline Buffer::Type
Buffer::type_of(const int8_t&) { return BYTE; }
template<> inline Buffer::Type
Buffer::type_of(const uint8_t&) { return BYTE; }

template<> inline Buffer::Type
Buffer::type_of(const int16_t&) { return INT16; }
template<> inline Buffer::Type
Buffer::type_of(const uint16_t&) { return UINT16; }

template<> inline Buffer::Type
Buffer::type_of(const int32_t&) { return INT32; }
template<> inline Buffer::Type
Buffer::type_of(const uint32_t&) { return UINT32; }

template<> inline Buffer::Type
Buffer::type_of(const int64_t&) { return INT64; }
template<> inline Buffer::Type
Buffer::type_of(const uint64_t&) { return UINT64; }

template<> inline Buffer::Type
Buffer::type_of(const float&) { return FLOAT; }
template<> inline Buffer::Type
Buffer::type_of(const double&) { return DOUBLE; }

template<> inline Buffer::Type
Buffer::type_of(const Buffer::cstr_t&) { return STRING; }
template<> inline Buffer::Type
Buffer::type_of(const std::string&) { return STRING; }

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

inline void Buffer::cp(const byte_t* _data,size_t _n) {
  assert(_n>4);
  m_buf.resize(_n);
  std::memcpy(&*m_buf.begin(),_data,_n);
  assert(is_valid_type(type()));
}

//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_BUFFER_HH defined
