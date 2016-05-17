//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================

#ifndef __VC_BASE_ENDIAN_HH
#define __VC_BASE_ENDIAN_HH

#include <algorithm> // std:swap

#if defined(__VC_USE_BOOST_STDINT)
# include <boost/cstdint.hpp>
using boost::int8_t;
using boost::uint8_t;
using boost::int16_t;
using boost::uint16_t;
using boost::int32_t;
using boost::uint32_t;
using boost::int64_t;
using boost::uint64_t;
#else
# include <cstdint>
#endif

#if defined(_WIN32) || defined(_WIN64)
# include <cstdlib>
# include <intrin.h>
#else
// __GNUC__, clang++
# include <byteswap.h>
#endif

namespace VC {
namespace base {
//=============================================================================

/// Query and convert endianess \ingroup vc_endian
namespace endian {

/** \defgroup vc_base_endian Query and convert endianess
    \ingroup vc_base_memory
 */

# if defined(DOXYGEN_SKIP)

  /// constants describing byte order \ingroup vc_base_endian
  enum ByteOrder {
    Little,        //!< little endian (e.g., Intel)
    Big,           //!< big endian (e.g., SPARC, MIPS)
    Native=Little  //!< native system byte order (equals Little or Big)
  };

  /** @name Swap bytes of operand (converts Little <-> Big endian)
      \ingroup vc_base_endian
      @{
   */

  inline uint16_t swab(uint16_t _value);  //!< swap byte order
  inline uiint32_t swab(uint32_t _value); //!< swap byte order
  inline uint64_t swab(uint64_t _value);  //!< swap byte order

# else

#  if defined(__GLIBC__)

#   include <endian.h>

#   if (__BYTE_ORDER != __LITTLE_ENDIAN && __BYTE_ORDER != __BIG_ENDIAN)
#    error Unknown machine endianness detected.
#   endif

  enum ByteOrder {
    Little = __LITTLE_ENDIAN,
    Big = __BIG_ENDIAN,
    Native = __BYTE_ORDER
  };

#  elif defined(_WIN32) || defined(_WIN64)

  enum ByteOrder {
    Little = 0x1234,
    Big = 0x4321,
    Native = Little
  };
#  elif defined(__MACH__)

#    include <machine/endian.h>

  enum ByteOrder {
    Little = LITTLE_ENDIAN,
    Big = BIG_ENDIAN,
    Native = BYTE_ORDER
  };

#  else
#   error Cannot detect endianess.
#  endif
#  if defined(__GNUC__) && !(defined(__MINGW32__) || defined(__MINGW64__))

#   if defined(__clang__)
#     pragma clang diagnostic push
#     pragma clang diagnostic ignored "-Wdeprecated-register"
#   endif // __clang__

  inline uint16_t swab(uint16_t _value) {
    return bswap_16(_value);
  }
  inline uint32_t swab(uint32_t _value) {
    return bswap_32(_value);
    //return __builtin_bswap32(_value);
  }
  inline uint64_t swab(uint64_t _value) {
    return bswap_64(_value);
    //return __builtin_bswap64(_value);
  }

#   if defined(__clang__)
#     pragma clang diagnostic pop
#   endif // __clang__

#  elif defined _MSC_VER
  inline uint16_t swab(uint16_t _value) {
    return _byteswap_ushort(_value);
  }
  inline uint32_t swab(uint32_t _value) {
    return _byteswap_ulong(_value);
  }
  inline uint64_t swab(uint64_t _value) {
    return _byteswap_uint64(_value);
  }
#  else
  inline uint16_t swab(uint16_t _value) {
    return ((_value >> 8) & 0xff) | ((_value & 0xff) << 8);
  }
  inline uint32_t swab(uint32_t _value) {
    return
      ((_value & 0xff000000) >> 24) | ((_value & 0x00ff0000) >>  8) |
      ((_value & 0x0000ff00) <<  8) | ((_value & 0x000000ff) << 24);
  }
  inline uint64_t swab(uint64_t _value) {
    union { unsigned char c[8]; uint16_t v; } rv;
    rv.v=_value;
    std::swap(rv.c[0],rv.c[7]);
    std::swap(rv.c[1],rv.c[6]);
    std::swap(rv.c[2],rv.c[5]);
    std::swap(rv.c[3],rv.c[4]); return rv.v;
  }
#  endif

# endif // DOXYGEN_SKIP

  inline uint8_t swab(uint8_t _value) {
    return _value;
  }
  inline int8_t swab(int8_t _value) {
    return _value;
  }
  inline int16_t swab(int16_t _value) {
    return int16_t(swab((uint16_t) _value));
  }
  inline int32_t swab(int32_t _value) {
    return int32_t(swab((uint32_t) _value));
  }
  inline int64_t swab(int64_t _value) {
    return int64_t(swab((uint64_t) _value));
  }

  inline float swab(float _value) {
    union { uint32_t i; float v; } rv; rv.v=_value;
    rv.i=swab(rv.i);
    return rv.v;
  }
  inline double swab(double _value) {
    union { uint64_t i; double v; } rv; rv.v=_value;
    rv.i=swab(rv.i);
    return rv.v;
  }

  /// @}

  /** @name Swap byte order for multiple values.
      @{
  */

  /// swab() _n 2-byte values at _ptr
  inline void swab16(void* _ptr,size_t _n) {
    uint16_t* p=(uint16_t*) _ptr;
    for (size_t i=0;i<_n;++i) p[i]=swab(p[i]);
  }
  /// swab() _n 4-byte values at _ptr
  inline void swab32(void* _ptr,size_t _n) {
    uint32_t* p=(uint32_t*) _ptr;
    for (size_t i=0;i<_n;++i) p[i]=swab(p[i]);
  }
  /// swab() _n 8-byte values at _ptr
  inline void swab64(void* _ptr,size_t _n) {
    uint64_t* p=(uint64_t*) _ptr;
    for (size_t i=0;i<_n;++i) p[i]=swab(p[i]);
  }

  /// @}

  /** @name Convert from/to native byte order.
      \ingroup vc_base_endian
      @{
  */

  /// convert from/to native() or from/to() specified byte order
  template <ByteOrder B>
  struct convert {
    template <typename T> static T native(T _value) { return swab(_value); }
    template <typename T,ByteOrder FROM> static T from(T _value) {
      return (FROM==B) ? _value : swab(_value);
    }
  };

# ifndef DOXYGEN_SKIP
  template <>
  struct convert<Native> {
    template <typename T> static T native(T _value) { return _value; }
    template <typename T,ByteOrder FROM> static T from(T _value) {
      return (FROM==Native) ? _value : swab(_value);
    }
  };
# endif

  /// @}

//-----------------------------------------------------------------------------

# ifndef DOXYGEN_SKIP
extern int32_t NATIVE_BYTEORDER; // check if ByteOrder is set correctly
# endif

//=============================================================================
} // namespace endian
} // namespace base
} // namespace VC

# endif // __VC_BASE_ENDIAN_HH
