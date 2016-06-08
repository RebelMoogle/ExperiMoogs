//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_PRIMTYPE_HH
#define VC_BASE_PRIMTYPE_HH

//== INCLUDES =================================================================

#include <cassert>

#include "../system/platform.hh"  // __RESTRICT

#include "endian.hh"
#include "exception.hh"

//== CLASS DEFINITION =========================================================

namespace std {
class type_info;
}

//-----------------------------------------------------------------------------

namespace VC {
namespace base {

/** \defgroup vc_base_primtype Copy and convert primitive types
    \ingroup vc_base_memory

    Copying and conversion of primitive types such as `char`, `int`,
     etc. (prefer `stdint` types `int8_t`, `int32_t`, etc.).
 */

/// Primitive types \ingroup vc_base_primtype
namespace primtype {

//-----------------------------------------------------------------------------

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4290)
#endif

/** @name Constants for primitive types
    @{
*/

//
// IMPORTANT IMPLEMENTATION NOTE:
//  Several tables in primtype.cc depend on the particular values
//  (their order) of constants in enum ID !!!
//

/** Constants describing primitive data types (e.g. \c int).
    \ingroup vc_base_primtype
 */
enum ID {
  NoPrimType=
      -1,  //!< not a primitive type (we associate \c void with thins one)
  Int8= 0,
  Char= Int8,
  UInt8,
  UChar= UInt8,
  Int16,
  Short= Int16,
  UInt16,
  UShort= UInt16,
  Int32,
  Int= Int32,
  UInt32,
  UInt= UInt32,
  Int64,
  Long= Int64,
  UInt64,
  ULong= UInt64,
  Float,
  Single= Float,
  Float32= Float,
  Double,
  Float64= Double,
  VoidPtr,
  Pointer= VoidPtr,
  NumPrimTypes  //!< total number of primitive types defined
};

/// @}

/** @name Mapping between prmitive type ID and C++ types
 */

/** map ID to C++ type
    \ingroup vc_base_primtype
*/
template <ID P>
struct cxx {
#ifdef DOXYGEN_SKIP
  typedef NATIVE_TYPE value_type;           //!< C++ type
  static ID id() { __ERROR_INVALID_ID(); }  //!< ID
 private:
  __ERROR_INVALID_ID();
#endif
};

#ifndef DOXYGEN_SKIP
template <>
struct cxx<NoPrimType> {
  typedef void value_type;
  static ID id() { return NoPrimType; }
};

template <>
struct cxx<Int8> {
  typedef std::int8_t value_type;
  static ID id() { return Int8; }
};
template <>
struct cxx<UInt8> {
  typedef std::uint8_t value_type;
  static ID id() { return UInt8; }
};
template <>
struct cxx<Int16> {
  typedef std::int16_t value_type;
  static ID id() { return Int16; }
};
template <>
struct cxx<UInt16> {
  typedef std::uint16_t value_type;
  static ID id() { return UInt16; }
};
template <>
struct cxx<Int32> {
  typedef std::int32_t value_type;
  static ID id() { return Int32; }
};
template <>
struct cxx<UInt32> {
  typedef std::uint32_t value_type;
  static ID id() { return UInt32; }
};
template <>
struct cxx<Int64> {
  typedef std::int64_t value_type;
  static ID id() { return Int64; }
};
template <>
struct cxx<UInt64> {
  typedef std::uint64_t value_type;
  static ID id() { return UInt64; }
};
template <>
struct cxx<Float> {
  typedef float value_type;
  static ID id() { return Float; }
};
template <>
struct cxx<Double> {
  typedef double value_type;
  static ID id() { return Double; }
};
template <>
struct cxx<VoidPtr> {
  typedef void* value_type;
  static ID id() { return VoidPtr; }
};
#endif

/** map C++ type to ID
    \ingroup vc_base_primtype
*/
template <typename T>
struct info {
  /// primitive type ID
  static ID id() { return NoPrimType; }
  /// C++ type
  typedef T value_type;
};

#ifndef DOXYGEN_SKIP
template <>
struct info<void> {
  static ID id() { return NoPrimType; }
  typedef void value_type;
};
// int8_t is defined as signed char which is treated != char
template <>
struct info<char> {
  static ID id() { return Int8; }
  typedef std::int8_t value_type;
};

template <>
struct info<std::int8_t> {
  static ID id() { return Int8; }
  typedef std::int8_t value_type;
};
template <>
struct info<std::uint8_t> {
  static ID id() { return UInt8; }
  typedef std::uint8_t value_type;
};
template <>
struct info<std::int16_t> {
  static ID id() { return Int16; }
  typedef std::int16_t value_type;
};
template <>
struct info<std::uint16_t> {
  static ID id() { return UInt16; }
  typedef std::uint16_t value_type;
};
template <>
struct info<std::int32_t> {
  static ID id() { return Int32; }
  typedef std::int32_t value_type;
};
template <>
struct info<std::uint32_t> {
  static ID id() { return UInt32; }
  typedef std::uint32_t value_type;
};
template <>
struct info<std::int64_t> {
  static ID id() { return Int64; }
  typedef std::int64_t value_type;
};
template <>
struct info<std::uint64_t> {
  static ID id() { return UInt64; }
  typedef std::uint64_t value_type;
};
template <>
struct info<float> {
  static ID id() { return Float; }
  typedef float value_type;
};
template <>
struct info<double> {
  static ID id() { return Double; }
  typedef double value_type;
};
template <>
struct info<void*> {
  static ID id() { return VoidPtr; }
  typedef void* value_type;
};
template <>
struct info<const void*> {
  static ID id() { return VoidPtr; }
  typedef const void* value_type;
};
#endif

/// @}

/** @name Descriptions of primitive types
    @{
*/

/** Check for NoPrimType.
    Throw exception if `_primType==NoPrimType` (includes *assertion* on valid
   range).
    \return _primType
 */
inline ID ensure_valid(ID _primType) throw(vc_runtime_error) {
  if (_primType == NoPrimType) throw VC_RUNTIME_ERROR("no primitive type");
  assert(0 <= int(_primType) && int(_primType) < NumPrimTypes);
  return _primType;
}

/// get VC::base::primtype::info from value
template <typename T>
info<T> get_info(const T&) {
  return info<T>();
}

/// get size of _primType in bytes \ingroup vc_base_primtype
size_t get_size_of(ID _primType);

/** Get name for _primType.
    \ingroup vc_base_primtype
    Name is
    \arg human readable and
    \arg a valid C++ type,
    we prefer to return types as defined by `stdint`, e.g., `"int8_t"`.
    \param _primType primitive type
    \return name or 0 (invalid ID or NoPrimType)
*/
const char* name(ID _primType);

/** Get ID from name().
    \ingroup vc_base_primtype
    \param _name as from name() and alternative names including char_code()s
    \return ID or NoPrimType on failure
*/
ID by_name(const char* _name);

/** Get 1-character code describing _primType.
    \ingroup vc_base_primtype
    \param _primType
    \return character in `[a-zA-Z]` or `-1` on failure
*/
int char_code(ID _primType);

/// get ID from character code (see char_code()) or NoPrimType \ingroup
/// vc_base_primtype
ID by_char_code(int _ccode);

/// get id from `_type`
ID by_type(const std::type_info& _type);

/// @}

/** @name Copying and converting primitive types.

    - May apply byte swap (VC::base::endian::swab()), e.g., for reading
      stored data.
    - Converts (see VC::base::primtype::copy_convert_t) without any range
      checks (plain casts).
    -  No conversion is available for VC::base::primtype::VoidPtr.
    @{
*/

/** Function handle for copy conversion.
    \ingroup vc_base_primtype

    Copy _n data items from _src to _dst, may
    VC::base::endian::swab() source data, see copy_convert().
    \param _dst destination
    \param _src source
    \param _n number of items to copy
    **Note:** source and destination **must not overlap**!

    @todo Do we want to convert VoidPtr? (take offsets)
*/
typedef void (*copy_convert_t)(void* __RESTRICT _dst,
                               const void* __RESTRICT _src, size_t _n);

/// argument to copy_convert() \ingroup vc_base_primtype
enum SwaB {
  NO_SWAB= 0,  //!< no byte swap
  SWAB_SRC,    //! swap byte order for source
  SWAB_DST     //!< swap byte order for destination
};

/** Get copy_convert_t function.
    \ingroup vc_base_primtype

    \param _dst destination type
    \param _src source type
    \param _swab apply VC::base::endian::swab()s on `_src` or on _dst or don't
    swap byte order
    \return function handle (see copy_convert_t) or `0` if no conversion is
    available or arguments are invalid
 */
copy_convert_t copy_convert(ID _dst, ID _src, SwaB _swab= NO_SWAB);

/// @}

/** @name Binary reading and writing primitive types
    @{
 */

/** Read `_n` data items from `_in` to `_dst`.
    \ingroup vc_base_primtype

    \param _in input stream
    \param _dst_type type of destination buffer
    \param _dst pointer to destination buffer
    \param _n read `_n` data items (not bytes !)
    \param _src_type type of source data (read from `_in`)
    \param _src_endian byte order of source data
    \return `_in`
    \throw vc_runtime_error on failure (IO error, checking istream::gcount(),
    or invalid types)
 */
std::istream& read(std::istream& _in, ID _dst_type, char* _dst, size_t _n,
                   ID _src_type, endian::ByteOrder _src_endian=
                                     endian::Native) throw(vc_runtime_error);

/** Read `_n` data items from `_in` to `_dst`.
    \ingroup vc_base_primtype

    \param _in input stream
    \param _dst pointer to destination buffer
    \param _n read `_n` data items (not bytes !)
    \param _src_type type of source data (read from `_in`)
    \param _src_endian byte order of source data
    \return `_in`
    \throw vc_runtime_error on failure (IO error, checking `istream::gcount()`,
    or invalid types)
 */
template <typename T>
std::istream& read(std::istream& _in, T* _dst, size_t _n, ID _src_type,
                   endian::ByteOrder _src_endian=
                       endian::Native) throw(vc_runtime_error) {
  return read(_in, info<T>::id(), (char*)_dst, _n, _src_type, _src_endian);
}

/** Write `_n` data items from `_src` to `_out`.
    \ingroup vc_base_primtype

    \param _out output stream
    \param _src_type type of source buffer
    \param _src pointer to destination buffer
    \param _n read `_n` data items (not bytes !)
    \param _out_type type of output data (as written to `_out`)
    \param _out_endian byte order of output data
    \return `_ou`t
    \throw vc_runtime_error on failure (IO error or invalid types)
 */
std::ostream& write(std::ostream& _out, ID _src_type, const char* _src,
                    size_t _n, ID _out_type,
                    endian::ByteOrder _out_endian=
                        endian::Native) throw(vc_runtime_error);

/** Write `_n` data items from `_src` to `_out`.
    \ingroup vc_base_primtype

    \param _out output stream
    \param _src_type type of source buffer
    \param _src pointer to destination buffer
    \param _n read `_n` data items (not bytes !)
    \param _out_type type of output data (as written to `_out`)
    \param _out_endian byte order of output data
    \return `_out`
    \throw vc_runtime_error on failure (IO error or invalid types)
 */
template <typename T>
std::ostream& write(std::ostream& _out, const T* _src, size_t _n, ID _out_type,
                    endian::ByteOrder _out_endian=
                        endian::Native) throw(vc_runtime_error) {
  return write(_out, info<T>::id(), (const char*)_src, _n, _out_type,
               _out_endian);
}

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

/// @}

//=============================================================================
}  // namespace primtype
}  // namespace base
}  // namespace VC
//=============================================================================
#endif  // VC_BASE_PRIMTYPE_HH defined
