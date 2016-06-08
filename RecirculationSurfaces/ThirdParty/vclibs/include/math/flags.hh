//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_BASE_FLAGS
#define VC_BASE_FLAGS

//== INCLUDES =================================================================

#include <cassert>
#include <string>
#include <iostream>
#include <typeinfo>
#include <cstdint>
#include <utility>
#include <vector>
#include <initializer_list>

//== CLASS DEFINITION =========================================================

namespace VC {
namespace base {

# ifndef DOXYGEN_SKIP

std::string
_map_flags_to_string(const std::type_info& _type,std::int64_t _value,
                     const std::string& _delimiter);

std::int64_t
_map_string_to_flags_int(const std::type_info& _type,const std::string& _str,
                         bool _ignore_errors=false);

void
_add_flag_to_string_mapping(const std::type_info& _type,
                            std::int64_t _value,const std::string& _str);

std::int64_t
_compute_range_mask(const std::type_info& _type);

std::vector<std::string> _get_symbols(const std::type_info& _type);

# endif

// TODO: remove this ASAP !!!
# if defined(_MSC_VER)
#  define constexpr
# endif

//-----------------------------------------------------------------------------

/** Type-safe manipulation of bit flags defined as an `enum`.

    Supports

    - casts to `enum` and integer type
    - operators `&,|,^,~`
    - mapping to_string() and parse_string() ("to flags")

    The user is responsible that bit operations generate valid
    results, in particular mind `operator~()`.

    Use define_symbol() or define_symbols() to setup mappings. The
    given pairs `(flag,string)` will be matched *in order* and every
    bit is matched *once*, i.e., string maps may include combination
    masks as "shortcuts" (define them "early"). The value `0` (no bit
    set) is treated as a special symbol.

    If there are bits left without symbols defined, a number is output
    instead. (Numbers and symbols may be mixed.)

    If no mappings exist, numbers are used instead of symbols. Numbers
    are parsed by `std::stol()` with base selected by prefix.

    The characters `|,;:+` and `space` are used as delimiters to
    separate OR-red symbols. The standard delimiter for *output* if
    `|`. Empty strings or empty substrings evaluate to `0`.

    \ingroup vc_base
 */
template <typename ENUM>
class Flags {
public:
  typedef ENUM enum_t;
  typedef typename std::underlying_type<ENUM>::type int_t;

  /// same as method lsb_index() on constants
  static constexpr int lsb_index(ENUM _value) {
    return lsb_index(static_cast<int_t>(_value),0);
  }

  /// empty flags (no bit set)
  Flags() {}
  /// flags from `_value`
  Flags(enum_t _value) : m_value(_value) {}
  /// flags from integer
  explicit Flags(int_t _value) : m_value(static_cast<enum_t>(_value)) {}
  /// set all flags in list
  Flags(std::initializer_list<ENUM> _list) : Flags() {
    for (auto f : _list) *this|=f;
  }
  Flags(const Flags&) = default;
  Flags& operator=(const Flags&) = default;

  /// cast to enum
  operator enum_t() const { return m_value; }
  /// test on empty flags: `to_int()!=0`
  operator bool() const { return to_int()!=0; }
  /// get integer
  int_t to_int() const { return static_cast<int_t>(m_value); }

  /// return index of least significant bit or `-1` if all zero
  int lsb_index() const {
    int n=-1,v=to_int();
    while (v!=0) {
      ++n;
      if ((v&1)!=0)
        return n;
      else
        v>>=1;
    }
    return n;
  }

  /// logical *or*
  Flags operator|(Flags _other) const {
    return Flags(to_int()|_other.to_int());
  }
  /// logical *or*
  Flags& operator|=(Flags _other) { return *this=(*this|_other); }

  /// logical *and*
  Flags operator&(Flags _other) const {
    return Flags(to_int()&_other.to_int());
  }
  /// logical *and*
  Flags& operator&=(Flags _other) { return *this=(*this&_other); }

  /// logical *xor*
  Flags operator^(Flags _other) const {
    return Flags(to_int()^_other.to_int());
  }
  /// logical *xor*
  Flags& operator^=(Flags _other) { return *this=(*this^_other); }

  /// negation (mind that result is not necessarily combination of or-ed values)
  Flags operator~() const {
    return Flags(~to_int());
  }

  /// get string representation
  std::string to_string(const std::string& _delimiter="|") const {
    return _map_flags_to_string(typeid(*this),to_int(),_delimiter);
  }
  /** Convert string to flags.
      \param _string input
      \param _ignore_errors don't throw on errors (fix or/and ignore)
      \return valid flag
      \throw `std::invalid_argument` and `std::range_error`
   */
  Flags& parse_string(const std::string& _string,
                      bool _ignore_errors=false) {
    return *this=Flags
      (int_t(_map_string_to_flags_int(typeid(*this),_string,_ignore_errors)));
  }

  /// define string representation for `this` value
  void define_symbol(const std::string& _str) const {
    _add_flag_to_string_mapping(typeid(*this),to_int(),_str);
  }
  /// define string representations (order may matter)
  static size_t
  define_symbols(std::initializer_list<std::pair<enum_t,std::string>> _list) {
    size_t n=0;
    for (auto const& f : _list) {
      Flags(f.first).define_symbol(f.second);
      ++n;
    }
    return n;
  }
  /// get flags with all bits set for which there are symbols defined
  static Flags defined_range() {
    return Flags(int_t(_compute_range_mask(typeid(Flags<ENUM>))));
  }
  /// get all defined symbols
  static std::vector<std::string> get_symbols() {
    return _get_symbols(typeid(Flags<ENUM>));
  }

private:

  static constexpr int lsb_index(int_t _value,int _count) {
    return (_value!=0 ?
            (_value&1)!=0 ? _count : lsb_index(_value>>1,_count+1) : -1);
  }

  enum_t m_value = static_cast<enum_t>(0);
};

//-----------------------------------------------------------------------------

/** convert `ENUM` to flag
    \ingroup vc_base
    \sa VC::base::Flags
 */
template <typename ENUM>
Flags<ENUM> to_flag(ENUM _f) { return Flags<ENUM>(_f); }

//-----------------------------------------------------------------------------

/** output Flags `_f` in string representation (Flags::to_string())
    \ingroup vc_base
    \sa VC::base::Flags
 */
template <typename ENUM>
std::ostream& operator<<(std::ostream& _out,Flags<ENUM> _f) {
  return _out << _f.to_string();
}

//-----------------------------------------------------------------------------

# if defined(_MSC_VER)
#  undef constexpr // TODO: remove this ASAP !!!
# endif

//=============================================================================
} //namespace base
} //namespace VC
//=============================================================================
#endif // VC_BASE_FLAGS defined
