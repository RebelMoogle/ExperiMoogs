//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_CONV_HH
#define VC_APPDB_CONV_HH


#include <boost/lexical_cast.hpp>
#include <boost/cast.hpp>

#include <vclibs/base/printf.hh>

#include "exceptions.hh"

//== INCLUDES =================================================================


//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

//-----------------------------------------------------------------------------

/** \defgroup vc_appdb_conv Data conversion (AppDB)
    \ingroup vc_appdb

    \b Note: You probably don't have to deal with these classed directly!

    \todo Consider use of `{std,boost}::is_{enum,integral,floating_point}()`
    \todo Consider use of `std::to_string` (in favor of `boost::lexical_cast`?!)
 */

/// Conversion from data type A to B \ingroup vc_appdb_conv
template <typename A,typename B>
class Conv {
public:
  /// convert `_src` to `_dst`, return success
  bool operator()(const A& _src,B& _dst) {
    try   { _dst=boost::lexical_cast<B>(_src);  }
    catch (boost::bad_lexical_cast&) { return false; }
    return true;
  }
};

#ifndef DOXYGEN_SKIP

template <typename A>
class Conv<A,A> {
public:
  bool operator()(const A& _src,A& _dst) { _dst=_src; return true; }
};

template<>
class Conv<int,std::string> {
public:
  bool operator()(const int& _src,std::string& _dst) {
    _dst=VC::base::formatf("%d",_src); return true;
  }
};
template<>
class Conv<float,std::string> {
public:
  bool operator()(const float& _src,std::string& _dst) {
    _dst=VC::base::formatf("%.8g",_src); return true;
  }
};
template<>
class Conv<double,std::string> {
public:
  bool operator()(const double& _src,std::string& _dst) {
    _dst=VC::base::formatf("%.16g",_src); return true;
  }
};
template<>
class Conv<char*,std::string> {
public:
  bool operator()(const char* _src,std::string& _dst) {
    _dst=std::string(_src); return true;
  }
};

template<>
class Conv<int,float> {
public:
  bool operator()(const int& _src,float& _dst) {
    try { _dst=boost::numeric_cast<float>(_src); }
    catch (...) { return false; }
    return true;
  }
};
template<>
class Conv<int,double> {
public:
  bool operator()(const int& _src,double& _dst) {
    try { _dst=boost::numeric_cast<float>(_src); }
    catch (...) { return false; }
    return true;
  }
};
template<>
class Conv<float,double> {
public:
  bool operator()(const float& _src,double& _dst) {
    _dst=_src; return true;
  }
};
template<>
class Conv<double,float> {
public:
  bool operator()(const double& _src,float& _dst) {
    _dst=float(_src); return true; // ignore overflow or underflow! inf->inf!
  }
};

template<>
class Conv<int,unsigned> {
public:
  bool operator()(const int& _src,unsigned& _dst) {
    _dst=boost::numeric_cast<unsigned>(_src); return true;
  }
};
template<>
class Conv<int,char> {
public:
  bool operator()(const int& _src,char& _dst) {
    _dst=boost::numeric_cast<char>(_src); return true;
  }
};
template<>
class Conv<int,unsigned char> {
public:
  bool operator()(const int& _src,unsigned char& _dst) {
    _dst=boost::numeric_cast<unsigned char>(_src); return true;
  }
};
template<>
class Conv<int,long> {
public:
  bool operator()(const int& _src,long& _dst) {
    _dst=_src; return true;
  }
};
template<>
class Conv<int,unsigned long> {
public:
  bool operator()(const int& _src,unsigned long& _dst) {
    _dst=boost::numeric_cast<unsigned long>(_src); return true;
  }
};


template<>
class Conv<unsigned,int> {
public:
  bool operator()(const unsigned& _src,int& _dst) {
    _dst=boost::numeric_cast<int>(_src); return true;
  }
};
template<>
class Conv<char,int> {
public:
  bool operator()(const char& _src,int& _dst) {
    _dst=_src; return true;
  }
};
template<>
class Conv<unsigned char,int> {
public:
  bool operator()(const unsigned char& _src,int& _dst) {
    _dst=_src; return true;
  }
};
template<>
class Conv<long,int> {
public:
  bool operator()(const long& _src,int& _dst) {
    _dst=boost::numeric_cast<int>(_src); return true;
  }
};
template<>
class Conv<unsigned long,int> {
public:
  bool operator()(const unsigned long& _src,int& _dst) {
    _dst=boost::numeric_cast<int>(_src); return true;
  }
};

template<>
class Conv<const char*,int> {
public:
  bool operator()(const char* _src,int& _dst) {
    char* end;
    long value=strtol(_src,&end,0);
    if (end==_src) return false;
    try { _dst=boost::numeric_cast<int>(value); }
    catch (...) { return false; }
    return true;
  }
};
template<>
class Conv<const char*,double> {
public:
  bool operator()(const char* _src,double& _dst) {
    char* end;
    _dst=strtod(_src,&end);
    if (end==_src) return false;
    return true;
  }
};
template<>
class Conv<const char*,float> {
public:
  bool operator()(const char* _src,float& _dst) {
    double value;
    if (Conv<const char*,double>()(_src,value)) {
      try { _dst=boost::numeric_cast<float>(value); return true; }
      catch (...) {};
    }
    return false;
  }
};

template<>
class Conv<std::string,int> {
public:
  bool operator()(const std::string& _src,int& _dst) {
    return Conv<const char*,int>()(_src.c_str(),_dst);
  }
};
template<>
class Conv<std::string,float> {
public:
  bool operator()(const std::string& _src,float& _dst) {
    return Conv<const char*,float>()(_src.c_str(),_dst);
  }
};
template<>
class Conv<std::string,double> {
public:
  bool operator()(const std::string& _src,double& _dst) {
    return Conv<const char*,double>()(_src.c_str(),_dst);
  }
};

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

/// Same as Conv but throw BadConversionException on failure \ingroup vc_appdb_conv
template <typename A,typename B>
class ConvThrow {
public:
  /// convert _src to _dst, throw on failure
  void operator()(const A& _src,B& _dst) throw(BadConversionException) {
    if (!Conv<A,B>()(_src,_dst))
      throw BadConversionException(" (Conv<>)",
                                   typeid(_src),typeid(_dst),VC_DBG_LOCATION());
  }
};

//-----------------------------------------------------------------------------

/// static definition of best type for Value to store T \ingroup vc_appdb_conv
template <typename T> struct BestType {
  typedef std::string type; //!< best type to store T
};
# ifndef DOXYGEN_SKIP
template <> struct BestType<int> { typedef int type; };
template <> struct BestType<double> { typedef double type; };

template <> struct BestType<bool> { typedef int type; };
template <> struct BestType<char> { typedef int type; };
template <> struct BestType<unsigned char> { typedef int type; };
template <> struct BestType<long> { typedef int type; }; // ??
template <> struct BestType<unsigned long> { typedef int type; };
template <> struct BestType<float> { typedef double type; };
# endif // DOXYGEN_SKIP


//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------

//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_CONV_HH defined
