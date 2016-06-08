//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_APPDB_EXCEPTIONS_HH
#define VC_APPDB_EXCEPTIONS_HH

#include <vclibs/base/exception.hh>

//== INCLUDES =================================================================


//== CLASS DEFINITION =========================================================

namespace VC {
namespace appdb {

//-----------------------------------------------------------------------------

class AbstractDatabase;
class Value;
class Node;
class ValueMap;

//-----------------------------------------------------------------------------

/** \defgroup vc_appdb_exc Exceptions (AppDB)
    \ingroup vc_appdb
*/

/// exception base class \ingroup vc_appdb_exc
class Exception :  public VC::base::vc_runtime_error {
public:
  Exception(const std::string& _msg,const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  Exception(const std::string& _msg,const VC::base::DbgLocation&) throw();
# endif
  virtual ~Exception() throw();
};

/// base class for exceptions related to a database \ingroup vc_appdb_exc
class DatabaseException :  public Exception {
public:
  AbstractDatabase* db;

  DatabaseException(AbstractDatabase* _db,
                    const std::string& _msg,
                    const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  DatabaseException(AbstractDatabase* _db,const std::string& _msg,
                    const VC::base::DbgLocation&) throw();
# endif
  virtual ~DatabaseException() throw();
};


/// node exception base class \ingroup vc_appdb_exc
class NodeException :  public DatabaseException {
public:
  const std::string key;

  NodeException(AbstractDatabase* _db,const std::string& _key,
            const std::string& _msg,const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  NodeException(AbstractDatabase* _db,const std::string& _key,
                const std::string& _msg,const VC::base::DbgLocation&) throw();
# endif
  virtual ~NodeException() throw();
};

/// invalid key exception \ingroup vc_appdb_exc
class InvalidKeyException : public NodeException {
public:
  InvalidKeyException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  InvalidKeyException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::DbgLocation& _loc) throw();
# endif
  virtual ~InvalidKeyException() throw();
};

/// node exists: AbstractDatabase::add_node() fails \ingroup vc_appdb_exc
class NodeExistsException : public NodeException {
public:
  NodeExistsException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  NodeExistsException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::DbgLocation&) throw();
# endif
  virtual ~NodeExistsException() throw();
};

/// cannot find node, AbstractDatabase::find_node() fails \ingroup vc_appdb_exc
class NoSuchNodeException : public NodeException {
public:
  NoSuchNodeException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  NoSuchNodeException(AbstractDatabase* _db,const std::string& _key,
                      const VC::base::DbgLocation&) throw();
# endif
  virtual ~NoSuchNodeException() throw();
};

/// cannot find node, AbstractDatabase::find_node() fails \ingroup vc_appdb_exc
class NoValueException : public NodeException {
public:
  NoValueException(AbstractDatabase* _db,const std::string& _key,
                   const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  NoValueException(AbstractDatabase* _db,const std::string& _key,
                   const VC::base::DbgLocation&) throw();
# endif
  virtual ~NoValueException() throw();
};

/// value exception base class \ingroup vc_appdb_exc
class ValueException :  public Exception {
public:
  Value* value;

  ValueException(Value* _value,const std::string& _msg,
                 const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  ValueException(Value* _value,const std::string& _msg,
                 const VC::base::DbgLocation&) throw();
# endif
  virtual ~ValueException() throw();
};

/// invalid index exception (e.g., exceed size) \ingroup vc_appdb_exc
class InvalidIndexException : public ValueException {
public:
  unsigned index;

  InvalidIndexException(Value* _value,unsigned _i,
                        const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  InvalidIndexException(Value* _value,unsigned _i,
                        const VC::base::DbgLocation&) throw();
# endif
  virtual ~InvalidIndexException() throw();
};

/// undefined symbol in SymbolTable \ingroup vc_appdb_exc
class UndefinedSymbolException : public Exception {
public:
  std::string map, symbol, value;

  UndefinedSymbolException(const std::string& _map,
                           const std::string& _symbol,const std::string& _value,
                           const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  UndefinedSymbolException(const std::string& _map,
                           const std::string& _symbol,const std::string& _value,
                           const VC::base::DbgLocation&) throw();
# endif
  virtual ~UndefinedSymbolException() throw();
};

/// cannot convert data \ingroup vc_appdb_exc
class BadConversionException : public Exception {
public:
  const std::type_info& src;
  const std::type_info& dst;
  BadConversionException(const std::string& _msg,
                         const std::type_info& _src,const std::type_info& _dst,
                         const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  BadConversionException(const std::string& _msg,
                         const std::type_info& _src,const std::type_info& _dst,
                         const VC::base::DbgLocation&) throw();
# endif
  virtual ~BadConversionException() throw();
};

/// exception thrown by VC::appdb::import_args() with user friendly message
class ParseArgsException : public NodeException {
public:
  ParseArgsException(AbstractDatabase* _db,const std::string& _key,
                      const std::string& _msg,
                      const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  ParseArgsException(AbstractDatabase* _db,const std::string& _key,
                     const std::string& _msg,
                     const VC::base::DbgLocation&) throw();
# endif
  virtual ~ParseArgsException() throw();
};


/// invalid type exception \ingroup vc_appdb_exc
class NoSuchTypeException :  public Exception {
public:
  NoSuchTypeException(const std::string& _type,
                      const VC::base::_DbgErrorLocation& _loc) throw();
# ifdef NDEBUG
  NoSuchTypeException(const std::string& _type,
                      const VC::base::DbgLocation&) throw();
# endif
  virtual ~NoSuchTypeException() throw();
};
//-----------------------------------------------------------------------------

//=============================================================================
} // namespace appdb
} // namespace VC
//=============================================================================
#endif // VC_APPDB_EXCEPTIONS_HH defined
