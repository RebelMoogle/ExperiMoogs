//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_BASE_MEMORY_HH
#define VC_BASE_MEMORY_HH


//== INCLUDES =================================================================

#include <stddef.h>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <vector>

//== CLASS DEFINITION =========================================================
namespace VC {
namespace base {

# ifndef DOXYGEN_SKIP
# ifndef NDEBUG
#  define COUNT_REQUEST(idx) ++sm_n_requests[idx];
#  define REPORT_CALL_TO_MALLOC(idx) ++sm_n_malloc_calls[idx];
  struct _ObjectLocalStorage_DebugReport {
    // counters for ObjectLocalStorage,ObjectSingleLocalStorage
    static size_t  sm_n_malloc_calls[2];
    static size_t  sm_n_requests[2];
    static void report();
  };
  struct _ObjectStorageBase : _ObjectLocalStorage_DebugReport {};
# else
  struct _ObjectStorageBase {};
#  define COUNT_REQUEST(n)
#  define REPORT_CALL_TO_MALLOC(n)
# endif // NDEBUG
# endif // DOXYGEN_SKIP

# ifdef DOXYGEN_SKIP
  /// (internal use only)
  struct _ObjectStorageBase {};
#  define COUNT_REQUEST(n)
#  define REPORT_CALL_TO_MALLOC(n)
# endif

/** \defgroup vc_base_memory Memory management
    \ingroup vc_base
 */

/** \class ObjectLocalStorage memory.hh
    \ingroup vc_base_memory
    \brief Provide `malloc()` and `free()` on pre-allocated storage block.

    **Basic idea:**

    - Preallocate storage block of Capacity bytes within object.
    - Requests to malloc() return pointer to this local storage if
      capacity is available otherwise `std::malloc()` will be called.
    - free() calls `std::free()` if `std::malloc()` had been called, i.e.,
      if pointer does not refer to local storage.
    - There is *no* `realloc()`.

    **Properties and assumptions:**

    - This is not a pool allocation scheme such as `boost::pool`.
    - The mechanism is lightweight and useful if generally only a small
      amount of memory is required, but the maximum amount is unknown.
    - malloc() and `ree() are assumed to be called \a once per data item.
      Essentially local storage is *never* freed until the ObjectLocalStorage
      is destroyed.
    - The user is responsible to call free()! There is no automatic
      tracking of pointers!

    \sa ObjectSingleLocalStorage
 */
template <unsigned _Capacity,unsigned _Alignment=8>
class ObjectLocalStorage : _ObjectStorageBase {
public:
  enum { Capacity=_Capacity };   //!< local storage capacity in bytes
  enum { Alignment=_Alignment }; //!< alignment in bytes

  ObjectLocalStorage() {
    m_idx=0;
    std::memset(m_data,0,Capacity);
  }

  /// Does _p point to local storage?
  bool is_local(void* _p) const {
    return (m_data<=_p && _p<m_data+Capacity);
  }

  /// allocate memory
  void* malloc(std::size_t _sz) {
    COUNT_REQUEST(0);
    if (m_idx+_sz<Capacity) {
      void* p=m_data+m_idx;
      m_idx+=unsigned((_sz/Alignment+1)*Alignment);
      return p;
    }
    else {
      REPORT_CALL_TO_MALLOC(0);
      return std::malloc(_sz);
    }
  }
  /// free memory
  void free(void* _p) {
    if (!is_local(_p))
      std::free(_p);
  }

private:
  unsigned        m_idx;
  unsigned char   m_data[Capacity];
};

//-----------------------------------------------------------------------------

/** \class ObjectSingleLocalStorage memory.hh
    \ingroup vc_base_memory
    \brief Provide `malloc()` and `free()` on preallocated storage block.

    Similar to ObjectLocalStorage but only a \a single request to malloc()
    is assumed. Subsequent requests (other than realloc()) are forwarded to
    `std::malloc()`. Enables realloc().

    \sa ObjectLocalStorage
 */
template <unsigned _Capacity>
class ObjectSingleLocalStorage : _ObjectStorageBase {
public:

  enum { Capacity=_Capacity };   //!< local storage capacity in bytes

  ObjectSingleLocalStorage() : m_size(0) {}

  /// Does _p point to local storage?
  bool is_local(void* _p) const {
    return (m_data<=_p && _p<m_data+Capacity);
  }

  /// allocate memory
  void* malloc(std::size_t _sz) {
    COUNT_REQUEST(1);
    if (m_size==0 && _sz<=Capacity) {
      m_size=_sz;
      return m_data;
    }
    else {
      REPORT_CALL_TO_MALLOC(1);
      return std::malloc(_sz);
    }
  }
  /// free memory
  void free(void* _p) {
    if (!is_local(_p))
      std::free(_p);
    else {
      assert(_p==(void*) m_data);
      m_size=0;
    }
  }
  /// reallocate memory
  void* realloc(void* _p,std::size_t _sz) {
    COUNT_REQUEST(1);
    if (!is_local()) {
      REPORT_CALL_TO_MALLOC(1);
      return std::realloc(_p,_sz);
    }
    else {
      assert(_p==(void*) m_data);
      if (_sz<=Capacity)
        return m_data;
      else {
        void* p=std::malloc(_sz);
        std::memcpy(p,m_data,m_size);
        m_size=0;
        return p;
      }
    }
  }

private:
  char        m_data[Capacity];
  std::size_t m_size;
};

# undef REPORT_CALL_TO_MALLOC
# undef COUNT_REQUEST

//-----------------------------------------------------------------------------

/** \class small_vector_allocator memory.hh
    \brief (internal use as allocator for small_vector)
    \ingroup vc_base_memory
    \internal
    - **DO NOT USE THIS AS A GENERAL ALLOCATOR!**
    - Use small_vector

    \sa small_vector, bounded_small_vector_allocator
 */
template<typename _Tp,std::size_t _Capacity>
class small_vector_allocator {
# ifndef DOXYGEN_SKIP
public:
  enum { Capacity=_Capacity };

  typedef std::size_t    size_type;
  typedef std::ptrdiff_t difference_type;
  typedef _Tp*           pointer;
  typedef const _Tp*     const_pointer;
  typedef _Tp&           reference;
  typedef const _Tp&     const_reference;
  typedef _Tp            value_type;

  template<typename _Tp1>
  struct rebind
  { typedef small_vector_allocator<_Tp1,_Capacity> other; };

  small_vector_allocator() throw() { }

  small_vector_allocator(const small_vector_allocator&) throw() { }

  template<typename _Tp1>
  small_vector_allocator(const small_vector_allocator<_Tp1,_Capacity>&) throw() { }

  ~small_vector_allocator() throw() { }

  pointer
  address(reference __x) const { return &__x; }

  const_pointer
  address(const_reference __x) const { return &__x; }

  // NB: __n is permitted to be 0.  The C++ standard says nothing
  // about what the return value is when __n == 0.
  pointer
  allocate(size_type __n, const void* = 0)
  {
# ifdef __GNUG__
    if (__builtin_expect(__n > this->max_size(), false))
      std::__throw_bad_alloc();
# else
        if (__n>this->max_size())
          throw std::bad_alloc();
# endif

    return static_cast<_Tp*>(m_storage.malloc(__n*sizeof(_Tp)));
  }

  // __p is not permitted to be a null pointer.
  void
  deallocate(pointer __p, size_type)
  { m_storage.free((void*) __p);
  }

  size_type
  max_size() const throw()
  { return std::size_t(-1) / sizeof(_Tp); }

  // _GLIBCXX_RESOLVE_LIB_DEFECTS
  // 402. wrong new expression in [some_] allocator::construct
  void
  construct(pointer __p, const _Tp& __val)
  { ::new((void *)__p) _Tp(__val); }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
  template<typename... _Args>
  void
  construct(pointer __p, _Args&&... __args)
  { ::new((void *)__p) _Tp(std::forward<_Args>(__args)...); }
#endif

  void
  destroy(pointer __p) { __p->~_Tp(); }
private:
  ObjectSingleLocalStorage <_Capacity*sizeof(value_type)> m_storage;

# endif // DOXYGEN_SKIP
};

# ifndef DOXYGEN_SKIP

template<typename _Tp,std::size_t _Capacity>
inline bool
operator==(const small_vector_allocator<_Tp,_Capacity>&,
           const small_vector_allocator<_Tp,_Capacity>&)
{ return true; }

template<typename _Tp,std::size_t _Capacity>
inline bool
operator!=(const small_vector_allocator<_Tp,_Capacity>&,
           const small_vector_allocator<_Tp,_Capacity>&)
{ return false; }

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

# if __cplusplus <= 199711L
/** \class small_vector memory.hh
    \ingroup vc_base_memory
    \brief Replaces `std::vector` for a small number (`_Capacity`) of entries.

    small_vector uses ObjectSingleLocalStorage (within small_vector_allocator)
    to provide storage for _Capacity elements.

    \tparam _Tp type as for `std::allocator`
    \tparam _Capacity number of `_Tp` objects stored in place

    *When to use?* If you are using lots of non-empty vectors with
    relatively few elements and the number of elements is rarely above
    `_Capacity`, then you might gain a substantial advantage.

    - A small_vector *always* consumes at least `_Capacity*sizeof(_Tp)` bytes
      of storage.
    - Memory is preallocated within the object, memory management is
      fast up to a maxmium of `_Capacity` entries. After construction
      `small_vector::capacity()>=_Capacity`.
    - If more entries are required, standard memory allocation is used (see
      ObjectSingleLocalStorage), and there is still `_Capacity*sizeof(Tp)` bytes
      overhead.

    \sa small_vector_allocator, bounded_small_vector
 */
template <typename Tp,std::size_t _Capacity=64>
class small_vector
  : public std::vector<Tp,small_vector_allocator<Tp,_Capacity> > {
public:
  typedef std::vector<Tp,small_vector_allocator<Tp,_Capacity> > vector_type;
  typedef typename vector_type::value_type                      value_type;
  typedef typename vector_type::pointer                         pointer;
  typedef typename vector_type::const_pointer                   const_pointer;
  typedef typename vector_type::reference                       reference;
  typedef typename vector_type::const_reference                 const_reference;
  typedef typename vector_type::iterator                        iterator;
  typedef typename vector_type::const_iterator                  const_iterator;
  typedef typename vector_type::const_reverse_iterator          const_reverse_iterator;
  typedef typename vector_type::reverse_iterator                reverse_iterator;
  typedef typename vector_type::size_type                       size_type;
  typedef typename vector_type::difference_type                 difference_type;
  typedef typename vector_type::allocator_type                  allocator_type;

  small_vector() { vector_type::reserve(_Capacity); }

  explicit
  small_vector(size_type __n, const value_type& _value = value_type()) {
    vector_type::reserve(_Capacity);
    vector_type::resize(__n,_value);
  }
  small_vector(const small_vector& __x) {
    vector_type::reserve(_Capacity);
    vector_type::insert(vector_type::begin(),__x.begin(),__x.end());
  }
  template<typename _InputIterator>
  small_vector(_InputIterator __first, _InputIterator __last) {
    vector_type::reserve(_Capacity);
    vector_type::insert(vector_type::begin(),__first,__last);
  }

};
# else // C++11
template <typename Tp, std::size_t _Capacity=64>
using small_vector=
    std::vector < Tp, small_vector_allocator<Tp, _Capacity> >;
#endif

//-----------------------------------------------------------------------------

/** \class bounded_small_vector_allocator memory.hh
    \brief (internal use as allocator for bounded_small_vector)
    \ingroup vc_base_memory
    \internal
    - <b>DO NOT USE THIS AS A GENERAL ALLOCATOR!</b>
    - Use bounded_small_vector

    \sa bounded_small_vector, small_vector_allocator
 */
template<typename _Tp,std::size_t _Capacity>
class bounded_small_vector_allocator {
# ifndef DOXYGEN_SKIP
public:
  enum { Capacity=_Capacity };

  typedef std::size_t    size_type;
  typedef std::ptrdiff_t difference_type;
  typedef _Tp*           pointer;
  typedef const _Tp*     const_pointer;
  typedef _Tp&           reference;
  typedef const _Tp&     const_reference;
  typedef _Tp            value_type;

  template<typename _Tp1>
  struct rebind
  { typedef bounded_small_vector_allocator<_Tp1,_Capacity> other; };

  bounded_small_vector_allocator() throw() { }

  bounded_small_vector_allocator(const bounded_small_vector_allocator&) throw() { }

  template<typename _Tp1>
  bounded_small_vector_allocator
  (const bounded_small_vector_allocator<_Tp1,_Capacity>&) throw() { }

  ~bounded_small_vector_allocator() throw() { }

  pointer
  address(reference __x) const { return &__x; }

  const_pointer
  address(const_reference __x) const { return &__x; }

  // NB: __n is permitted to be 0.  The C++ standard says nothing
  // about what the return value is when __n == 0.
  pointer
  allocate(size_type __n, const void* = 0)
  {
# ifdef __GNUG__
    if (__builtin_expect(__n > this->max_size(), false))
      std::__throw_bad_alloc();
# else
        if (__n>this->max_size())
          throw std::bad_alloc();
# endif

    return m_storage;
  }

  // __p is not permitted to be a null pointer.
  void
  deallocate(pointer /*__p*/, size_type)
  {
  }

  size_type
  max_size() const throw() { return size_t(_Capacity); }

  // _GLIBCXX_RESOLVE_LIB_DEFECTS
  // 402. wrong new expression in [some_] allocator::construct
  void
  construct(pointer __p, const _Tp& __val)
  { ::new((void *)__p) _Tp(__val); }

#ifdef __GXX_EXPERIMENTAL_CXX0X__
  template<typename... _Args>
  void
  construct(pointer __p, _Args&&... __args)
  { ::new((void *)__p) _Tp(std::forward<_Args>(__args)...); }
#endif

  void
  destroy(pointer __p) { __p->~_Tp(); }
private:
  value_type  m_storage[_Capacity];

# endif // DOXYGEN_SKIP
};

# ifndef DOXYGEN_SKIP

template<typename _Tp,std::size_t _Capacity>
inline bool
operator==(const bounded_small_vector_allocator<_Tp,_Capacity>&,
           const bounded_small_vector_allocator<_Tp,_Capacity>&)
{ return true; }

template<typename _Tp,std::size_t _Capacity>
inline bool
operator!=(const bounded_small_vector_allocator<_Tp,_Capacity>&,
           const bounded_small_vector_allocator<_Tp,_Capacity>&)
{ return false; }

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

# if __cplusplus <= 199711L
/** \class bounded_small_vector memory.hh
    \ingroup vc_base_memory
    \brief Same as small_vector but with maximum _Capacity.

    Provides storage within `vector` object just as small_vector, but is
    not able to grow beyond `_Capacity` elements!

    A `std::bad_alloc` exception is thrown if the size exceeds this limit.
 */
template <typename Tp,std::size_t _Capacity>
class bounded_small_vector
  : public std::vector<Tp,bounded_small_vector_allocator<Tp,_Capacity> > {
public:
  typedef std::vector<Tp,bounded_small_vector_allocator<Tp,_Capacity> > vector_type;
  typedef typename vector_type::value_type                      value_type;
  typedef typename vector_type::pointer                         pointer;
  typedef typename vector_type::const_pointer                   const_pointer;
  typedef typename vector_type::reference                       reference;
  typedef typename vector_type::const_reference                 const_reference;
  typedef typename vector_type::iterator                        iterator;
  typedef typename vector_type::const_iterator                  const_iterator;
  typedef typename vector_type::const_reverse_iterator          const_reverse_iterator;
  typedef typename vector_type::reverse_iterator                reverse_iterator;
  typedef typename vector_type::size_type                       size_type;
  typedef typename vector_type::difference_type                 difference_type;
  typedef typename vector_type::allocator_type                  allocator_type;

  bounded_small_vector() { vector_type::reserve(_Capacity); }

  explicit
  bounded_small_vector(size_type __n, const value_type& _value = value_type()) {
    vector_type::reserve(_Capacity);
    vector_type::resize(__n,_value);
  }
  bounded_small_vector(const bounded_small_vector& __x) {
    vector_type::reserve(_Capacity);
    vector_type::insert(vector_type::begin(),__x.begin(),__x.end());
  }
  template<typename _InputIterator>
  bounded_small_vector(_InputIterator __first, _InputIterator __last) {
    vector_type::reserve(_Capacity);
    vector_type::insert(vector_type::begin(),__first,__last);
  }
};
# else // C++11
template <typename Tp, std::size_t _Capacity>
using bounded_small_vector=
    std::vector < Tp, bounded_small_vector_allocator<Tp, _Capacity> >;
#endif

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace base
} // namespace VC
//=============================================================================
#endif // VC_BASE_MEMORY_HH defined
