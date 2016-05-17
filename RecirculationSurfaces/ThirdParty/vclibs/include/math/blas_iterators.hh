//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
//
//=============================================================================

#ifndef VC_MATH_BLAS_ITERATORS_HH
#define VC_MATH_BLAS_ITERATORS_HH

# if (defined(_WIN32) || defined(_WIN64)) && !defined(VC_PLATFORM_WINDOWS_MINGW)
# pragma warning(push)
# pragma warning(disable: 4512)
# endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** \defgroup vc_blas_iter BLAS matrix iterators.
    \ingroup vc_blas
    Iterators for matrix access or initialization.

    \arg Initialization from comma-separated lists does \a not use \c template
    list construct.
    \arg Assignment to \c *ii can take special values [\ref vc_blas_iter_special].
    \arg Iterators can be used with matrix_reference_t, interfaces exist, but
    they are generally independent of [\ref vc_blas].
    \arg Focus is on ease of handling, not on efficiency.
    \arg The main purpose of iterators is initialization: there are currently no
    "const" versions.
    \arg Prefer VC::math::blas::Foreach [\ref vc_blas_mat] if order of iteration
    does not matter (always fastest, e.g., column order for GE matrices) and no
    iterator state is required.

    See GE_row_iterator and [\ref vc_blas_iter_special] for examples.

    \todo TR_iterator (DIAG), SY_iterator; TP_Iterator (DIAG), SP_iterator

    \sa VC::math::blas::Foreach
 */

# ifndef DOXYGEN_SKIP

template <typename iterator>
struct _ListInitializer {

  _ListInitializer(iterator& _ii) : ii(_ii) {}

  template <typename val>
  _ListInitializer& operator,(const val& _value) {
    //VC_DBG_P(_value);
    (*ii).pass(_value);
    return *this;
  }

  iterator& ii;
};
# endif

/** \defgroup vc_blas_iter_special Special values for iterator assignment.
    \ingroup vc_blas_iter
    Values references by iterators may be assigned special symbols/values
    (e.g., i_skip_row) which invoke methods of iterator (e.g.,
    GE_row_iterator::skip_row()) which also advance the iterator.
    \code
    GE_row_iterator ii=...;
    GE_column_iterator jj=...;

    *jj=i_skip_column;        // equivalent to GE_column_iterator::skip_column()
    *ii=i_skip_row;           // equivalent to GE_row_iterator::skip_row()
    *ii=i_skip(n);            // equivalent to GE_row_iterator::skip(n)

    *ii=i_end;                // equivalent to GE_row_iterator::assert_end()

    *ii=i_fill(value,n);      // equivalent to GE_row_iterator::fill(_value,n)
    *ii=i_fill_row(value);    // equivalent to GE_row_iterator::fill_row(value)
    *jj=i_fill_column(value); // equivalent to GE_column_iterator::fill_colum(value)

    *ii=i_range(begin,end);   // equivalent to GE_row_iterator::copy(begin,end)

    // Values can be used in combination with comma separated lists:

    float data[]={7,8,9};

    ii << 1,2,3,i_skip_r,
          i_fill_r(4),
          5,6,i_range(data,data+3),i_fill_end(10), i_end;
    \endcode
    \sa _ListInitializer<GE_row_iterator<T> > operator<<(GE_row_iterator<T>& _ii,const value& _value)
 */

/// matrix iterator: fill value \ingroup vc_blas_iter_special
struct i_skip_t {
  i_skip_t(unsigned _n) : n(_n) {}
  unsigned n;     //!< skip count
};
/// matrix iterator: insert _value _n times \ingroup vc_blas_iter_special
inline i_skip_t i_skip(unsigned _n) { return i_skip_t(_n); }

/// matrix iterator: skip current row \ingroup vc_blas_iter_special
enum i_skip_row_t {
  i_skip_row=0, i_skip_r=0 //!< *ii=i_skip_row equals ii.skip_row()
};
/// matrix iterator: skip current row \ingroup vc_blas_iter_special
enum i_skip_column_t {
  i_skip_column=0, i_skip_c //!< *ii=i_skip_column equals ii.skip_column()
};
/// matrix iterator: skip current row \ingroup vc_blas_iter_special
enum i_end_t {
  i_end=0 //!< *ii=i_end equals ii.assert_end()
};

/// matrix iterator: fill value \ingroup vc_blas_iter_special
template <typename T>
struct i_fill_t {
  i_fill_t(const T& _value,unsigned _n) : value(_value), n(_n) {}
  T        value; //!< value
  unsigned n;     //!< fill count
};
/// matrix iterator: insert _value _n times \ingroup vc_blas_iter_special
template <typename T> i_fill_t<T> i_fill(const T& _value,unsigned _n) {
  return i_fill_t<T>(_value,_n);
}

/// matrix iterator: fill current row with value \ingroup vc_blas_iter_special
template <typename T>
struct i_fill_row_t {
  i_fill_row_t(const T& _value) : value(_value) {}
  T value; //!< value
};
/// matrix iterator: fill current row with _value \ingroup vc_blas_iter_special
template <typename T> i_fill_row_t<T> i_fill_row(const T& _value) {
  return i_fill_row_t<T>(_value);
}
/// matrix iterator: synonym for i_fill_row() \ingroup vc_blas_iter_special
template <typename T> i_fill_row_t<T> i_fill_r(const T& _value) {
  return i_fill_row_t<T>(_value);
}

/// matrix iterator: fill current column with value \ingroup vc_blas_iter_special
template <typename T>
struct i_fill_column_t {
  i_fill_column_t(const T& _value) : value(_value) {}
  T value; //!< value
};
/// matrix iterator: fill current column with _value \ingroup vc_blas_iter_special
template <typename T> i_fill_column_t<T> i_fill_column(const T& _value) {
  return i_fill_column_t<T>(_value);
}
/// matrix iterator: synonym for i_fill_column() \ingroup vc_blas_iter_special
template <typename T> i_fill_column_t<T> i_fill_c(const T& _value) {
  return i_fill_column_t<T>(_value);
}

/// matrix iterator: fill until end with value \ingroup vc_blas_iter_special
template <typename T>
struct i_fill_end_t {
  i_fill_end_t(const T& _value) : value(_value) {}
  T value; //!< value
};
/// matrix iterator: fill until end with _value \ingroup vc_blas_iter_special
template <typename T> i_fill_end_t<T> i_fill_end(const T& _value) {
  return i_fill_end_t<T>(_value);
}
/// matrix iterator: synonym for i_fill_end() \ingroup vc_blas_iter_special
template <typename T> i_fill_end_t<T> i_fill(const T& _value) {
  return i_fill_end_t<T>(_value);
}

/// matrix iterator: copy range  \ingroup vc_blas_iter_special
template <typename const_iterator>
struct i_range_t {
  i_range_t(const_iterator _begin,const_iterator _end)
    : begin(_begin), end(_end) {}
  const_iterator begin, end;
};
/// matrix iterator: copy range \ingroup vc_blas_iter_special
template <typename const_iterator>
i_range_t<const_iterator> i_range(const_iterator _begin,const_iterator _end) {
  return i_range_t<const_iterator>(_begin,_end);
}

//-----------------------------------------------------------------------------

/** Proxy to current element of iterator.
    \ingroup vc_blas_iter
    Implements special assignment operations which advance the iterator.
*/
template <typename iterator>
struct iterator_element_proxy {

  typedef iterator iterator_t;
  typedef typename iterator_t::value_type value_type;
  typedef iterator_element_proxy<iterator_t> self_t;

  iterator_element_proxy(iterator_t& _ii) : ii(_ii)  {}

  /// read access to element (via cast)
  operator value_type() const {
    assert(!ii.passed_end() && "exceeded matrix size");
    return ii.current();
  }
  /// write access to element
  value_type& operator=(const value_type& _value) {
    assert(!ii.passed_end() && "exceeded matrix size");
    return ii.current()=_value;
  }

  /** @name pass() \a always advances the iterator
      (This interface is used by _ListInitializer.)
      @{
   */

  template <typename T>
  self_t& pass(const T& _v) {
    assert(!ii.passed_end() && "exceeded matrix size");
    ii.current()=typename iterator::value_type(_v); ++ii;
    return *this;
  }
  self_t& pass(const i_skip_t& _v) { ii.skip(_v.n); return *this; }
  self_t& pass(const i_skip_row_t&) { ii.skip_row(); return *this; }
  self_t& pass(const i_skip_column_t&) { ii.skip_column(); return *this; }
  self_t& pass(const i_end_t&) { ii.assert_end(); return *this; }
  template <typename S> self_t& pass(const i_fill_t<S>& _v) {
    ii.fill(_v.value,_v.n); return *this;
  }
  template <typename S> self_t& pass(const i_fill_row_t<S>& _v) {
    ii.fill_row(_v.value); return *this;
  }
  template <typename S> self_t& pass(const i_fill_column_t<S>& _v) {
    ii.fill_column(_v.value); return *this;
  }
  template <typename S> self_t& pass(const i_fill_end_t<S>& _v) {
    ii.fill_end(_v.value); return *this;
  }
  template <typename const_iterator>
  self_t& pass(const i_range_t<const_iterator>& _v) {
    ii.copy(_v.begin,_v.end); return *this;
  }
  /// @}

  /** @name assign special values
      \sa [\ref vc_blas_iter_special]
      @{
  */

  self_t& operator=(const i_skip_t& _v) {
    ii.skip(_v.n); return *this;
  }
  self_t& operator=(const i_skip_row_t&) {
    ii.skip_row(); return *this;
  }
  self_t& operator=(const i_skip_column_t&) {
    ii.skip_column(); return *this;
  }
  self_t& operator=(const i_end_t&) {
    ii.assert_end(); return *this;
  }

  template <typename T> self_t& operator=(const i_fill_t<T>& _v) {
    ii.fill(_v.value,_v.n); return *this;
  }
  template <typename T> self_t& operator=(const i_fill_row_t<T>& _v) {
    ii.fill_row(_v.value); return *this;
  }
  template <typename T> self_t& operator=(const i_fill_column_t<T>& _v) {
    ii.fill_colum(_v.value); return *this;
  }
  template <typename T> self_t& operator=(const i_fill_end_t<T>& _v) {
    ii.fill_end(_v.value); return *this;
  }

  template <typename const_iterator> self_t&
  operator=(const i_range_t<const_iterator>& _v) {
    ii.copy(_v.begin._v.end); return *this;
  }

  /// @}

  iterator_t& ii; //!< iterator
};

//-----------------------------------------------------------------------------

/** Iterates row-wise over GE matrix.
    \ingroup vc_blas_iter
    Usage:
    \code
    GE_row_iterator ii(ptr,m,n,ld);

    *ii=1; ++ii; *ii=2; ++ii;   // forward iterator

    ii << 4,5,6;                // using lists of values

    *ii=i_fill_row(5);          // using special values

    float data[]={8,9,10};

    ii << 6,i_fill(7,3),        // using special values in list
          i_range(data,data+3),i_fill_end;

    \endcode
    \sa [\ref vc_blas_iter_special], VC::math::blas::ge_mat,
    _ListInitializer<GE_row_iterator<T> > operator<<(GE_row_iterator<T>& _ii,const value& _value)
 */
template <typename T>
class GE_row_iterator {
public:

  typedef GE_row_iterator<T> iterator_t;
  typedef T value_type;

  GE_row_iterator(T* _a,int _m,int _n,int _ld)
    : base(_a), m(_m), n(_n), ld(_ld), i(0), j(0) {}

  /** @name forward iterator interface
      @{
  */

  /// access current element (via iterator_element_proxy)
  iterator_element_proxy<iterator_t> operator*() {
    return iterator_element_proxy<iterator_t>(*this);
  }

  /// advance iterator
  GE_row_iterator& operator++() {
    if (++j==n) {
      j=0; ++base; i=0; --m;
    }

    else
      i+=ld;

    return *this;
  }

  /// @}

  /** @name extra information
      @{
  */

  /// iterator passed end of matrix
  bool passed_end() const { return m<=0; }
  /// return current element
  const T& current() const { return base[i]; }
  /// return current element
  T& current() { return base[i]; }

  /// @}

  /** @name insert and advance iterator
      @{
  */

  /// copy range (_begin,_end) to iterator
  template <typename const_iter>
  GE_row_iterator& copy(const_iter _begin,const_iter _end) {
    for (const_iter ii=_begin;ii!=_end;++ii,++(*this))
      *(*this)=*ii;

    return *this;
  }

  /// repeat _n times: set current element to _value and advance iterator
  GE_row_iterator& fill(const T& _value,int _n=1) {
    for (unsigned k=0;k<_n;++k,++(*this))
      *(*this)=_value;
    return *this;
  }

  /// fill the current row with _value (advances iterator)
  GE_row_iterator& fill_row(const T& _value) {
    int m0=m;
    do {
      *(*this)=_value; ++(*this);
    } while (m0==m);
    return *this;
  }
  /// skip _n entries (advances iterator)
  GE_row_iterator& skip(unsigned _n) {
    for (unsigned i=0;i<+_n;++i)
      ++(*this);
    return *this;
  }
  /// skip (remainder of) current row (advances iterator)
  GE_row_iterator& skip_row() {
    i=j=0; --m; ++base;
    return *this;
  }

  /// fill the matrix until end with _value (advances iterator to end)
  GE_row_iterator& fill_end(const T& _value) {
    do {
      *(*this)=_value; ++(*this);
    } while (m>0);
    return *this;
  }
  /// assert on iterator has reached end ("matrix is filled")
  GE_row_iterator& assert_end() {
    assert(m==0 && j==0);
    return *this;
  }

  /// @}

private:

  friend struct iterator_element_proxy<iterator_t>;

  T*  base;       //!< current row
  int m,n,ld,i,j;
};

/** Iterates over vector.
    \ingroup vc_blas_iter
    Note: no specialization for INC=1 <=> could use GE_column_iterator.
    \sa VC::math::blas::vec,
    _ListInitializer<GE_column_iterator<T> > operator<<(GE_row_iterator<T>& _ii,const value& _value)
 */
template <typename T>
struct vec_iterator : GE_row_iterator<T> {
  typedef GE_row_iterator<T> iterator_t;
  typedef T value_type;

  vec_iterator(T* _a,int _n,int _inc)
    : GE_row_iterator<T>(_a,1,_n,_inc) {}
};

//-----------------------------------------------------------------------------

/** Iterates column-wise over GE matrix.
    \ingroup vc_blas_iter
    \sa VC::math::blas::ge_mat,
    _ListInitializer<GE_column_iterator<T> > operator<<(GE_column_iterator<T>& _ii,const value& _value)
 */
template <typename T>
class GE_column_iterator {
public:

  typedef GE_column_iterator<T> iterator_t;
  typedef T value_type;

  GE_column_iterator(T* _a,int _m,int _n,int _ld)
    : base(_a), m(_m), n(_n), ld(_ld), i(0) {}

  /** @name forward iterator interface
      @{
  */

  /// access current element (via iterator_element_proxy)
  iterator_element_proxy<iterator_t> operator*() {
    return iterator_element_proxy<iterator_t>(*this);
  }

  /// advance iterator
  GE_column_iterator& operator++() {
    if (++i==m) {
      base+=ld; --n; i=0;
    }
    return *this;
  }

  /// @}

  /** @name extra information
      @{
  */

  /// iterator passed end of matrix
  bool passed_end() const { return m<=0; }
  /// return current element
  const T& current() const { return base[i]; }
  /// return current element
  T& current() { return base[i]; }

  /// @}

  /** @name insert and advance iterator
      @{
  */

  /// copy range (_begin,_end) to iterator
  template <typename const_iter>
  GE_column_iterator& copy(const_iter _begin,const_iter _end) {
    for (const_iter ii=_begin;ii!=_end;++ii,++(*this))
      *(*this)=*ii;

    return *this;
  }

  /// repeat _n times: set current element to _value and advance iterator
  GE_column_iterator& fill(const T& _value,int _n=1) {
    for (int k=0;k<_n;++k,++(*this))
      *(*this)=_value;
    return *this;
  }

  /// fill the current column with _value (advances iterator)
  GE_column_iterator& fill_column(const T& _value) {
    int n0=n;
    do {
      *(*this)=_value; ++(*this);
    } while (n0==n);
    return *this;
  }
  /// skip _n entries (advances iterator)
  GE_column_iterator& skip(unsigned _n) {
    for (unsigned i=0;i<+_n;++i)
      ++(*this);
    return *this;
  }
  /// skip (remainder of) current column (advances iterator)
  GE_column_iterator& skip_column() {
    i=0; --n; base+=ld;
    return *this;
  }

  /// fill the matrix until end with _value (advances iterator to end)
  GE_column_iterator& fill_end(const T& _value) {
    do {
      *(*this)=_value; ++(*this);
    } while (n>0);
    return *this;
  }
  /// assert on iterator has reached end ("matrix is filled")
  GE_column_iterator& assert_end() {
    assert(n==0 && i==0);
    return *this;
  }

  /// @}

private:

  friend struct iterator_element_proxy<iterator_t>;

  T*  base;       //!< current column
  int m,n,ld,i;
};

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Store values or list of values in iterator.
    \ingroup vc_blas_iter
    \tparam T value type for matrix iterator
    \tparam value value type, e.g., scalar or a special value [\ref vc_blas_iter_special]
    \param _ii iterator
    \param _value value
    Example
    \code
    GE_row_iterator ii=...;
    ii << 1,2,3,4; // same as *ii=1; ++ii; *ii=2; ++ii; ...
    \endcode
 */
template <typename T,typename value>
inline _ListInitializer<GE_row_iterator<T> >
operator<<(GE_row_iterator<T>& _ii,const value& _value) {
  return _ListInitializer<GE_row_iterator<T> >(_ii),_value;
}

/** Store values or list of values in iterator.
    \ingroup vc_blas_iter
    \tparam T value type for matrix iterator
    \tparam value value type, e.g., scalar or a special value [\ref vc_blas_iter_special]
    \param _ii iterator
    \param _value value
    Example
    \code
    GE_column_iterator ii=...;
    ii << 1,2,3,4; // same as *ii=1; ++ii; *ii=2; ++ii; ...
    \endcode
 */
template <typename T,typename value>
inline _ListInitializer<GE_column_iterator<T> >
operator<<(GE_column_iterator<T>& _ii,const value& _value) {
  return _ListInitializer<GE_column_iterator<T> >(_ii),_value;
}


//-----------------------------------------------------------------------------


//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

# if (defined(_WIN32) || defined(_WIN64)) && !defined(VC_PLATFORM_WINDOWS_MINGW)
# pragma warning(pop)
# endif


#endif // VC_MATH_BLAS_ITERATORS_HH
