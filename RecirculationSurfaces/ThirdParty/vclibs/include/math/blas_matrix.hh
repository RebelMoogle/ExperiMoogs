//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================

//
// have a "detail" namespace / doxygen group
//


#ifndef VC_MATH_BLAS_MATRIX_HH
#define VC_MATH_BLAS_MATRIX_HH

#include <cassert>
#include <cmath>
#include <string>
#include <algorithm> // std::min
#include <iostream>  // operator<<
#include <limits>
#include <limits.h>  // require constant INT_MAX

#include "blas.hh"

//=============================================================================

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable:4100) // unreferenced formal parameter
# pragma warning(disable:4127) // conditional expression is constant
#endif

#include "blas_iterators.hh" // mostly independent, we provide interfaces
#include "mat4.hh"           // completely independent, we provide interfaces

# include "blas_matrix_inc.hh"
# include "blas_matrix_error_inc.hh"
# include "blas_matrix_idx_range_inc.hh"
# include "blas_matrix_foreach_inc.hh"
# include "blas_matrix_detail_inc.hh"
# include "blas_matrix_expr_inc.hh"
# include "blas_matrix_io_inc.hh"
# include "blas_matrix_at_inc.hh"
# include "blas_matrix_functions_inc.hh"

#ifdef _MSC_VER
# pragma warning(pop)
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** (used internally: base class of vector_const_reference_t)
    \ingroup vc_blas
 */
template <typename T,typename V>
struct vector_const_reference_base_t {
  typedef V vector_t;
  typedef T value_t;
  vector_const_reference_base_t(const V& _v,const T* _data) : m_v(_v), m_data(_data) {}

  const V& vector() const { return m_v; }
  const T* data() const { return m_data; }
  operator const T*() const { return m_data; }

protected:
  V        m_v;
  const T* m_data;
};

/** (used internally: base class of vector_reference_t)
    \ingroup vc_blas
 */
template <typename T,typename V>
struct vector_reference_base_t {
  typedef V vector_t;
  typedef T value_t;
  vector_reference_base_t(const V& _v,T* _data) : m_v(_v), m_data(_data) {}
  const V& vector() const { return m_v; }
  T* data() const { return m_data; }
  operator T*() const { return m_data; }

protected:
  V  m_v;
  T* m_data;
};

/** (used internally: mixin to vector_const_reference_t and vector_reference_t)
    \ingroup vc_blas
    Mixin common ("const") functionality.
 */
template <typename REFBASE>
struct vector_const_reference_mixin
  : public REFBASE {

  typedef typename REFBASE::value_t  value_t;
  typedef typename REFBASE::vector_t vector_t;
  typedef vector_const_reference_t<value_t,vector_t> const_reference_t;

  vector_const_reference_mixin(const vector_t& _v,value_t* _data)
    : REFBASE(_v,_data) {}
  vector_const_reference_mixin(const vector_t& _v,const value_t* _data)
    : REFBASE(_v,_data) {}


  /** @name vec properties
      @{
  */

  int n() const { return this->m_v.n(); }      //!< number of elements
  int inc() const { return this->m_v.inc(); }  //!< increment

  /// @}

  /// get const_reference_t
  const_reference_t const_ref() const {
    return const_reference_t(this->m_v,this->m_data);
  }
  /// get element
  const value_t& operator()(int _i) const { return at(const_ref(),_i); }
   /// get element
  const value_t& operator[](int _i) const { return at(const_ref(),_i); }
  /// get element
  const value_t& operator()(const idx_end& _i) const { return at(const_ref(),this->n()-_i.offset-1); }
   /// get element
  const value_t& operator[](const idx_end& _i) const { return at(const_ref(),this->n()-_i.offset-1); }

  /// get reference to block (C++11)
  template <bool END1,bool END2>
  auto operator()(const idx_range<END1,END2>& _ri) const
    -> decltype(block(*this,_ri))
  {
    return block(*this,_ri);
  }

  /** Dot product with _x.
      \f[ dot \leftarrow x^T y\f]
      \sa VC::math::blas::dot()
   */
  template <typename W>
  value_t dot(const vector_const_reference_t<value_t,W>& _x) const;

  template <typename W>
  value_t dot(const vector_reference_t<value_t,W>& _x) const {
    return this->dot(_x.const_ref());
  }

  /** Double precision dot product with _x.
      \f[ dot \leftarrow x^T y\f]
      \sa VC::math::blas::ddot()
   */
  template <typename W>
  double ddot(const vector_const_reference_t<value_t,W>& _x) const;

  template <typename W>
  double ddot(const vector_reference_t<value_t,W>& _x) const {
    return this->ddot(_x.const_ref());
  }

  /** Get 2-norm.
      \f[ nrm2 \leftarrow ||x||_2\f]
      \sa VC::math::blas::nrm2()
   */
  value_t nrm2() const;

  /** Get 1-norm.
      \f[ asum \leftarrow ||x||_1\f]
      \sa VC::math::blas::asum()
   */
  value_t asum() const;
};

/** Reference to a constant vector.
    \ingroup vc_blas
    \sa vector_reference_t
 */
template <typename T,typename V>
struct vector_const_reference_t
  : public vector_const_reference_mixin<vector_const_reference_base_t<T,V> > {

  typedef vector_reference_t<T,V> self_t;
  typedef T value_t;
  typedef V vector_t;

  vector_const_reference_t(const V& _v,const T* _data)
    : vector_const_reference_mixin<vector_const_reference_base_t<T,V> > (_v,_data) {}

  vector_const_reference_t(const vector_reference_t<T,V>& _v)
    : vector_const_reference_mixin<vector_const_reference_base_t<T,V> >
      (_v.vector(),_v.data()) {}

  /// get same reference with variable n(), inc()
  vector_const_reference_t<T,vec<VarInt,VarInt> > to_var() const {
    return vector_const_reference_t<T,vec<VarInt,VarInt> >
      (vec<VarInt,VarInt>(this->n(),this->inc()),this->data());
  }

  /** @name set properties
      @{
  */

  /// set dimension (only if N=VarInt)
  self_t& set_n(int _n) { this->m_v.set_n(_n); return *this; }
  /// set dimension (only if INC=VarInt)
  self_t& set_inc(int _inc) { this->m_v.set_inc(_inc); return *this; }
  /// set data pointer
  self_t& set_data(const value_t* _data) { this->m_data=_data; return *this; }

  /// @}
};


/** Reference to a vector.
    \ingroup vc_blas

    \arg All operations that yield a vector as result modify \c *this
    vector!

    \arg Some operations are defined only for specific matrices
    (triangular, non-triangular, see mv(), sv())!

    \arg Operands must \a not refer to destination data (this->data())!
    For example, `x=x+y` ins \a invalid (run-time
    assertion fails) -- use `x+=y` (`x=y+x` is
    fine, too, in this case.)

    \arg The operations below can be combined to sums, e.g., `x=2*y+A*x`

    The following operators are defined

    \code
    value_t                    alpha;
    vector_reference_t<>       y;
    vector_const_reference_t<> x;
    matrix_const_reference_t<> A;

    x(i);          // access element (VC::math::blas::at(int))

    y=zeros;       // load zero vector (ld_all())
    y=all(alpha);  // load vector with constant elements (ld_all())
    y=ones;        // load vector with all elements equal to 1 (ld_all())
    y=unit(i,alpha);  // load (unit) vector with only i-th element !=0 (ld_unit())
    y=random_values;  // load random values in [0,1] (ld_rand())

    y+=alpha;      // element-wise add (adds())
    y-=alpha;      // element-wise add (adds())

    y=x;           // copy()
    y=alpha*x;     // copy_scal()

    y*=alpha;      // scal()
    y/=alpha;      // scal()

    y=-x;          // axpy()
    y=y*alpha;     // axpy(), multiplication by scalar from right

    y=A*x;         // mv(), NOT for triangular matrices (TR,TB,SP)
    y=A*x*alpha;   // mv(), NOT for triangular matrices (TR,TB,SP)

    y*=A;          // mv(), ONLY for triangular matrices (TR,TB,SP)
    y/=A;          // triangular solve sv(), ONLY for triangular matrices (TR,TB,SP)

    y+=x;          // axpy()
    y-=x;          // axpy()
    y+=x*alpha;    // axpy()
    y-=x*alpha;    // axpy()
    y+=A*x*alpha;  // mv()
    y-=A*x*alpha;  // mv()

    x=A.column(j); // cp_column()
    x=A.row(i);    // cp_row()
    x=A.diag();    // cp_diag()

    x=column(A,j); // vector reference (VarInt), no copy; GE only
    x=row(A,i);    // vector reference (VarInt), no copy; GE only
    x=diag(A);     // vector reference (VarInt), no copy; GE only
    x=block(y,i,j);// get vector reference (VarInt) to sub vector

    std::cout << x
              << pp(x)          // pretty print
              << matlab(x);     // print as Matlab vector

    write_mat4(cout,trans(trans(v))); // write into MAT-file (V4)

    \endcode

    \sa vector_const_reference_t, matrix_reference_t
*/
template <typename T,typename V>
struct vector_reference_t
  :  public vector_const_reference_mixin<vector_reference_base_t<T,V> >  {

  typedef vector_reference_t<T,V> self_t;
  typedef T value_t;
  typedef V vector_t;

  vector_reference_t(const V& _v,T* _data)
    : vector_const_reference_mixin<vector_reference_base_t<T,V> >(_v,_data) {}
  template <typename W>
  vector_reference_t(const vector_reference_t<T,W>& _w)
    : vector_const_reference_mixin<vector_reference_base_t<T,V> >
      (_w.vector(),_w.data()) {}

  /// get same reference with variable n(), inc()
  vector_reference_t<T,vec<VarInt,VarInt> > to_var() const {
    return vector_reference_t<T,vec<VarInt,VarInt> >
      (vec<VarInt,VarInt>(this->n(),this->inc()),this->data());
  }

  /** @name set properties
      @{
  */

  /// set dimension (only if N=VarInt)
  self_t& set_n(int _n) { this->m_v.set_n(_n); return *this; }
  /// set dimension (only if INC=VarInt)
  self_t& set_inc(int _inc) { this->m_v.set_inc(_inc); return *this; }
  /// set data pointer
  self_t& set_data(value_t* _data) { this->m_data=_data; return *this; }

  /// @}


  /** @name BLAS functions.
      \arg Method names map to respective BLAS functions.
      \arg The destination of operations (\f$lhs \leftarrow\f$) is always \c *this!
      @{
  */

  /** Swap contents with _x.
      \f[ x \leftrightarrow y\f]
      \sa VC::math::blas::swap()
   */
  template <typename W>
  const self_t& swap(const vector_reference_t<T,W>& _x) const;

  /** Scale *this.
      \f[ x \leftarrow \alpha x\f]
      \sa VC::math::blas::scal()
   */
  const self_t& scal(const T&  _alpha) const;

  /** Copy contents from _x.
      \f[ y \leftrightarrow x\f]
      \sa VC::math::blas::copy()
   */
  template <typename W>
  const self_t& copy(const vector_const_reference_t<T,W>& _x) const;

  /** Load _alpha*_x.
      \f[ y \leftrightarrow \alpha x\f]
      \sa VC::math::blas::copy_scal()
   */
  template <typename W>
  const self_t& copy_scal(const T& _alpha,const vector_const_reference_t<T,W>& _x) const;

  /** VC::math::blas::axpy():
      \f[ y \leftarrow \alpha x + y\f]
      \sa VC::math::blas::copy()
   */
  template <typename W>
  const self_t& axpy(const T& _alpha,const vector_const_reference_t<T,W>& _x) const;


  /** Matrix vector multiplication for \b general and \b symmetric matrices.
      \f[ y \leftarrow \alpha A x + \beta y \f]
      \sa VC::math::blas::mv()
   */
  template <typename A,typename W>
  const self_t& mv(const T& _alpha,
                   const matrix_const_reference_t<T,A>& _A,
                   const vector_const_reference_t<T,W>& _x,const T& _beta) const;

  /** Matrix vector multiplication for \b triangular matrices.
      \f[ x \leftarrow A x \f]
      \sa VC::math::blas::mv()
   */
  template <typename A>
  const self_t& mv(const matrix_const_reference_t<T,A>& _A) const;


  /** Solve \b triangular system.
      \f[ x \leftarrow A^{-1} x \f]
      \sa VC::math::blas::sv()
   */
  template <typename A>
  const self_t& sv(const matrix_const_reference_t<T,A>& _A) const;

  /// @}

  /** @name BLAS extensions.
      \arg The destination of operations (\f$lhs \leftarrow\f$) is always \c *this!
      @{
  */

  /// load _a to all elements
  const self_t& ld_all(const T& _a) const;
  /// load zeros (ld_all(0))
  const self_t& ld_zero() const { ld_all(T(0)); return *this; }
  /// load zero vector and set x[_i]=_alpha
  const self_t& ld_unit(int _i,const T& _a=T(1)) const;
  /// load random values uniformly distributed in in [0,1]
  const self_t& ld_rand() const;

  /// add _alpha to all elements
  const self_t& adds(const T& _alpha) const;

  /// @}

  /// access element _i
  T& operator()(int _i) const { return at(*this,_i); }
  /// access element _i
  T& operator[](int _i) const { return at(*this,_i); }
  /// access element
  T& operator()(const idx_end& _i) const { return at(*this,this->n()-_i.offset-1); }
  /// access element
  T& operator[](const idx_end& _i) const { return at(*this,this->n()-_i.offset-1); }

  /// get reference to block (C++11)
  template <bool END1,bool END2>
  auto operator()(const idx_range<END1,END2>& _ri) const
    -> decltype(block(*this,_ri))
  {
    return block(*this,_ri);
  }


  /** @name evaluate expressions
      The Expression classes define functionality for operator=() and
      operator+=().
      @{
   */

  template <typename EXPR,typename ARG1,typename ARG2>
  const self_t& operator=(const Expression<Type_vector,EXPR,T,ARG1,ARG2>& _expr) const {
    _expr.assign(*this);
    return *this;
  }
  template <typename EXPR,typename ARG1,typename ARG2>
  const self_t& operator+=(const Expression<Type_vector,EXPR,T,ARG1,ARG2>& _expr) const {
    _expr.plus_assign(*this);
    return *this;
  }

  template <typename EXPR,typename ARG1,typename ARG2,typename ARG3>
  const self_t& operator=(const Expression<Type_initializer,EXPR,ARG1,ARG2,ARG3>& _expr) const {
    _expr.assign(*this);
    return *this;
  }
  template <typename EXPR,typename ARG1,typename ARG2,typename ARG3>
  const self_t& operator+=(const Expression<Type_initializer,EXPR,ARG1,ARG2,ARG3>& _expr) const {
    _expr.plus_assign(*this);
    return *this;
  }

  template <typename EXPR>
  const self_t& operator-=(const EXPR& _expr) const {
    return *this+=(-_expr); // DELEGATE ALL
  }

  /// @}

  /// assign y=x (copy()) -- never assign reference!
  const self_t& operator=(const vector_const_reference_t<value_t,vector_t>& _x) const {
    copy(_x);
    return *this;
  }
  /// assign y=x (copy()) -- never assign reference!
  const self_t& operator=(const vector_reference_t<value_t,vector_t>& _x) const {
    copy(_x.const_ref());
    return *this;
  }
  /// assign y=x (copy())
  template <typename VX>
  const self_t& operator=(const vector_const_reference_t<value_t,VX>& _x) const {
    copy(_x);
    return *this;
  }
  /// assign y=_data (cp(), no checks on _data bounds)
  template <typename S>
  const self_t& operator=(const S* _data) const {
    cp(_data,*this);
    return *this;
  }

  /// @}

  /// assign y=x (copy())
  template <typename VX>
  const self_t& operator=(const vector_reference_t<value_t,VX>& _x) const {
    copy(_x.const_ref());
    return *this;
  }
  /// assign y+=x (axpy())
  template <typename VX>
  const self_t& operator+=(const vector_const_reference_t<value_t,VX>& _x) const {
    axpy(value_t(1),_x);
    return *this;
  }
  /// assign y+=x (axpy())
  template <typename VX>
  const self_t& operator+=(const vector_reference_t<value_t,VX>& _x) const {
    axpy(value_t(1),_x.const_ref());
    return *this;
  }
  /// assign y-=x (axpy())
  template <typename VX>
  const self_t& operator-=(const vector_const_reference_t<value_t,VX>& _x) const {
    axpy(value_t(-1),_x);
    return *this;
  }
  /// assign y-=x (axpy())
  template <typename VX>
  const self_t& operator-=(const vector_reference_t<value_t,VX>& _x) const {
    axpy(value_t(-1),_x.const_ref());
    return *this;
  }


  template <int NX,int INCX>
  const self_t& operator+=(const vector_const_reference_t<T,vec<NX,INCX> >& _x) const {
    return axpy(T(1),_x);
  }
  template <int NX,int INCX>
  const self_t& operator-=(const vector_const_reference_t<T,vec<NX,INCX> >& _x) const {
    return axpy(T(-1),_x);
  }


  /// y*=a (scal())
  const self_t& operator*=(const T& _a) const {
    scal(_a);
    return *this;
  }
  /// y/=a (scal())
  const self_t& operator/=(const T& _a) const {
    scal(value_t(1)/_a);
    return *this;
  }

  /// y*=x (elementwise multiplication emul())
  template <typename VX>
  const self_t& operator*=(const vector_const_reference_t<T,VX>& _x) const {
    emul(T(1),_x,*this);
    return *this;
  }
  /// y*=x (elementwise multiplication emul())
  template <typename VX>
  const self_t& operator*=(const vector_reference_t<T,VX>& _x) const {
    emul(T(1),_x.const_ref(),*this);
    return *this;
  }
  /// y*=a*x (elementwise multiplication emul())
  template <typename VX>
  const self_t& operator*=(const Expression<Type_vector,Expr_ax,T,VX>& _expr) const {
    emul(_expr.a,_expr.x,*this);
    return *this;
  }
  /// y/=a (elementwise division ediv())
  template <typename VX>
  const self_t& operator/=(const vector_const_reference_t<T,VX>& _x) const {
    ediv(T(1),_x,*this);
    return *this;
  }
  /// y/=a (elementwise division ediv())
  template <typename VX>
  const self_t& operator/=(const vector_reference_t<T,VX>& _x) const {
    ediv(T(1),_x.const_ref(),*this);
    return *this;
  }
  /// y/=a*x (elementwise division ediv())
  template <typename VX>
  const self_t& operator/=(const Expression<Type_vector,Expr_ax,T,VX>& _expr) const {
    ediv(_expr.a,_expr.x,*this);
    return *this;
  }

  /// *left* multiplication y=A*y \a triangular matrices only! (mv())
  template <typename MA>
  const self_t& operator*=(const matrix_const_reference_t<value_t,MA>& _A) const {
    static_assert(std::is_same<typename MA::matrix_id_t,TR_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TB_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TP_Mat>::value,
                  "self-multiplication x=A*x requires for TR,TB,TP matrix A");
    mv(_A);
    return *this;
  }
  /// y=A*y \a triangular matrices only! (mv())
  template <typename MA>
  const self_t& operator*=(const matrix_reference_t<value_t,MA>& _A) const {
    static_assert(std::is_same<typename MA::matrix_id_t,TR_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TB_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TP_Mat>::value,
                  "self-multiplication x=A*x requires for TR,TB,TP matrix A");
    return *this*=_A.const_ref();
  }
  /// solve \a triangular system Ax=b ("x=A\x") (sv())
  template <typename MA>
  const self_t& operator/=(const matrix_const_reference_t<value_t,MA>& _A) const {
    static_assert(std::is_same<typename MA::matrix_id_t,TR_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TB_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TP_Mat>::value,
                  "triangular solve x/=A requires for TR,TB,TP matrix A");
    sv(_A);
    return *this;
  }
  /// solve \a triangular system Ax=b ("x=A\x") (sv())
  template <typename MA>
  const self_t& operator/=(const matrix_reference_t<value_t,MA>& _A) const {
    static_assert(std::is_same<typename MA::matrix_id_t,TR_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TB_Mat>::value ||
                  std::is_same<typename MA::matrix_id_t,TP_Mat>::value,
                  "triangular solve x/=A requires for TR,TB,TP matrix A");
    return *this/=_A.const_ref();
  }

  /// @}

};

//-----------------------------------------------------------------------------

/** \class matrix_const_reference_base_t
    (used internally: base class of matrix_const_reference_t)
    \ingroup vc_blas
    \tparam T `float` or `double`
    \tparam M matrix type, e.g., VC::math::blas::ge_mat
 */
template <typename T,typename M>
struct matrix_const_reference_base_t  {

  typedef T value_t;
  typedef M matrix_t;

  matrix_const_reference_base_t(const M& _m,const T* _data)
  : m_m(_m), m_data(_data) {}

  const M& matrix() const { return m_m; }
  const T* data() const { return m_data; }
  operator const T*() const { return m_data; }

protected:
  M        m_m;
  const T* m_data;
};

/** \class matrix_reference_base_t
    (used internally: base class of matrix_reference_t)
    \ingroup vc_blas
    \tparam T `float` or `double`
    \tparam M matrix type, e.g., VC::math::blas::ge_mat
 */
template <typename T,typename M>
struct matrix_reference_base_t {

  typedef T value_t;
  typedef M matrix_t;

  matrix_reference_base_t(const M& _m,T* _data) : m_m(_m), m_data(_data) {}

  const M& matrix() const { return m_m; }
  T* data() const { return m_data; }
  operator T*() const { return m_data; }

protected:
  M    m_m;
  T*   m_data;
};


/** \class matrix_const_reference_mixin
    (used internally: mixin to matrix_const_reference_t and matrix_reference_t)
    \ingroup vc_blas
    Mixin common ("const") functionality.
 */
template <typename REFBASE>
struct matrix_const_reference_mixin
  : public REFBASE {

  typedef typename REFBASE::value_t  value_t;
  typedef typename REFBASE::matrix_t matrix_t;
  typedef matrix_const_reference_t<value_t,matrix_t> const_reference_t;

  matrix_const_reference_mixin(const matrix_t& _m,const value_t* _data)
    : REFBASE(_m,_data) {}
  matrix_const_reference_mixin(const matrix_t& _m,value_t* _data)
    : REFBASE(_m,_data) {}

  /// BLAS type (GE_Mat,...)
  int matrix_type() const { return this->m_m.matrix_type(); }

  /** @name matrix properties (as for mat_prop_base)
      @{
  */

  /// number of rows (trans() flag is \a not considered)
  int m() const { return this->m_m.m(); }
  /// number of columns (trans() flag is \a not considered)
  int n() const { return this->m_m.n(); }
  /// upper/lower flag (triangular/symmetric matrices TR,TB,TP,SY,SB,SP only)
  UpperLowerFlag uplo() const { return this->m_m.uplo(); }
  /// transpose flag
  TransposeFlag trans() const { return this->m_m.trans(); }
  /// diagonal flag (triangular TR,TB,TP matrices only)
  DiagonalFlag unit_diag() const { return this->m_m.diag(); }
  /// leading dimension (\a not for packed storage SP,TP)
  int ld() const { return this->m_m.ld(); }
  /// number of sub-diagonals (\a only for general banded GB matrices)
  int kl() const { return this->m_m.kl(); }
  /// number of super-diagonals (\a only for general banded GB matrices)
  int ku() const { return this->m_m.ku(); }
  /// number of sub-/super-diagonals (\a only for symetric banded SB matrices)
  int k() const { return this->m_m.k(); }

  /// check if any dimension equals 1
  bool is_vector() const { return this->m_m.m()==1 || this->m_m.n()==1; }

  /// @}


  /** @name derived  from attributes
      @{
  */

  /// number of rows (respects trans() flag)
  int m_rows() const { return this->trans()==NoT ? this->m() : this->n(); }
  /// number of rows (respects trans() flag)
  int n_cols() const { return this->trans()==NoT ? this->n() : this->m(); }

  /// size required to store this matrix with current settings (e.g., ld()*n())
  int size() const { return this->m_m.size(); }

  /// check if m_rows()==1
  bool is_row_vector() const { return this->m_rows()==1; }
  /// check if m_rows()==1
  bool is_column_vector() const { return this->n_cols()==1; }

  /// @}


  /** @name read access
      @{
   */

  /// get const reference to this matrix
  const_reference_t const_ref() const {
    return const_reference_t(this->m_m,this->m_data);
  }

  /// return element (_i,_j)
  value_t operator()(int _i,int _j) const { return at(const_ref(),_i,_j); }
  /// return element (_i,_j)
  value_t operator()(int _i,const idx_end& _j) const {
    return at(const_ref(),_i,this->n_cols()-1-_j.offset);
  }
  /// return element (_i,_j)
  value_t operator()(const idx_end& _i,int _j) const {
    return at(const_ref(),this->m_rows()-1-_i.offset,_j);
  }
  /// return element (_i,_j)
  value_t operator()(const idx_end& _i,const idx_end& _j) const {
    return at(const_ref(),this->m_rows()-1-_i.offset,this->n_cols()-1-_j.offset);
  }
  
  /** get/copy column as vector_reference_t
      \arg same as `column(*this,_j)` for GE matrices (see column())
      \arg Evaluates to cp_column() for any other matrices
   */
  typename CP_Column<value_t,matrix_t>::op_t column(int _j) const {
    return CP_Column<value_t,matrix_t>::get(this->const_ref(),_j);
  }
 
  /** get/copy row as vector_reference_t
      \arg same as `row(*this,_j)` for GE matrices (see row())
      \arg Evaluates to cp_row() for any other matrices
   */
  typename CP_Row<value_t,matrix_t>::op_t row(int _i) const {
    return CP_Row<value_t,matrix_t>::get(this->const_ref(),_i);
  }

  /** get/copy diagonal as vector_reference_t
      \arg same as `diag(*this)` for GE matrices (see diag())
      \arg Evaluates to cp_diag() for any other matrices
   */
  typename CP_Diag<value_t,matrix_t>::op_t diag() const {
    return CP_Diag<value_t,matrix_t>::get(this->const_ref());
  }

  /** support `A=B(i1:i2,j1:j2)`
      \arg Same as `block(*this,_i1,_i2,_j1,_j2)` for **GE** matrices (see block()).
      \arg Evaluates to cp_block() for any other matrices.
   */
  typename CP_Block<value_t,matrix_t>::op_t
  operator()(int _i1,int _i2,int _j1,int _j2) const {
    return CP_Block<value_t,matrix_t>::get(this->const_ref(),_i1,_i2,_j1,_j2);
  }

  /** Support `A=B(i1:i2,j1:j2)` using ranges.

      \arg Same as `block(A,_ri,_rj)` for **GE** matrices (see block()).
      \arg Evaluates to cp_block() for any other matrices.

      The following syntax is used.

      C++                                  | Matlab
      -------------------------------------|--------------
      `B(_(i1,i2),_(j1,j2))`               | `B(i1:i2,j1:j2)`
      `B($,_(j1,j2))` or `B(_(),_(j1,j2))` | `B(:,j1:j2)`
      `B($,_(j1,end))`                     | `B(:,j1:end)`
      `B($,_(j1,end-k))`                   | `B(:,j1:end-k)`
      `B(_(k),$)`                          | `B(k,:)`

      Note that `"_"`, `"$"`, and `"end"` are symbols in the namespace
      VC::math::blas, i.e., VC::math::blas::_ , VC::math::blas::$ ,
      VC::math::blas::end .
  */
  template <bool END1,bool END2,bool END3,bool END4>
  typename CP_Block<value_t,matrix_t>::op_t 
  operator()
  (const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) const {
    return CP_Block<value_t,matrix_t>::get(this->const_ref(),_ri,_rj);
  }

  /// @}

};

//-----------------------------------------------------------------------------

/** \class matrix_const_reference_t
    Reference to a constant matrix.
    \ingroup vc_blas

    Use VC::math::blas::trans() to get a reference to the transposed.

    \sa matrix_reference_t, trans()
 */
template <typename T,typename M>
struct matrix_const_reference_t :
    public matrix_const_reference_mixin<matrix_const_reference_base_t<T,M> > {

  typedef matrix_const_reference_t<T,M> self_t;
  typedef T value_t;
  typedef M matrix_t;

  matrix_const_reference_t(const M& _m,const T* _data)
    : matrix_const_reference_mixin<matrix_const_reference_base_t<T,M> >(_m,_data) {}
  matrix_const_reference_t(const matrix_reference_t<T,M>& _m)
    : matrix_const_reference_mixin<matrix_const_reference_base_t<T,M> >(_m.matrix(),_m.data()) {}

  /** @name set properties
      @{
  */

  /// set dimension (only if M=VarInt), calls set_ld(_m) if ld()<_m
  self_t& set_m(int _m) {
    this->m_m.set_m(_m);
    if (this->m_m.ld()<_m)
      this->m_m.set_ld(_m);
    return *this;
  }
  /// set dimension (only if N=VarInt)
  self_t& set_n(int _n) { this->m_m.set_n(_n); return *this; }
  /// set dimension
  self_t& set_size(int _m,int _n) { return set_m(_m).set_n(_n); }
  /// set dimension
  self_t& set_size(int _m,int _n,int _ld) { return set_m(_m).set_n(_n).set_ld(_ld); }
  /// set leading dimension (only if LD=VarInt)
  self_t& set_ld(int _ld) { assert(_ld>=this->m()); this->m_m.set_ld(_ld); return *this; }
  /// set data pointer
  self_t& set_data(const value_t* _data) { this->m_data=_data; return *this; }

  /// @}
};


/** \class matrix_reference_t
    Reference to a matrix.
    \ingroup vc_blas

    \arg All operations that yield a matrix as result modify \c *this
    matrix!

    \arg Some operations are defined only for specific matrices
    (triangular, non-triangular, see mm(), sm())!

    \arg Operands must \a not refer to destination data
    (this->data())! (With exception of element-wise operations, e.g., madd()).
    For example, `A=A+B` is invalid (run-time assertion fails), use
    `A+=B` instead (in this case, `A=B+A` is fine, too.)

    The following operators are defined

    \code
    value_t                    alpha;
    vector_const_reference_t<> x,y;
    matrix_reference_t<>       A;
    matrix_const_reference_t<> B,C;

    A(i,j);                    // access element (VC::math::blas::at(int,int)
                               //  (r/w accessed depends on matrix type)

    A=B(_(i1,i2),_(j1,j2))              // B(i1:i2,j1:j2) -> cp_block()
    A=B($,_(j1,j2)) == B(_(),_(j1,j2))  // B(:,j1:j2)
    A=B($,_(j1,end))                    // B(:,j1:end)
    A=B($,_(j1,end-k))                  // B(:,j1:end-k)

    block(A,i1,i2,j1,j2)                // submatrix as matrix_reference_t (GE only)
                                        //  same indexing with ranges _(i1,i2)
                                        //  as above applies

    A=B;                       // GE matrix A (cp())
    A=trans(B);
    trans(A)=B;

    A=alpha*B;

    A=B(i1,i2,j1,j2);          // GE matrix A=B(i1:i2,j1:j2) (cp_block())
    A=trans(B(i1,i2,j1,j2));
    trans(A)=B(i1,i2,j1,j2);
    A=block(B,i1,i2,j1,j2);    // GE matrix B, GE reference (VarInt) A (no copy)

    A=zeros;                   // yields zero matrix (ld_zero())
    A=all(alpha);              // set all element to _alpha (ld_all(), GE,SY,SP only)
    A=ones;                    // same as A=all(1)
    A=eye;                     // load identity matrix (ld_eye())
    A=eye*alpha;               // load identity matrix (ld_eye(alpha))
    A=random_values;           // load random values in [0,1] (ld_rand())


    A=diag(x);                 // load diagonal matrix defined by vector x (ld_diag())
                               //  (not if DiagonalFlag==UnitDiag)
    A+=diag(x);                // only for GE,GB matrices where diag(A) is defined
    A-=diag(x);                //  (same as above)

    A*=alpha;                  // mscal()
    A/=alpha;                  // mscal()
    A+=alpha;                  // add constant to all elements (adds(), only GE,SY,SP)
    A-=alpha;                  // A+=(-alpha)

    A+=B*alpha;                // add matrices (madd(), not TR,TB,SP)
    A-=B;

    A.add_trans();             // A+=A' (quadratic GE matrix only)
    A.add_trans(alpha);        // A=(A+A')*alpha (quadratic GE matrix only)

                               // Note: A+=trans(A) fails at run-time!

    A*=B;                      // ELEMENTWISE multiplication (emul(), GE only, NoTranspose)
    A/=B;                      // ELEMENTWISE division (ediv(), GE only, NoTranspose)
                               //  Note that neither mm() nor sm() for TR matrices
                               //  are mapped to any operator!

    A*=diag(x);                // right-multiply by diagonal matrix: scale columns
                               //  (scal_cols(), not SY,SB,SP)
    trans(A)*=diag(x);         // left-multiply by diagonal matrix: scale rows
                               //  (scal_cols(), not GB,TB,SY,SB,SP)

    A=B*C;                     // matrix-matrix multiplication (mm())
    A=B*C*alpha;               //  (only GE,SY)
    A+=B*C;
    A-=B*C;                    //  Note that there is a special form of mm() for TR
    A+=B*C*alpha;              //  matrices. This one is NOT overloaded to A*=B (which is
    A-=B*C*alpha;              //  ELEMENTWISE multiplication).

    A=outer_prod(x);           // rank-1 update using outer product xx^T (r())
    A=outer_prod(x)*alpha;     //  (only SY,SP)
    A+=outer_prod(x);
    A+=outer_prod(x)*alpha;
    A-=outer_prod(x);
    A-=outer_prod(x)*alpha;

    A=ata(B)*alpha;            // rank-k update A=alpha*B'*B
    A=ata(B);                  //   (only SY,SP)
    A=ata(B)*alpha;
    A+=ata(B)*alpha;
    A-=ata(B)*alpha;

    A=ata(trans(B))*alpha;     // A=B*B'*alpha
    // etc.                    //   (only SY,SP)

    swap(A,B);                 // swap data of A,B (requires same matrix types)

    map_f(f,A);                // element-wise A(i,j)=f(A(i,j))
                               //   available for all matrices: processes only
                               //   elements that are stored explicitly (no implicit
                               //   zeros, ones, e.g., TR)
    map_f2(f,A,B);             // element-wise A=f(B), see above; requires same matrix
                               //  types
    map_f12(f,A,B);            // element-wise A=f(A,B)
    map_f23(f,A,B,C);          // element-wise A=f(B,C)

    std::cout << A
              << pp(A)         // pretty print
              << matlab(A)     // print as Matlab vector
              << pdata(A);     // print data elements

    write_mat4(cout,A);        // write into MAT-file (V4)

    \endcode

    \sa matrix_const_reference_t, trans(), mat
*/
template <typename T,typename M>
struct matrix_reference_t
  : matrix_const_reference_mixin<matrix_reference_base_t<T,M> > {
  typedef matrix_reference_t<T,M> self_t;
  typedef T value_t;
  typedef M matrix_t;

  matrix_reference_t(const M& _m,T* _data)
    : matrix_const_reference_mixin<matrix_reference_base_t<T,M> >(_m,_data) {}

  /// get reference to this matrix
  self_t ref() const { return *this; }

  /** @name set properties
      @{
  */

  /// set dimension (only if M=VarInt), calls set_ld(_m) if ld()<_m
  self_t& set_m(int _m) {
    this->m_m.set_m(_m);
    if (this->m_m.ld()<_m)
      this->m_m.set_ld(_m);
    return *this;
  }
  /// set dimension (only if N=VarInt)
  self_t& set_n(int _n) { this->m_m.set_n(_n); return *this; }
  /// set dimension
  self_t& set_size(int _m,int _n) { return set_m(_m).set_n(_n); }
  /// set dimension
  self_t& set_size(int _m,int _n,int _ld) { return set_m(_m).set_n(_n).set_ld(_ld); }
  /// set leading dimension (only if LD=VarInt)
  self_t& set_ld(int _ld) { assert(_ld>=this->m()); this->m_m.set_ld(_ld); return *this; }
  /// set data pointer
  self_t& set_data(const value_t* _data) { this->m_data=_data; return *this; }

  /// @}

  /** @name BLAS functions.
      \arg Method names map to respective BLAS functions.
      \arg The destination of operations (\f$lhd \leftarrow\f$) is always \c *this!
      @{
   */

  /** Matrix-matrix multiplication for \b general matrices and \b symmetric _a.
      \f[ C \leftarrow \alpha A B + \beta C \f]
      \sa VC::math::blas::mm()
   */
  template <typename A,typename B>
  const self_t& mm(const T& _alpha,
                   const matrix_const_reference_t<T,A>& _a,
                   const matrix_const_reference_t<T,B>& _b,
                   const T& _beta) const;

  /** Matrix-matrix multiplication with symmetric matrix _a.
      \f[ C \leftarrow \alpha A B + \beta C \f] for _side==Left or
      \f[ C \leftarrow \alpha B A + \beta C \f] for _side==Right
      \sa VC::math::blas::mm()
  */
  template <typename A,typename B>
  const self_t& mm(SideFlag _side,
                   const T& _alpha,
                   const matrix_const_reference_t<T,A>& _a,
                   const matrix_const_reference_t<T,B>& _b,
                   const T& _beta) const;

  /** Matrix-matrix multiplication by a symmetric matrix _a.
      \f[ B \leftarrow \alpha A B \f]
      \sa VC::math::blas::mm()
  */
  template <typename A>
  const self_t& mm(const T& _alpha,
             const matrix_const_reference_t<T,A>& _a) const;

  /** Rank-1 update for \b symmetric matrix.
      \f[ A \leftarrow \alpha x x^{T} + A \f]
      \sa VC::math::blas::r()
   */
  template <typename V>
  const self_t& r(const T& _alpha,const vector_const_reference_t<T,V>& _x) const;

  /** Rank-2 update for \b symmetric matrix.
      \f[ A \leftarrow \alpha x y^{T} + \alpha y x^{T} + A \f]
      \sa VC::math::blas::r2()
   */
  template <typename V,typename W>
  const self_t& r2(const T& _alpha,
                   const vector_const_reference_t<T,V>& _x,
                   const vector_const_reference_t<T,W>& _y) const;


  /** Rank-1 update for \b symmetric matrix.
      \f[ C \leftarrow \alpha A A^{T} + \beta C \f]
      \sa VC::math::blas::rk()
  */
  template <typename A>
  const self_t& rk(const T& _alpha,
                   const matrix_const_reference_t<T,A>& _a,const T& _beta) const;

  /** Rank-2 update for \b symmetric matrix.
      \f[ C \leftarrow \alpha A B^{T} + \alpha B A^{T} + \beta C \f]
      \sa VC::math::blas::r2k()
  */
  template <typename A,typename B>
  const self_t& r2k(const T& _alpha,
                    const matrix_const_reference_t<T,A>& _a,
                    const matrix_const_reference_t<T,B>& _b, const T& _beta) const;

  /** Solve \b triangular system.
      \f[ B \leftarrow \alpha A^{-1} B \f]
      \sa VC::math::blas::sm()
  */
  template <typename A>
  const self_t& sm(const T& _alpha,
                   const matrix_const_reference_t<T,A>& _a) const;
  /// @}


  /** @name BLAS extensions.
      \arg The destination of operations (\f$lhd \leftarrow\f$) is always \c *this!
      @{
   */

  /// \f$A\leftarrow 0\f$
  const self_t& ld_zero() const;

  /// \f$A\leftarrow 0\f$
  const self_t& ld_all(const T& _a) const;

  /// \f$A\leftarrow I\f$
  const self_t& ld_eye(const T& _a=T(1)) const;

  /// load random values uniformly distributed in in [0,1]
  const self_t& ld_rand() const;

  /// load diagonal matrix (scaled by _alpha) from _d
  template <typename V>
  const self_t& ld_diag(const vector_const_reference_t<T,V>& _d,
                        const T& _alpha=1.0) const;
  /// load diagonal matrix  (scaled by _alpha) from _d
  template <typename V>
  const self_t& ld_diag(const vector_reference_t<T,V>& _d,
                        const T& _alpha=1.0) const {
    return ld_diag(_d.const_ref(),_alpha);
  }

  /// \f$A\leftarrow \alpha A\f$
  const self_t& mscal(const T& _alpha) const;
  /// add _alpha to all elements
  const self_t& adds(const T& _alpha) const;

  /// element-wise multiplication this *= _b
  template <typename B>
  const self_t& emul(const matrix_const_reference_t<T,B>& _b,const T& _alpha=T(1)) const;
  /// element-wise division this/= _b
  template <typename B>
  const self_t& ediv(const matrix_const_reference_t<T,B>& _b,const T& _alpha=T(1)) const;

  /// matrix addition \f$A\leftarrow \alpha B + A\f$
  template <typename B>
  const self_t& madd(const T& _alpha,const matrix_const_reference_t<T,B>& _b) const;
  /// Same as madd()
  template <typename B>
  const self_t& add(const matrix_const_reference_t<T,B>& _b,const T& _alpha=T(1)) const {
    return madd(_alpha,_b);
  }
  /// matrix subtraction based (madd())
  template <typename B>
  const self_t& sub(const matrix_const_reference_t<T,B>& _b,const T& _alpha=T(1)) const {
    return madd(-_alpha,_b);
  }

  /// right multiplication by diagonal matrix
  template <typename V>
  const self_t& scal_cols(const vector_const_reference_t<T,V>& _d,const T& _alpha=T(1)) const;

  /// add transposed and scale \f$A\leftarrow \alpha (A + A')\f$ (quadratic GE matrix)
  const self_t& add_trans(const T& _alpha=T(1)) const;

  /// @}

  /// access element (_i,_j)
  T& operator()(int _i,int _j) const {  // const-ness is inconsistent
    return at(*this,_i,_j);
  }
  /// access element (_i,_j)
  T& operator()(int _i,const idx_end& _j) const {
    return at(*this,_i,this->n_cols()-1-_j.offset);
  }
  /// access element (_i,_j)
  T& operator()(const idx_end& _i,int _j) const {
    return at(*this,this->m_rows()-1-_i.offset,_j);
  }
  /// access element (_i,_j)
  T& operator()(const idx_end& _i,const idx_end& _j) const {
    return at(*this,this->m_rows()-1-_i.offset,this->n_cols()-1-_j.offset);
  }


  /// assign matrix (cp()) -- never assign reference!
  const self_t& operator=(const matrix_const_reference_t<T,matrix_t>& _a) const {
    cp(_a,*this);
    return *this;
  }
  /// assign matrix (cp()) -- never assign reference!
  const self_t& operator=(const matrix_reference_t<T,matrix_t>& _a) const {
    cp(_a.const_ref(),*this);
    return *this;
  }
  /// assign matrix (cp())
  template <typename A>
  const self_t& operator=(const matrix_const_reference_t<T,A>& _a) const {
    cp(_a,*this);
    return *this;
  }

  /// assign matrix (cp())
  template <typename A>
  const self_t& operator=(const matrix_reference_t<T,A>& _a) const {
    cp(_a.const_ref(),*this);
    return *this;
  }

  /// assign A=_data (cp(), no checks on _data bounds)
  template <typename S>
  const self_t& operator=(const S* _data) const {
    cp(_data,*this);
    return *this;
  }

  /** @name evaluate expressions
      The Expression classes define functionality for operator=() and
      operator+=().
      @{
   */

  template <typename EXPR,typename ARG1,typename ARG2>
  const self_t& operator=(const Expression<Type_matrix,EXPR,T,ARG1,ARG2>& _expr) const {
    _expr.assign(*this);
    return *this;
  }
  template <typename EXPR,typename ARG1,typename ARG2>
  const self_t& operator+=(const Expression<Type_matrix,EXPR,T,ARG1,ARG2>& _expr) const {
    _expr.plus_assign(*this);
    return *this;
  }
  template <typename EXPR,typename ARG1,typename ARG2,typename ARG3>
  const self_t& operator=(const Expression<Type_initializer,EXPR,ARG1,ARG2,ARG3>& _expr) const {
    _expr.assign(*this);
    return *this;
  }
  template <typename EXPR,typename ARG1,typename ARG2,typename ARG3>
  const self_t& operator+=(const Expression<Type_initializer,EXPR,ARG1,ARG2,ARG3>& _expr) const {
    _expr.plus_assign(*this);
    return *this;
  }

  template <typename EXPR>
  const self_t& operator-=(const EXPR& _expr) const {
    return *this+=(-_expr); // DELEGATE ALL
  }

   /// right-multiply by diagonal matrix (scal_cols())
  template <typename V>
  const self_t& operator*=(const Expression<Type_matrix,Expr_diag_mat,T,V>& _d) const {
    return this->scal_cols(_d.x,_d.a);
  }

  /// @}

  /** @name Assignment operators for \b elementwise mutliplication and division.
      @{
  */

  /// multiply by scalar (mscal())
  const self_t& operator*=(const T& _a) const {
    mscal(_a);
    return *this;
  }
  /// divide by scalar (mscal())
  const self_t& operator/=(const T& _a) const {
    mscal(T(1)/T(_a));
    return *this;
  }


  /// \b elementwise multiplication (emul())
  template <typename B>
  const self_t& operator*=(const matrix_const_reference_t<T,B>& _b) const {
    emul(_b,T(1));
    return *this;
  }
  /// \b elementwise multiplication (emul())
  template <typename B>
  const self_t& operator*=(const matrix_reference_t<T,B>& _b) const {
    emul(_b.const_ref(),T(1));
    return *this;
  }
  /// \b elementwise multiplication (emul())
  template <typename MA>
  const self_t& operator*=(const Expression<Type_matrix,Expr_aA,T,MA>& _expr) const {
    emul(_expr.A,_expr.a);
    return *this;
  }
  /// \b elementwise division (ediv())
  template <typename B>
  const self_t& operator/=(const matrix_const_reference_t<T,B>& _b) const {
    ediv(_b,T(1));
    return *this;
  }
  /// \b elementwise division (ediv())
  template <typename B>
  const self_t& operator/=(const matrix_reference_t<T,B>& _b) const {
    ediv(_b.const_ref(),T(1));
    return *this;
  }
  /// \b elementwise division (ediv())
  template <typename MA>
  const self_t& operator/=(const Expression<Type_matrix,Expr_aA,T,MA>& _expr) const {
    ediv(_expr.A,_expr.a);
    return *this;
  }

  /// @}

  /// matrix addition (madd())
  template <typename B>
  const self_t& operator+=(const matrix_const_reference_t<T,B>& _b) const {
    madd(T(1),_b);
    return *this;
  }
  /// matrix addition (madd())
  template <typename B>
  const self_t& operator+=(const matrix_reference_t<T,B>& _b) const {
    madd(T(1),_b.const_ref());
    return *this;
  }
  /// matrix subtraction (madd())
  template <typename B>
  const self_t& operator-=(const matrix_const_reference_t<T,B>& _b) const {
    madd(T(-1),_b);
    return *this;
  }
  /// matrix subtraction (madd())
  template <typename B>
  const self_t& operator-=(const matrix_reference_t<T,B>& _b) const {
    madd(T(-1),_b.const_ref());
    return *this;
  }

  /** get/copy column as vector_reference_t
      \arg same as `column(*this,_j)` for GE matrices (see column())
      \arg Evaluates to cp_column() for any other matrices
   */
  typename CP_Column<value_t,matrix_t>::op_t column(int _j) const {
    return CP_Column<value_t,matrix_t>::get(*this,_j);
  }
 
  /** get/copy row as vector_reference_t
      \arg same as `row(*this,_j)` for GE matrices (see row())
      \arg Evaluates to cp_row() for any other matrices
   */
  typename CP_Row<value_t,matrix_t>::op_t row(int _i) const {
    return CP_Row<value_t,matrix_t>::get(*this,_i);
  }

  /** get/copy diagonal as vector_reference_t
      \arg same as `diag(*this)` for GE matrices (see diag())
      \arg Evaluates to cp_diag() for any other matrices
   */
  typename CP_Diag<value_t,matrix_t>::op_t diag() const {
    return CP_Diag<value_t,matrix_t>::get(*this);
  }

  /** support `A=B(i1:i2,j1:j2)`
      \arg Same as `block(A,_i1,_i2,_j1,_j2)` for **GE** matrices (see block()).
      \arg Evaluates to cp_block() for any other matrices.
   */
  typename CP_Block<value_t,matrix_t>::op_t
  operator()(int _i1,int _i2,int _j1,int _j2) const {
    return CP_Block<value_t,matrix_t>::get(this->ref(),_i1,_i2,_j1,_j2);
  }

  /** Support `A=B(i1:i2,j1:j2)` using ranges.

      \arg Same as `block(A,_ri,_rj)` for **GE** matrices (see block()).
      \arg Evaluates to cp_block() for any other matrices.

      The following syntax is used.

      C++                                  | Matlab
      -------------------------------------|--------------
      `B(_(i1,i2),_(j1,j2))`               | `B(i1:i2,j1:j2)`
      `B($,_(j1,j2))` or `B(_(),_(j1,j2))` | `B(:,j1:j2)`
      `B($,_(j1,end))`                     | `B(:,j1:end)`
      `B($,_(j1,end-k))`                   | `B(:,j1:end-k)`
      `B(_(k),$)`                          | `B(k,:)`

      Note that `"_"`, `"$"`, and `"end"` are symbols in the namespace
      VC::math::blas, i.e., VC::math::blas::_ , VC::math::blas::$ ,
      VC::math::blas::end .
  */
  template <bool END1,bool END2,bool END3,bool END4>
  typename CP_Block<value_t,matrix_t>::op_t 
  operator()
  (const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) const {
    return CP_Block<value_t,matrix_t>::get(this->ref(),_ri,_rj);
  }  

  // MISSING: *= TR_only  xxx
};

//-----------------------------------------------------------------------------

/// get iterator on vector _x
template <typename T,typename V>
vec_iterator<T>
iterate(const vector_reference_t<T,V>& _x)  {
  return vec_iterator<T>(_x.data(),_x.n(),_x.inc());
}

/// get row iterator on ge_mat
template <typename T,int M,int N,int LD>
GE_row_iterator<T>
iterate_rows(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _A) {
   return GE_row_iterator<T>(_A.data(),_A.m(),_A.n(),_A.ld());
}

/// get column iterator on ge_mat
template <typename T,int M,int N,int LD>
GE_column_iterator<T>
iterate_columns(const matrix_reference_t<T,ge_mat<NoT,M,N,LD> >& _A) {
   return GE_column_iterator<T>(_A.data(),_A.m(),_A.n(),_A.ld());
}

/// get row iterator on ge_mat
template <typename T,int M,int N,int LD>
GE_column_iterator<T>
iterate_rows(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _A) {
   return GE_column_iterator<T>(_A.data(),_A.m(),_A.n(),_A.ld());
}

/// get column iterator on ge_mat
template <typename T,int M,int N,int LD>
GE_row_iterator<T>
iterate_columns(const matrix_reference_t<T,ge_mat<Transpose,M,N,LD> >& _A) {
   return GE_row_iterator<T>(_A.data(),_A.m(),_A.n(),_A.ld());
}

//-----------------------------------------------------------------------------

} // namespace blas
} // namespace math
} // namespace VC

// doing the namespace thing for DOXYGEN

# include "blas_matrix_mat_inc.hh"
# include "blas_matrix_glue_inc.hh"
# include "blas_matrix_methods_inc.hh"
# include "blas_matrix_trans_inc.hh"
# include "blas_matrix_storage_inc.hh"
# include "blas_matrix_algorithms_inc.hh"

# include "blas_matrix_doc_inc.hh"

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

namespace VC {
namespace math {
namespace blas {
//=============================================================================

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

#endif // VC_MATH_BLAS_MATRIX_HH
