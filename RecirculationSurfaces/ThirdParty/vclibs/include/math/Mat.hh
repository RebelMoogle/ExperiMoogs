//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef __VC_MATH_MAT_HH__
#define __VC_MATH_MAT_HH__

#include <utility> // std::move

#include "VecN.hh"
#include "lu.hh"
#include "blas_matrix.hh"
#include "Mat2x2.hh" // eig2x2
#include "roots.hh" // real_roots_3

#include "../base/debug.hh"

namespace VC {
namespace math {


//-----------------------------------------------------------------------------

/** \file Mat.hh
    small matrices, [\ref vc_math_lav]
*/

/** \defgroup vc_math_lam Linear algebra: fixed-size (small) matrices
    \ingroup vc_math

    Data types, arithmetic and utility functions for matrix-vector and
    matrix-matrix operations. Targeted at convenient use of **small**
    matrices, see [here](@ref vc_math_mat). For efficient processing
    of larger (or special types of) matrices see [\ref vc_blas].

    \arg Mat<T,M,N> defines an MxN matrix data type.
    \arg The implementation builds on VecN.
    \arg There are few specializations for small quadratic matrices
    (see Mat_op<T,N,N>). In general, we rely on the compiler
    optimization. Some functions are implemented using [\ref
    vc_vcblas].
    \arg blas::matrix_reference_t defines an interface to using BLAS
    operations on matrices, see Mat::ref().
    \arg ex_trans_t<T,M,N> and ex_diag_t<T,M,N> define expressions
    referring to a Mat<T,M,N> instance, which in most sitautions can
    be treated like matrix instances. They are typically not evaluated
    directly (instead allow use of optimized operators).

    \sa VecN, [\ref vc_vcblas], [\ref vc_math_mat]
*/

template <typename T,int M,int N> class Mat;

template <typename T,int M,int K,int N>
Mat<T,M,K+N> h_cat(const Mat<T,M,K>& _a,const Mat<T,M,N>& _b);
template <typename T,int M,int N>
Mat<T,M,1+N> h_cat(const VecN<T,unsigned(M)>& _a,const Mat<T,M,N>& _b);
template <typename T,int M,int N>
Mat<T,M,N+1> h_cat(const Mat<T,M,N>& _b,const VecN<T,unsigned(M)>& _a);
template <typename T,unsigned M>
Mat<T,M,2> h_cat(const VecN<T,M>& _a,const VecN<T,M>& _b);

template <typename T,int M,int K,int N>
Mat<T,M+K,N> v_cat(const Mat<T,M,K>& _a,const Mat<T,M,N>& _b);
template <typename T,int M,int N>
Mat<T,1+M,N> v_cat(const VecN<T,unsigned(M)>& _a,const Mat<T,M,N>& _b);
template <typename T,int M,int N>
Mat<T,M+1,N> v_cat(const Mat<T,M,N>& _b,const VecN<T,unsigned(M)>& _a);
template <typename T,unsigned N>
Mat<T,2,N> v_cat(const VecN<T,N>& _a,const VecN<T,N>& _b);

template <int R,int C> struct repmat;

template <typename T,int M,int N> struct inverse_t;
template <typename T,int M,int N> struct lu_factor_t;
template <typename T,int M,int N> struct det_t;
template <typename T,int M,int N,bool FORCE_LAPACK=false> struct sy_eig_t;
template <typename T,int M,int N,bool FORCE_LAPACK=false> struct eig_t;

//-----------------------------------------------------------------------------

/** expression for transposed matrix `A'` (used by Mat)
    \ingroup vc_math_lam
    \see VC::math::Mat<T,M,N>::trans()
*/
template <typename T,int M,int N>
struct ex_trans_t {
  ex_trans_t(const Mat<T,M,N>& _mat) : matrix(_mat) {}
  const Mat<T,M,N>& matrix;
  /// evaluate transpose of `matrix`
  Mat<T,N,M> operator*() const {
    Mat<T,N,M> c;
    for (int j=0;j<N;++j)
      for (int i=0;i<M;++i)
        c.data()[j+i*N]=matrix.data()[i+j*M];
    return c;
  }
  /// evaluate matrix-vector product `A'*x`
  typename Mat<T,N,M>::col_t
  operator*(const typename Mat<T,N,M>::row_t& _x) const {
    typename Mat<T,N,M>::col_t y;
    for (int j=0;j<N;++j)
      y[j]=(_x|matrix.col(j));
    return y;
  }
  /// transpose
  const Mat<T,M,N>& trans() const { return matrix; }
  /// trace
  T trace() const { return matrix.trace(); }
private:
  ex_trans_t() {}
};

//-----------------------------------------------------------------------------

/** expression for diagonal matrix `diag(A)` (used by Mat)
    \ingroup vc_math_lam
    \see VC::math::diag()
*/
template <typename T,int M,int N>
struct ex_diag_t {
  typedef typename Mat<T,M,N>::diag_t diag_t;
  ex_diag_t(const diag_t& _d) : diagonal(_d) {
    static_assert(M==N,"available only for square matrices");
  }
  const diag_t diagonal;
  /// evaluate diagonal matrix
  Mat<T,M,N> operator*() const {
    Mat<T,M,N> a; a.load_diag_matrix(diagonal); return a;
  }
  /// transpose (identity)
  const ex_diag_t& trans() const { return *this; }
  /// trace (sum)
  T trace() const { return diagonal.sum(); }
};

//-----------------------------------------------------------------------------

/// argument to CTOR of VC::math::Mat \ingroup vc_math_lam
enum RowMajor_t { RowMajor };

/** \class VC::math::Mat VecN.hh
    \brief Small matrices.
    \ingroup vc_math_lam

    - Operators `A*=B` and `A/=B` work **element-wise**!

    - Use elts() for missing **element-wise** operations (relies on
       VecN).

    - The **transpose matrix** trans() yields an **expression** that
    "acts as" transposed. For instance, `C=A'\*B` (or `C=A.trans()*B`)
    does not evaluate `A'` but uses a different (more efficient)
    multiplication. Evaluation of the transposed matrix can be
    enforced by ex_trans_t<T,M,N>::operator*().

    - The use of Mat::trans() is similar to (but simpler and more
    restrictive than) VC::math::blas::trans().

    - Similarly, there are expressions to have a vector "act as" as
    **diagonal matrix**, see VC::math::diag(). For instance,
    `A=diag(x)*A` scales the columns of `A` (without constructing a
    matrix `diag(x)`). Evaluation of a diagonal matrix can be enforced
    by `ex_diag_t<T,M,N>::operator*()`.

    - Note that the inherited methods diag() and set_diag() are only
    used to extract or manipulate the matrix diagonal.

    - Support for diagonal matrices is currently restricted to *square
    matrices*.


    \sa [\ref vc_math_lam_mops], [\ref vc_math_lam]
*/

# define m_raw m_elts.data() // mscv13

template <typename T,int M,int N>
class Mat {

  // union { // mscv13
  //   VecN<T,M>   m_cols[N];  //!< matrix as column vectors
    VecN<T,M*N> m_elts;     //!< matrix as stacked columns
  //   T           m_raw[M*N]; //!M matrix as C-array
  // };      // mscv13

public:

  typedef VecN<T,M> col_t;  //!< column vector type
  typedef VecN<T,N> row_t;  //!< row vector type
  typedef VecN<T,M> diag_t; //!< diagonal vector type

  typedef T value_type;            //!< scalar type
  typedef Mat<T,M,N> self_t;       //!< matrix type

  enum { ROWS=M, //!< number of rows (n_rows())
         COLS=N  //!< number of rows (n_cols())
  };

  /** @name matrix and vector reference types
      @{
  */

  /// general matrix type
  typedef blas::ge_mat<blas::NoT,M,N,M> matrix_t;
  /// matrix reference returned by const_ref()
  typedef blas::matrix_const_reference_t<T,matrix_t> const_ref_t;
  /// matrix reference returned by ref()
  typedef blas::matrix_reference_t<T,matrix_t> ref_t;\

  /// column vector type
  typedef blas::vec<N,1> col_vec_t;
  /// column reference returned by const_col_ref()
  typedef blas::vector_const_reference_t<T,col_vec_t> const_col_ref_t;
  /// column reference returned by col_ref()
  typedef blas::vector_reference_t<T,col_vec_t> col_ref_t;

  /// row vector type
  typedef blas::vec<M,N> row_vec_t;
  /// row reference returned by const_row_ref()
  typedef blas::vector_const_reference_t<T,row_vec_t> const_row_ref_t;
  /// row reference returned by row_ref()
  typedef blas::vector_reference_t<T,row_vec_t> row_ref_t;

  /// @}


  /// create matrix w/o any initialization
  Mat() {}
  /// create matrix from other matrix `_a`
  template <typename S> Mat(const Mat<S,M,N>& _a) : m_elts(_a.m_elts) {}
  /// copy CTOR
  Mat(const Mat<T,M,N>& _a) : m_elts(_a.m_elts) {}
  /// create constant matrix with all elements equal to `_x`
  Mat(T _x) { load_all(_x); }
  /// create matrix and load elements from **column major** array `_x`
  Mat(const T* _x) { for (int i=0;i<M*N;++i) this->m_raw[i]=_x[i]; }
  /// create matrix and load elements from **column major** array `_x`
  template <typename S>
  Mat(const S* _x) { for (int i=0;i<M*N;++i) this->m_raw[i]=T(_x[i]); }
  /// create matrix and load elements from **row major** array `_x`
  Mat(const T* _x,RowMajor_t) {
    for (int j=0;j<N;++j)
      for (int i=0;i<M;++i)
        this->m_raw[i+j*M]=_x[j+i*N];
  }
  /// create matrix and load elements from **row major** array `_x`
  template <typename S>
  Mat(const S* _x,RowMajor_t) {
    for (int j=0;j<N;++j)
      for (int i=0;i<M;++i)
        this->m_raw[i+j*M]=_x[j+i*N];
  }
  /// create `Mx2` matrix from column vectors
  Mat(const VecN<T,M>& _a0,const VecN<T,M>& _a1) {
    static_assert(N==2,"dimension mismatch");
    this->col(0)=_a0; this->col(1)=_a1;
  }
  /// create `Mx3` matrix from column vectors
  Mat(const VecN<T,M>& _a0,const VecN<T,M>& _a1,const VecN<T,M>& _a2) {
    static_assert(N==3,"dimension mismatch");
    this->col(0)=_a0; this->col(1)=_a1; this->col(2)=_a2;
  }
  /// create `Mx4` matrix from column vectors
  Mat(const VecN<T,M>& _a0,const VecN<T,M>& _a1,const VecN<T,M>& _a2,
      const VecN<T,M>& _a3) {
    static_assert(N==4,"dimension mismatch");
    this->col(0)=_a0; this->col(1)=_a1; this->col(2)=_a2; this->col(3)=_a3;
  }
  /// create matrix from transposed matrix (load transposed)
  Mat(const ex_trans_t<T,N,M>& _other_trans) { *this=*_other_trans; }
  /// create matrix from diagonal matrix (load_diag_matrix())
  Mat(const ex_diag_t<T,M,N>& _diag) {load_diag_matrix(_diag); }

  /// define `A=blas::zeros` (load_zeros())
  Mat(const blas::Expression<blas::Type_initializer,blas::Expr_zeros>&) { load_zeros(); }
  /// define `A=blas::random_values` (load_rand())
  Mat(const blas::Expression<blas::Type_initializer,blas::Expr_random_values>&) { load_rand(); }
  /// define `A=blas::ones` (load_all(1))
  Mat(const blas::Expression<blas::Type_initializer,blas::Expr_ones>&) { load_all(T(1)); }
  /// define `A=blas::all(x)` (load_all(x))
  Mat(const blas::Expression<blas::Type_initializer,blas::Expr_all,T>& _all) { load_all(_all.a); }
  /// define `A=blas::eye` (load_eye())
  Mat(blas::Expression<blas::Type_initializer,blas::Expr_eye>&) { load_eye(); }

  ~Mat() {}


  /** @name dimensions
      @{
    */

    unsigned n_rows() const { return M; } //!< number of rows (`=M`)
    unsigned n_cols() const { return N; } //!< number of columns (`=N`)
    unsigned size() const { return M*N; }   //!< number of elements (`=M*N`)

  /// @}

  /** @name access and assignment
      @{
  */

  /// get raw data
  T* data() { return m_raw; }
  /// get raw data
  const T* data() const { return m_raw; }

  /// alias for raw()
  T* begin() { return m_raw; }
  /// alias for raw()
  const T* begin() const { return m_raw; }

  /// alias for raw()+size()
  T* end() { return m_raw+M*N; }
  /// alias for raw()+size()
  const T* end() const { return m_raw+M*N; }

  /// access element `A(i,j)` (zero based, respect VC_VECN_CHECK_BOUNDS)
  T& operator()(int _i,int _j) {
# ifdef VC_VECN_CHECK_BOUNDS
    assert(0<=_i && _i<M);
    assert(0<=_j && _j<N);
#endif
    return m_raw[_i+_j*M];
  }

  /// access element `A(i,j)` (zero based, respect VC_VECN_CHECK_BOUNDS)
  const T& operator()(int _i,int _j) const {
# ifdef VC_VECN_CHECK_BOUNDS
    assert(0<=_i && _i<M);
    assert(0<=_j && _j<N);
#endif
    return m_raw[_i+_j*M];
  }

  /// cast **column-major** C-array (alias VC::math::to_m<M,N>())
  static self_t&
  to_m(T* _p) { return (self_t&) *_p; }

  /// cast **column-major** `const` C-array (alias VC::math::to_cm<M,N>())
  static const self_t&
  to_cm(const T* _p) { return *((const self_t*) _p); }

  /// external definition of Mat<T,M,N>::to_m()
  template <int MM,int NN,typename U>
  Mat<U,MM,NN>& to_m(U* _v); // for documentation only
  /// external definition of Mat<T,M,N>::to_cm()
  template <int MM,int NN,typename U>
  const Mat<U,MM,NN>& to_cm(const U* _v); // for documentation only


  /// access all elements as (stacked column) vector
  VecN<T,M*N>& elts() { return m_elts; }
  /// access all elements as (stacked column) vector
  const VecN<T,M*N>& elts() const { return m_elts; }

  /// access column `_j`
  // col_t& col(int _j) { assert(0<=_j && _j<N); return m_cols[_j]; } // mscv13
  col_t& col(int _j) {
    assert(0<=_j && _j<N);
    col_t* cols=(col_t*) m_raw;
    return cols[_j];
  }
  /// access column `_j`
  // const col_t& col(int _j) const { assert(0<=_j && _j<N); return m_cols[_j]; } // mscv13
  const col_t& col(int _j) const {
    assert(0<=_j && _j<N);
    const col_t* cols=(col_t*) m_raw;
    return cols[_j];
  }

  /// get (copy of) row `_i`
  row_t row(int _i) const {
    assert(0<=_i && _i<M);
    row_t r;
    for (int j=0;j<N;++j) r[j]=m_raw[_i+j*M];
    return r;
  }
  /// set row `_i`
  void set_row(int _i,const row_t& _r) {
    assert(0<=_i && _i<M);
    for (int j=0;j<N;++j) m_raw[_i+j*M]=_r[j];
  }

  /// get row reference
  //  blas::

  /// swap columns `_i` and `_j)
  self_t& swap_cols(int _i,int _j) {
    if (_i!=_j)
      std::swap(col(_i),col(_j));
    return *this;
  }
  /// swap rows `_i` and `_j)
  self_t& swap_rows(int _i,int _j) {
    if (_i!=_j) {
      row_t ri=row(_i), rj=row(_j);
      set_row(_i,rj); set_row(_j,ri);
    }
    return *this;
  }

  /// get (copy of) diagonal
  diag_t diag() const {
    static_assert(M==N,"available only for square matrices");
    diag_t d;
    for (int i=0;i<M;++i) d[i]=m_raw[i+i*M];
    return d;
  }
  /// set matrix diagonal
  void set_diag(const diag_t& _d) {
    static_assert(M==N,"available only for square matrices");
    for (int i=0;i<M;++i) m_raw[i+i*M]=_d[i];
  }

  /// alias for `VC::math::v_cat(*this,a)`
  template <int K>
  Mat<T,M+K,N>
  v_cat(const Mat<T,K,N>& _a) const { return VC::math::v_cat(*this,_a); }

  /// alias for `VC::math::h_cat(*this,a)`
  template <int K>
  Mat<T,M,N+K>
  h_cat(const Mat<T,M,K>& _a) const { return VC::math::h_cat(*this,_a); }

  /// alias for `VC::math::repmat<R,C>::get(*this)`
  template <int R,int C>
  Mat<T,R*M,C*N> repmat() const { return VC::math::repmat<R,C>::get(*this); }

  /// @}


  /** @name access via references
      \sa vc_blas
      @{
  */

  /// access matrix via \ref vc_blas
  const_ref_t const_ref() const {
    return const_ref_t(matrix_t(M,N,M),m_raw);
  }
  /// access matrix via \ref vc_blas
  ref_t ref() {
    return ref_t(matrix_t(M,N,M),m_raw);
  }

  /// alias for const_ref()
  const_ref_t operator()() const { return const_ref(); }
  /// alias for ref()
  ref_t operator()() { return ref(); }


  /// access column via \ref vc_blas
  const_col_ref_t const_col_ref(int _j) const {
    assert(0<=_j && _j<N);
    return const_col_ref_t(col_vec_t(N,1),m_raw[_j*M]);
  }
  /// access column via \ref vc_blas
  col_ref_t col_ref(int _j) {
    assert(0<=_j && _j<N);
    return col_ref_t(col_vec_t(N,1),m_raw[_j*M]);
  }

  /// access row via \ref vc_blas
  const_row_ref_t const_row_ref(int _i) const {
    assert(0<=_i && _i<M);
    return const_row_ref_t(row_vec_t(M,N),m_raw[_i]);
  }
  /// access row via \ref vc_blas
  row_ref_t row_ref(int _i) {
    assert(0<=_i && _i<M);
    return row_ref_t(row_vec_t(N,N),m_raw[_i]);
  }

  /// @}


  /** @name initialization
      @{
  */
  /// set all element to `_x`
  self_t& load_all(T _x) { for (int i=0;i<M*N;++i) m_raw[i]=_x; return *this; }
  /// set all elements to `0` (load_all(0))
  self_t& load_zeros() { load_all(T(0)); return *this; }
  /// load identity matrix
  self_t& load_eye() {
    for (int j=0;j<N;++j)
      for (int i=0;i<M;++i)
        this->m_raw[i+j*M]=(i==j) ? T(1) : T(0);
    return *this;
  }
  /// load random values in t`[0,`]`
  self_t& load_rand() { m_elts.rand(); return *this; }
  /// load a diagonal matrix with diagonal `_d`
  self_t& load_diag_matrix(const col_t& _d) {
    static_assert(M==N,"available only for square matrices");
    for (int j=0;j<N;++j)
      for (int i=0;i<M;++i)
        this->m_raw[i+j*M]=(i==j) ? _d[i] : T(0);
    return *this;
  }
  /// load rank-1 matrix `x*y'`
  self_t& load_rank1_matrix(const VecN<T,M>& _x,const VecN<T,N>& _y) {
    for (int j=0;j<N;++j)
      this->col(j)=_x*_y[j];
    return *this;
  }

  /// @}


  /** @name operators
      @{
  */

  // TODO: improve operatorX=(const ex_trans_t)

  /// assign transposed matrix from expression (trans())
  self_t& operator=(const ex_trans_t<T,N,M>& _a) {
    return (*this=*_a);
  }
  /// assign diagonal matrix from expression (VA::math::diag())
  self_t& operator=(const ex_diag_t<T,M,N>& _d) { return (*this=*_d); }
  /// assignment
  self_t& operator=(const self_t& _a) { this->m_elts=_a.m_elts; return *this; }

  /// assign w/ element-wise cast
  template <typename S>
  self_t&
  operator=(const Mat<S,M,N>& _a) { this->m_elts=_a.m_elts; return *this; }

  /// define `A=blas::zeros` (load_zeros())
  self_t& operator=(const blas::Expression<blas::Type_initializer,blas::Expr_zeros>&) {
    load_zeros(); return *this;
  }
  /// define `A=blas::random_values` (load_rand())
  self_t& operator=(const blas::Expression<blas::Type_initializer,blas::Expr_random_values>&) {
    load_rand(); return *this;
  }
  /// define `A=blas::ones` (load_all(1))
  self_t& operator=(const blas::Expression<blas::Type_initializer,blas::Expr_ones>&) {
    load_all(T(1)); return *this;
  }
  /// define `A=blas::all(x)` (load_all(x))
  self_t& operator=(const blas::Expression<blas::Type_initializer,blas::Expr_all,T>& _all) {
    load_all(_all.a); return *this;
  }
  /// define `A=blas::eye` (load_eye(1))
  self_t& operator=(blas::Expression<blas::Type_initializer,blas::Expr_eye>&) {
    load_eye(); return *this;
  }

  /// get transposed matrix as a simple exession
  ex_trans_t<T,M,N> trans() const { return ex_trans_t<T,M,N>(*this); }

  /// add `A+=B`
  self_t&
  operator+=(const self_t& _a) { this->m_elts+=_a.m_elts; return *this; }
  /// subtract `A-=B`
  self_t&
  operator-=(const self_t& _a) { this->m_elts-=_a.m_elts; return *this; }
  /// **element-wise** multiplication `A=A.*B`
  self_t&
  operator*=(const self_t& _a) { this->m_elts*=_a.m_elts; return *this; }
  /// **element-wise** division `A=A./B`
  self_t&
  operator/=(const self_t& _a) { this->m_elts/=_a.m_elts; return *this; }

  /// add transposed `A+=B'`
  self_t&
  operator+=(const ex_trans_t<T,N,M>& _a) { return *this+=*_a; }
  /// subtract transposed `A+=B'`
  self_t&
  operator-=(const ex_trans_t<T,N,M>& _a) { return *this-=*_a; }
  /// **element-wise** multiplication w/ transposed `A=A.*B'`
  self_t&
  operator*=(const ex_trans_t<T,N,M>& _a) { return *this*=*_a; }
  /// **element-wise** division w/ transposed `A=A./B'`
  self_t&
  operator/=(const ex_trans_t<T,N,M>& _a) { return *this/=*_a; }

  /// add diagonal `A+=diag(x)`
  self_t& operator+=(const ex_diag_t<T,M,N>& _d) {
    for (int i=0;i<M;++i) this->m_raw[i+i*M]+=_d.diagonal[i]; return *this;
  }
  /// subtract diagonal `A+=diag(x)`
  self_t& operator-=(const ex_diag_t<T,M,N>& _d) {
    for (int i=0;i<M;++i) this->m_raw[i+i*M]-=_d.diagonal[i]; return *this;
  }

  /// add scalar element-wise `A+=_x`
  self_t& operator+=(T _x) { this->m_elts+=_x; return *this; }
  /// subtract scalar element-wise `A-=_x`
  self_t& operator-=(T _x) { this->m_elts-=_x; return *this; }
  /// multiply by scalar `A*=_x
  self_t& operator*=(T _x) { this->m_elts*=_x; return *this; }
  /// divide by scalar `A/=_x`
  self_t& operator/=(T _x) { this->m_elts/=_x; return *this; }

  /// (right) multiply by scalar `A*_x`
  self_t operator*(T _x) const { self_t c=*this; c*=_x; return c; } // TODO: std::move
  /// divide by scalar `A*_x`
  self_t operator/(T _x) const { self_t c=*this; c/=_x; return c; }

  /// matrix-vector product `A*x` (ex_trans_t defines `A'*x`)
  col_t operator*(const row_t& _x) const {
    col_t y=this->col(0)*_x[0];
    for (int j=1;j<N;++j)
      y+=this->col(j)*_x[j];
    return y;
  }

  /// @}

  /** @name Operations (some may not be available depending on matrix type)
      @{
  */

  /// get trace
  T trace() const {
    static_assert(M==N,"require square matrix");
    T s=m_raw[0];
    for (int i=1;i<N;++i)
      s+=m_raw[i+i*N];
    return s;
  }

  /// get Frobenius norm
  T normf() const { return m_elts.norm(); }

  /// get maximum norm
  T normmax() const { return m_elts.norminf(); }

  /// get infinity norm (largest row sum of the absolute values)
  T norminf() const {
    col_t s=col(0).abs();
    for (int j=1;j<N;++j)
      s+=col(j).abs();
    return s.norminf();
  }

  /** Load inverse.
      \return `false` if matrix is singular (in this case the matrix is *not*
      modified)
      \sa inverse_t
  */
  bool invert() { return inverse_t<T,M,N>::get(*this); }

  /** Same as invert() but force use of lu() to load inverse.
      \return `false` if matrix is singular (in this case the matrix is *not*
      modified)
      \sa lu_factor_t
   */
  bool lu_invert() {
    static_assert(M==N,"require square matrix");
    VecN<int,N> pinv;
    Mat<T,M,N> tmp(*this);
    if (!tmp.lu(pinv))
      return false;
    this->load_eye();
    tmp.lu_subs(pinv,*this);
    return true;
  }

  /** Compute factorization `A=L*U` with row permutation.
      Load factors `L` and `U` as lower an upper triangle of `*this`.
      \param[out] _pinv stores row permutation
      \return `false` if matrix is singular, the contents of `*this` are
      **destroyed** in this case
      \sa lu_factor_t
   */
  bool lu(VecN<int,M>& _pinv) { return lu_factor_t<T,N,N>::get(*this,_pinv); }

  /** Solve linear system `A*x=b` from lu() factorization `A=L*U`.
      Requires successful call to lu() immediately before solving.
      \param _pinv permutation as computed by lu()
      \param[in,out] _x on **input**: `_nrhs` columns of right hand side,
      on **output**: solution
      _param _nrhs number of right hand sides
      \sa lu_factor_t
   */
  void lu_subs(const VecN<int,M>& _pinv,VecN<T,M>* _x,int _nrhs) {
    lu_factor_t<T,M,N>::apply(*this,_pinv,_x,_nrhs);
  }
  /// same as above for `_nrhs==1`
  void lu_subs(const VecN<int,M>& _pinv,VecN<T,M>& _x) {
    lu_factor_t<T,M,N>::apply(*this,_pinv,&_x,1);
  }
  /// same as above for `_nrhs==NRHS`
  template <int NRHS>
  void lu_subs(const VecN<int,M>& _pinv,Mat<T,M,NRHS>& _x) {
    lu_factor_t<T,M,N>::apply(*this,_pinv,&_x.col(0),NRHS);
  }

  /// get determinant \sa det_t
  T det() const { return det_t<T,M,N>::get(*this); }

  /** Get eigenvalues of symmetric matrix.
      Only the **upper** triangle of the matrix is used.
      Symmetry is **not** tested!
      \return vector of real eigenvalues, on error (e.g.,
      VC::math::lapack::syev() failed to converge) test VecN::is_finite() fails
      \sa sy_eig_t
   */
  col_t sy_eig() const { return sy_eig_t<T,M,N>::get(*this); }

  /** Get eigenvalues and eigenvectors of symmetric matrix.
      Only the **upper** triangle of the matrix is used.
      Symmetry is **not** tested!

      **Note:** As for `g++-4.8.1` math optimization seems more
      aggressive (and unfortunately less correct), we currently
      recommend use of sy_eig_lapack() for matrices with `M>2` to
      avoid numerical errors. (Although this is not expected,
      non-optimized and optimized code may behave very differently in
      terms or errors for specializations.)

      \param _v eigenvectors corresponding to returned eigenvalues
      \return vector of real eigenvalues, on error (e.g.,
      VC::math::lapack::syev() failed to converge) test VecN::is_finite() fails
      \sa sy_eig_t

      \bug Aggressive math optimization of new gcc versions may have
      large impact on results for specialized versions. Check this!
   */
  col_t sy_eig(self_t& _v) const { return sy_eig_t<T,M,N>::get(*this,_v); }

  /// same as above but force using VC::math::lapack::syev()
  col_t sy_eig_lapack() const { return sy_eig_t<T,M,N,true>::get(*this); }
  /// same as above but force using VC::math::lapack::syev()
  col_t sy_eig_lapack(self_t& _v) const {
    return sy_eig_t<T,M,N,true>::get(*this,_v);
  }

  /** Get eigenvalues of square matrix.
      \return real (column 0) and imaginary (column 1) parts of
      eigenvalues.  On error (e.g., VC::math::lapack::geev() may fail
      to converge) the test `is_finite()` fails for the first column
      of the returned matrix).
   */
  Mat<T,N,2> eig() const { return eig_t<T,M,N>::get(*this); }

  /** Get eigenvalues and (right) eigenvectors of square matrix.
      \param[out] _v (complex conjugate pairs of) eigenvectors,
      see see VC::math::lapack::geev().
      \return real (column 0) and imaginary (column 1) parts of
      eigenvalues On error (e.g., VC::math::lapack::geev() may fail to
      converge) the test `is_finite()` fails for the first column of
      the returned matrix).
   */
  Mat<T,N,2> eig(self_t& _v) const { return eig_t<T,M,N>::get(*this,_v); }

  /// same as above but force using VC::math::lapack::geev()
  Mat<T,N,2> eig_lapack() const { return eig_t<T,M,N,true>::get(*this); }
  /// same as above but force using VC::math::lapack::geev()
  Mat<T,N,2> eig_lapack(self_t& _v) const {
    return eig_t<T,M,N,true>::get(*this,_v);
  }

  /// @}
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** \defgroup vc_math_lam_mops Linear algebra: externally defined operations on Mat
    \ingroup vc_math_lam
    \sa Mat
*/

/** @name Externally defined operations on Mat
    @{
*/

/// external definition of Mat<T,M,N>::to_m() \ingroup vc_math_lam_mops
template <int M,int N,typename T>
Mat<T,M,N>& to_m(T* _m) {
  return (Mat<T,M,N>&) *_m;
}
/// external definition of Mat<T,M,N>::to_cm() \ingroup vc_math_lam_mops
template <int M,int N,typename T>
const Mat<T,M,N>& to_cm(const T* _m) {
  return *((const Mat<T,M,N>*) _m);
}

/// cast `M`-vector to `Mx1` matrix
template <typename T,int M>
Mat<T,M,1>& to_m(VecN<T,M>& _v) {
  return Mat<T,M,1>::to_m(_v.data());
}
/// cast `M`-vector to `Mx1` matrix
template <typename T,int M>
const Mat<T,M,1>& to_cm(const VecN<T,M>& _v) {
  return Mat<T,M,1>::to_cm(_v.data());
}

/// cast row `N`-vector to `1xN` matrix
template <typename T,int N>
Mat<T,1,N>& to_m(VecN<T,N>& _v,RowMajor_t) {
  return Mat<T,1,N>::to_m(_v.data());
}
/// cast row `N`-vector to `1xN` matrix
template <typename T,int N>
const Mat<T,1,N>& to_cm(const VecN<T,N>& _v,RowMajor_t) {
  return Mat<T,1,N>::to_cm(_v.data());
}

// NEED TO MAKE binop from left w/ SCALAR EXPLICIT !!!

// Make sure there is no confusion with VecN operators !!!


template <typename T,int M,int N>
Mat<T,M,N> operator*(T _b,const Mat<T,M,N>& _a) { return _a*_b; }

template <typename T,int M,int N>
Mat<T,M,N> operator/(T _b,const Mat<T,M,N>& _a) { return _a*(T(1)/_b); }


/// matrix-vector product `(x'*A_' = A'*x'` \ingroup vc_math_lam_mops
template <typename T,int M,int N>
typename Mat<T,M,N>::row_t
operator*(const typename Mat<T,M,N>::col_t& _y,const Mat<T,M,N>& _a) {
  return _a.trans()*_y;
}
/// matrix-vector product with transposed `A'*x` \ingroup vc_math_lam_mops
template <typename T,int M,int N>
typename Mat<T,M,N>::row_t
operator*(const ex_trans_t<T,M,N>& _a,const typename Mat<T,M,N>::col_t& _y) {
  typename Mat<T,M,N>::row_t x=_a.row(0)*_y[0];
  for (int i=1;i<M;++i)
    x+=_a.row(i)*_y[i];
  return x;
}

// >> candidates for optimization

template <typename T,int M,int K,int N>
Mat<T,M,N> operator*(const Mat<T,M,K>& _a,const Mat<T,K,N>& _b) {
  Mat<T,M,N> c;
  blas::vc_gemm(blas::NoT,blas::NoT, M,N,K,
                T(1),_a.data(),M,_b.data(),K,T(0),c.data(),M);
  return c;
}
template <typename T,int M,int K,int N>
Mat<T,M,N> operator*(const ex_trans_t<T,K,M>& _a,const Mat<T,K,N>& _b) {
  Mat<T,M,N> c;
  blas::vc_gemm(blas::Transpose,blas::NoT, M,N,K,
                T(1),_a.matrix.data(),M,_b.data(),K,T(0),c.data(),M);
  return c;
}
template <typename T,int M,int K,int N>
Mat<T,M,N> operator*(const Mat<T,M,K>& _a,const ex_trans_t<T,N,K>& _b) {
  Mat<T,M,N> c;
  blas::vc_gemm(blas::NoT,blas::Transpose, M,N,K,
                T(1),_a.data(),M,_b.matrix.data(),K,T(0),c.data(),M);
  return c;
}
template <typename T,int M,int K,int N>
Mat<T,M,N> operator*(const ex_trans_t<T,K,M>& _a,const ex_trans_t<T,N,K>& _b) {
  Mat<T,M,N> c;
  blas::vc_gemm(blas::Transpose,blas::Transpose, M,N,K,
                T(1),_a.matrix.data(),M,_b.matrix.data(),K,T(0),c.data(),M);
  return c;
}

// <<


template <typename T,int M,int N>
Mat<T,M,N> operator+(const Mat<T,M,N>& _a,const Mat<T,M,N>& _b) {
  Mat<T,M,N> c=_a; c+=_b;
  return c;
}
template <typename T,int M,int N>
Mat<T,M,N> operator+(Mat<T,M,N>&& _a,const Mat<T,M,N>& _b) {
  _a+=_b; return std::move(_a);
}
template <typename T,int M,int N>
Mat<T,M,N> operator+(const Mat<T,M,N>& _a,Mat<T,M,N>&& _b) {
  _b+=_a; return std::move(_b);
}

template <typename T,int M,int N>
Mat<T,M,N> operator-(const Mat<T,M,N>& _a,const Mat<T,M,N>& _b) {
  Mat<T,M,N> c=_a; c-=_b;
  return c;
}
template <typename T,int M,int N>
Mat<T,M,N> operator-(Mat<T,M,N>&& _a,const Mat<T,M,N>& _b) {
  _a-=_b; return std::move(_a);
}

template <typename T,int M,int N>
Mat<T,M,N> operator+(const Mat<T,M,N>& _a,const ex_trans_t<T,M,N>& _b) {
  return _a+*_b;
}
template <typename T,int M,int N>
Mat<T,M,N> operator+(const ex_trans_t<T,M,N>& _a,const Mat<T,M,N>& _b) {
  return *_a+_b;
}

template <typename T,int M,int N>
Mat<T,M,N> operator-(const Mat<T,M,N>& _a,const ex_trans_t<T,M,N>& _b) {
  return _a-*_b;
}
template <typename T,int M,int N>
Mat<T,M,N> operator-(const ex_trans_t<T,M,N>& _a,const Mat<T,M,N>& _b) {
  return *_a-_b;
}

template <typename T,int M,int N>
const Mat<T,M,N>& operator+(const Mat<T,M,N>& _a) { return _a; }

template <typename T,int M,int N>
Mat<T,M,N> operator-(const Mat<T,M,N>& _a) {
  Mat<T,M,N> c=_a; c*=T(-1); return c;
}

template <typename T,int M,int N>
Mat<T,M,N> operator-(Mat<T,M,N>&& _a) {
  _a*=T(-1); return std::move(_a);
}

//-----------------------------------------------------------------------------

/// expression for diagonal matrix from vector \ingroup vc_math_lam_mops
template <typename T,unsigned N>
ex_diag_t<T,N,N> diag(const VecN<T,N>& _d) {
  return ex_diag_t<T,N,N>(_d);
}

template <typename T,int M,int N>
Mat<T,M,N>
operator*(const Mat<T,M,N>& _a,const ex_diag_t<T,M,N>& _d) {
  Mat<T,M,N> c;
  for (int j=0;j<N;++j)
    c.col(j)=_a.col(j)*_d.diagonal; // scale rows w/ element-wise c.*d
  return c;
}
template <typename T,int M,int N>
Mat<T,M,N>
operator*(const Mat<T,M,N>&& _a,const ex_diag_t<T,M,N>& _d) {
  for (int j=0;j<N;++j)
    _a.col(j)*=_d.diagonal; // scale rows w/ element-wise c.*d
  return std::move(_a);
}

template <typename T,int M,int N>
Mat<T,M,N>
operator/(const Mat<T,M,N>& _a,const ex_diag_t<T,M,N>& _d) {
  Mat<T,M,N> c;
  for (int j=0;j<N;++j)
    c.col(j)=_a.col(j)/_d.diagonal; // scale rows w/ element-wise c./d
  return c;
}
template <typename T,int M,int N>
Mat<T,M,N>
operator/(const Mat<T,M,N>&& _a,const ex_diag_t<T,M,N>& _d) {
  for (int j=0;j<N;++j)
    _a.col(j)/=_d.diagonal; // scale rows w/ element-wise c./d
  return std::move(_a);
}

template <typename T,int M,int N>
Mat<T,M,N>
operator*(const ex_diag_t<T,M,N>& _d,const Mat<T,M,N>& _a) {
  Mat<T,M,N> c;
  for (int j=0;j<N;++j)
    c.col(j)=_a.col(j)*_d.diagonal[j]; // scale columns
  return c;
}
template <typename T,int M,int N>
Mat<T,M,N>
operator*(const ex_diag_t<T,M,N>& _d,const Mat<T,M,N>&& _a) {
  for (int j=0;j<N;++j)
    _a.col(j)*=_d.diagonal[j]; // scale columns
  return std::move(_a);
}

template <typename T,int M,int N>
Mat<T,M,N>
operator/(const ex_diag_t<T,M,N>& _d,const Mat<T,M,N>& _a) {
  Mat<T,M,N> c;
  for (int j=0;j<N;++j)
    c.col(j)=_a.col(j)/_d.diagonal[j]; // scale columns
}
template <typename T,int M,int N>
Mat<T,M,N>
operator/(const ex_diag_t<T,M,N>& _d,const Mat<T,M,N>&& _a) {
  for (int j=0;j<N;++j)
    _a.col(j)/=_d.diagonal[j]; // scale columns
  return std::move(_a);
}

//@}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/// stack horizontally \ingroup vc_blas_lam
template <typename T,int M,int K,int N>
Mat<T,M,K+N> h_cat(const Mat<T,M,K>& _a,const Mat<T,M,N>& _b) {
  Mat<T,M,K+N> c;
  Mat<T,M,K>::to_m(c.data())=_a;
  Mat<T,M,N>::to_m(c.data()+M*K)=_b;
  return c;
}

/// stack horizontally \ingroup vc_blas_lam
template <typename T,int M,int N>
Mat<T,M,1+N> h_cat(const VecN<T,unsigned(M)>& _a,const Mat<T,M,N>& _b) {
  Mat<T,M,1+N> c;
  c.col(0)=_a;
  Mat<T,M,N>::to_m(c.data()+M)=_b;
  return c;
}

/// stack horizontally \ingroup vc_blas_lam
template <typename T,int M,int N>
Mat<T,M,N+1> h_cat(const Mat<T,M,N>& _b,const VecN<T,unsigned(M)>& _a) {
  Mat<T,M,N+1> c;
  Mat<T,M,N>::to_m(c.data())=_b;
  c.col(N)=_a;
  return c;
}

/// stack horizontally \ingroup vc_blas_lam
template <typename T,unsigned M>
Mat<T,M,2> h_cat(const VecN<T,M>& _a,const VecN<T,M>& _b) {
  Mat<T,M,2> c;
  c.col(0)=_a;
  c.col(1)=_b;
  return c;
}

//-----------------------------------------------------------------------------

/// stack vertically \ingroup vc_blas_lam
template <typename T,int M,int K,int N>
Mat<T,M+K,N> v_cat(const Mat<T,M,K>& _a,const Mat<T,M,N>& _b) {
  Mat<T,M+K,N> c;
  c()(blas::_(0,M-1),blas::$)=_a();
  c()(blas::_(M,blas::end),blas::$)=_b();
  return c;
}

/// stack vertically \ingroup vc_blas_lam
template <typename T,int M,int N>
Mat<T,1+M,N> v_cat(const VecN<T,unsigned(M)>& _a,const Mat<T,M,N>& _b) {
  Mat<T,1+M,N> c;
  c.set_row(0)=_a;
  c()(blas::_(1,blas::end),blas::$)=_b();
  return c;
}

/// stack vertically \ingroup vc_blas_lam
template <typename T,int M,int N>
Mat<T,M+1,N> v_cat(const Mat<T,M,N>& _b,const VecN<T,unsigned(M)>& _a) {
  Mat<T,M+1,1> c;
  c()(blas::_(0,blas::end-1),blas::$)=_b();
  c.set_row(M)=_a;
  return c;
}

/// stack vertically \ingroup vc_blas_lam
template <typename T,unsigned N>
Mat<T,2,N> v_cat(const VecN<T,N>& _a,const VecN<T,N>& _b) {
  Mat<T,2,N> c;
  c.set_row(0,_a);
  c.set_row(1,_b);
  return c;
}

//-----------------------------------------------------------------------------

/// Kronecker product of matrices \ingroup vc_blas_lam
template <typename T,int M,int N,int P,int R>
Mat<T,M*P,N*R> kron(const Mat<T,M,N>& _a,const Mat<T,P,R>& _b) {
  Mat<T,M*P,N*R> c;
  using blas::_;
  for (int j=0;j<N;++j)
    for (int i=0;i<M;++i) {
      c()(_(i*P,(i+1)*P-1),_(j*R,(j+1)*R-1))=_b()*_a(i,j);
    }
  return c;
}

//-----------------------------------------------------------------------------

/// repeat matrix block \ingroup vc_math_lam
template <int R,int C>
struct repmat {
  /// `repmat<r,c>::get(A)` is `repmat(A,r,c)` in Matlab
  template <typename T,int M,int N>
  static Mat<T,M*R,N*C> get(const Mat<T,M,N>& _a) {
    Mat<T,M*R,N*C> c;

    for (int j=0;j<N;++j)
      for (int i=0;i<R;++i)
        VecN<T,M>::to_v(c.data()+j*M*R+i*M)=_a.col(j);

    for (int jj=1;jj<C;++jj)
      for (int j=0;j<N;++j)
        c.col(jj*N+j)=c.col(j);
    return c;
  }
};

//-----------------------------------------------------------------------------

/// get column sums \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,N> sum(const Mat<T,M,N>& _a) {
  VecN<T,N> s;
  for (int j=0;j<N;++j)
    s[j]=_a.col(j).sum();
  return s;
}
/// get row sums (column sums of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,M> sum(const ex_trans_t<T,M,N>& _a) {
  VecN<T,M> s;
  for (int i=0;i<M;++i)
    s[i]=_a.matrix.row(i).sum();
  return s;
}

/// get mean of columns \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,N> mean(const Mat<T,M,N>& _a) { return sum(_a)/T(M); }
/// get mean of rows (column mean of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,M> mean(const ex_trans_t<T,M,N>& _a) { return sum(_a)/T(N); }

/// get column minima \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,N> min(const Mat<T,M,N>& _a) {
  VecN<T,N> s;
  for (int j=0;j<N;++j)
    s[j]=_a.col(j).min();
  return s;
}
/// get row minima (column minima of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,M> min(const ex_trans_t<T,M,N>& _a) {
  VecN<T,M> s;
  for (int i=0;i<M;++i)
    s[i]=_a.matrix.row(i).min();
  return s;
}

/// get column minimum indices \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<unsigned,N> imin(const Mat<T,M,N>& _a) {
  VecN<unsigned,N> s;
  for (int j=0;j<N;++j)
    s[j]=_a.col(j).imin();
  return s;
}
/// get row minimum indices (column indices of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<unsigned,M> imin(const ex_trans_t<T,M,N>& _a) {
  VecN<unsigned,M> s;
  for (int i=0;i<M;++i)
    s[i]=_a.matrix.row(i).imin();
  return s;
}

/// get column maxima \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,N> max(const Mat<T,M,N>& _a) {
  VecN<T,N> s;
  for (int j=0;j<N;++j)
    s[j]=_a.col(j).max();
  return s;
}
/// get row maxima (column maxima of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<T,M> max(const ex_trans_t<T,M,N>& _a) {
  VecN<T,M> s;
  for (int i=0;i<M;++i)
    s[i]=_a.matrix.row(i).max();
  return s;
}

/// get column maximum indices \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<unsigned,N> imax(const Mat<T,M,N>& _a) {
  VecN<unsigned,N> s;
  for (int j=0;j<N;++j)
    s[j]=_a.col(j).imax();
  return s;
}
/// get row maximum indices (column indices of transposed) \ingroup vc_math_lam
template <typename T,int M,int N>
VecN<unsigned,M> imax(const ex_trans_t<T,M,N>& _a) {
  VecN<unsigned,M> s;
  for (int i=0;i<M;++i)
    s[i]=_a.matrix.row(i).imax();
  return s;
}

/// get Frobenius norm of `_a` \ingroup vc_math_lam
template <typename T,int M,int N>
T normf(const Mat<T,M,N>& _a) { return _a.normf(); }

/// get maximum norm of `_a` \ingroup vc_math_lam
template <typename T,int M,int N>
T normmax(const Mat<T,M,N>& _a) { return _a.norminf(); }

/** get infinity norm of `_a` (largest row sum of the absolute values)
    \ingroup vc_math_lam
*/
template <typename T,int M,int N>
T norminf(const Mat<T,M,N>& _a) { return _a.norminf(); }


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Stream output of VC::math::Mat.
    \ingroup vc_math_lam
    Outputs Mat::elts(), i.e., elements column-wise w/o any formatting.
    Use Mat::const_ref() and, e.g., VC::math::blas::pmatlab(), for formatted
    output.
*/
template <typename T,int M,int N>
std::ostream& operator<<(std::ostream& _s,const ::VC::math::Mat<T,M,N>& _a) {
  return _s << _a.elts();
}
/** Stream output of VC::math::Mat.
    \ingroup vc_math_lam
    Reads Mat::elts(), i.e., elements column-wise.
*/
template <typename T,int M,int N>
std::istream& operator>>(std::istream& _s,::VC::math::Mat<T,M,N>& _a) {
  return _s >> _a.elts();
}

/** Read matrix from MATLAB V4 file.
    \ingroup vc_math_lam
    Read matrix from Level 4 MAT-file, see VC::math::blas::write_mat4().
    In contrast to VC::math::mat4::read_matrix() this function does *not*
    throw exceptions.
    \param _in input stream
    \param _name[out] matrix name *or* error message
    \return success, on failure an error message is returned in `_name`
    (e.g., dimension mismatch)
    \sa write_mat4(), [\ref vc_mat4], VC::math::blas::GE_Matrix::read_mat4()
*/
template <typename T,int M,int N>
bool read_mat4(std::istream& _in,std::string& _name,Mat<T,M,N>& _a) {
  mat4::MatrixInfo mat4;
  try {
    mat4.read(_in);
    _name=mat4.name();

    if (mat4.n_imag()>0) {
      _name="cannot read complex matrix";
      return false;
    }
    if (mat4.mrows()!=M || mat4.ncols()!=N) {
      _name=VC::base::formatf("dimension mismatch: expected %dx%d, got %dx%d",
                              M,N,mat4.mrows(),mat4.ncols());
      return false;
    }
    mat4.read_matrix(_in,M,_a.data());
  } catch (base::vc_runtime_error& e) {
    _name=e.what();
    return false;
  }
  return true;
}

/** Write matrix to MATLAB V4 file.
    \ingroup vc_math_lam
    Short for VC::math::blas::write_mat4(_out,_a->const_ref(),_name).
*/
template <typename T,int M,int N>
std::ostream& write_mat4(std::ostream& _out,
                         const std::string& _name,const Mat<T,M,N>& _a) {
  return ::VC::math::blas::write_mat4(_out,_a->const_ref(),_name);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Interpret square matrix as transformation of homogeneous coordinates.
    \ingroup vc_math_lam
    This class requires specialization!
 */
template <typename T,int N>
struct transformation_t {
  Mat<T,N,N>& m; //!< the "wrapped" matrix
  transformation_t(Mat<T,N,N>&) {
    static_assert(N!=N,"no specialization available");
  }
  /// get matrix
  const Mat<T,N,N> get() const { return m; }
};

/// Wrap `_m` as a homogeneous transformation \ingroup vc_math_lam
template <typename T,int N>
transformation_t<T,N>
transformation(Mat<T,N,N>& _m) { return transformation_t<T,N>(_m); }

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Compute and apply LU factorization of square matrix.
    \ingroup vc_math_lam
    The standard implementation uses lu() and lu_subs(). There may be
    specializations.
*/
template <typename T,int M,int N>
struct lu_factor_t {
  /** Compute factor.
    \param _a square matrix
    \param _pinv row permutation
    \return success
  */
  static bool get(Mat<T,M,N>& _a,VecN<int,M>& _pinv) {
    static_assert(M==N,"require square matrix");
    return lu<T,N>(_a.data(),_pinv.data());
  }
  /** Apply factor.
    \param _a factorized matrix (from lu())
    \param _pinv row permutation (from lu())
    \param[in,out] _x `_nrhs` col;umns of RHS
    \param _nrhs number of RHS
  */
  static void apply(const Mat<T,M,N>& _a,const VecN<int,M>& _pinv,
                    VecN<T,M>* _x,int _nrhs) {
    static_assert(M==N,"require square matrix");
    lu_subs<T,N>(_a.data(),_pinv.data(),_x->data(),_nrhs);
  }
};

/** Compute inverse of square matrix.
    \ingroup vc_math_lam
    The standard implementation uses lu() and lu_subs(). There may be
    specializations.
*/
template <typename T,int M,int N>
struct inverse_t {
  /** Compute inverse of square matrix.
    \param _a square matrix
    \return success
  */
  static bool get(Mat<T,M,N>& _a) {
    static_assert(M==N,"require square matrix");
    VecN<int,N> pinv;
    Mat<T,M,N> tmp(_a);
    if (!lu<T,N>(tmp.data(),pinv.data()))
      return false;
    _a.load_eye();
    lu_subs<T,N>(tmp.data(),pinv.data(),_a.data(),N);
    return true;
  }
};

/** Compute determinant of square matrix.
    \ingroup vc_math_lam
    This function is currently implemented **only** in specializations!
*/
template <typename T,int M,int N>
struct det_t {
  /** Get determinant.
      \param _a square matrix
      \return determinant
   */
  static T get(const Mat<T,M,N>& _a) {
    use_nowarn(_a);
    static_assert(M==M,"not implemented"); return T(0);
  }
};

/** Compute eigenvalues of symmetric matrix.
    \ingroup vc_math_lam
    \arg Only the **upper** triangle of the matrix is used.
    \arg Symmetry is **not** tested!
    \arg The standard implementation uses lapack::syev(). There may be
    specializations.
*/
template <typename T,int M,int N,bool FORCE_LAPACK>
struct sy_eig_t {

  /** Compute eigenvalues.
      \param _a symmetric matrix
      \return eigenvalues `w` in **ascending** order, on failure (convergence
      failed) `VecN::is_finite(w)==false`
  */
  static typename Mat<T,M,N>::col_t get(const Mat<T,M,N>& _a) {
    static_assert(M==M,"require square matrix");
    static_assert(std::is_same<T,float>::value ||
                  std::is_same<T,double>::value,"require float or double");
    static const int LWORK=lapack::syev_lwork(N);
    T* work=(T*) alloca(LWORK*sizeof(T));
    Mat<T,M,N> a=_a;
    typename Mat<T,M,N>::col_t w;
    int info=lapack::syev('N',blas::Upper,N,a.data(),N,w.data(),work,LWORK);
    assert(info>=0 && "invalid argument"!=0);
    if (info!=0)
      w[0]=std::numeric_limits<T>::quiet_NaN();
    return w;
  }

  /** Compute eigenvalues and eigenvectors.
      \param _a symmetric matrix
      \param[out] _v eigenvectors corresponding to returned eigenvalues
      \return eigenvalues `w` in **ascending** order, on failure (convergence
      failed) `VecN::is_finite(w)==false`
  */
  static typename Mat<T,M,N>::col_t get(const Mat<T,M,N>& _a, Mat<T,M,N>& _v) {
    static_assert(M==M,"require square matrix");
    static_assert(std::is_same<T,float>::value ||
                  std::is_same<T,double>::value,"require float or double");
    static const int LWORK=lapack::syev_lwork(N);
    T* work=(T*) alloca(LWORK*sizeof(T));
    _v=_a;
    typename Mat<T,M,N>::col_t w;
    int info=lapack::syev('V',blas::Upper,N,_v.data(),N,w.data(),work,LWORK);
    assert(info>=0 && "invalid argument"!=0);
    if (info!=0)
      w[0]=std::numeric_limits<T>::quiet_NaN();
    return w;
  }
};

/** Compute eigenvalues of matrix.
    \ingroup vc_math_lam
    \arg The standard implementation uses VC::math::lapack::geev().
    There may be specializations.
*/
template <typename T,int M,int N,bool FORCE_LAPACK>
struct eig_t {

  /** Compute eigenvalues.
      \param _a square matrix
      \return `w` with real (`w.col(0)`) and imaginary (`w.col(1)`)
      parts eigenvalues, on failure (convergence failed)
      `VecN::is_finite(w.col(0))==false`
  */
  static Mat<T,M,2> get(const Mat<T,M,N>& _a) {
    static_assert(M==M,"require square matrix");
    static_assert(std::is_same<T,float>::value ||
                  std::is_same<T,double>::value,"require float or double");
    static const int LWORK=lapack::geev_lwork('N','N',N);
    T* work=(T*) alloca(LWORK*sizeof(T));
    Mat<T,M,N> a=_a;
    Mat<T,M,2> w;
    int info=lapack::geev('N','N',N,a.data(),N,w.data(),w.data()+M,
                          nullptr,1,nullptr,1,work,LWORK);
    assert(info>=0 && "invalid argument"!=0);
    if (info!=0)
      w.data()[0]=std::numeric_limits<T>::quiet_NaN();
    return w;
  }

  /** Compute eigenvalues and (right) eigenvectors.
      \param _a square matrix
      \param[out] _v (complex conjugate pairs of) eigenvectors,
      see VC::math::lapack::geev()
      \return `w` with real (`w.col(0)`) and imaginary (`w.col(1)`)
      parts eigenvalues, on failure (convergence failed)
      `VecN::is_finite(w.col(0))==false`
  */
  static Mat<T,M,2> get(const Mat<T,M,N>& _a,Mat<T,M,N>& _v) {
    static_assert(M==M,"require square matrix");
    static_assert(std::is_same<T,float>::value ||
                  std::is_same<T,double>::value,"require float or double");
    static const int LWORK=lapack::geev_lwork('N','V',N);
    T* work=(T*) alloca(LWORK*sizeof(T));
    Mat<T,M,N> a=_a;
    Mat<T,M,2> w;
    int info=lapack::geev('N','V',N,a.data(),N,w.data(),w.data()+M,
                          nullptr,1,_v.data(),N,work,LWORK);
    assert(info>=0 && "invalid argument"!=0);
    if (info!=0)
      w.data()[0]=std::numeric_limits<T>::quiet_NaN();
    return w;
  }

};

//-----------------------------------------------------------------------------

//# ifndef DOXYGEN_SKIP

//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------


//# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
} // namespace math
} // namespace VC
//-----------------------------------------------------------------------------

/** \def VC_MAT_NO_SPECIALIZATIONS
    \ingroup vc_math_lam
    Don't use specializations for certain matrix types or dimensions [\ref vc_math_lam].
*/

# ifdef DOXYGEN_SKIP
#  define VC_MAT_NO_SPECIALIZATIONS "this is an external switch"
#  error "doxygen only"
# endif

// # ifndef DOXYGEN_SKIP

/** \defgroup vc_math_lam_special Specializations for VC::math::Mat
    \ingroup vc_math_lam
 */

# include "mat2_special.hh"
# include "mat3_special.hh"
# include "mat4_special.hh"

// # endif // DOXYGEN_SKIP

# undef m_raw

//-----------------------------------------------------------------------------

#endif // __VC_MATH_MAT_HH__
