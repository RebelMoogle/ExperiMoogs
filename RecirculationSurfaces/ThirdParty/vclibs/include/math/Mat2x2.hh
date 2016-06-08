//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: VecN.hh 318 2010-04-29 14:21:51Z roessl $
// $Revision$
//
//=============================================================================

#ifndef __VC_MATH_MAT2X2_HH__
#define __VC_MATH_MAT2X2_HH__

#include "VecN.hh"
#include "blas_matrix.hh"

namespace VC {
namespace math {

//=============================================================================

/** Solve linear system A*x=b for 2x2 matrix A.
    \ingroup vc_math_lam
    \param[in] _a 2x2 matrix (column storage, leading dimension is 2)
    \param[in,out] _b right hand side b on input, solution x on output
    (dimension 2x_n)
    \param _n number of right hand sides
    \param _ld leading dimension of _b
    \return true on success, false if matrix _a is singular
    \sa Mat2x2
    \deprecated prefer VC::math::Mat 
 */
template <typename T>
bool solve2x2(const T* _a,T* _b,int _n=1,int _ld=2) {
  double a11=_a[0], a21=_a[1], a12=_a[2], a22=_a[3];

  double d=a11*a22-a21*a12;

  if (fabs(d)<std::numeric_limits<T>::epsilon())
    return false;

  for (int i=0;i<_n;++i) {
      double u=_b[_ld*i], v=_b[_ld*i+1];

      _b[_ld*i  ]=(u*a22-v*a12)/d;
      _b[_ld*i+1]=(a11*v-a21*u)/d;
  }

  return true;
}

/** Compute eigenvalues and eigenvectors of symmetric 2x2 matrix.
    \ingroup vc_math_lam
    Computes roots (real_roots_2()) of characteristic polynomial
    and solves for null space.
    \param[in] _a 2x2 \b symmetric matrix, symmetry is \a not checked,
    only \a upper triangle is referenced
    \param[out] _w store 2 eigenvalues (ascending values)
    \param[out] _e store 2x2 matrix with eigenvectors to _w as columns
    (for _e==0, no eigenvectors are computed)
    \deprecated prefer VC::math::Mat 
 */
template <typename T>
void eig2x2(const T* _a,T* _w,T* _e=0) {   

  // diagonal matrix
  if (fabs(_a[2])==T(0)) {
    _w[0]=std::min(_a[0],_a[3]);
    _w[1]=std::max(_a[0],_a[3]);
    if (_e!=0) {
      _e[0]=T(1); _e[1]=T(0);
      _e[2]=T(0); _e[3]=T(1);      
    }
    return;
  }

  double a11=_a[0], a12=_a[2], a22=_a[3];  


  double tr=a11+a22, det=a11*a22-a12*a12;
  double discr=tr*tr-4.0*det;
  
  assert(discr>=0.0);

  discr=sqrt(discr);

  double x[2];
  // Numerical Recipes (quadratic and cubic equations) 
  //  (instead 0f x=(tr +/- discr)*0.5;)
  double q=-0.5*(-tr+(-tr>0.0 ? +1.0 : -1.0)*discr);
  x[0]=q;
  x[1]=det/q;

  if (x[0]>x[1])
    std::swap(x[0],x[1]);

  _w[0]=x[0];
  _w[1]=x[1];


  // see also (I use opposite sign [a11-x,x12; a12,a22-x])
  // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html

  if (_e==0)
    return;
  
  T e[2],len;
  e[0]=a22-x[0];
  e[1]=-a12;
  
  len=sqrt(e[0]*e[0]+e[1]*e[1]);
  e[0]/=len;
  e[1]/=len;
  
  _e[0]=e[0];
  _e[1]=e[1];

  _e[2]= _e[1]; // rotate by -pi/2 (lapack chooses clockwise)
  _e[3]=-_e[0];
}

//-----------------------------------------------------------------------------

/** \defgroup vc_math_lam_old Linear algebra: matrices
    \ingroup vc_math
    Data types, arithmetic and utility functions for
    \arg Mat2x2 defines a 2x2 matrix data type, (see also
    VecN <T,N> and other small, quadratic matrix types).
    \arg blas::matrix_reference_t defines an interface to using BLAS
    operations on various types of matrices.
    \arg blas::GE_Matrix <T,alloc> defines a variable size matrix
    based on matrix_reference_t (see also blas::SP_Matrix,
    blas::TP_Matrix).
    \arg blas::GE_Mat_MxN defines a fixed size MxN matrix and
    matrix_reference_t types (see also blas::SP_Mat,
    blas::TP_Mat).

    \sa  [\ref vc_blas], [\ref vc_vcblas], [\ref vc_math_lav]
    \deprecated prefer VC::math::Mat 
*/

/** \class VC::math::Mat2x2 Mat2x2.hh
    \brief 2x2 matrices.

    \ingroup vc_math_lam

    \todo DOCUMENTATION

    \sa VC::math::VecN, [\ref vc_math_lam], [\ref vc_math_lav], [\ref vc_blas]
*/
template <typename T>
class Mat2x2 {
public:
  typedef T value_type;       //!< scalar type

  typedef Mat2x2<T> Self;     //!< synonym for VecN<T,K>

  enum Dimensions {
    ROWS=2,  //!< number of rows (=2=n_rows())
    COLS=2,  //!< number of columns (=2=n_cols())
    M=2,     //!  M=rows, synonym for ROWS
    N=2,     //!< N=COLS, synonym for COLS
    SIZE=4   //!< number of elements
  };

  /// matrix type (for [\ref vc_blas])
  typedef VC::math::blas::ge_mat<VC::math::blas::NoT,2,2,2> matrix_t;
  /// matrix reference
  typedef VC::math::blas::matrix_const_reference_t<T,matrix_t> const_reference_t;
  /// matrix reference
  typedef VC::math::blas::matrix_reference_t<T,matrix_t> reference_t;

  /** @name initialization
        @{
    */

  Mat2x2() {} //!< uninitialized matrix

  Mat2x2(const Self& _v) { memcpy(m_m,_v.m_m,sizeof(*this)); } //!< copy constructor

  /// elementwise copy with type cast
  template <typename S>
  Mat2x2(const Mat2x2<S>& _m) {
    m_m[0]=_m.m_m[0]; m_m[1]=_m.m_m[1]; m_m[2]=_m.m_m[2]; m_m[3]=_m.m_m[3];
  }
  /// set all elements to _s
  explicit Mat2x2(T _s) { for (unsigned i=0;i<SIZE;++i) m_m=_s; }

  /// initialize from pointer (4 values stored column-wise)
  template <typename S>
  explicit Mat2x2(const S* _m) {
    m_m[0]=_m[0]; m_m[1]=_m[1]; m_m[2]=_m[2]; m_m[3]=_m[3];
  }
  /// initialize from column vectors
  template <typename U,typename W>
  Mat2x2(const VecN<U,2>& _a,const VecN<W,2>& _b) {
    m_m[0]=_a[0]; m_m[1]=_a[1];
    m_m[2]=_b[0]; m_m[3]=_b[1];
  }
  /// initialize diagonal matrix
  template <typename S>
  explicit Mat2x2(const VecN<S,2>& _diag) {
    m_m[0]=_diag[0]; m_m[1]=m_m[2]=T(0); m_m[3]=_diag[1];
  }

  /// (*this)=_v
  Self& operator=(const Self& _m) {
    memcpy(m_m,_m.m_m,sizeof(*this)); return *this;
  }

  /// (*this)=_v
  template <typename S>
  Self& operator=(const Mat2x2<S>& _m) {
    for (unsigned i=0;i<SIZE;++i) m_m[i]=_m.data()[i];
    return *this;
  }

  /// assign matrix_const_reference_t
  template <typename A>
  Self& operator=(const blas::matrix_const_reference_t<T,A>& _m) {
    ref()=_m;
    return *this;
  }
  /// assign matrix_reference_t
  template <typename A>
  Self& operator=(const blas::matrix_reference_t<T,A>& _m) {
    ref()=_m;
    return *this;
  }

  /// @}

  /** @name dimensions
        @{
    */

  unsigned m() const { return M; }        //!< number of rows (=M)
  unsigned n() const { return N; }        //!< number of columns (=N)
  unsigned n_rows() const { return M; }   //!< number of rows (=M)
  unsigned n_cols() const { return N; }   //!< number of columns (=N)
  unsigned size() const { return SIZE; }  //!< number of elements (=SIZE)

  /// @}

  /** @name access and assignment
      If \c VC_MATH_CHECK_BOUNDS is defined then
      VecN::operator[]() and function nth() include an assertion on index
      within bounds.
      @{
  */

  T* data() { return m_m; }              //!< get raw data
  const T* data() const { return m_m; }  //!< get raw data

  T* begin() { return m_m; }              //!< synonym for data()
  const T* begin() const { return m_m; }  //!< synonym for data()

  VecN<T,N>* columns() { return (VecN<T,N>*) m_m; } //!< get data() as VecN
  const VecN<T,N>* columns() const { return (const VecN<T,N>*)m_m; } //!< get data() as VecN

  T* end() { return m_m+SIZE; }              //!< synonym for data()+SIZE
  const T* end() const { return m_m+SIZE; }  //!< synonym for data()+SIZE

# ifdef VC_MATH_CHECK_BOUNDS
  /// access element
  T& operator()(int _i,int _j)  {
    assert(0<=i && i<M && 0<=j && _j<N); return m_m[_j*2+_i];
  }
  /// access element
  const T& operator()(int _i,int _j) const {
    assert(0<=i && i<M && 0<=j && _j<N); return m_m[_j*2+_i];
  }

  /// get row vector
  VecN<T,M> row(int _i) const {
    assert(0<=_i && _i<M); return VecN<T,M>(m_m[_i],m_m[_i+M]);
  }
  /// access column vector
  VecN<T,N>& column(int _j) {
    assert(0<=_j && _j<N); return *(m_m+2*_j);
  }
  /// access column vector
  const VecN<T,N>& column(int _j) const {
    assert(0<=_j && _j<N); return *(m_m+2*_j);
  }

# else
  T& operator()(int _i,int _j)  { return m_m[_j*2+_i]; }
  const T& operator()(int _i,int _j) const { return m_m[_j*2+_i]; }
  VecN<T,M> row(int _i) const { return VecN<T,M>(m_m[_i],m_m[_i+M]); }
  VecN<T,N>& column(int _j) { return *((VecN<T,N>*) (m_m+2*_j)); }
  const VecN<T,N>& column(int _j) const { return *((VecN<T,N>*) m_m+2*_j); }
#endif

  /** @name matrix_reference_t and pointer casts
      \sa [\ref vc_blas]
      @{
  */

  /// get matrix reference
  reference_t ref() { return reference_t(matrix_t(),m_m); }
  /// get matrix reference
  const_reference_t const_ref() const { return const_reference_t(matrix_t(),m_m); }

  /// cast C-array (also defined externally as to_m(p))
  static Self&
  to_2x2(T* _p) { return (Self&) *_p; }

  /// cast const C-array (also defined externally as to_cm(p))
  static const Self&
  to_c2x2(const T* _p) { return *((const Self*) _p); }

  // @}

  /** @name special initializations
      @{
  */

  /// return transposed
  Self trans() const { T t[]={m_m[0],m_m[2],m_m[1],m_m[3]}; return Self(t); }
  /// transpose *this
  Self& transpose() { std::swap(m_m[1],m_m[2]); return *this; }
  /// add transposed
  Self& add_trans() {
    T x=m_m[1]; m_m[1]+=m_m[2]; m_m[2]+=x;
    m_m[0]*=T(2); m_m[3]*=T(2);
    return *this;
  }

  /// load outer product x'*x
  template <typename S>
  Self& outer(const VecN<S,M>& _x) {
    m_m[0]=_x[0]*_x[0]; m_m[1]=m_m[2]=_x[0]*_x[1]; m_m[3]=_x[1]*_x[1];
    return *this;
  }
  /// add outer product x'*x
  template <typename S>
  Self& add_outer(const VecN<S,M>& _x) {
    return add_outer(T(1),_x);
  }
  /// add outer product x'*x*_alpha
  template <typename S>
  Self& add_outer(T _alpha,const VecN<S,M>& _x) {
    m_m[0]+=_x[0]*_x[0]*_alpha;
    m_m[1]+=_x[1]*_x[0]*_alpha;
    m_m[2]+=_x[0]*_x[1]*_alpha;
    m_m[3]+=_x[1]*_x[1]*_alpha;
    return *this;
  }
  /// add outer product x'*y
  template <typename S>
  Self& add_outer(const VecN<S,M>& _x,const VecN<S,M>& _y) {
    return add_outer(T(1),_x,_y);
  }
  /// add outer product x'*x*_alpha
  template <typename S>
  Self& add_outer(T _alpha,const VecN<S,M>& _x,const VecN<S,M>& _y) {
    m_m[0]+=_x[0]*_y[0]*_alpha;
    m_m[1]+=_x[1]*_y[0]*_alpha;
    m_m[2]+=_x[0]*_y[1]*_alpha;
    m_m[3]+=_x[1]*_y[1]*_alpha;
    return *this;
  }

  /// load identity matrix
  Self& eye() const {
    m_m[0]=m_m[3]=T(1); m_m[1]=m_m[2]=T(0); return *this;
  }
  /// load zeros
  Self& zeros() { return all(T(0)); }
  // load ones
  Self& ones() { return all(T(1)); }
  /// load all values with _a
  Self& all(T _a) {
    m_m[0]=m_m[1]=m_m[2]=m_m[3]=_a; return *this;
  }
  /// load random values
  Self& rand() { rand(m_m,SIZE); return *this; }

  /// load rotation matrix [cos(phi),-sin(phi);sin(phi),cos(phi)[
  Self& rot(T _phi) {
    m_m[0]=cos(_phi); m_m[1]=sin(_phi);
    m_m[2]=-m_m[1]; m_m[3]=m_m[0];
    return *this;
  }

  /// @}

   /** @name arithmetic: evaluate functions elementwise
       @{
    */

  /// replace each element a(i,j) by _op(a(i,j))
  template <typename Op>
  Self& map(Op& _op) {
    for (unsigned i=0;i<SIZE;++i)
      m_m[i]=_op(m_m[i]);
    return *this;
  }
  /// replace each element (*this)(i,j) by _op(_m(i,j))
  template <typename Op,typename S>
  Self& map(Op& _op,const Mat2x2<S>& _m) {
    for (unsigned i=0;i<SIZE;++i)
      m_m[i]=_op(_m.data()[i]);
    return *this;
  }
  /// replace each element (*this)(i,j) by _op(_a(i,j),_b(i,j))
  template <typename Op,typename A,typename B>
  Self& map(Op& _op,const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
    for (unsigned i=0;i<SIZE;++i)
      m_m[i]=_op(_a.data()[i],_b.data()[i]);
    return *this;
  }

  //@}

  /** @name sums, extrema, and permutation.
      @{
  */

  /// synonym for sum_cols()
  VecN<T,N> sum() const { return VecN<T,N>(m_m[0]+m_m[1],m_m[2]+m_m[3]); }
  /// return sum of columns
  VecN<T,N> sum_cols() const { return sum(); }
  /// return sum of rows
  VecN<T,M> sum_rows() const { return VecN<T,M>(m_m[0]+m_m[2],m_m[1]+m_m[3]); }

  /// synonym for mean_rows()
  VecN<T,N> mean() const { return sum()/T(M); }
  /// return means of columns
  VecN<T,N> mean_cols() const { return sum()/T(M); }
  /// return means of rows
  VecN<T,M> mean_rows() const { return sum_rows()/T(N); }

  /// synomym for min_cols()
  VecN<T,N> min() const {
    return VecN<T,N>(std::min(m_m[0],m_m[1]),std::min(m_m[2],m_m[3]));
  }
  /// get minima of columns
  VecN<T,N> min_cols() const { return min(); }
  /// get minima of rows
  VecN<T,M> min_rows() const {
    return VecN<T,M>(std::min(m_m[0],m_m[2]),std::min(m_m[1],m_m[3]));
  }

  /// synomym for max_cols()
  VecN<T,N> max() const {
    return VecN<T,N>(std::max(m_m[0],m_m[1]),std::max(m_m[2],m_m[3]));
  }
  /// get maxima of columns
  VecN<T,N> max_cols() const { return max(); }
  /// get maxima of rows
  VecN<T,M> max_rows() const {
    return VecN<T,M>(std::max(m_m[0],m_m[2]),std::max(m_m[1],m_m[3]));
  }

  /// synomym for imin_cols()
  VecN<unsigned,N> imin() const {
    return VecN<unsigned,N>(m_m[0]<=m_m[1] ? 0 : 1,m_m[2]<=m_m[3] ? 0 : 1);
  }
  /// get indices of minima of columns
  VecN<unsigned,N> imin_cols() const { return imin(); }
  /// get indices of minima of rows
  VecN<unsigned,M> imin_rows() const {
    return VecN<unsigned,N>(m_m[0]<=m_m[2] ? 0 : 1,m_m[1]<=m_m[3] ? 0 : 1);
  }

  /// synomym for imax_cols()
  VecN<unsigned,N> imax() const {
    return VecN<unsigned,N>(m_m[0]>=m_m[1] ? 0 : 1,m_m[2]>=m_m[3] ? 0 : 1);
  }
  /// get indices of maxima of columns
  VecN<unsigned,N> imax_cols() const { return imax(); }
  /// get indices of maxima of rows
  VecN<unsigned,M> imax_rows() const {
    return VecN<unsigned,N>(m_m[0]>=m_m[2] ? 0 : 1,m_m[1]>=m_m[3] ? 0 : 1);
  }

  /// select from index vector _idx: a(idx,:)
  template <typename IDX>
  VecN<T,N> select_in_cols(const VecN<IDX,N>& _idx) const {
    return VecN<T,N>((*this)(_idx[0],0),(*this)(_idx[1]),1);
  }
  /// select from index vector _idx: a(:,idx)
  template <typename IDX>
  VecN<T,M> select_in_rows(const VecN<IDX,M>& _idx) const {
    return VecN<T,M>((*this)(0,_idx[0]),(*this)(0,_idx[1]));
  }

  // @}

  /** @name matrix-scalar operations
      @{
   */

  /// multiply by scalar _alpha
  Self& operator*=(T _alpha) {
    for (int i=0;i<SIZE;++i) m_m[i]*=_alpha; return *this;
  }
  /// divide by scalar _alpha
  Self& operator/=(T _alpha) {
    for (int i=0;i<SIZE;++i) m_m[i]/=_alpha; return *this;
  }
  /// add _alpha to all elements
  Self& operator+=(T _alpha) {
    for (int i=0;i<SIZE;++i) m_m[i]+=_alpha; return *this;
  }
  /// subtract _alpha from all elements
  Self& operator-=(T _alpha) {
    for (int i=0;i<SIZE;++i) m_m[i]-=_alpha; return *this;
  }
  /// elementwise a(i,j) -> _alpha/a(_i,_j)
  Self& rdiv(T _alpha=T(1)) {
    for (int i=0;i<SIZE;++i) m_m[i]=_alpha/m_m[i]; return *this;
  }

  /// @}

  /** @name matrix-vector operations
      @{
   */

  /// matrix-vector multiplication

  template <typename S>
  VecN<typename Result<S,T>::type,M>
  operator*(const VecN<S,M>& _v) const {
    return VecN<typename Result<S,T>::type,M>
      (m_m[0]*_v.data()[0]+m_m[2]*_v.data()[1],
       m_m[1]*_v.data()[0]+m_m[3]*_v.data()[1]);
  }

  /// @}


  /** @name matrix-matrix operations
      \ref vc_math_lam_mops
      @{
  */

  /// \a element-wise multiplication
  template <typename S>
  Self& operator*=(const Mat2x2<S>& _m) {
    for (int i=0;i<SIZE;++i) m_m[i]*=_m.data()[i]; return *this;
  }
  /// \a element-wise division
  template <typename S>
  Self& operator/=(const Mat2x2<S>& _m) {
    for (int i=0;i<SIZE;++i) m_m[i]/=_m.data()[i]; return *this;
  }
  /// matrix addition
  template <typename S>
  Self& operator+=(const Mat2x2<S>& _m) {
    for (int i=0;i<SIZE;++i) m_m[i]+=_m.data()[i]; return *this;
  }
  /// matrix subtraction
  template <typename S>
  Self& operator-=(const Mat2x2<S>& _m) {
    for (int i=0;i<SIZE;++i) m_m[i]-=_m.data()[i]; return *this;
  }

  /// square all elements
  Self& sqr()  { for (int i=0;i<SIZE;++i) m_m[i]*=m_m[i]; return *this;  }

  /// matrix multiplication (returns *this|_m)
  template <typename S>
  Self mult(const Mat2x2<S>& _m) const {
    const S* a=_m.data();
    T m[SIZE]={m_m[0]*a[0]+m_m[2]*a[1],
               m_m[1]*a[0]+m_m[3]*a[1],
               m_m[0]*a[2]+m_m[2]*a[3],
               m_m[1]*a[2]+m_m[3]*a[3]};
    return Self(m);
  }

  /// swap contents of (*this) and _m
  template <typename S>
  inline Self& swap(Mat2x2<S>& _m)  {
    FixedAryOps<SIZE>::swap(m_m,_m.data()); return *this;
  }

  /// return trans(*this)|*this
  Self ata() const {
    T a=m_m[0]*m_m[0]*m_m[1]*m_m[1], b=m_m[0]*m_m[2]+m_m[1]*m_m[3], c=m_m[2]*m_m[2]+m_m[3]*m_m[3];
    return Self(a,b,b,c);
  }

  /// @}

  /** \defgroup vc_math_lam_mops_old Linear algebra: externally defined operations on Mat2x2
      \ingroup vc_math_lam
      \sa Mat2x2
   */

  /** @name Externally defined operations on Mat2x2
      @{
  */

  /// matrix + scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator+(const Mat2x2<A>& _a,const double& _b);

  /// matrix - scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator-(const Mat2x2<A>& _a,const double& _b);

  /// matrix + scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator+(const double& _b,const Mat2x2<A>& _a);

  /// matrix - scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator-(const double& _b,const Mat2x2<A>& _a);

  /// matrix * scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator*(const Mat2x2<A>& _a,const double& _b);

  /// matrix / scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator/(const Mat2x2<A>& _a,const double& _b);

  /// matrix * scalar
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator*(const double& _b,const Mat2x2<A>& _a);

  /// scalar/ matrix
  template <typename A>
  friend Mat2x2<typename Result<A,double>::type>
  operator/(const double& _b,const Mat2x2<A>& _a);

  /// matrix + matrix
  template <typename A,typename B>
  friend Mat2x2<typename Result<A,B>::type>
  operator+(const Mat2x2<A>& _a,const Mat2x2<B>& _b);

  /// matrix - matrix
  template <typename A,typename B>
  friend Mat2x2<typename Result<A,B>::type>
  operator-(const Mat2x2<A>& _a,const Mat2x2<B>& _b);

  /// \a element-wise matrix * matrix
  template <typename A,typename B>
  friend Mat2x2<typename Result<A,B>::type>
  operator*(const Mat2x2<A>& _a,const Mat2x2<B>& _b);

  /// \a element-wise matrix / matrix
  template <typename A,typename B>
  friend Mat2x2<typename Result<A,B>::type>
  operator/(const Mat2x2<A>& _a,const Mat2x2<B>& _b);

  /// matrix-matrix multiplication
  template <typename A,typename B>
  friend Mat2x2<typename Result<A,B>::type>
  operator|(const Mat2x2<A>& _a,const Mat2x2<B>& _b);

  /// (row) vector-matrix multiplication
  template <typename A,typename B>
  friend VecN<typename Result<A,B>::type,2>
  operator*(const VecN<A,2>& _a,const Mat2x2<B>& _b);

  /// alternative to std::swap(_a,_b)
  template <typename A>
  friend void swap(Mat2x2<A>& _a,Mat2x2<A>& _b);

  /// matrix from row vectors _r0,r1
  template <typename T0,typename T1>
  friend Mat2x2<typename Result<T0,T1>::type>
  v_cat(const VecN<T0,2>& _r0,const VecN<T1,2>& _r1);

  /// matrix from column vectors _c0,c1
  template <typename T0,typename T1>
  friend Mat2x2<typename Result<T0,T1>::type>
  h_cat(const VecN<T0,2>& _c0,const VecN<T1,2>& _c1);

  /// cast C-array (also defined externally as to_m(p))
  template <typename S>
  friend Mat2x2<S>& to_2x2(T* _p);

  /// cast const C-array (also defined externally as to_cm(p))
  template <typename S>
  friend const Mat2x2<S>& to_c2x2(const T* _p);

  /// solve linear system (see solve2x2())
  template <typename S>
  friend bool solve(const Mat2x2<S>& _A,VecN<S,2>& _x,int _n);

  /// Output matrix row-wise
  template <typename S>
  friend std::ostream& operator<<(std::ostream& _out,const Mat2x2<S>& _m);

  /// Input matrix row-wise
  template <typename S>
  friend std::istream& operator>>(std::istream& _in,Mat2x2<S>& _m);


  /// @}

  /** @name Matrix properties
      Note: some properties (norm2(),cond2(),rank()) use VC::math::blas::svd()
      which yields precise results but is rather slow. (Taking eigenvalues
      of A'*A produces unacceptable rounding errors.)
      @{
   */

  /// get determinant
  T det() const { return m_m[0]*m_m[3]-m_m[1]*m_m[2]; }
  /// get determinant as \c double
  double ddet() const { return double(m_m[0])*m_m[3]-double(m_m[1])*m_m[2]; }

  /// get trace
  T trace() const { return m_m[0]+m_m[3]; }

  /// synomym for sum().is_finite() \sa VecN::is_finite()
  bool is_finite() const {
    return sum().is_finite();
  }

  /// Frobenius norm
  T normf() const {
    double a=0.0;
    for (int i=0;i<SIZE;++i) a+=m_m[i]*m_m[i];
    return sqrt(a);
  }
  /// maximum norm
  T norminf() const {
    T a=fabs(m_m[0]);
    for (int i=1;i<SIZE;++i) if (fabs(m_m[i])>a) a=fabs(m_m[i]);
    return a;
  }
  /// 2-norm (calls VC::math::blas::norm2() using SVD)
  T norm2() const { return blas::norm2(const_ref()); }
  /// condition w.r.t. to norm2() (uses SVD)
  T cond2() const { return blas::cond2(const_ref()); }
  /// get effective rank (calls VC::math::blas::rank() using SVD)
  int rank() const { return blas::rank(const_ref()); }

  /// compute eigenvalues _w of \a symmetric matrix (see eig2x2())
  void eig(VecN<T,2>& _w) const {
    eig2x2(m_m,_w.data());
  }
  /// compute eigenvalues _w and eigenvectors _e of \a symmetric matrix (see eig2x2())
  void eig(Self& _e,VecN<T,2>& _w) const {
    eig2x2(m_m,_w.data(),_e.m_m);
  }

  /// @}

  /** Solve linear system (*this)*x=n.
      Same as VC::math::solve(const Mat2x2<S>& _A,VecN<S,2>& _x,int _n).
  */
  bool solve(VecN<T,2>& _x,int _n=1) {
    return solve2x2(m_m,_x.data(),_n);
  }
  

private:
  union {
    T               m_m[4];        //!< data fields
    VecN_align<T,4> m_align;       //!< force alignment
  };
};

  // list initialization

  //  | 2xn matrix

  // qr (yields angle)



//-----------------------------------------------------------------------------

/// matrix + scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator+(const Mat2x2<A>& _a,const double& _b) {
  Mat2x2<typename Result<A,double>::type> c(_a); c+=_b;
  return c;
}

/// matrix + scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A> 
Mat2x2<typename Result<A,float>::type>
operator+(const Mat2x2<A>& _a,const float& _b) {
  Mat2x2<typename Result<A,float>::type> c(_a); c+=_b;
  return c;
}

/// matrix - scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator-(const Mat2x2<A>& _a,const B& _b) {
  Mat2x2<typename Result<A,B>::type> c(_a); c-=_b;
  return c;
}

/// matrix + scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator+(const B& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,B>::type> c(_a); c+=_b;
  return c;
}

/// matrix - scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator-(const double& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,double>::type> c(_a); c-=_b;
  return c;
}

/// matrix - scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat 
template <typename A>
Mat2x2<typename Result<A,float>::type>
operator-(const float& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,float>::type> c(_a); c-=_b;
  return c;
}


/// matrix * scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator*(const Mat2x2<A>& _a,const double& _b) {
  Mat2x2<typename Result<A,double>::type> c(_a); c*=_b;
  return c;
}
/// matrix * scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,float>::type>
operator*(const Mat2x2<A>& _a,const float& _b) {
  Mat2x2<typename Result<A,float>::type> c(_a); c*=_b;
  return c;
}

/// matrix / scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator/(const Mat2x2<A>& _a,const double& _b) {
  Mat2x2<typename Result<A,double>::type> c(_a); c/=_b;
  return c;
}

/// matrix / scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,float>::type>
operator/(const Mat2x2<A>& _a,const float& _b) {
  Mat2x2<typename Result<A,float>::type> c(_a); c/=_b;
  return c;
}

/// matrix * scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator*(const double& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,double>::type> c(_a); c*=_b;
  return c;
}

/// matrix * scalar \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,float>::type>
operator*(const float& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,float>::type> c(_a); c*=_b;
  return c;
}

/// scalar/ matrix  \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,double>::type>
operator/(const double& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,double>::type> c(_a); c/=_b;
  return c;
}

/// scalar/ matrix  \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A>
Mat2x2<typename Result<A,float>::type>
operator/(const float& _b,const Mat2x2<A>& _a) {
  Mat2x2<typename Result<A,float>::type> c(_a); c/=_b;
  return c;
}

/// matrix + matrix \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator+(const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
  Mat2x2<typename Result<A,B>::type> c(_a); c+=_b;
  return c;
}

/// matrix - matrix \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator-(const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
  Mat2x2<typename Result<A,B>::type> c(_a); c-=_b;
  return c;
}

/// \a element-wise matrix * matrix \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator*(const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
  Mat2x2<typename Result<A,B>::type> c(_a); c*=_b;
  return c;
}

/// \a element-wise matrix / matrix \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator/(const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
  Mat2x2<typename Result<A,B>::type> c(_a); c/=_b;
  return c;
}


/// matrix-matrix multiplication \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
Mat2x2<typename Result<A,B>::type>
operator|(const Mat2x2<A>& _a,const Mat2x2<B>& _b) {
  return _a.mult(_b);
}

/// (row) vector-matrix multiplication \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename A,typename B>
VecN<typename Result<A,B>::type,2>
operator*(const VecN<A,2>& _a,const Mat2x2<B>& _b) {
  return VecN<typename Result<A,B>::type,2>
    (_a.data()[0]*_b.data()[0]+_a.data()[1]*_b.data()[1],
     _a.data()[0]*_b.data()[2]+_a.data()[1]*_b.data()[3]);
}

/// alternative to std::swap(_a,_b) \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T>
void swap(Mat2x2<T>& _a,Mat2x2<T>& _b) {
  FixedAryOps<Mat2x2<T>::SIZE>::swap(_a.m_m,_b.m_m);
}

/*

  interferes with new Mat<T,M,N>

/// load matrix from row vectors _r0,r1  \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T0,typename T1>
Mat2x2<typename Result<T0,T1>::type>
v_cat(const VecN<T0,2>& _r0,const VecN<T1,2>& _r1) {
  return Mat2x2<typename Result<T0,T1>::type>
    (_r0.data()[0],_r1.data()[0],_r0.data()[1],_r1.data()[1]);
}

/// load matrix from column vectors _c0,c1  \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T0,typename T1>
Mat2x2<typename Result<T0,T1>::type>
v_cat(const VecN<T0,2>& _c0,const VecN<T1,2>& _c1) {
  return Mat2x2<typename Result<T0,T1>::type>
    (_c0.data()[0],_c0.data()[1],_c1.data()[0],_c1.data()[1]);
}
*/

/// cast C-array (also defined externally as to_m(p))  \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T>
Mat2x2<T>& to_2x2(T* _p) { return *((Mat2x2<T>*) _p); }

/// cast const C-array (also defined externally as to_cm(p)) \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T>
const Mat2x2<T>& to_c2x2(const T* _p) { return *((const Mat2x2<T>*) _p); }

/** Solve linear system A*x=b.
    \ingroup vc_math_lam_mops
    \param[in] _A matrix
    \param[in,out] _x right hand side (input) and solution (output)
    \param _n number of right hand sides
    \return false if _A is singular, true else
    \sa solve2x2()
    \deprecated prefer VC::math::Mat
 */
template <typename S>
bool solve(const Mat2x2<S>& _A,VecN<S,2>& _x,int _n=1) {
  return solve2x2(_A.data(),_x.data(),_n);
}

/// Output matrix row-wise \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T>
std::ostream& operator<<(std::ostream& _out,const Mat2x2<T>& _m) {
  return _out << _m(0,0) << ' '  << _m(0,1) << ' '
              << _m(1,0) << ' '  << _m(1,1);
}

/// Input matrix row-wise \ingroup vc_math_lam_mops \deprecated prefer VC::math::Mat
template <typename T>
std::istream& operator>>(std::istream& _in,Mat2x2<T>& _m) {
  return _in >> _m(0,0) >> _m(0,1) >>_m(1,0) >> _m(1,1);
}

//=============================================================================

} // namespace math
} // namespace VC

#endif // __VC_MATH_MAT2X2_HH__
