//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id: blas.hh 105 2009-10-14 18:18:57Z roessl $
// $Revision$
//
//=============================================================================

/** \file

    This file is included by blas_matrix.hh.

    DO NOT INCLUDE IT DIRECTLY!

    \arg Provide classes GE_Matrix, etc. (derived from Matrix_alloc
    and matrix_reference_t) which allocate storage for variable size
    matrices.

    \internal

    \todo mixins: access, load(), assign self_t -- copyable, loadable,
    vector_accessible, matrix_accessible, ge_matrix_accessible,...
    \todo delegates for VecN, MatMxN
    \todo evaluate to temporary, e.g., cout << pmatlab(eval(A+B))
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** \class Matrix_alloc
    (used internally: allocate storage for variable size matrices)
    \ingroup vc_blas_matrices
 */
template <typename allocator>
class Matrix_alloc {
public:
  Matrix_alloc() : m_capacity(0) {}

  void reserve(unsigned _capacity,
               typename allocator::pointer& _data) {
    if (m_capacity<_capacity) {
      if (_data!=0) m_alloc.deallocate(_data,m_capacity);
      m_capacity=_capacity;
      _data=m_alloc.allocate(_capacity);
      assert(_data!=0);
    }
  }
  void free(typename allocator::pointer& _data) {
    if (_data!=0) {
      m_alloc.deallocate(_data,m_capacity);
      _data=0;
      m_capacity=0;
    }
  }

  /// query capacity (number of elements for which storage is allocated)
  unsigned capacity() const { return m_capacity; }

protected:
  allocator   m_alloc;
  unsigned    m_capacity;
private:
  Matrix_alloc(const Matrix_alloc&); //!< non_copyable
};

//-----------------------------------------------------------------------------

/** \class Vector
    Vector (see vector_reference_t and vector_const_reference_t).
    \ingroup vc_blas_matrices
 */
template <typename T,typename allocator=std::allocator<T> >
class Vector : public vector_reference_t<T,vec<VarInt,1> >,
               public Matrix_alloc<allocator> {
public:
  typedef Vector<T,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef vec<VarInt,1> vector_t;
  typedef vector_reference_t<value_t,vector_t> reference_t;
  typedef allocator allocator_t;

  Vector()
    : reference_t(vector_t(0,1),0) {}
  Vector(int _n,unsigned _capacity=0)
    : reference_t(vector_t(0,1),0) {
    reserve(_capacity);
    resize(_n);
  }
  /// copy ctor calls load()
  Vector(const self_t& _x) : reference_t(vector_t(0,1),0) {
    this->load(_x.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~Vector() { this->free(this->m_data); /* not global std::free()!!!*/  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by Vector"); return *this; }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by Vector"); return *this; }
   /// forbidden
  const reference_t& set_inc(int); // { assert(!"data is managed by Vector"); return *this; }

public:

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize length
  void resize(int _n) {
    this->m_v=vector_t(_n,1);
    this->reserve(this->n());
  }
  /// resize to same length as _other
  template <typename V>
  void resize(const vector_const_reference_t<T,V>& _other) {
    resize(_other.n());
  }
  /// load _other (resize() and vector_reference_t::copy())
  template <typename V>
  reference_t& load(const vector_const_reference_t<T,V>& _other) {
    resize(_other);
    if (_other.data()!=0)
      ((const reference_t*) this)->copy(_other);
    return *this;
  }

  /// get reference to n() x 1 matrix
  matrix_reference_t<T,ge_mat<NoT,VarInt,VarInt,VarInt> > as_matrix() {
    return mat<GE>::ref(this->m_data,this->n(),1,this->n());
  }
  /// get reference to n() x 1 matrix
  matrix_const_reference_t<T,ge_mat<NoT,VarInt,VarInt,VarInt> > as_matrix() const {
    return mat<GE>::ref(this->m_data,this->n(),1,this->n());
  }

  /// access element _i
  T& operator()(int _i) { return at(*this,_i); }
  /// access element _i
  T& operator[](int _i) { return at(*this,_i); }
  /// access element
  T& operator()(const idx_end& _i) { return at(*this,this->n()-_i.offset-1); }
  /// access element
  T& operator[](const idx_end& _i) { return at(*this,this->n()-_i.offset-1); }

  /// access element _i
  const T& operator()(int _i) const { return at(*this,_i); }
  /// access element _i
  const T& operator[](int _i) const { return at(*this,_i); }
  /// access element
  const T& operator()(const idx_end& _i) const { return at(*this,this->n()-_i.offset-1); }
  /// access element
  const T& operator[](const idx_end& _i) const { return at(*this,this->n()-_i.offset-1); }



  /// get reference to subvector
  vector_reference_t<T,vector_t> // VarInt
  operator()(int _i1,int _i2) const {
    return block(*this,_i1,_i2);
  }
  /// get reference to subvector
  template <bool END1,bool END2>
  vector_reference_t<T,vector_t> // VarInt
  operator()(const idx_range<END1,END2>& _ri) const {
    return block(*this,_ri);
  }

  /// delegate any assignment to vector_reference_t
  template <typename A>
  const reference_t& operator=(const A& _a) const {
    return reference_t::operator=(_a);
  }
};

//-----------------------------------------------------------------------------

/** \class GE_Matrix
    General matrix (see matrix_reference_t and matrix_const_reference_t).
    \ingroup vc_blas_matrices

    Use GE_matrix to allocate storage for tr_mat, sy_mat.
 */
template <typename T,typename allocator=std::allocator<T> >
class GE_Matrix : public matrix_reference_t<T,ge_mat<NoT,VarInt,VarInt,VarInt> >,
                  public Matrix_alloc<allocator> {
public:
  typedef GE_Matrix<T,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef ge_mat<NoT,VarInt,VarInt,VarInt> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef allocator allocator_t;

  typedef sy_mat<Upper,VarInt,VarInt> sy_upper_t; //!< sy_mat
  typedef sy_mat<Lower,VarInt,VarInt> sy_lower_t; //!< sy_mat
  typedef tr_mat<Upper,NoT,NoU,VarInt,VarInt> tr_upper_t; //!< tr_mat
  typedef tr_mat<Lower,NoT,NoU,VarInt,VarInt> tr_lower_t; //!< tr_mat
  typedef tr_mat<Upper,NoT,UnitDiag,VarInt,VarInt> tr_upper_u_t; //!< tr_mat (UnitDiag)
  typedef tr_mat<Lower,NoT,UnitDiag,VarInt,VarInt> tr_lower_u_t; //!< tr_mat (UnitDiag)

  GE_Matrix()
    : reference_t(matrix_t(0,0,0),0) {}
  GE_Matrix(int _m,int _n,int _ld=0,unsigned _capacity=0)
    : reference_t(matrix_t(0,0,0),0) {
    if (_capacity>0) reserve(_capacity);
    resize(_m,_n,_ld);
  }
  /// copy ctor calls load()
  GE_Matrix(const self_t& _a) : reference_t(matrix_t(0,0,0),0) {
    this->load(_a.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~GE_Matrix() { this->free(this->m_data); /* not global std::free()!!!*/  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by GE_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_m(int);   // { assert(!"data is managed by GE_Matrix"); return *this;  }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by GE_Matrix"); return *this;  }
  /// forbidden
  const reference_t& set_ld(int);  // { assert(!"data is managed by GE_Matrix"); return *this;  }

public:

  /** @name Access submatrix/block.
      For GE_Matrix operator()(i1,i2,j1,j2)/operator()(idx_range,idx_range)
      is equal to block(*this,...), i.e., always provides a reference to
      a submatrix.
      @{
  */

  /// access element
  T& operator()(int _i,int _j) const {
    return at(*this,_i,_j);
  }
  /// access element
  T& operator()(const idx_end& _i,int _j) const {
    return at(*this,this->m()-_i.offset-1,_j);
  }
  /// access element
  T& operator()(int _i,const idx_end& _j) const {
    return at(*this,_i,this->n()-_j.offset-1);
  }
  /// access element
  T& operator()(const idx_end& _i,const idx_end& _j) const {
    return at(*this,this->m()-_i.offset-1,this->n()-_j.offset-1);
  }

  /// get reference to submatrix
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(int _i1,int _i2,int _j1,int _j2) const {
    return block(*this,_i1,_i2,_j1,_j2);
  }
  /// get reference to submatrix A(i1:i2,j1:j2)
  template <bool END1,bool END2,bool END3,bool END4>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) const {
    return block(*this,_ri,_rj);
  }
  /// get reference to submatrix A(i1:i2,j) (IDX={int,idx_end})
  template <bool END1,bool END2,typename IDX>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const idx_range<END1,END2>& _ri,const IDX& _j) const {
    return block(*this,_ri,_j);
  }
  /// get reference to submatrix A(i,j1:j2) (IDX={int,idx_end})
  template <bool END1,bool END2,typename IDX>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const IDX& _i,const idx_range<END1,END2>& _rj) const {
    return block(*this,_i,_rj);
  }

  /// @}

  /** @name iterators
      @{
   */

  typedef GE_row_iterator<T> row_iterator_t; //!< iterate rows
  typedef GE_column_iterator<T> column_iterator_t; //!< iterate columns

  /// iterate rows
  row_iterator_t row_iter() {
    return GE_row_iterator<T>(this->data(),this->m(),this->n(),this->ld());
  }
  /// iterate columns
  column_iterator_t column_iterator() {
    return GE_column_iterator<T>(this->data(),this->m(),this->n(),this->ld());
  }

  /// @}

  /** @name Reference upper and lower triangle (as sy_mat, tr_mat).
      @{
  */

  /// get reference to lower triangle (as sy_mat)
  matrix_reference_t<T,sy_lower_t> sy_lower() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,sy_lower_t>(sy_lower_t(this->n(),this->ld()),this->m_data);
  }
  /// get reference to upper triangle (as sy_mat)
  matrix_reference_t<T,sy_upper_t> sy_upper() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,sy_upper_t>(sy_upper_t(this->n(),this->ld()),this->m_data);
  }
  /// get reference to lower triangle (as tr_mat)
  matrix_reference_t<T,tr_lower_t> tr_lower() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_lower_t>(tr_lower_t(this->n(),this->ld()),this->m_data);
  }
  /// get reference to upper triangle (as tr_mat)
  matrix_reference_t<T,tr_upper_t> tr_upper() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_upper_t>(tr_upper_t(this->n(),this->ld()),this->m_data);
  }
  /// get reference to lower triangle (as tr_mat w/ UnitDiag)
  matrix_reference_t<T,tr_lower_u_t> tr_lower_u() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_lower_u_t>(tr_lower_u_t(this->n(),this->ld()),this->m_data);
  }
  /// get reference to upper triangle (as tr_mat w/ UnitDiag)
  matrix_reference_t<T,tr_upper_t> tr_upper_u() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_upper_u_t>(tr_upper_u_t(this->n(),this->ld()),this->m_data);
  }

  /// @}

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize matrix dimensions and LD (default _ld=_m)
  void resize(int _m,int _n,int _ld=0) {
    if (_ld==0) _ld=_m;
    this->m_m=matrix_t(_m,_n,_ld);
    reserve(this->size());
  }
  /// resize to same dimensions as _other (_other.mrows(), _other.mcols(), sets LD=M)
  template <typename M>
  void resize(const matrix_const_reference_t<T,M>& _other) {
    resize(_other.m_rows(),_other.n_cols());
  }
  /// resize to same dimensions as _other (_other.mrows(), _other.mcols(), sets LD=M)
  template <typename M>
  void resize(const matrix_reference_t<T,M>& _other) { resize(_other.const_ref()); }
  /// load _other (resize() and cp())
  template <typename M>
  reference_t& load(const matrix_const_reference_t<T,M>& _other) {
    resize(_other);
    if (_other.data()!=0)
      cp(_other,*this);
    return *this;
  }
  template <typename M>
  reference_t& load(const matrix_reference_t<T,M>& _other) {
    return load(_other.const_ref());
  }


  /** Read data from MATLAB V4 file.
      \ingroup vc_blas_io
      Read matrix from Level 4 MAT-file, see also VC::math::blas::write_mat4().
      On failure the matrix size is set to m()=n()=0.
      \param _in input stream
      \return matrix name (or failure description).
      Note: in contrast to, e.g., mat4::read_matrix() this function does
      \a not throw an exception but signal failure by the return value.
      \sa [\ref vc_blas_io], write_mat4(), [\ref vc_mat4]
   */
  std::string read_mat4(std::istream& _in);

  /** Short for VC::math::blas::write_mat4(_out,*this,_name).
      \ingroup vc_blas_io
  */
  std::ostream& write_mat4(std::ostream& _out,const std::string& _name) const {
    return ::VC::math::blas::write_mat4(_out,this->const_ref(),_name);
  }

  /// delegate any assignment to matrix_reference_t
  template <typename A>
  const reference_t& operator=(const A& _a) { return reference_t::operator=(_a); }

};

//-----------------------------------------------------------------------------

/** \class GB_Matrix
    General banded matrix (see matrix_reference_t and matrix_const_reference_t).
    \ingroup vc_blas_matrices
 */
template <typename T,typename allocator=std::allocator<T> >
class GB_Matrix : public matrix_reference_t<T,gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> >,
                  public Matrix_alloc<allocator> {
public:
  typedef GB_Matrix<T,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef allocator allocator_t;

  GB_Matrix()
    : reference_t(matrix_t(0,0,0,0,0),0) {}
  GB_Matrix(int _m,int _n,int _kl,int _ku,int _ld=0,unsigned _capacity=0)
    : reference_t(matrix_t(0,0,0,0,0),0) {
    if (_capacity>0)
      reserve(_capacity);
    resize(_m,_n,_kl,_ku,_ld);
  }
  /// copy ctor calls load()
  GB_Matrix(const self_t& _a)
    : reference_t(matrix_t(0,0,0),0) {
    this->load(_a.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~GB_Matrix() { this->free(this->m_data);  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by GB_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_m(int);   // { assert(!"data is managed by GB_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by GB_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_ld(int);  // { assert(!"data is managed by GB_Matrix"); return *this; }

public:

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize matrix dimensions and KL,KU,LD (default _ld=_kl+_ku+1)
  void resize(int _m,int _n,int _kl,int _ku,int _ld=0) {
    if (_ld==0) _ld=_kl+_ku+1;
    this->m_m=matrix_t(_m,_n,_kl,_ku,_ld);
    reserve(this->size());
  }
  /// load _other (resize() and cp())
  template <typename M>
  reference_t& load(const matrix_const_reference_t<T,M>& _other) {
    resize(_other);
    if (_other.data()!=0)
      cp(_other,*this);
    return *this;
  }

  /// delegate any assignment to matrix_reference_t
  template <typename A>
  reference_t& operator=(const A& _a) { return reference_t::operator=(_a); }

};

//-----------------------------------------------------------------------------

/** \class SB_Matrix
    Symmetric banded matrix (see matrix_reference_t and matrix_const_reference_t).
    \ingroup vc_blas_matrices
 */
template <typename T,typename allocator=std::allocator<T> >
class SB_Matrix : public matrix_reference_t<T,gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> >,
                  public Matrix_alloc<allocator> {
public:
  typedef SB_Matrix<T,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef allocator allocator_t;

  SB_Matrix()
    : reference_t(matrix_t(0,0,0,0,0),0) {}
  SB_Matrix(int _m,int _n,int _k,int _ld=0,unsigned _capacity=0)
    : reference_t(matrix_t(0,0,0,0),0) {
    if (_capacity>0)
      reserve(_capacity);
    resize(_m,_n,_k,_ld);
  }
  /// copy ctor calls load()
  SB_Matrix(const self_t& _a)
    : reference_t(matrix_t(0,0,0),0) {
    this->load(_a.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~SB_Matrix() { this->free(this->m_data);  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by SB_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by SB_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_ld(int);  // { assert(!"data is managed by SB_Matrix"); return *this; }

public:

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize matrix dimensions and KL,KU,LD (default _ld=_kl+_ku+1)
  void resize(int _m,int _n,int _k,int _ld=0) {
    if (_ld==0) _ld=_k+1;
    this->m_m=matrix_t(_m,_n,_k,_ld);
    reserve(this->size());
  }
  /// load _other (resize() and cp())
  template <typename M>
  reference_t& load(const matrix_const_reference_t<T,M>& _other) {
    resize(_other);
    if (_other.data()!=0)
      cp(_other,*this);
    return *this;
  }

  /// delegate any assignment to matrix_reference_t
  template <typename A>
  reference_t& operator=(const A& _a) { return reference_t::operator=(_a); }

};

//-----------------------------------------------------------------------------

/** \class SP_Matrix
    Symmetric packed matrix (see matrix_reference_t and matrix_const_reference_t).
    \ingroup vc_blas_matrices
 */
template <typename T,UpperLowerFlag UPLO,
          typename allocator=std::allocator<T> >
class SP_Matrix : public matrix_reference_t<T,sp_mat<UPLO,VarInt> >,
                  public Matrix_alloc<allocator> {
public:
  typedef SP_Matrix<T,UPLO,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef sp_mat<UPLO,VarInt> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef allocator allocator_t;

  SP_Matrix()
    : reference_t(matrix_t(0),0) {}
  SP_Matrix(int _n,unsigned _capacity=0)
    : reference_t(matrix_t(0),0) {
    if (_capacity>0)
      reserve(_capacity);
    resize(_n);
  }
  /// copy ctor calls load()
  SP_Matrix(const self_t& _a)
    : reference_t(matrix_t(0,0,0),0) {
    this->load(_a.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~SP_Matrix() { this->free(this->m_data);  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by SP_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by SP_Matrix"); return *this; }

public:

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize matrix dimensions
  void resize(int _n) {
    this->m_m=matrix_t(_n);
    reserve(this->size());
  }
  /// resize to same dimensions as _other (assumes quadratic matrix)
  template <typename M>
  void resize(const matrix_const_reference_t<T,M>& _other) {
    assert(_other.m()==_other.n());
    resize(_other.n());
  }
  /// load _other (resize() and cp())
  template <typename M>
  reference_t& load(const matrix_const_reference_t<T,M>& _other) {
    resize(_other);
    if (_other.data()!=0)
      cp(_other,*this);
    return *this;
  }

  /// delegate any assignment to matrix_reference_t
  template <typename A>
  const reference_t& operator=(const A& _a) { return reference_t::operator=(_a); }

};

//-----------------------------------------------------------------------------

/** \class TP_Matrix
    Packed triangular matrix (see matrix_reference_t and matrix_const_reference_t).
    \ingroup vc_blas_matrices
 */
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG=NoU,
          typename allocator=std::allocator<T> >
class TP_Matrix : public matrix_reference_t<T,tp_mat<UPLO,TRANS,DIAG,VarInt> >,
                  public Matrix_alloc<allocator> {
public:
  typedef TP_Matrix<T,UPLO,TRANS,DIAG,allocator> self_t;
  typedef T value_t;
  typedef T value_type;
  typedef tp_mat<UPLO,TRANS,DIAG,VarInt> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef allocator allocator_t;

  TP_Matrix()
    : reference_t(matrix_t(0),0) {}
  TP_Matrix(int _n,unsigned _capacity=0)
    : reference_t(matrix_t(0),0) {
    if (_capacity>0)
      reserve(_capacity);
    resize(_n);
  }
  /// copy ctor calls load()
  TP_Matrix(const self_t& _a)
    : reference_t(matrix_t(0,0,0),0) {
    this->load(_a.const_ref());
  }
  /// assignment calls load()
  self_t& operator=(const self_t& _x) {
    this->load(_x.const_ref());
    return *this;
  }

  ~TP_Matrix() { this->free(this->m_data);  }

private:
  /// forbidden
  const reference_t& set_data(T*); // { assert(!"data is managed by TP_Matrix"); return *this; }
  /// forbidden
  const reference_t& set_n(int);   // { assert(!"data is managed by TP_Matrix"); return *this; }

public:

  /// reserve storage for _capacity elements (called by resize())
  void reserve(unsigned _capacity) {
    Matrix_alloc<allocator>::reserve(_capacity,this->m_data);
  }
  /// resize matrix dimensions
  void resize(int _n) {
    this->m_m=matrix_t(_n);
    reserve(this->size());
  }
  /// resize to same dimensions as _other (assumes quadratic matrix)
  template <typename M>
  void resize(const matrix_const_reference_t<T,M>& _other) {
    assert(_other.m()==_other.n());
    resize(_other.n());
  }
  /// load _other (resize() and cp())
  template <typename M>
  reference_t& load(const matrix_const_reference_t<T,M>& _other) {
    resize(_other);
    if (_other.data()!=0)
      cp(_other,*this);
    return *this;
  }

  /// delegate any assignment to matrix_reference_t
  template <typename A>
  reference_t& operator=(const A& _a) { return reference_t::operator=(_a); }

};

//-----------------------------------------------------------------------------

template <typename T,typename allocator>
std::string GE_Matrix<T,allocator>::read_mat4(std::istream& _in) {
  resize(0,0); // indicates failure
  mat4::MatrixInfo mat4;
  try {
    mat4.read(_in);

    if (mat4.n_imag()>0)
      return "cannot read complex matrix";

    // we ignore mat4::Text and mat4::Sparse and just read the data

    resize(mat4.mrows(),mat4.ncols());
    mat4.read_matrix(_in,this->ld(),this->data());

   } catch (base::vc_runtime_error& _e) {
    resize(0,0); // indicates failure
    return _e.what();
  }
  return mat4.name();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** \class Vec_N
    Vector with fixed size.
    \ingroup vc_blas_matrices
    This class provides only storage and defines types.
 */
template <typename T,int N>
class Vec_N : public vec<N,1> {
public:
  typedef T value_t;
  typedef T value_type;
  typedef vec<N,1> vector_t;
  typedef vector_reference_t<value_t,vector_t> reference_t;
  typedef vector_const_reference_t<value_t,vector_t> const_reference_t;

  Vec_N() : vec<N,1>(N,1) {}

  /** @name matrix_reference_t interface
      @{
   */

  /// get matrix_const_reference_t
  const_reference_t const_ref() const {
    return const_reference_t(vector_t(N,1),m_v);
  }
  /// get matrix_reference_t
  reference_t ref() {
    return reference_t(vector_t(N,1),m_v);
  }
  /// same as const_ref()
  const_reference_t operator()() const {
    return const_reference_t(vector_t(N,1),m_v);
  }
  /// same as ref()
  reference_t operator()() {
    return reference_t(vector_t(N,1),m_v);
  }

  /// @}

  /// get element
  T& operator()(int _i) const { return at(reference_t(vector_t(N,1),m_v),_i); }
  /// get element
  T& operator[](int _i) const { return at(reference_t(vector_t(N,1),m_v),_i); }

private:
  T m_v[N]; //!< data
};

//-----------------------------------------------------------------------------

/** \class GE_Mat_MxN
    General matrix with fixed size.
    \ingroup vc_blas_matrices
    This class provides only storage and defines types.
 */
template <typename T,int M,int N>
class GE_Mat_MxN : public ge_mat<NoT,M,N,M> {
public:
  typedef T value_t;
  typedef T value_type;
  typedef ge_mat<NoT,M,N,M> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef matrix_const_reference_t<value_t,matrix_t> const_reference_t;

  typedef sy_mat<Upper,N,N> sy_upper_t; //!< sy_mat
  typedef sy_mat<Lower,N,N> sy_lower_t; //!< sy_mat
  typedef tr_mat<Upper,NoT,NoU,N,N> tr_upper_t; //!< tr_mat
  typedef tr_mat<Lower,NoT,NoU,N,N> tr_lower_t; //!< tr_mat
  typedef tr_mat<Upper,NoT,UnitDiag,N,N> tr_upper_u_t; //!< tr_mat (UnitDiag)
  typedef tr_mat<Lower,NoT,UnitDiag,N,N> tr_lower_u_t; //!< tr_mat (UnitDiag)

  typedef vec<M,1> column_vec_t; //!< column vector
  typedef vec<N,M> row_vec_t; //!< row vector

  /// reference to row vector
  typedef vector_const_reference_t<row_vec_t,T> row_const_reference_t;
  /// reference to row vector
  typedef vector_reference_t<row_vec_t,T> row_reference_t;
  /// reference to column vector
  typedef vector_const_reference_t<column_vec_t,T> column_const_reference_t;
  /// reference to column vector
  typedef vector_reference_t<column_vec_t,T> column_reference_t;


  GE_Mat_MxN() : ge_mat<NoT,M,N,M>(M,N,M) {}

  /** @name matrix_reference_t interface
      @{
   */

  /// get matrix_const_reference_t
  const_reference_t const_ref() const {
    return const_reference_t(matrix_t(M,N,M),m_m);
  }
  /// get matrix_reference_t
  reference_t ref() {
    return reference_t(matrix_t(M,N,M),m_m);
  }
  /// same as const_ref()
  const_reference_t operator()() const {
    return const_reference_t(matrix_t(M,N,M),m_m);
  }
  /// same as ref()
  reference_t operator()() {
    return reference_t(matrix_t(M,N,M),m_m);
  }  

  /// @}

  /** @name iterators
      @{
   */

  typedef GE_row_iterator<T> row_iterator_t; //!< iterate rows
  typedef GE_column_iterator<T> column_iterator_t; //!< iterate columns

  /// iterate rows
  row_iterator_t row_iter() {
    return GE_row_iterator<T>(this->data(),this->m(),this->n(),this->ld());
  }
  /// iterate columns
  column_iterator_t column_iterator() {
    return GE_column_iterator<T>(this->data(),this->m(),this->n(),this->ld());
  }

  /// @}

  /** @name Access submatrix/block.
      For GE_Mat_MxN operator()(i1,i2,j1,j2)/operator()(idx_range,idx_range)
      is equal to block(*this,...), i.e., always provides a reference to
      a submatrix.
      @{
  */

  /// access element
  T& operator()(int _i,int _j) const {
    return at(reference_t(matrix_t(M,N,M),m_m),_i,_j);
  }
  /// get reference to submatrix
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(int _i1,int _i2,int _j1,int _j2) const {
    return block(reference_t(matrix_t(M,N,M),m_m),_i1,_i2,_j1,_j2);
  }
  /// get reference to submatrix A(i1:i2,j1:j2)
  template <bool END1,bool END2,bool END3,bool END4>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const idx_range<END1,END2>& _ri,const idx_range<END3,END4>& _rj) const {
    return block(reference_t(matrix_t(M,N,M),m_m),_ri,_rj);
  }
  /// get reference to submatrix A(i1:i2,j) (IDX={int,idx_end})
  template <bool END1,bool END2,typename IDX>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const idx_range<END1,END2>& _ri,const IDX& _j) const {
    return block(reference_t(matrix_t(M,N,M),m_m),_ri,_j);
  }
  /// get reference to submatrix A(i,j1:j2) (IDX={int,idx_end})
  template <bool END1,bool END2,typename IDX>
  matrix_reference_t<T,matrix_t> // VarInt
  operator()(const IDX& _i,const idx_range<END1,END2>& _rj) const {
    return block(reference_t(matrix_t(M,N,M),m_m),_i,_rj);
  }

  /// get column vector
  column_reference_t column(int _j) const {
    assert(0<=_j && _j<N);
    return column_reference_t(column_vec_t(M,1),m_m+_j*M);
  }
  /// get row vector
  row_reference_t row(int _i) const {
    assert(0<=_i && _i<M);
    return row_reference_t(row_vec_t(N,M),m_m+_i);
  }


  /// @}

  /** @name Reference upper and lower triangle (as sy_mat, tr_mat).
      @{
  */

  /// get reference to lower triangle (as sy_mat)
  matrix_reference_t<T,sy_lower_t> sy_lower() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,sy_lower_t>(sy_lower_t(this->n(),this->ld()),(T*) this->m_m);
  }
  /// get reference to upper triangle (as sy_mat)
  matrix_reference_t<T,sy_upper_t> sy_upper() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,sy_upper_t>(sy_upper_t(this->n(),this->ld()),(T*) this->m_m);
  }
  /// get reference to lower triangle (as tr_mat)
  matrix_reference_t<T,tr_lower_t> tr_lower() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_lower_t>(tr_lower_t(this->n(),this->ld()),(T*) this->m_m);
  }
  /// get reference to upper triangle (as tr_mat)
  matrix_reference_t<T,tr_upper_t> tr_upper() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_upper_t>(tr_upper_t(this->n(),this->ld()),(T*) this->m_m);
  }
  /// get reference to lower triangle (as tr_mat w/ UnitDiag)
  matrix_reference_t<T,tr_lower_u_t> tr_lower_u() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_lower_u_t>(tr_lower_u_t(this->n(),this->ld()),(T*) this->m_m);
  }
  /// get reference to upper triangle (as tr_mat w/ UnitDiag)
  matrix_reference_t<T,tr_upper_t> tr_upper_u() const {
    assert(this->m()==this->n());
    return matrix_reference_t<T,tr_upper_u_t>(tr_upper_u_t(this->n(),this->ld()),(T*) this->m_m);
  }

  /// @}

private:
  T m_m[M*N]; //!< data
};

//-----------------------------------------------------------------------------

/** \class SP_Mat_NxN
    Symmetric packed matrix with fixed size.
    \ingroup vc_blas_matrices
    This class provides only storage and defines types.
 */
template <typename T,UpperLowerFlag UPLO,int N >
class SP_Mat_NxN : public sp_mat<UPLO,N> {
public:
  typedef T value_t;
  typedef T value_type;
  typedef sp_mat<UPLO,N> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef matrix_const_reference_t<value_t,matrix_t> const_reference_t;

  SP_Mat_NxN() : sp_mat<UPLO,N>(N) {}

  /** @name matrix_reference_t interface
      @{
   */

  /// get matrix_const_reference_t
  const_reference_t const_ref() const {
    return const_reference_t(matrix_t(N),m_m);
  }
  /// get matrix_reference_t
  reference_t ref() {
    return reference_t(matrix_t(N),m_m);
  }
  /// same as const_ref()
  const_reference_t operator()() const {
    return const_reference_t(matrix_t(N),m_m);
  }
  /// same as ref()
  reference_t operator()() {
    return reference_t(matrix_t(N),m_m);
  }

  /// @}

private:
  T m_m[N*(N-1)/2]; //!< data
};

//-----------------------------------------------------------------------------

/** \class TP_Mat_NxN
    Packed triangular matrix with fixed size.
    \ingroup vc_blas_matrices
    This class provides only storage and defines types.
 */
template <typename T,UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N>
class TP_Mat_NxN : public tp_mat<UPLO,TRANS,DIAG,N> {
public:
  typedef T value_t;
  typedef T value_type;
  typedef tp_mat<UPLO,TRANS,DIAG,N> matrix_t;
  typedef matrix_reference_t<value_t,matrix_t> reference_t;
  typedef matrix_const_reference_t<value_t,matrix_t> const_reference_t;

  TP_Mat_NxN() : tp_mat<UPLO,TRANS,DIAG,N>(N) {}

  /** @name matrix_reference_t interface
      @{
   */

  /// get matrix_const_reference_t
  const_reference_t const_ref() const {
    return const_reference_t(matrix_t(N),m_m);
  }
  /// get matrix_reference_t
  reference_t ref() {
    return reference_t(matrix_t(N),m_m);
  }
  /// same as const_ref()
  const_reference_t operator()() const {
    return const_reference_t(matrix_t(N),m_m);
  }
  /// same as ref()
  reference_t operator()() {
    return reference_t(matrix_t(N),m_m);
  }

  /// @}

private:
  T m_m[N*(N-1)/2]; //!< data
};

//-----------------------------------------------------------------------------


//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
