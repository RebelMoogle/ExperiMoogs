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

    \arg Define BLAS matrix specifications.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================


/// run-time variable integer property \ingroup vc_blas_mat
enum VarInt_t { VarInt = INT_MAX };
/// undefined integer property \ingroup vc_blas_mat
enum UndfInt_t { UndefInt = VarInt-1 };

  //
  // VarInt ==> int_max !!! inc may be negative, k may be 0 !!!
  //

# ifndef DOXYGEN_SKIP

# define define_mat_prop(type,name)                  \
  template <type N>                                  \
  struct mat_prop_##name {                           \
    mat_prop_##name<N>(type _i=N) { assert(_i==N); } \
    static type name() { return N; }                   \
    static void set_##name(type /*_arg*/) {            \
      assert( ! #name " is not variable at run-time"); \
    }                                                  \
  }

# define define_mat_prop2(type,name)            \
  define_mat_prop(type,name);                   \
  template <> struct mat_prop_##name<VarInt> {                          \
    mat_prop_##name<VarInt>(type _i) : m_##name(_i) { assert(_i!=VarInt); } \
    int name() const { return m_##name; }                               \
    void set_##name(type _arg) { m_##name=_arg; }                       \
  protected: type m_##name; }

define_mat_prop(UpperLowerFlag,uplo);
define_mat_prop(TransposeFlag,trans);
define_mat_prop(DiagonalFlag,diag);
define_mat_prop2(int,m);
define_mat_prop2(int,n);
define_mat_prop2(int,ld);
define_mat_prop2(int,k);
define_mat_prop2(int,kl);
define_mat_prop2(int,ku);
define_mat_prop2(int,inc);

# endif

//-----------------------------------------------------------------------------

template <typename T,typename N> struct vector_reference_t;
template <typename T,typename N> struct vector_const_reference_t;

/** (Internally used as base of vec.) \ingroup vc_blas_mat
    \ingroup vc_blas_mat
 */
template <int N,int INC>
struct vec_base
  : public mat_prop_n<N>,
    public mat_prop_inc<INC> {

  vec_base(int _n=N,int _inc=INC)
    : mat_prop_n<N>(_n), mat_prop_inc<INC>(_inc) {}
};

/** Specify BLAS vector.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <int N=VarInt,int INC=1>
struct vec : public vec_base<N,INC> {
  vec(int _n=N,int _inc=INC) : vec_base<N,INC>(_n,_inc) {}
  typedef vec<N,INC> self_t;

  // get vector_const_reference_t to _data
  template <typename T> vector_const_reference_t<T,self_t>
  const_reference(const T* _data) const {
    return vector_const_reference_t<T,self_t>(*this,_data);
  }
  /// get vector_const_reference_t to _data
  template <typename T>  vector_reference_t<T,self_t>
  reference(T* _data) const {
    return vector_reference_t<T,self_t>(*this,_data);
  }

  /** @name get reference types from scalar type (as for mat)
      @{
   */
  template <typename T>
  struct ref_t {
    /// matrix_reference_t
    typedef vector_reference_t<T,self_t> type;
  };
  template <typename T>
  struct const_ref_t {
    /// matrix_const_reference_t
    typedef vector_const_reference_t<T,self_t> type;
  };

  /// @}


  /// get vector_const_reference() (as in mat)
  template <typename T>
  static  vector_const_reference_t<T,self_t>
  const_ref(const T* _data,int _n=N,int _inc=INC) {
    return vector_const_reference_t<T,self_t>(self_t(_n,_inc),_data);
  }
  /// get vector_reference() (as in mat)
  template <typename T>
  static vector_reference_t<T,self_t>
  ref(T* _data,int _n=N,int _inc=INC)  {
    return vector_reference_t<T,self_t>(self_t(_n,_inc),_data);
  }
};

//-----------------------------------------------------------------------------


/// BLAS GEneral matrix id \ingroup vc_blas_mat
enum GE_Mat { GE=0 };
/// BLAS General Banded matrix id \ingroup vc_blas_mat
enum GB_Mat { GB=1 };
/// BLAS SYmmetric matrix id \ingroup vc_blas_mat
enum SY_Mat { SY=2 };
/// BLAS Symmetric Banded matrix id \ingroup vc_blas_mat
enum SB_Mat { SB=3 };
/// BLAS Symmetric Packed matrix id \ingroup vc_blas_mat
enum SP_Mat { SP=4 };
 /// BLAS TRiangular matrix id \ingroup vc_blas_mat
enum TR_Mat { TR=5 };
/// BLAS Triangular Banded matrix id \ingroup vc_blas_mat
enum TB_Mat { TB=6 };
/// BLAS Triangular Packed matrix id \ingroup vc_blas_mat
enum TP_Mat { TP=7 };

/** (Internally used as base class to matrix classes.)
    \ingroup vc_blas_mat
 */
struct mat_prop_base {
  static bool is_symmetric() { return false; }
  static bool is_triangular() { return false; }
  static bool is_dense() { return true; }
  static bool is_banded() { return false; }
  static bool is_packed() { return false; }
};

  //
  // Sorry, I felt this was a pain in the ass doing the following
  // "correctly...
  //

  template <typename T,typename M> struct matrix_reference_t;
  template <typename T,typename M> struct matrix_const_reference_t;

# define _MIXIN_MATRIX_REFERENCE \
  template <typename T>  matrix_const_reference_t<T,self_t>\
  const_reference(const T* _data) const { \
    return matrix_const_reference_t<T,self_t>(*this,_data); \
  } \
  template <typename T> matrix_reference_t<T,self_t>    \
  reference(T* _data) const { \
    return matrix_reference_t<T,self_t>(*this,_data);   \
  }

  //
  //
  //

/** (Internally used as base class to ge_mat.)
    \ingroup vc_blas_mat
 */
template <TransposeFlag TRANS,int M,int N,int LD>
struct ge_mat_base
  : public mat_prop_base,
    public mat_prop_trans<TRANS>,
    public mat_prop_m<M>,
    public mat_prop_n<N>,
    public mat_prop_ld<LD> {
  typedef GE_Mat matrix_id_t;
  static int matrix_type() { return GE; }
  ge_mat_base(int _m=M,int _n=N,int _ld=LD)
    : mat_prop_m<M>(_m), mat_prop_n<N>(_n), mat_prop_ld<LD>(_ld) {
    assert(_ld>=_m);
  }

  unsigned size() const { return this->ld()*this->n(); }
};

/** Specify BLAS GEneral matrix.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <TransposeFlag TRANS,int M,int N,int LD=M>
struct ge_mat : public ge_mat_base<TRANS,M,N,LD> {
  ge_mat(int _m=M,int _n=N,int _ld=LD)
    : ge_mat_base<TRANS,M,N,LD>(_m,_n,_ld) {}
  typedef ge_mat<TRANS,M,N,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to gb_mat.)
    \ingroup vc_blas_mat
 */
template <TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
struct gb_mat_base
  : public mat_prop_base,
    public mat_prop_trans<TRANS>,
    public mat_prop_m<M>,
    public mat_prop_n<N>,
    public mat_prop_kl<KL>,
    public mat_prop_ku<KU>,
    public mat_prop_ld<LD> {
  typedef GB_Mat matrix_id_t;
  static int matrix_type() { return GB; }
  gb_mat_base(int _m=M,int _n=N,int _kl=KL,int _ku=KU,int _ld=LD)
    : mat_prop_m<M>(_m), mat_prop_n<N>(_n), mat_prop_kl<M>(_kl), mat_prop_ku<N>(_ku),
      mat_prop_ld<LD>(_ld) {
    assert(this->ld()>=this->kl()+this->ku()+1);
  }
  static bool is_banded() { return true; }

  unsigned size() const { return this->ld()*this->n(); }
};

/** Specify BLAS General Banded matrix.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
struct gb_mat : public gb_mat_base<TRANS,M,N,KL,KU,LD> {
  gb_mat(int _m=M,int _n=N,int _kl=KL,int _ku=KU,int _ld=LD)
    : gb_mat_base<TRANS,M,N,KL,KU,LD>(_m,_n,_kl,_ku,_ld) {}
  typedef gb_mat<TRANS,M,N,KL,KU,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to sy_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,int N,int LD>
struct sy_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_n<N>,
    public mat_prop_ld<LD> {
  typedef SY_Mat matrix_id_t;
  static int matrix_type() { return SY; }
  sy_mat_base(int _n=N,int _ld=LD)
    : mat_prop_n<N>(_n), mat_prop_ld<LD>(_ld) { assert(_ld>=_n); }
  static bool is_symmetric() { return true; }

  int m() const { return this->n(); }
  TransposeFlag trans() const { return NoT; }
  unsigned size() const { return this->ld()*this->n(); }
};

/** Specify BLAS SYmmetric matrix.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <UpperLowerFlag UPLO,int N,int LD=N>
struct sy_mat : public sy_mat_base<UPLO,N,LD> {
  sy_mat(int _n=N,int _ld=LD) : sy_mat_base<UPLO,N,LD>(_n,_ld) {}
  typedef sy_mat<UPLO,N,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to sb_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,int N,int K,int LD>
struct sb_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_n<N>,
    public mat_prop_k<K>,
    public mat_prop_ld<LD> {
  typedef SB_Mat matrix_id_t;
  static int matrix_type() { return SB ; }
  sb_mat_base(int _n=N,int _k=K,int _ld=LD)
    : mat_prop_n<N>(_n), mat_prop_k<K>(_k), mat_prop_ld<LD>(_ld) {
    assert(_ld>=_k+1);
  }

  static bool is_symmetric() { return true; }
  static bool is_banded() { return true; }

  int m() const { return this->n(); }
  TransposeFlag trans() const { return NoT; }
  unsigned size() const { return this->ld()*this->n(); }
};

  // Better no default LD for sb mat! (K+1 would make sense)

/** Specify BLAS Symmetric Banded matrix.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <UpperLowerFlag UPLO,int N,int K,int LD>
struct sb_mat : public sb_mat_base<UPLO,N,K,LD> {
  typedef SB_Mat matrix_id_t;
  sb_mat(int _n=N,int _k=K,int _ld=LD) : sb_mat_base<UPLO,N,K,LD>(_n,_k,_ld) {}
  typedef sb_mat<UPLO,N,K,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to sp_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,int N>
struct sp_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_n<N> {
  typedef SP_Mat matrix_id_t;
  static int matrix_type() { return SP; }
  sp_mat_base(int _n=N)
    : mat_prop_n<N>(_n) {}

  static bool is_symmetric() { return true; }
  static bool is_packed() { return true; }

  int m() const { return this->n(); }
  TransposeFlag trans() const { return NoT; }
  unsigned size() const {  return this->n()*(this->n()+1)/2; }
};

/** Specify BLAS Symmetric Packed matrix.
    \ingroup vc_blas_mat
    \sa matrix_reference_t
 */
template <UpperLowerFlag UPLO,int N>
struct sp_mat : public sp_mat_base<UPLO,N> {
  sp_mat(int _n=N) : sp_mat_base<UPLO,N>(_n) {}
  typedef sp_mat<UPLO,N> self_t;
  _MIXIN_MATRIX_REFERENCE
};

/** (Internally used as base class to tr_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int LD>
struct tr_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_trans<TRANS>,
    public mat_prop_diag<DIAG>,
    public mat_prop_n<N>,
    public mat_prop_ld<LD> {
  typedef TR_Mat matrix_id_t;
  static int matrix_type() { return TR; }
  tr_mat_base(int _n=N,int _ld=LD)
    : mat_prop_n<N>(_n), mat_prop_ld<LD>(_ld) { assert(_ld>=_n); }

  static bool is_triangular() { return true; }

  int m() const { return this->n(); }
  unsigned size() const { return this->ld()*this->n(); }
};

template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int LD=N>
struct tr_mat : public tr_mat_base<UPLO,TRANS,DIAG,N,LD> {
  tr_mat(int _n=N,int _ld=LD) : tr_mat_base<UPLO,TRANS,DIAG,N,LD>(_n,_ld) {}
  typedef tr_mat<UPLO,TRANS,DIAG,N,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to tb_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int K,int LD>
struct tb_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_trans<TRANS>,
    public mat_prop_diag<DIAG>,
    public mat_prop_n<N>,
    public mat_prop_k<K>,
    public mat_prop_ld<LD> {
  typedef TB_Mat matrix_id_t;
  static int matrix_type() { return TB; }
  tb_mat_base(int _n=N,int _k=K,int _ld=LD)
    : mat_prop_n<N>(_n), mat_prop_k<K>(_k), mat_prop_ld<LD>(_ld) {
    assert(_ld>=_k+1);
  }

  static bool is_triangular() { return true; }
  static bool is_banded() { return true; }

  int m() const { return this->n(); }
  unsigned size() const { return this->ld()*this->n(); }
};


template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int K,int LD>
struct tb_mat : public tb_mat_base<UPLO,TRANS,DIAG,N,K,LD> {
  tb_mat(int _n=N,int _k=K,int _ld=LD) : tb_mat_base<UPLO,TRANS,DIAG,N,K,LD>(_n,_k,_ld) {}
  typedef tb_mat<UPLO,TRANS,DIAG,N,K,LD> self_t;
  _MIXIN_MATRIX_REFERENCE
};


/** (Internally used as base class to tp_mat.)
    \ingroup vc_blas_mat
 */
template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N>
struct tp_mat_base
  : public mat_prop_base,
    public mat_prop_uplo<UPLO>,
    public mat_prop_trans<TRANS>,
    public mat_prop_diag<DIAG>,
    public mat_prop_n<N> {
  typedef TP_Mat matrix_id_t;
  static int matrix_type() { return TP; }
  tp_mat_base(int _n=N)
    : mat_prop_n<N>(_n)  {}

  static bool is_triangular() { return true; }
  static bool is_packed() { return true; }

  int m() const { return this->n(); }
  unsigned size() const {  return this->n()*(this->n()+1)/2; }
};

template <UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N>
struct tp_mat : public tp_mat_base<UPLO,TRANS,DIAG,N> {
  tp_mat(int _n=N) : tp_mat_base<UPLO,TRANS,DIAG,N>(_n) {}
  typedef tp_mat<UPLO,TRANS,DIAG,N> self_t;
  _MIXIN_MATRIX_REFERENCE
};

template <UpperLowerFlag UPLO,DiagonalFlag DIAG,int N>
struct tp_mat<UPLO,NoT,DIAG,N> : public tp_mat_base<UPLO,NoT,DIAG,N> {
  tp_mat(int _n=N) : tp_mat_base<UPLO,NoT,DIAG,N>(_n) {}
  typedef tp_mat<UPLO,NoT,DIAG,N> self_t;
  _MIXIN_MATRIX_REFERENCE
};

# undef _MIXIN_MATRIX_REFERENCE

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
