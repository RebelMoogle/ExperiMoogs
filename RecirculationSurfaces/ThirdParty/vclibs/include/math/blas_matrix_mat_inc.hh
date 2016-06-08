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

    \arg Provide mapping < T,matrix_t > -> matrix_reference_t:
    simplify construction of matrix_reference_t

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

/** (Internally used: provides types ref_t and const_ref_t to mat.
    \ingroup vc_blas
    \internal
    Base class providing ref_t and const_ref_t as matrix_reference_t
    and matrix_const_reference_t, respectively.
*/
template <typename A>
struct mat_base {
  template <typename T>
  struct ref_t {
    /// matrix_reference_t
    typedef matrix_reference_t<T,A> type;
  };
  template <typename T>
  struct const_ref_t {
    /// matrix_const_reference_t
    typedef matrix_const_reference_t<T,A> type;
  };
};

/** Provide matrix_reference_t and matrix_const_reference_t objects.
    \ingroup vc_blas

    Specializations apply to ge_mat, etc. The first template parameter
    specifies the matrix structure, e.g., VC::math::blas::GE,
    subsequent parameters depend on matrix structure, e.g., M,N,LD for
    ge_mat. The VC::math::blas::TransposeFlag is \a always omitted!

    Specializations are available for all dimensions/parameters fixed
    and all parameters run-time variable (as for using
    VC::math::blas::VarInt).

    Example
    \code
    mat<GE,3,4>::const_ref_t<double>::type A=mat<GE,3,4>::const_ref(ptr); // fixed
    mat<GE>::const_ref_t<double>::type B=mat<GE>::const_ref(ptr,3,4);     // variable
    \endcode

    The default \a leading \a dimension parameter is either m() (or
    n() for quadratic matrices, or k()+1 for symmetric banded matrices
    (sy_mat) or kl()+ku()+1 for general banded matrices (gb_mat()).

    Note that the vec structure which specifies a vector (n() and
    inc()) provides the same interface:
    \code
    vec<3>::ref_t<double>::type x=vec<3>::ref(dx); // fixed
    vec<>::ref_t<double>::type y=vec<>::ref(dx,3); // variable
    \endcode
*/
template <int TYPE,
          int ARG1=UndefInt,int ARG2=UndefInt,int ARG3=UndefInt,
          int ARG4=UndefInt,int ARG5=UndefInt,int ARG6=UndefInt>
struct mat {
# ifdef DOXYGEN_SKIP
  typedef XX_mat<> matrix_t; //!< matrix type, e.g., ge_mat

  /** @name get reference types from scalar type (mat_base mixin)
      @{
   */
  template <typename T>
  struct ref_t {
    /// matrix_reference_t
    typedef matrix_reference_t<T,A> type;
  };
  template <typename T>
  struct const_ref_t {
    /// matrix_const_reference_t
    typedef matrix_const_reference_t<T,A> type;
  };

  /// @}

  /// get matrix_reference_t from pointer (and variable parameters "...")
  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data,...) {
    return matrix_reference_t<T,matrix_t>(matrix_t(M,N,LD),_data);
  }
  /// get matrix_const_reference from pointer (and variable parameters "...")
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data,...) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(M,N,LD),_data);
  }
# else
 
  static_assert(ARG1!=ARG1,"invalid specialization");

# endif
};

// GE

template <int M,int N,int LD>
struct mat<GE,M,N,LD> : public mat_base<ge_mat<NoT,M,N,LD> > {
  typedef ge_mat<NoT,M,N,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(M,N,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(M,N,LD),_data);
  }
};

template <int M,int N>
struct mat<GE,M,N> : public mat_base<ge_mat<NoT,M,N,M> > {
  typedef ge_mat<NoT,M,N,M> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(M,N,M),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(M,N,M),_data);
  }
};

template <>
struct mat<GE> : public mat_base<ge_mat<NoT,VarInt,VarInt,VarInt> > {
  typedef ge_mat<NoT,VarInt,VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t>
  ref(T* _data,int _m,int _n,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_m,_n,_ld!=UndefInt ? _ld : _m),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _m,int _n,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_m,_n,_ld!=UndefInt ? _ld : _m),_data);
  }
};

// GB

template <int M,int N,int KL,int KU,int LD>
struct mat<GB,M,N,KL,KU,LD> : public mat_base<gb_mat<NoT,M,N,KL,KU,LD> > {
  typedef gb_mat<NoT,M,N,KL,KU,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(M,N,KL,KU,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(M,N,KL,KU,LD),_data);
  }
};

template <int M,int N,int KL,int KU>
struct mat<GB,M,N,KL,KU> : public mat_base<gb_mat<NoT,M,N,KL,KU,KL+KU+1> > {
  typedef gb_mat<NoT,M,N,KL,KU,KL+KU+1> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(M,N,KL,KU,KL+KU+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(M,N,KL,KU,KL+KU+1),_data);
  }
};

template <>
struct mat<GB> : public mat_base<gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> > {
  typedef gb_mat<NoT,VarInt,VarInt,VarInt,VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t>
  ref(T* _data,int _m,int _n,int _kl,int _ku,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_m,_n,_kl,_ku,_ld!=UndefInt ? _ld : _kl+_ku+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _m,int _n,int _kl,int _ku,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_m,_n,_kl,_ku,_ld!=UndefInt ? _ld : _kl+_ku+1),_data);
  }
};

// SY

template <int UPLO,int N,int LD>
struct mat<SY,UPLO,N,LD> : public mat_base<sy_mat<UpperLowerFlag(UPLO),N,LD> > {
  typedef sy_mat<UpperLowerFlag(UPLO),N,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,LD),_data);
  }
};

template <int UPLO,int N>
struct mat<SY,UPLO,N> : public mat_base<sy_mat<UpperLowerFlag(UPLO),N,N> > {
  typedef sy_mat<UpperLowerFlag(UPLO),N,N> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,N),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,N),_data);
  }
};

template <int UPLO>
struct mat<SY,UPLO> : public mat_base<sy_mat<UpperLowerFlag(UPLO),VarInt,VarInt> > {
  typedef sy_mat<UpperLowerFlag(UPLO),VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data,int _n,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_n,_ld!=UndefInt ? _ld : _n),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _n,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_n,_ld!=UndefInt ? _ld : _n),_data);
  }
};

// SB

template <int UPLO,int N,int K,int LD>
struct mat<SB,UPLO,N,K,LD> : public mat_base<sb_mat<UpperLowerFlag(UPLO),N,K,LD> > {
  typedef sb_mat<UpperLowerFlag(UPLO),N,K,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,K,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,K,LD),_data);
  }
};

template <int UPLO,int N,int K>
struct mat<SB,UPLO,N,K,N> : public mat_base<sb_mat<UpperLowerFlag(UPLO),N,K,K+1> > {
  typedef sb_mat<UpperLowerFlag(UPLO),N,K,K+1> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,K,K+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,K,K+1),_data);
  }
};

template <int UPLO>
struct mat<SB,UPLO>
  : public mat_base<sb_mat<UpperLowerFlag(UPLO),VarInt,VarInt,VarInt> > {
  typedef sb_mat<UpperLowerFlag(UPLO),VarInt,VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t>
  ref(T* _data,int _n,int _k,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_n,_k,_ld!=UndefInt ? _ld : _k+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _n,int _k,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_n,_ld!=UndefInt ? _ld : _k+1),_data);
  }
};


// SP

template <int UPLO,int N>
struct mat<SP,UPLO,N> : public mat_base<sp_mat<UpperLowerFlag(UPLO),N> > {
  typedef sp_mat<UpperLowerFlag(UPLO),N> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
};

template <int UPLO>
struct mat<SP,UPLO> : public mat_base<sp_mat<UpperLowerFlag(UPLO),VarInt> > {
  typedef sp_mat<UpperLowerFlag(UPLO),VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data,int _n) {
    return matrix_reference_t<T,matrix_t>(matrix_t(_n),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data,int _n) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(_n),_data);
  }
};

// TR

template <int UPLO,int DIAG,int N,int LD>
struct mat<TR,UPLO,DIAG,N,LD>
  : public mat_base<tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,LD> > {
  typedef tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,LD),_data);
  }
};

template <int UPLO,int DIAG,int N>
struct mat<TR,UPLO,DIAG,N>
  : public mat_base<tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,N> > {
  typedef tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,N> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
};

template <int UPLO,int DIAG>
struct mat<TR,UPLO,DIAG>
  : public mat_base<tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),VarInt,VarInt> > {
  typedef tr_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data,int _n,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_n,_ld!=UndefInt ? _ld : _n),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _n,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_n,_ld!=UndefInt ? _ld : _n),_data);
  }
};

// TB

template <int UPLO,int DIAG,int N,int K,int LD>
struct mat<TB,UPLO,DIAG,N,K,LD>
  : public mat_base<tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,K,LD> > {
  typedef tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,K,LD> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,K,LD),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,K,LD),_data);
  }
};

template <int UPLO,int DIAG,int N,int K>
struct mat<TB,UPLO,DIAG,N,K>
  : public mat_base<tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,K,K+1> > {
  typedef tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N,K,K+1> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N,K,K+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N,K,K+1),_data);
  }
};

template <int UPLO,int DIAG>
struct mat<TB,UPLO,DIAG>
  : public mat_base<tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),
                           VarInt,VarInt,VarInt> > {
  typedef tb_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),
                 VarInt,VarInt,VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t>
  ref(T* _data,int _n,int _k,int _ld=UndefInt) {
    return matrix_reference_t<T,matrix_t>
      (matrix_t(_n,_k,_ld!=UndefInt ? _ld : _k+1),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t>
  const_ref(const T* _data,int _n,int _k,int _ld=UndefInt) {
    return matrix_const_reference_t<T,matrix_t>
      (matrix_t(_n,_k,_ld!=UndefInt ? _ld : _k+1),_data);
  }
};

// TP

template <int UPLO,int DIAG,int N>
struct mat<TP,UPLO,DIAG,N>
  : public mat_base<tp_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N> > {
  typedef tp_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),N> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data) {
    return matrix_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(N),_data);
  }
};

template <int UPLO,int DIAG>
struct mat<TP,UPLO,DIAG>
  : public mat_base<tp_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),VarInt> > {
  typedef tp_mat<UpperLowerFlag(UPLO),NoT,DiagonalFlag(DIAG),VarInt> matrix_t ;

  template <typename T>
  static matrix_reference_t<T,matrix_t> ref(T* _data,int _n) {
    return matrix_reference_t<T,matrix_t>(matrix_t(_n),_data);
  }
  template <typename T>
  static matrix_const_reference_t<T,matrix_t> const_ref(const T* _data,int _n) {
    return matrix_const_reference_t<T,matrix_t>(matrix_t(_n),_data);
  }
};

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
