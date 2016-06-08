//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_MATH_LAPACK_HH
#define VC_MATH_LAPACK_HH

#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "blas.hh"
#include "lapack_prototypes.hh"
#include "../system/platform.hh" // THREAD_LOCAL, use_nowarn

/** \file lapack.hh Math: C++-wrappers for some LAPACK routines
 */

// debugging aid: define to track calls to LAPACK routines (wrappers)
# ifndef _VC_DBG_LAPACK_SCOPE
#  define _VC_DBG_LAPACK_SCOPE()
//#  define _VC_DBG_LAPACK_SCOPE() VC_DBG_SCOPE("lapack") // debug: show what's going on
# endif

/** Customized LAPACK error handler.
    Calls VC::math::lapack::xerbla_handler
    \ingroup vc_lapack
*/
extern "C" void XERBLA(char*,VC::math::lapack::INTEGER*);


namespace VC {
namespace math {
/// LAPACK wrappers, helpers, prototypes \sa [\ref vc_lapack] \ingroup vc_lapack
namespace lapack {

/** \defgroup vc_lapack LAPACK: C++ interface
    \ingroup vc_math
    \sa VC::math::lapack
 */

/// use [\ref vc_blas] types \ingroup vc_lapack
typedef blas::TransposeFlag TransposeFlag;
/// use [\ref vc_blas] types \ingroup vc_lapack
typedef blas::UpperLowerFlag UpperLowerFlag;


//=============================================================================


/** User defined xerbla error handler.
    \ingroup vc_lapack
*/
typedef void (*xerbla_handler_t)(const char* _routine,int _argi);

/** User defined xerbla error handler.
    \ingroup vc_lapack
    Default: xerbla_print_to_cerr()
*/

extern xerbla_handler_t xerbla_handler;

/** Standard xerbla_handler.
    \ingroup vc_lapack
 */
void xerbla_print_to_cerr(const char* _routine,int _argi);

//-----------------------------------------------------------------------------

/** Allocates temporary data.
    \ingroup vc_lapack
    \arg Use for `work` and `iwork` parameters to minimize calls to `malloc()`.
    \arg The idea is to keep a *single* instance of Workspace (per thread) such
    that allocated temporary storage is reused between calls.
    \arg Temporary storage is acquired by a *single* call to `Workspace::alloc()` 
    immediately before calling the LAPACK function.
    \arg There is no `"free()"`: ideally, `alloc()` always returns a pointer to
    the same memory block which is reused.
    \arg Data size only grows. Data is overwritten/reused with next call to
    `alloc()`. Only `free_all()` releases the whole data block.\
    \arg There is *no initialization* of storage!

    **Notes:**

    \arg The previous implementation used a global but *thread-local*
    Workspace instance. Initialization lead to problems with `OpenMP`, and 
    using C++11's `thread_local` is not yet supported by `clang++`. -- So I
    gave up on this pass the Workspace explicitly to functions.<br>
    You can  use the instance() method as a replacement for the previous
    global().
    
    \arg We try to *align* memory: reserve() allocates page-aligned
    memory on Unix systems, pointers from alloc() should be
    "well-aligned".   
 */
class Workspace {
public:

  /// empty workspace
  Workspace() : m_work(0), m_size(0) {}
  /// calls free_all()
  ~Workspace() { free_all(); }

  /// reserve _b bytes
  void reserve(size_t _n) {
    if (m_size<_n) {
      free(m_work); // no copy (realloc)
# if (_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600)
      static const size_t PAGESIZE=getpagesize();
      m_work=0;
      int rv=posix_memalign((void**) &m_work,PAGESIZE,_n); use_nowarn(rv);
      assert(rv==0);
# else
      m_work=malloc(_n);
# endif
      if (m_work==0) throw std::bad_alloc();
      m_size=_n;
    }
  }

  /// allocate _n instaces of T
  template<typename T>
  void alloc(T** _ptr,size_t _n) {
    size_t n=_n*sizeof(T);
    reserve(n);
    *_ptr=(T*) m_work;
  }
  /// allocate _m instances of TM and _n instaces of TN
  template<typename TM,typename TN>
  void alloc(TM** _pm,size_t _m,TN** _pn,size_t _n) {
    size_t m=align(_m*sizeof(TM)), n=align(_n*sizeof(TN));
    reserve(m+n);
    *_pm=(TM*) m_work;
    *_pn=(TN*) (((char*) m_work)+m);
  }
  /// allocate _m instances of TM and _n instaces of TN, and _o instances of TO
  template<typename TM,typename TN,typename TO>
  void alloc(TM** _pm,size_t _m,TN** _pn,size_t _n,TO** _po,size_t _o) {
    size_t m=align(_m*sizeof(TM)), n=align(_n*sizeof(TN)),
      o=align(_o*sizeof(TO));
    reserve(m+n+o);
    *_pm=(TM*) m_work;
    *_pn=(TN*) (((char*) m_work)+m);
    *_po=(TO*) (((char*) m_work)+m+n);
  }

  /// query current maxmium size in bytes
  size_t max_size() const { return m_size; }

  /// free allocated memory block
  void free_all() { std::free(m_work); m_work=0; m_size=0; }

  /// Get (per thread) global instance.       
# if defined(__GNUC__) && __GNUC_PREREQ(4,8)
  static Workspace& instance() {
    static thread_local Workspace ws;
    return ws;
  }
# else
  static Workspace& instance() {
    static THREAD_LOCAL Workspace* ws=0;
    if (ws==0)
      ws=new Workspace; // never freed!
    return *ws;
  }
# endif


  /// \deprecated
  static Workspace& global();
  /// \deprecated
  static void global_init();
  /// \deprecated
  static void global_destroy();
 

private:
  Workspace(const Workspace&);

  enum { ALIGN2=5 /**!< alignment in `2^ALIGN2` bytes */ };

  size_t align(size_t _n) { return ((_n+(1<<ALIGN2)-1)>>ALIGN2)<<ALIGN2; }

  void*   m_work;
  size_t  m_size;
};

//-----------------------------------------------------------------------------

/** Copy _m x _n source matrix _c from C-style to FORTRAN-style array _f.
    \ingroup vc_lapack
    \arg _c 2-dimensional C-style array (rows first)
    \arg _f 2-dimensional FORTRAN style array (columns first) of
    dimension M x _n, where M<=_ldf, use _ldf=_m if _ldf<_m
    \arg _m, _n dimensions of _c
    \arg _ldf leading dimension of _f
    Note: _c and _f must \a not overlap!
*/
template <typename C,typename F>
inline void
cp_c2f(const C* _c,F* _f,unsigned _m,unsigned _n,unsigned _ldf=0) {

  _ldf=std::max(_ldf,_m);

  if ((_m==1 && _ldf==1) || _n==1) {
    for (unsigned i=0;i<_m*_n;++i)
      _f[i]=_c[i];
  }
  else {
    F* f=_f;
    for (unsigned j=0;j<_n;++j) {
      const C* c=_c+j;
      for (unsigned i=0;i<_m;++i,c+=_n)
        f[i]=*c;
      f+=_ldf;
    }
  }
}

/** Copy _m x _n source matrix _f from FORTRAN-style to C-style array _c.
    \ingroup vc_lapack
    \arg _c 2-dimensional C-style array (rows first)
    \arg _f 2-dimensional FORTRAN style array (columns first) of
    dimension M x _n, where M<=_ldf, use _ldf=M if _ldf<M
    \arg _m, _n dimensions of _c and (submatrix of) _f (depending on _ldf)
    \arg _ldf leading dimension of _f,where M<=_ldf, use _ldf=M if _ldf<M
    Note: _c and _f must \a not overlap!
*/
template <typename F,typename C>
inline void
cp_f2c(const F* _f,C* _c,unsigned _m,unsigned _n,unsigned _ldf=0) {
  _ldf=std::max(_ldf,_m);

  if ((_m==1 && _ldf==1) || _n==1) {
    for (unsigned i=0;i<_m*_n;++i)
      _c[i]=_f[i];
  }
  else {
    const F* f=_f;
    for (unsigned j=0;j<_n;++j) {
      C* c=_c+j;
      for (unsigned i=0;i<_m;++i,c+=_n)
        *c=f[i];
      f+=_ldf;
    }
  }
}

//-----------------------------------------------------------------------------

/// wrap <tt>::ILAENV</tt> \ingroup vc_lapack
int ilaenv(int _spec,const char* _name,const char* _opts,
           int _n1=-1,int _n2=-1,int _n3=-1,int _n4=-1);

  //
  // NOTE: The current implementation is ** POTENTIALLY NOT PORTABLE ** !!!
  //       (Passing of variable length string parameters.)
  //

//-----------------------------------------------------------------------------

/// wrap ::DGELSD \ingroup vc_lapack
inline
int gelsd(int _m,int _n, int _nrhs,
          const double* _a,int _lda,
          double* _b,int _ldb, double* _s,
          double _rcond,int& _rank,
          double* _work,int _lwork,INTEGER* _iwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,rank,lwork=_lwork,info;

  DGELSD(&m,&n,&nrhs,_a,&lda,_b,&ldb,_s,&_rcond,&rank,_work,&lwork,_iwork,&info);

  _rank=rank;

  return info;
}

/// wrap ::SGELSD \ingroup lapack
inline
int gelsd(int _m,int _n, int _nrhs,
          const float* _a,int _lda,
          float* _b,int _ldb, float* _s,
          float _rcond,int& _rank,
          float* _work,int _lwork,INTEGER* _iwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,rank,lwork=_lwork,info;

  SGELSD(&m,&n,&nrhs,_a,&lda,_b,&ldb,_s,&_rcond,&rank,_work,&lwork,_iwork,&info);

  _rank=rank;

  return info;
}

/// query storage required by gelsd(): returns lwork
inline
int gelsd_lwork(int _m,int _n,int _nrhs) {
  double w;
  int irank, liwork;
  int info=lapack::gelsd(_m,_n,_nrhs, 0,_m, 0,_m, 0, 0.0, irank,&w,-1,&liwork);
  assert(info==0); use_nowarn(info);
  return int(w);
}

/// query storage required by gelsd(): size of iwork
inline
int gelsd_liwork(int _m,int _n,int /*_nrhs*/) {
  unsigned minmn=std::min(_m,_n);
  unsigned nlvl=unsigned(log2(float(minmn)))+1;
  return 3*minmn*nlvl+11*minmn;
}


/// wrap ::DGELSY \ingroup vc_lapack
inline
int gelsy(int _m,int _n, int _nrhs,
          double* _a,int _lda,
          double* _b,int _ldb,
          INTEGER* _jpvt,
          double _rcond,int& _rank,
          double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,rank,lwork=_lwork,info;

  DGELSY(&m,&n,&nrhs,_a,&lda,_b,&ldb,_jpvt,&_rcond,&rank,_work,&lwork,&info);

  _rank=rank;

  return info;
}

/// wrap ::SGELSY \ingroup vc_lapack
inline
int gelsy(int _m,int _n, int _nrhs,
          float* _a,int _lda,
          float* _b,int _ldb,
          INTEGER* _jpvt,
          float _rcond,int& _rank,
          float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,rank,lwork=_lwork,info;

  SGELSY(&m,&n,&nrhs,_a,&lda,_b,&ldb,_jpvt,&_rcond,&rank,_work,&lwork,&info);

  _rank=rank;

  return info;
}

/// query storage required by gelsy(): returns lwork
inline
int gelsy_lwork(int _m,int _n,int _nrhs) {
  double w;
  int irank;
  int info=lapack::gelsy(_m,_n,_nrhs, 0,_m, 0,_m, 0, 0.0, irank,&w,-1);
  assert(info==0); use_nowarn(info);
  return int(w);
}

//-----------------------------------------------------------------------------

/// wrap LU decomposition ::DGETRF \ingroup vc_lapack
inline
int getrf(int _m,int _n,double *_a,int _lda,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,info;
  DGETRF(&m,&n,_a,&lda,_ipiv,&info);
  return info;
}
/// wrap LU decomposition ::SGETRF \ingroup vc_lapack
inline
int getrf(int _m,int _n,float *_a,int _lda,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,info;
  SGETRF(&m,&n,_a,&lda,_ipiv,&info);
  return info;
}

/// wrap LU decomposition ::DGBTRF \ingroup vc_lapack
inline
int gbtrf(int _m,int _n,int _kl,int _ku,double *_ab,int _ldab,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,kl=_kl,ku=_ku,ldab=_ldab,info;
  DGBTRF(&m,&n,&kl,&ku,_ab,&ldab,_ipiv,&info);
  return info;
}
/// wrap LU decomposition ::DGBTRF \ingroup vc_lapack
inline
int gbtrf(int _m,int _n,int _kl,int _ku,float *_ab,int _ldab,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,kl=_kl,ku=_ku,ldab=_ldab,info;
  SGBTRF(&m,&n,&kl,&ku,_ab,&ldab,_ipiv,&info);
  return info;
}

/// wrap LU decomposition ::DGTTRF \ingroup vc_lapack
inline
int gttrf(int _n,double* _dl,double* _d,double* _du,double* _du2,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,info;
  DGTTRF(&n,_dl,_d,_du,_du2,_ipiv,&info);
  return info;
}
/// wrap LU decomposition ::DGTTRF \ingroup vc_lapack
inline
int gttrf(int _n,float* _dl,float* _d,float* _du,float* _du2,INTEGER* _ipiv) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,info;
  SGTTRF(&n,_dl,_d,_du,_du2,_ipiv,&info);
  return info;
}

/// wrap LU substitution ::DGETRS \ingroup vc_lapack
inline
int getrs(blas::TransposeFlag _trans,
          int _n,int _nrhs,const double *_a,int _lda,const INTEGER* _ipiv,
          double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  DGETRS(&trans,&n,&nrhs,_a,&lda,_ipiv,_b,&ldb,&info);
  return info;
}
/// wrap LU substitution ::SGETRS \ingroup vc_lapack
inline
int getrs(TransposeFlag _trans,
          int _n,int _nrhs,const float *_a,int _lda,const INTEGER* _ipiv,
          float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  SGETRS(&trans,&n,&nrhs,_a,&lda,_ipiv,_b,&ldb,&info);
  return info;
}


/// wrap LU substitution ::DGBTRS \ingroup vc_lapack
inline
int gbtrs(blas::TransposeFlag _trans,
          int _n,int _kl,int _ku,int _nrhs,const double *_a,int _lda,const INTEGER* _ipiv,
          double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,kl=_kl,ku=_ku,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  DGBTRS(&trans,&n,&kl,&ku,&nrhs,_a,&lda,_ipiv,_b,&ldb,&info);
  return info;
}
/// wrap LU substitution ::SGBTRS \ingroup vc_lapack
inline
int gbtrs(blas::TransposeFlag _trans,
          int _n,int _kl,int _ku,int _nrhs,const float *_a,int _lda,const INTEGER* _ipiv,
          float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,kl=_kl,ku=_ku,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  SGBTRS(&trans,&n,&kl,&ku,&nrhs,_a,&lda,_ipiv,_b,&ldb,&info);
  return info;
}

/// wrap LU substitution ::DGTTRS \ingroup vc_lapack
inline
int gttrs(blas::TransposeFlag _trans,
          int _n,int _nrhs,
          const double *_dl,const double* _d,const double* _du,const double* _du2,
          const INTEGER* _ipiv, double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  DGTTRS(&trans,&n,&nrhs,_dl,_d,_du,_du2,_ipiv,_b,&ldb,&info);
  return info;
}
/// wrap LU substitution ::SGTTRS \ingroup vc_lapack
inline
int gttrs(blas::TransposeFlag _trans,
          int _n,int _nrhs,
          const float *_dl,const float* _d,const float* _du,const float* _du2,
          const INTEGER* _ipiv, float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char trans=char(_trans);
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  SGTTRS(&trans,&n,&nrhs,_dl,_d,_du,_du2,_ipiv,_b,&ldb,&info);
  return info;
}


/// wrap LU solve tridiagonal ::DGTSV \ingroup vc_lapack
inline
int gtsv(int _n,int _nrhs,double *_dl,double* _d,double* _du,double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  DGTSV(&n,&nrhs,_dl,_d,_du,_b,&ldb,&info);
  return info;
}
/// wrap LU solve  tridiagonal ::SGTSV \ingroup vc_lapack
inline
int gtsv(int _n,int _nrhs,float *_dl,float* _d,float* _du,float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  SGTSV(&n,&nrhs,_dl,_d,_du,_b,&ldb,&info);
  return info;
}

/// wrap LU inverse ::DGETRI \ingroup vc_lapack
inline
int getri(int _n,double* _a,int _lda,INTEGER* _ipiv,double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,lda=_lda,lwork=_lwork,info;
  DGETRI(&n,_a,&lda,_ipiv,_work,&lwork,&info);
  return info;
}
/// wrap LU inverse ::SGETRI \ingroup vc_lapack
inline
int getri(int _n,float* _a,int _lda,INTEGER* _ipiv,float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,lda=_lda,lwork=_lwork,info;
  SGETRI(&n,_a,&lda,_ipiv,_work,&lwork,&info);
  return info;
}

//-----------------------------------------------------------------------------

/// wrap Cholesky decomposition ::DPOTRF \ingroup vc_lapack
inline
int potrf(UpperLowerFlag _uplo,int _n,double *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,info;
  DPOTRF(&uplo,&n,_a,&lda,&info);
  return info;
}
/// wrap Cholesky decomposition ::SPOTRF \ingroup vc_lapack
inline
int potrf(UpperLowerFlag _uplo,int _n,float *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,info;
  SPOTRF(&uplo,&n,_a,&lda,&info);
  return info;
}

/// wrap Cholesky decomposition ::DPBTRF \ingroup vc_lapack
inline
int pbtrf(UpperLowerFlag _uplo,int _n,int _k,double *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,k=_k,lda=_lda,info;
  DPBTRF(&uplo,&n,&k,_a,&lda,&info);
  return info;
}
/// wrap Cholesky decomposition ::SPBTRF \ingroup vc_lapack
inline
int pbtrf(UpperLowerFlag _uplo,int _n,int _k,float *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,k=_k,lda=_lda,info;
  SPBTRF(&uplo,&n,&k,_a,&lda,&info);
  return info;
}

/// wrap tridiagonal Cholesky decomposition ::DPTTRF \ingroup vc_lapack
inline
int pttrf(int _n,double *_d,double* _e) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,info;
  DPTTRF(&n,_d,_e,&info);
  return info;
}
/// wrap tridiagonal Cholesky decomposition ::SPTTRF \ingroup vc_lapack
inline
int pttrf(int _n,float *_d,float* _e) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,info;
  SPTTRF(&n,_d,_e,&info);
  return info;
}


/// wrap Cholesky decomposition ::DPPTRF \ingroup vc_lapack
inline
int pptrf(UpperLowerFlag _uplo,int _n,double *_ap) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,info;
  DPPTRF(&uplo,&n,_ap,&info);
  return info;
}
/// wrap Cholesky decomposition ::SPPTRF \ingroup vc_lapack
inline
int pptrf(UpperLowerFlag _uplo,int _n,float *_ap) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,info;
  SPPTRF(&uplo,&n,_ap,&info);
  return info;
}

/// wrap Cholesky substitution ::DPOTRS \ingroup vc_lapack
inline
int potrs(UpperLowerFlag _uplo,int _n,int _nrhs,const double *_a,int _lda,
          double* _b,int _ldb) {
 _VC_DBG_LAPACK_SCOPE()
   char uplo=char(_uplo);
  INTEGER n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  DPOTRS(&uplo,&n,&nrhs,_a,&lda,_b,&ldb,&info);
  return info;
}
/// wrap Cholesky substitution ::SPOTRS \ingroup vc_lapack
inline
int potrs(UpperLowerFlag _uplo,int _n,int _nrhs,const float *_a,int _lda,
          float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  SPOTRS(&uplo,&n,&nrhs,_a,&lda,_b,&ldb,&info);
  return info;
}

/// wrap Cholesky substitution ::DPBTRS \ingroup vc_lapack
inline
int pbtrs(UpperLowerFlag _uplo,int _n,int _k,int _nrhs,const double *_a,int _lda,
 double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,k=_k,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  DPBTRS(&uplo,&n,&k,&nrhs,_a,&lda,_b,&ldb,&info);
  return info;
}
/// wrap Cholesky substitution ::SPBTRS \ingroup vc_lapack
inline
int pbtrs(UpperLowerFlag _uplo,int _n,int _k,int _nrhs,const float *_a,int _lda,
          float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,k=_k,nrhs=_nrhs,lda=_lda,ldb=_ldb,info;
  SPBTRS(&uplo,&n,&k,&nrhs,_a,&lda,_b,&ldb,&info);
  return info;
}

/// wrap Cholesky substitution ::DPPTRS \ingroup vc_lapack
inline
int pptrs(UpperLowerFlag _uplo,int _n,int _nrhs,const double *_ap,
  double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  DPPTRS(&uplo,&n,&nrhs,_ap,_b,&ldb,&info);
  return info;
}
/// wrap Cholesky substitution ::SPPTRS \ingroup vc_lapack
inline
int pptrs(UpperLowerFlag _uplo,int _n,int _nrhs,const float *_ap,float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  SPPTRS(&uplo,&n,&nrhs,_ap,_b,&ldb,&info);
  return info;
}

/// wrap Cholesky inverse ::DPOTRI \ingroup vc_lapack
inline
int potri(UpperLowerFlag _uplo,int _n,double *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,info;
  DPOTRI(&uplo,&n,_a,&lda,&info);
  return info;
}
/// wrap Cholesky inverse ::SPOTRI \ingroup vc_lapack
inline
int potri(UpperLowerFlag _uplo,int _n,float *_a,int _lda) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,lda=_lda,info;
  SPOTRI(&uplo,&n,_a,&lda,&info);
  return info;
}

/// wrap Cholesky inverse ::DPPTRI \ingroup vc_lapack
inline
int pptri(UpperLowerFlag _uplo,int _n,double *_ap) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,info;
  DPPTRI(&uplo,&n,_ap,&info);
  return info;
}
/// wrap Cholesky inverse ::SPPTRI \ingroup vc_lapack
inline
int pptri(UpperLowerFlag _uplo,int _n,float *_ap) {
  _VC_DBG_LAPACK_SCOPE()
  char uplo=char(_uplo);
  INTEGER n=_n,info;
  SPPTRI(&uplo,&n,_ap,&info);
  return info;
}

/// wrap Cholesky tridiagonal solve ::DPTTRS \ingroup vc_lapack
inline
int pttrs(int _n,int _nrhs,const double* _d,const double* _e,double* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  DPTTRS(&n,&nrhs,_d,_e,_b,&ldb,&info);
  return info;
}
/// wrap Cholesky tridiagonal solve ::SPTTRS \ingroup vc_lapack
inline
int pttrs(int _n,int _nrhs,const float* _d,const float* _e,float* _b,int _ldb) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,nrhs=_nrhs,ldb=_ldb,info;
  SPTTRS(&n,&nrhs,_d,_e,_b,&ldb,&info);
  return info;
}

//-----------------------------------------------------------------------------

/// wrap QR decomposition ::DGEQRF \ingroup vc_lapack
inline
int geqrf(int _m,int _n,double* _a,int _lda,double* _tau,double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,lwork=_lwork,info;
  DGEQRF(&m,&n,_a,&lda,_tau,_work,&lwork,&info);
  return info;
}

/// wrap QR decomposition ::SGEQRF \ingroup vc_lapack
inline
int geqrf(int _m,int _n,float* _a,int _lda,float* _tau,float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,lwork=_lwork,info;
  SGEQRF(&m,&n,_a,&lda,_tau,_work,&lwork,&info);
  return info;
}

/// query storage required by geqrf(): returns lwork \ingroup vc_lapack
inline
int geqrf_lwork(int _m,int _n) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_m,lwork=-1,info;
  double w;
  DGEQRF(&m,&n,0,&lda,0,&w,&lwork,&info);
  assert(info==0);
  return int(w);
}

/// wrap non-blocked QR decomposition ::DGEQR2 \ingroup vc_lapack
inline
int geqr2(int _m,int _n,double* _a,int _lda,double* _tau,double* _work) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,info;
  DGEQR2(&m,&n,_a,&lda,_tau,_work,&info);
  return info;
}

/// wrap non-blocked QR decomposition ::SGEQR2 \ingroup vc_lapack
inline
int geqr2(int _m,int _n,float* _a,int _lda,float* _tau,float* _work) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER m=_m,n=_n,lda=_lda,info;
  SGEQR2(&m,&n,_a,&lda,_tau,_work,&info);
  return info;
}

//-----------------------------------------------------------------------------

/// wrap computation of eigenvalues ::DSYEV \ingroup vc_lapack
inline
int syev(int _jobz,UpperLowerFlag _uplo,int _n,double *_a,int _lda,double* _w,
         double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobz=char(_jobz), uplo=char(_uplo);
  assert(jobz=='N' || jobz=='V');
  INTEGER n=_n,lda=_lda,lwork=_lwork,info;
  DSYEV(&jobz,&uplo,&n,_a,&lda,_w,_work,&lwork,&info);
  return info;
}
/// wrap computation of eigenvalues ::SSYEV \ingroup vc_lapack
inline
int syev(int _jobz,UpperLowerFlag _uplo,int _n,float *_a,int _lda,float* _w,
         float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobz=char(_jobz), uplo=char(_uplo);
  assert(jobz=='N' || jobz=='V');
  INTEGER n=_n,lda=_lda,lwork=_lwork,info;
  SSYEV(&jobz,&uplo,&n,_a,&lda,_w,_work,&lwork,&info);
  return info;
}
/// query storage required by syev(): returns lwork \ingroup vc_lapack
inline
int syev_lwork(int _n) {
  double w;
  int info=lapack::syev('V',blas::Upper,_n,0,_n,0,&w,-1);
  assert(info==0); use_nowarn(info);
  return int(w);
}

/// wrap computation of eigenvalues ::DSYEVR \ingroup vc_lapack
inline
int syevr(int _jobz,int _range,UpperLowerFlag _uplo,int _n,double *_a,int _lda,
          double _vl,double _vu,int _il,int _iu,double _abstol,int& _m,
          double* _w,double* _z,int _ldz,INTEGER* _isuppz,
          double* _work,int _lwork,INTEGER* _iwork,int _liwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobz=char(_jobz), range=char(_range), uplo=char(_uplo);
  assert(jobz=='N' || jobz=='V');
  assert(range=='A' || range=='V' || range=='I');
  INTEGER n=_n,lda=_lda,il=_il,iu=_iu,m,ldz=_ldz,lwork=_lwork,liwork=_liwork,info;
  ::DSYEVR(&jobz,&range,&uplo,
           &n,_a,&lda,
           &_vl,&_vu,&il,&iu,&_abstol,
           &m,_w,_z,&ldz,_isuppz,_work,&lwork,_iwork,&liwork,&info);
  _m=m;
  return info;
}
/// wrap computation of eigenvalues ::SSYEVR \ingroup vc_lapack
inline
int syevr(int _jobz,int _range,UpperLowerFlag _uplo,int _n,float *_a,int _lda,
          float _vl,float _vu,int _il,int _iu,float _abstol,int& _m,
          float* _w,float* _z,int _ldz,INTEGER* _isuppz,
          float* _work,int _lwork,INTEGER* _iwork,int _liwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobz=char(_jobz), range=char(_range), uplo=char(_uplo);
  assert(jobz=='N' || jobz=='V');
  assert(range=='A' || range=='V' || range=='I');
  INTEGER n=_n,lda=_lda,il=_il,iu=_iu,m,ldz=_ldz,lwork=_lwork,liwork=_liwork,info;
  ::SSYEVR(&jobz,&range,&uplo,
           &n,_a,&lda,
           &_vl,&_vu,&il,&iu,&_abstol,
           &m,_w,_z,&ldz,_isuppz,_work,&lwork,_iwork,&liwork,&info);
  _m=m;
  return info;
}
/// query storage required by syev(): returns lwork, set_liwork \ingroup vc_lapack
inline
int syevr_lwork(int _n,int& _liwork) {
  double w;
  INTEGER iw;
  int m;
  int info=lapack::syevr('V','A',blas::Upper,_n,0,_n,
                         0.0,0.0,1,_n,0.0,m,
                         0,0,_n,0,
                         &w,-1,&iw,-1);
  assert(info==0); use_nowarn(info);
  _liwork=iw;
  return int(w);
}

//-----------------------------------------------------------------------------

/// wrap computation of eigenvalues ::DGEEV \ingroup vc_lapack
inline
int geev(int _jobvl,int _jobvr,int _n,double* _a,int _lda,double* _wr,double* _wi,
         double* _vl,int _ldvl,double* _vr,int _ldvr,double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobvl=_jobvl, jobvr=_jobvr;
  assert(jobvl=='N' || jobvl=='V');
  assert(jobvr=='N' || jobvr=='V');
  INTEGER n=_n,lda=_lda,ldvl=_ldvl,ldvr=_ldvr,lwork=_lwork,info;
  ::DGEEV(&jobvl,&jobvr,
          &n,_a,&lda,_wr,_wi,_vl,&ldvl,_vr,&ldvr,_work,&lwork,&info);
  return int(info);
}

/// wrap computation of eigenvalues ::DGEEV \ingroup vc_lapack
inline
int geev(int _jobvl,int _jobvr,int _n,float* _a,int _lda,float* _wr,float* _wi,
         float* _vl,int _ldvl,float* _vr,int _ldvr,float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobvl=_jobvl, jobvr=_jobvr;
  assert(jobvl=='N' || jobvl=='V');
  assert(jobvr=='N' || jobvr=='V');
  INTEGER n=_n,lda=_lda,ldvl=_ldvl,ldvr=_ldvr,lwork=_lwork,info;
  ::SGEEV(&jobvl,&jobvr,
          &n,_a,&lda,_wr,_wi,_vl,&ldvl,_vr,&ldvr,_work,&lwork,&info);
  return int(info);
}

/// query storage required by dgeev() \ingroup vc_lapack
inline
int geev_lwork(int _jobvl,int _jobvr,int _n) {
  double w;
  int info=lapack::geev(_jobvl,_jobvr,_n,0,_n,0,0,0,_n,0,_n,&w,-1);
  assert(info==0); use_nowarn(info);
  return int(w);
}

//-----------------------------------------------------------------------------

/// wrap SVD ::DGESVD \ingroup vc_lapack
inline
int gesvd(int _jobu,int _jobvt,int _m,int _n,double* _a,int _lda,
          double* _s,double* _u,int _ldu,double* _vt,int _ldvt,
          double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobu=char(_jobu), jobvt=char(_jobvt);
  INTEGER m=_m,n=_n,lda=_lda,ldu=_ldu,ldvt=_ldvt,lwork=_lwork,info;
  ::DGESVD(&jobu,&jobvt,&m,&n,_a,&lda,_s,_u,&ldu,_vt,&ldvt,_work,&lwork,&info);
  return info;
}
/// wrap SVD ::SGESVD \ingroup vc_lapack
inline
int gesvd(int _jobu,int _jobvt,int _m,int _n,float* _a,int _lda,
          float* _s,float* _u,int _ldu,float* _vt,int _ldvt,
          float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char jobu=char(_jobu), jobvt=char(_jobvt);
  INTEGER m=_m,n=_n,lda=_lda,ldu=_ldu,ldvt=_ldvt,lwork=_lwork,info;
  ::SGESVD(&jobu,&jobvt,&m,&n,_a,&lda,_s,_u,&ldu,_vt,&ldvt,_work,&lwork,&info);
  return info;
}

/// query storage required by gesvd(): returns lwork \ingroup vc_lapack
inline
int gesvd_lwork(int _m,int _n) {
  double w;
  int info=lapack::gesvd('A','A',_m,_n,0,_m,0,0,_m,0,_n,&w,-1);
  assert(info==0); use_nowarn(info);
  return int(w);
}

//-----------------------------------------------------------------------------

//
// missing: xLASR (sequence of Givens rotations)
//

//-----------------------------------------------------------------------------

/// wrap application of Householder reflector ::DLASR \ingroup vc_lapack
inline
void larf(blas::SideFlag _side,int _m,int _n,const double* _v,int _incv,double _tau,
          double* _c,int _ldc,double* _work) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side);
  INTEGER m=_m,n=_n,incv=_incv,ldc=_ldc;
  ::DLARF(&side,&m,&n,_v,&incv,&_tau,_c,&ldc,_work);
}

/// wrap application of Householder reflector ::SLASR \ingroup vc_lapack
inline
void larf(blas::SideFlag _side,int _m,int _n,const float* _v,int _incv,float _tau,
          float* _c,int _ldc,float* _work) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side);
  INTEGER m=_m,n=_n,incv=_incv,ldc=_ldc;
  ::SLARF(&side,&m,&n,_v,&incv,&_tau,_c,&ldc,_work);
}

//-----------------------------------------------------------------------------

/// wrap forming of block reflectors <tt>::DLARFT</tt> \ingroup vc_lapack
inline
void larft(int _direct,int _strorev,int _n,int _k,double* _v,int _ldv,
           const double* _tau,double* _t,int _ldt) {
  _VC_DBG_LAPACK_SCOPE()
  char direct=char(_direct), storev=char(_strorev);
  INTEGER n=_n,k=_k,ldv=_ldv,ldt=_ldt;
  ::DLARFT(&direct,&storev,&n,&k,_v,&ldv,_tau,_t,&ldt);
}

/// wrap forming of block reflectors ::SLARFT \ingroup vc_lapack
inline
void larft(int _direct,int _strorev,int _n,int _k,float* _v,int _ldv,
           const float* _tau,float* _t,int _ldt) {
  _VC_DBG_LAPACK_SCOPE()
  char direct=char(_direct), storev=char(_strorev);
  INTEGER n=_n,k=_k,ldv=_ldv,ldt=_ldt;
  ::SLARFT(&direct,&storev,&n,&k,_v,&ldv,_tau,_t,&ldt);
}


/// wrap application of block reflectors <tt>::DLARFB</tt> \ingroup vc_lapack
inline
void larfb(blas::SideFlag _side,blas::TransposeFlag _trans,int _direct,int _strorev,
           int _m,int _n,int _k,const double* _v,int _ldv,
           const double* _t,int _ldt,double* _c,int _ldc,
           double* _work,int _ldwork) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side), trans=char(_trans);
  char direct=char(_direct), storev=char(_strorev);
  INTEGER m=_m,n=_n,k=_k,ldv=_ldv,ldt=_ldt,ldc=_ldc,ldwork=_ldwork;
  ::DLARFB(&side,&trans,&direct,&storev,
           &m,&n,&k,_v,&ldv,_t,&ldt,_c,&ldc,_work,&ldwork);
}

/// wrap application of block reflectors ::SLARFB \ingroup vc_lapack
inline
void larfb(blas::SideFlag _side,blas::TransposeFlag _trans,int _direct,int _strorev,
           int _m,int _n,int _k,const float* _v,int _ldv,
           const float* _t,int _ldt,float* _c,int _ldc,
           float* _work,int _ldwork) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side), trans=char(_trans);
  char direct=char(_direct), storev=char(_strorev);
  INTEGER m=_m,n=_n,k=_k,ldv=_ldv,ldt=_ldt,ldc=_ldc,ldwork=_ldwork;
  ::SLARFB(&side,&trans,&direct,&storev,
           &m,&n,&k,_v,&ldv,_t,&ldt,_c,&ldc,_work,&ldwork);
}

//-----------------------------------------------------------------------------

/// wrap generation of Householder reflector ::DLARFG \ingroup vc_lapack
inline
void larfg(int _n,double& _alpha,double* _x,int _incx,double& _tau) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,incx=_incx;
  ::DLARFG(&n,&_alpha,_x,&incx,&_tau);
}

/// wrap generation of Householder reflector `::DLARFG` \ingroup vc_lapack
inline
void larfg(int _n,float& _alpha,float* _x,int _incx,float& _tau) {
  _VC_DBG_LAPACK_SCOPE()
  INTEGER n=_n,incx=_incx;
  ::SLARFG(&n,&_alpha,_x,&incx,&_tau);
}

//-----------------------------------------------------------------------------

/// wrap multiplication with Householder reflectors ::DORMQR \ingroup vc_lapack
inline
int ormqr(blas::SideFlag _side,blas::TransposeFlag _trans,
          int _m,int _n,int _k,const double* _a,int _lda,const double* _tau,
          double* _c,int _ldc,double* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side), trans=char(_trans);
  INTEGER m=_m,n=_n,k=_k,lda=_lda,ldc=_ldc,lwork=_lwork,info;
  ::DORMQR(&side,&trans,&m,&n,&k,_a,&lda,_tau,_c,&ldc,_work,&lwork,&info);
  return info;
}

/// wrap multiplication with Householder reflectors ::DORMQR \ingroup vc_lapack
inline
int ormqr(blas::SideFlag _side,blas::TransposeFlag _trans,
          int _m,int _n,int _k,const float* _a,int _lda,const float* _tau,
          float* _c,int _ldc,float* _work,int _lwork) {
  _VC_DBG_LAPACK_SCOPE()
  char side=char(_side), trans=char(_trans);
  INTEGER m=_m,n=_n,k=_k,lda=_lda,ldc=_ldc,lwork=_lwork,info;
  ::SORMQR(&side,&trans,&m,&n,&k,_a,&lda,_tau,_c,&ldc,_work,&lwork,&info);
  return info;
}

//-----------------------------------------------------------------------------

/// wrap matrix (1,fro,inf) norm <tt>::DLANGE</tt> \ingroup vc_lapack
inline
double lange(int _norm,int _m,int _n,const double* _a,int _lda,double* _work) {
  _VC_DBG_LAPACK_SCOPE()
  char norm=char(_norm);
  INTEGER m=_m,n=_n,lda=_lda;
  return ::DLANGE(&norm,&m,&n,_a,&lda,_work);
}

/// wrap matrix (1,fro,inf) norm ::SLANGE \ingroup vc_lapack
inline
float lange(int _norm,int _m,int _n,const float* _a,int _lda,float* _work) {
  _VC_DBG_LAPACK_SCOPE()
  char norm=char(_norm);
  INTEGER m=_m,n=_n,lda=_lda;
  return ::SLANGE(&norm,&m,&n,_a,&lda,_work);
}

//-----------------------------------------------------------------------------

} // namespace lapack
} // namespace math
} // namespace newlib

#endif // VC_MATH_LAPACK_HH
