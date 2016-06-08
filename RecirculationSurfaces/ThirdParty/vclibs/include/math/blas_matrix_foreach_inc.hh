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

    \arg Iterators over non-zero matrix elements in any order.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {

//=============================================================================

/** Iterate over stored (non-zero) elements.
    \internal \ingroup vc_blas_mat
    Call _op(_a) for each element.
*/
template <typename T,typename Op,typename A>
struct Foreach {
private:
  Foreach(const A& _ma,T* _a,Op _op);
  // undefined
};
/** Same as Foreach but respect TransposeFlag.
    \internal \ingroup vc_blas_mat
    Call _op(_a) for each element.
*/
template <typename T,typename Op,typename A>
struct ForeachT {
private:
  ForeachT(const A& _ma,T* _a,Op _op);
  // undefined
};
/** Iterate over stored (non-zero) elements of same matrix types (e.g., GE).
    \internal \ingroup vc_blas_mat
    Call _op(_a,_b) for each element.
*/
  template <typename T,typename Op,typename A,typename B>
struct Foreach2 {
private:
  Foreach2(const A& _ma,T* _a,const B& _mb,T* _b,Op _op);
  // undefined
};
/** Iterate over stored (non-zero) elements of same matrix types (e.g., GE).
    \internal \ingroup vc_blas_mat
    Call _op(_a,_b,_c) for each element.
*/
template <typename T,typename Op,typename A,typename B,typename C>
struct Foreach3 {
private:
  Foreach3(const A& _ma,T* _a,const B& _mb,T* _b,const C& _mc,T* _c,Op _op);
  // undefined
};

# ifndef DOXYGEN_SKIP

template <typename T,typename Op>
struct Assign11 {
  Assign11(Op _op) : op(_op) {}
  Op op;
  void operator()(T& _a) { _a=op(_a);  }
};

template <typename T>
struct Assign11Call {
  Assign11Call(T (*_op)(T)) : op(_op) {}
  T (*op)(T);
  void operator()(T& _a) { _a=op(_a);  }
};

template <typename T,typename Op>
struct Assign21 {
  Assign21(Op _op) : op(_op) {}
  Op op;
  void operator()(T& _a,const T& _b) { _a=op(_b);  }
};

template <typename T>
struct Assign21Call {
  Assign21Call(T (*_op)(T)) : op(_op) {}
  T (*op)(T);
  void operator()(T& _a,const T& _b) { _a=op(_b);  }
};

template <typename T,typename Op>
struct Assign22 {
  Assign22(Op _op) : op(_op) {}
  Op op;
  void operator()(T& _a,const T& _b) { _a=op(_a,_b);  }
};

template <typename T>
struct Assign22Call {
  Assign22Call(T (*_op)(T,T)) : op(_op) {}
  T (*op)(T,T);
  void operator()(T& _a,const T& _b) { _a=op(_a,_b);  }
};

template <typename T,typename Op>
struct Assign32 {
  Assign32(Op _op) : op(_op) {}
  Op op;
  void operator()(T& _a,const T& _b,const T& _c) { _a=op(_b,_c);  }
};

template <typename T>
struct Assign32Call {
  Assign32Call(T (*_op)(T,T)) : op(_op) {}
  T (*op)(T,T);
  void operator()(T& _a,const T& _b,const T& _c) { _a=op(_b,_c);  }
};

template <typename T>
struct Swap {
  void operator()(T& _a,T& _b) { T t=_a; _a=_b; _b=t;  }
};

template <typename S,typename T>
struct CopyFromArray {
  CopyFromArray(const S* _s) : s(_s) {}
  void operator()(T& _b) { _b=T(*s++); }
  const S* s;
};

template <typename S,typename T>
struct CopyToArray {
  CopyToArray(S* _s) : s(_s) {}
  void operator()(T& _b) { *s++=S(_b); }
  S* s;
};

# endif // DOXYGEN_SKIP

/** Assign array _ary elementwise to matrix.
    \ingroup vc_blas
 */
template <typename S,typename T,typename A>
void cp(const S* _ary,const matrix_reference_t<T,A>& _a) {
  CopyFromArray<S,T> cpop(_ary);
  ForeachT<T,CopyFromArray<S,T>,A>(_a.matrix(),_a.data(),cpop);
}
/** Assign matrix elements to array _ary.
    \ingroup vc_blas
 */
template <typename S,typename T,typename A>
void cp(const matrix_const_reference_t<T,A>& _a,const S* _ary) {
  CopyToArray<S,T> cpop(_ary);
  ForeachT<T,CopyToArray<S,T>,A>(_a.matrix(),(T*) _a.data(),cpop);
}

/** Assign each element _a(ij)=_op(_a(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat).
*/
template <typename T,typename Op,typename A>
void map_f(Op _op,const matrix_reference_t<T,A>& _a) {
  Assign11<T,Op> assign(_op);
  Foreach<T,Assign11<T,Op>,A>(_a.matrix(),_a.data(),assign);
}
template <typename T,typename A>
void map_f(T (*_op)(T),const matrix_reference_t<T,A>& _a) {
  Assign11Call<T> assign(_op);
  Foreach<T,Assign11Call<T>,A>(_a.matrix(),_a.data(),assign);
}

/** Assign each element _a(ij)=_op(_b(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f2(Op _op,
            const matrix_reference_t<T,A>& _a,
            const matrix_const_reference_t<T,B>& _b) {
  Assign21<T,Op> assign(_op);
  Foreach2<T,Assign21<T,Op>,A,B>(_a.matrix(),_a.data(),
                                 _b.matrix(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f2(T (*_op)(T),
            const matrix_reference_t<T,A>& _a,
            const matrix_const_reference_t<T,B>& _b) {
  Assign21Call<T> assign(_op);
  Foreach2<T,Assign21Call<T>,A,B>(_a.matrix(),_a.data(),
                                  _b.matrix(),(T*) _b.data(),assign);
}
/** Assign each element _a(ij)=_op(_b(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f2(Op _op,
            const matrix_reference_t<T,A>& _a,
            const matrix_reference_t<T,B>& _b) {
  Assign21<T,Op> assign(_op);
  Foreach2<T,Assign21<T,Op>,A,B>(_a.matrix(),_a.data(),
                                 _b.matrix(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f2(T (*_op)(T),
            const matrix_reference_t<T,A>& _a,
            const matrix_reference_t<T,B>& _b) {
  Assign21Call<T> assign(_op);
  Foreach2<T,Assign21Call<T>,A,B>(_a.matrix(),_a.data(),
                                  _b.matrix(),(T*) _b.data(),assign);
}

/** Assign each element _a(ij)=_op(_a(ij),_b(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f12(Op _op,
             const matrix_reference_t<T,A>& _a,
             const matrix_const_reference_t<T,B>& _b) {
  Assign22<T,Op> assign(_op);
  Foreach2<T,Assign22<T,Op>,A,B>(_a.matrix(),_a.data(),
                                 _b.matrix(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f12(T (*_op)(T,T),
             const matrix_reference_t<T,A>& _a,
             const matrix_const_reference_t<T,B>& _b) {
  Assign22Call<T> assign(_op);
  Foreach2<T,Assign22Call<T>,A,B>(_a.matrix(),_a.data(),
                                  _b.matrix(),(T*) _b.data(),assign);
}
/** Assign each element _a(ij)=_op(_a(ij),_b(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f12(Op _op,
             const matrix_reference_t<T,A>& _a,
             const matrix_reference_t<T,B>& _b) {
  Assign22<T,Op> assign(_op);
  Foreach2<T,Assign22<T,Op>,A,B>(_a.matrix(),_a.data(),
                                 _b.matrix(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f12(T (*_op)(T,T),
             const matrix_reference_t<T,A>& _a,
             const matrix_reference_t<T,B>& _b) {
  Assign22Call<T> assign(_op);
  Foreach2<T,Assign22Call<T>,A,B>(_a.matrix(),_a.data(),
                                  _b.matrix(),(T*) _b.data(),assign);
}

/** Assign each element _a(ij)=_op(_b(ij),_c(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B,typename C>
void map_f23(Op _op,
             const matrix_reference_t<T,A>& _a,
             const matrix_const_reference_t<T,B>& _b,
             const matrix_const_reference_t<T,C>& _c) {
  Assign32<T,Op> assign(_op);
  Foreach3<T,Assign32<T,Op>,A,B,C>(_a.matrix(),_a.data(),
                                   _b.matrix(),(T*) _b.data(),
                                   _c.matrix(),(T*) _c.data(),assign);
}
template <typename T,typename A,typename B,typename C>
void map_f23(T (*_op)(T,T),
             const matrix_reference_t<T,A>& _a,
             const matrix_const_reference_t<T,B>& _b,
             const matrix_const_reference_t<T,C>& _c) {
  Assign32Call<T> assign(_op);
  Foreach3<T,Assign32Call<T>,A,B,C>(_a.matrix(),_a.data(),
                                    _b.matrix(),(T*) _b.data(),
                                    _c.matrix(),(T*) _c.data(),assign);
}
/** Assign each element _a(ij)=_op(_b(ij),_c(ij)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B,typename C>
void map_f23(Op _op,
             const matrix_reference_t<T,A>& _a,
             const matrix_reference_t<T,B>& _b,
             const matrix_reference_t<T,C>& _c) {
  Assign32<T,Op> assign(_op);
  Foreach3<T,Assign32<T,Op>,A,B,C>(_a.matrix(),_a.data(),
                                   _b.matrix(),(T*) _b.data(),
                                   _c.matrix(),(T*) _c.data(),assign);
}
template <typename T,typename A,typename B,typename C>
void map_f23(T (*_op)(T,T),
             const matrix_reference_t<T,A>& _a,
             const matrix_reference_t<T,B>& _b,
             const matrix_reference_t<T,C>& _c) {
  Assign32Call<T> assign(_op);
  Foreach3<T,Assign32Call<T>,A,B,C>(_a.matrix(),_a.data(),
                                    _b.matrix(),(T*) _b.data(),
                                    _c.matrix(),(T*) _c.data(),assign);
}

/// swap contents of _a and _b (same matrix types and dimensions) \ingroup vc_blas
template <typename T,typename A,typename B>
void swap(const matrix_reference_t<T,A>& _a,
          const matrix_reference_t<T,B>& _b) {
  Foreach2<T,Swap<T>,A,B>(_a.matrix(),_a.data(),
                          _b.matrix(),_b.data(),Swap<T>());
}

//-----------------------------------------------------------------------------

/** Assign each element _a(i)=_op(_a(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat).
*/
template <typename T,typename Op,typename A>
void map_f(Op _op,const vector_reference_t<T,A>& _a) {
  Assign11<T,Op> assign(_op);
  Foreach<T,Assign11<T,Op>,A>(_a.vector(),_a.data(),assign);
}
template <typename T,typename A>
void map_f(T (*_op)(T),const vector_reference_t<T,A>& _a) {
  Assign11Call<T> assign(_op);
  Foreach<T,Assign11Call<T>,A>(_a.vector(),_a.data(),assign);
}


/** Assign each element _a(i)=_op(_b(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f2(Op _op,
            const vector_reference_t<T,A>& _a,
            const vector_const_reference_t<T,B>& _b) {
  Assign21<T,Op> assign(_op);
  Foreach2<T,Assign21<T,Op>,A,B>(_a.vector(),_a.data(),
                                 _b.vector(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f2(T (*_op)(T),
            const vector_reference_t<T,A>& _a,
            const vector_const_reference_t<T,B>& _b) {
  Assign21Call<T> assign(_op);
  Foreach2<T,Assign21Call<T>,A,B>(_a.vector(),_a.data(),
                                  _b.vector(),(T*) _b.data(),assign);
}
/** Assign each element _a(i)=_op(_b(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f2(Op _op,
            const vector_reference_t<T,A>& _a,
            const vector_reference_t<T,B>& _b) {
  Assign21<T,Op> assign(_op);
  Foreach2<T,Assign21<T,Op>,A,B>(_a.vector(),_a.data(),
                                 _b.vector(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f2(T (*_op)(T),
            const vector_reference_t<T,A>& _a,
            const vector_reference_t<T,B>& _b) {
  Assign21Call<T> assign(_op);
  Foreach2<T,Assign21Call<T>,A,B>(_a.vector(),_a.data(),
                                  _b.vector(),(T*) _b.data(),assign);
}

/** Assign each element _a(i)=_op(_a(i),_b(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f12(Op _op,
             const vector_reference_t<T,A>& _a,
             const vector_const_reference_t<T,B>& _b) {
  Assign22<T,Op> assign(_op);
  Foreach2<T,Assign22<T,Op>,A,B>(_a.vector(),_a.data(),
                                 _b.vector(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f12(T (*_op)(T,T),
             const vector_reference_t<T,A>& _a,
             const vector_const_reference_t<T,B>& _b) {
  Assign22Call<T> assign(_op);
  Foreach2<T,Assign22Call<T>,A,B>(_a.vector(),_a.data(),
                                  _b.vector(),(T*) _b.data(),assign);
}
/** Assign each element _a(i)=_op(_a(i),_b(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B>
void map_f12(Op _op,
             const vector_reference_t<T,A>& _a,
             const vector_reference_t<T,B>& _b) {
  Assign22<T,Op> assign(_op);
  Foreach2<T,Assign22<T,Op>,A,B>(_a.vector(),_a.data(),
                                 _b.vector(),(T*) _b.data(),assign);
}
template <typename T,typename A,typename B>
void map_f12(T (*_op)(T,T),
             const vector_reference_t<T,A>& _a,
             const vector_reference_t<T,B>& _b) {
  Assign22Call<T> assign(_op);
  Foreach2<T,Assign22Call<T>,A,B>(_a.vector(),_a.data(),
                                  _b.vector(),(T*) _b.data(),assign);
}

/** Assign each element _a(i)=_op(_b(i),_c(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B,typename C>
void map_f23(Op _op,
             const vector_reference_t<T,A>& _a,
             const vector_const_reference_t<T,B>& _b,
             const vector_const_reference_t<T,C>& _c) {
  Assign32<T,Op> assign(_op);
  Foreach3<T,Assign32<T,Op>,A,B,C>(_a.vector(),_a.data(),
                                   _b.vector(),(T*) _b.data(),
                                   _c.vector(),(T*) _c.data(),assign);
}
template <typename T,typename A,typename B,typename C>
void map_f23(T (*_op)(T,T),
             const vector_reference_t<T,A>& _a,
             const vector_const_reference_t<T,B>& _b,
             const vector_const_reference_t<T,C>& _c) {
  Assign32Call<T> assign(_op);
  Foreach3<T,Assign32Call<T>,A,B,C>(_a.vector(),_a.data(),
                                    _b.vector(),(T*) _b.data(),
                                    _c.vector(),(T*) _c.data(),assign);
}
/** Assign each element _a(i)=_op(_b(i),_c(i)).
    \ingroup vc_blas
    Processes only stored elements (no implicit zeros or ones, e.g.,
    in tr_mat). _b may be trans()posed for ge_mat.
*/
template <typename T,typename Op,typename A,typename B,typename C>
void map_f23(Op _op,
             const vector_reference_t<T,A>& _a,
             const vector_reference_t<T,B>& _b,
             const vector_reference_t<T,C>& _c) {
  Assign32<T,Op> assign(_op);
  Foreach3<T,Assign32<T,Op>,A,B,C>(_a.vector(),_a.data(),
                                   _b.vector(),(T*) _b.data(),
                                   _c.vector(),(T*) _c.data(),assign);
}
template <typename T,typename A,typename B,typename C>
void map_f23(T (*_op)(T,T),
             const vector_reference_t<T,A>& _a,
             const vector_reference_t<T,B>& _b,
             const vector_reference_t<T,C>& _c) {
  Assign32Call<T> assign(_op);
  Foreach3<T,Assign32Call<T>,A,B,C>(_a.vector(),_a.data(),
                                    _b.vector(),(T*) _b.data(),
                                    _c.vector(),(T*) _c.data(),assign);
}

//=============================================================================

# ifndef DOXYGEN_SKIP

template <typename T,typename Op,int N,int INC>
struct Foreach<T,Op,vec<N,INC> > {
  Foreach(const vec<N,INC>& _va,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();
    int k=_va.inc();
    for (int j=0;j<_va.n();++j)
      _op(_a[j*k]);
  }
};

template <typename T,typename Op,int N,int INC>
struct ForeachT<T,Op,vec<N,INC> > : Foreach<T,Op,vec<N,INC> > {
  ForeachT(const vec<N,INC>& _va,T* _a,Op _op)
    : Foreach<T,Op,vec<N,INC> >(_va,_a,_op) {}
};

template <typename T,typename Op,
          int NA,int INCA,int NB,int INCB>
struct Foreach2<T,Op,vec<NA,INCA>,vec<NB,INCB> > {
  Foreach2(const vec<NA,INCA>& _va,T* _a,
           const vec<NB,INCB>& _vb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();
    assert(_va.n()==_vb.n());
    int n=_va.n(), ka=_va.inc(), kb=_vb.inc();
    for (int j=0;j<n;++j)
      _op(_a[j*ka],_b[j*kb]);
  }
};

template <typename T,typename Op,
          int NA,int INCA,int NB,int INCB,int NC,int INCC>
struct Foreach3<T,Op,vec<NA,INCA>,vec<NB,INCB>,vec<NC,INCC> > {
  Foreach3(const vec<NA,INCA>& _va,T* _a,
           const vec<NB,INCB>& _vb,T* _b,
           const vec<NC,INCC>& _vc,T* _c,
           Op _op) {
    _VC_DBG_BLAS_SCOPE();
    assert(_va.n()==_vb.n());
    assert(_va.n()==_vb.n());
    int n=_va.n(), ka=_va.inc(), kb=_vb.inc(), kc=_vc.inc();
    for (int j=0;j<n;++j)
      _op(_a[j*ka],_b[j*kb],_c[j*kc]);
  }
};

//-----------------------------------------------------------------------------

template <typename T,typename Op,
          TransposeFlag TRANS,int M,int N,int LD>
struct Foreach<T,Op,ge_mat<TRANS,M,N,LD> > {
  Foreach(const ge_mat<TRANS,M,N,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<_ma.m();++i)
        _op(_a[j*_ma.ld()+i]);
  }
};

template <typename T,typename Op,
          TransposeFlag TRANS,int M,int N,int LD>
struct ForeachT<T,Op,ge_mat<TRANS,M,N,LD> > {
  ForeachT(const ge_mat<TRANS,M,N,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    if (TRANS==NoT)
      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<_ma.m();++i)
          _op(_a[j*_ma.ld()+i]);
    else
      for (int i=0;i<_ma.m();++i)
        for (int j=0;j<_ma.n();++j)
          _op(_a[j*_ma.ld()+i]);
  }
};

template <typename T,typename Op,
          int MA,int NA,int LDA,
          TransposeFlag TRANSB,int MB,int NB,int LDB>
struct Foreach2<T,Op,ge_mat<NoT,MA,NA,LDA>,ge_mat<TRANSB,MB,NB,LDB> > {
  Foreach2(const ge_mat<NoT,MA,NA,LDA>& _ma,T* _a,
           const ge_mat<TRANSB,MB,NB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    if (TRANSB==NoT) {
      assert(_ma.m()==_mb.m());
      assert(_ma.n()==_mb.n());

      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<_ma.m();++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
    }
    else {
      assert(_ma.n()==_mb.m());
      assert(_ma.m()==_mb.n());

      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<_ma.m();++i)
          _op(_a[j*_ma.ld()+i],_b[i*_mb.ld()+j]);
    }
  }
};

template <typename T,typename Op,
          int MA,int NA,int LDA,
          TransposeFlag TRANSB,int MB,int NB,int LDB,
          int MC,int NC,int LDC>
struct Foreach3<T,Op,
                ge_mat<NoT,MA,NA,LDA>,
                ge_mat<TRANSB,MB,NB,LDB>,
                ge_mat<NoT,MC,NC,LDC> > {
  Foreach3(const ge_mat<NoT,MA,NA,LDA>& _ma,T* _a,
           const ge_mat<TRANSB,MB,NB,LDB>& _mb,T* _b,
           const ge_mat<NoT,MC,NC,LDC>& _mc,T* _c,
           Op _op) {
    _VC_DBG_BLAS_SCOPE();
    assert(_ma.m()==_mc.m());
    assert(_ma.n()==_mc.n());

    if (TRANSB==NoT) {
      assert(_ma.m()==_mb.m());
      assert(_ma.n()==_mb.n());

      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<_ma.m();++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mc.ld()+i]);
    }
    else {
      assert(_ma.n()==_mb.m());
      assert(_ma.m()==_mb.n());

      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<_ma.m();++i)
          _op(_a[j*_ma.ld()+i],_b[i*_mb.ld()+j],_c[j*_mc.ld()+i]);
    }
  }
};

//-----------------------------------------------------------------------------

template <typename T,typename Op,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
struct Foreach<T,Op,gb_mat<TRANS,M,N,KL,KU,LD> > {
  Foreach(const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    int m=_ma.kl()+_ma.ku()+1; // ignore extra work
    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<m;++i)
        _op(_a[j*_ma.ld()+i]);
  }
};

template <typename T,typename Op,
          TransposeFlag TRANS,int M,int N,int KL,int KU,int LD>
struct ForeachT<T,Op,gb_mat<TRANS,M,N,KL,KU,LD> > {
  ForeachT(const gb_mat<TRANS,M,N,KL,KU,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    int m=_ma.kl()+_ma.ku()+1; // ignore extra work
    if (TRANS==NoT)
      for (int j=0;j<_ma.n();++j)
        for (int i=0;i<m;++i)
          _op(_a[j*_ma.ld()+i]);
    else
      for (int i=0;i<m;++i)
        for (int j=0;j<_ma.n();++j)
          _op(_a[j*_ma.ld()+i]);
  }
};

template <typename T,typename Op,TransposeFlag TRANS,
          int MA,int NA,int KLA,int KUA,int LDA,
          int MB,int NB,int KLB,int KUB,int LDB>
struct Foreach2<T,Op,gb_mat<TRANS,MA,NA,KLA,KUA,LDA>,gb_mat<TRANS,MB,NB,KLB,KUB,LDB> > {
  Foreach2(const gb_mat<TRANS,MA,NA,KLA,KUA,LDA>& _ma,T* _a,
           const gb_mat<TRANS,MB,NB,KLB,KUB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.ku()==_mb.ku());
    assert(_ma.kl()==_mb.kl());

    int m=_ma.kl()+_ma.ku()+1; // ignore extra work

    assert(_ma.m()==_mb.m());
    assert(_ma.n()==_mb.n());

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<m;++i)
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
  }
};

template <typename T,typename Op,TransposeFlag TRANS,
          int MA,int NA,int KLA,int KUA,int LDA,
          int MB,int NB,int KLB,int KUB,int LDB,
          int MC,int NC,int KLC,int KUC,int LDC>
struct Foreach3<T,Op,
                gb_mat<TRANS,MA,NA,KLA,KUA,LDA>,
                gb_mat<TRANS,MB,NB,KLB,KUB,LDB>,
                gb_mat<TRANS,MC,NC,KLC,KUC,LDC> > {
  Foreach3(const gb_mat<TRANS,MA,NA,KLA,KUA,LDA>& _ma,T* _a,
           const gb_mat<TRANS,MB,NB,KLB,KUB,LDB>& _mb,T* _b,
           const gb_mat<TRANS,MC,NC,KLC,KUC,LDC>& _mc,T* _c,
           Op _op) {
    _VC_DBG_BLAS_SCOPE();
    assert(_ma.m()==_mc.m());
    assert(_ma.n()==_mc.n());

    assert(_ma.m()==_mb.m());
    assert(_ma.n()==_mb.n());

    assert(_ma.ku()==_mb.ku());
    assert(_ma.kl()==_mb.kl());

    assert(_ma.ku()==_mc.ku());
    assert(_ma.kl()==_mc.kl());

    int m=_ma.kl()+_ma.ku()+1; // ignore extra work

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<m;++i)
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mc.ld()+i]);
  }
};

//-----------------------------------------------------------------------------

template <typename T,typename Op,UpperLowerFlag UPLO,int N,int LD>
struct Foreach<T,Op,sy_mat<UPLO,N,LD> > {
  Foreach(const sy_mat<UPLO,N,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          _op(_a[j*_ma.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
          _op(_a[j*_ma.ld()+i]);
    }
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,int N,int LD>
struct ForeachT<T,Op,sy_mat<UPLO,N,LD> > {
  ForeachT(const sy_mat<UPLO,N,LD>& _ma,T* _a,Op _op) {
    Foreach<T,Op,sy_mat<UPLO,N,LD> >(_ma,_a,_op);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,int LDA,
          int NB,int LDB>
struct Foreach2<T,Op,sy_mat<UPLO,NA,LDA>,sy_mat<UPLO,NB,LDB> > {
  Foreach2(const sy_mat<UPLO,NA,LDA>& _ma,T* _a,
           const sy_mat<UPLO,NB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
    }
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,int LDA,
          int NB,int LDB,
          int NC,int LDC>
struct Foreach3<T,Op,
                sy_mat<UPLO,NA,LDA>,
                sy_mat<UPLO,NB,LDB>,
                sy_mat<UPLO,NC,LDC> > {
  Foreach3(const sy_mat<UPLO,NA,LDA>& _ma,T* _a,
           const sy_mat<UPLO,NB,LDB>& _mb,T* _b,
           const sy_mat<UPLO,NC,LDC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
          _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
    }
  }
};


//-----------------------------------------------------------------------------

template <typename T,typename Op,UpperLowerFlag UPLO,int N,int K,int LD>
struct Foreach<T,Op,sb_mat<UPLO,N,K,LD> > {
  Foreach(const sb_mat<UPLO,N,K,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<=_ma.k();++i) // ignore the extra work
        _op(_a[j*_ma.ld()+i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,int N,int K,int LD>
struct ForeachT<T,Op,sb_mat<UPLO,N,K,LD> > {
  ForeachT(const sb_mat<UPLO,N,K,LD>& _ma,T* _a,Op _op) {
    Foreach<T,Op,sb_mat<UPLO,N,K,LD> >(_ma,_a,_op);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,int KA,int LDA,
          int NB,int KB,int LDB>
struct Foreach2<T,Op,sb_mat<UPLO,NA,KA,LDA>,sb_mat<UPLO,NB,KB,LDB> > {
  Foreach2(const sb_mat<UPLO,NA,KA,LDA>& _ma,T* _a,
           const sb_mat<UPLO,NB,KB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.k()==_mb.k());


    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<=_ma.k();++i) // ignore the extra work
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,int KA,int LDA,
          int NB,int KB,int LDB,
          int NC,int KC,int LDC>
struct Foreach3<T,Op,
                sb_mat<UPLO,NA,KA,LDA>,
                sb_mat<UPLO,NB,KB,LDB>,
                sb_mat<UPLO,NC,KC,LDC> > {
  Foreach3(const sb_mat<UPLO,NA,KA,LDA>& _ma,T* _a,
           const sb_mat<UPLO,NB,KB,LDB>& _mb,T* _b,
           const sb_mat<UPLO,NC,KC,LDC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());
    assert(_ma.k()==_mb.k());
    assert(_ma.k()==_mc.k());

    for (int j=0;j<_ma.n();++j)
      for (int i=0;i<=_ma.k();++i) // ignore the extra work
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
  }
};


//-----------------------------------------------------------------------------

template <typename T,typename Op,UpperLowerFlag UPLO,int N>
struct Foreach<T,Op,sp_mat<UPLO,N> > {
  Foreach(const sp_mat<UPLO,N>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,int N>
struct ForeachT<T,Op,sp_mat<UPLO,N> > {
  ForeachT(const sp_mat<UPLO,N>& _ma,T* _a,Op _op) {
    Foreach<T,Op,sp_mat<UPLO,N> >(_ma,_a,_op);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,
          int NB>
struct Foreach2<T,Op,sp_mat<UPLO,NA>,sp_mat<UPLO,NB> > {
  Foreach2(const sp_mat<UPLO,NA>& _ma,T* _a,
           const sp_mat<UPLO,NB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i],_b[i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO,
          int NA,
          int NB,
          int NC>
struct Foreach3<T,Op,
                sp_mat<UPLO,NA>,
                sp_mat<UPLO,NB>,
                sp_mat<UPLO,NC> > {
  Foreach3(const sp_mat<UPLO,NA>& _ma,T* _a,
           const sp_mat<UPLO,NB>& _mb,T* _b,
           const sp_mat<UPLO,NC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i],_b[i],_c[i]);
  }
};

//-----------------------------------------------------------------------------

template <typename T,typename Op,
          UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,int N,int LD>
struct Foreach<T,Op,tr_mat<UPLO,TRANS,DIAG,N,LD> > {
  Foreach(const tr_mat<UPLO,TRANS,DIAG,N,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
           if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i]);
    }
  }
};

  // MISSING: ForeachT

template <typename T,typename Op,
          UpperLowerFlag UPLO,TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int LDA,
          int NB,int LDB>
struct Foreach2<T,Op,tr_mat<UPLO,TRANS,DIAG,NA,LDA>,tr_mat<UPLO,TRANS,DIAG,NB,LDB> > {
  Foreach2(const tr_mat<UPLO,TRANS,DIAG,NA,LDA>& _ma,T* _a,
           const tr_mat<UPLO,TRANS,DIAG,NB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
          if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
    }
  }
};

template <typename T,typename Op,
          UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int LDA,
          int NB,int LDB,
          int NC,int LDC>
struct Foreach3<T,Op,
                tr_mat<UPLO,TRANS,DIAG,NA,LDA>,
                tr_mat<UPLO,TRANS,DIAG,NB,LDB>,
                tr_mat<UPLO,TRANS,DIAG,NC,LDC> > {
  Foreach3(const tr_mat<UPLO,TRANS,DIAG,NA,LDA>& _ma,T* _a,
           const tr_mat<UPLO,TRANS,DIAG,NB,LDB>& _mb,T* _b,
           const tr_mat<UPLO,TRANS,DIAG,NC,LDC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());

    if (UPLO==Lower) {
      int m0=0;
      for (int j=0;j<_ma.n();++j,++m0)
        for (int i=m0;i<_ma.n();++i)
          if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
    }
    else {
      int m=1;
      for (int j=0;j<_ma.n();++j,++m)
        for (int i=0;i<m;++i)
          if (DIAG==NoU)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
          else if (i!=j)
            _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
    }
  }
};


//-----------------------------------------------------------------------------

template <typename T,typename Op,
          UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,int N,int K,int LD>
struct Foreach<T,Op,tb_mat<UPLO,TRANS,DIAG,N,K,LD> > {
  Foreach(const tb_mat<UPLO,TRANS,DIAG,N,K,LD>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    int i0=(DIAG==UnitDiag && UPLO==Lower) ? 1 : 0;
    int m=_ma.k()+((DIAG==UnitDiag && UPLO==Lower) ? 0 : 1);

    for (int j=0;j<_ma.n();++j)
      for (int i=i0;i<m;++i) // ignore the extra work
          _op(_a[j*_ma.ld()+i]);
  }
};

  // MISSING: ForeachT

template <typename T,typename Op,UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int KA,int LDA,
          int NB,int KB,int LDB>
struct Foreach2<T,Op,tb_mat<UPLO,TRANS,DIAG,NA,KA,LDA>,tb_mat<UPLO,TRANS,DIAG,NB,KB,LDB> > {
  Foreach2(const tb_mat<UPLO,TRANS,DIAG,NA,KA,LDA>& _ma,T* _a,
           const tb_mat<UPLO,TRANS,DIAG,NB,KB,LDB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.k()==_mb.k());

    int i0=(DIAG==UnitDiag && UPLO==Lower) ? 1 : 0;
    int m=_ma.k()+((DIAG==UnitDiag && UPLO==Lower) ? 0 : 1);

    for (int j=0;j<_ma.n();++j)
      for (int i=i0;i<m;++i) // ignore the extra work
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,int KA,int LDA,
          int NB,int KB,int LDB,
          int NC,int KC,int LDC>
struct Foreach3<T,Op,
                tb_mat<UPLO,TRANS,DIAG,NA,KA,LDA>,
                tb_mat<UPLO,TRANS,DIAG,NB,KB,LDB>,
                tb_mat<UPLO,TRANS,DIAG,NC,KC,LDC> > {
  Foreach3(const tb_mat<UPLO,TRANS,DIAG,NA,KA,LDA>& _ma,T* _a,
           const tb_mat<UPLO,TRANS,DIAG,NB,KB,LDB>& _mb,T* _b,
           const tb_mat<UPLO,TRANS,DIAG,NC,KC,LDC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());
    assert(_ma.k()==_mb.k());
    assert(_ma.k()==_mc.k());

    int i0=(DIAG==UnitDiag && UPLO==Lower) ? 1 : 0;
    int m=_ma.k()+((DIAG==UnitDiag && UPLO==Lower) ? 0 : 1);

    for (int j=0;j<_ma.n();++j)
      for (int i=i0;i<m;++i) // ignore the extra work
        _op(_a[j*_ma.ld()+i],_b[j*_mb.ld()+i],_c[j*_mb.ld()+i]);
  }
};


//-----------------------------------------------------------------------------

template <typename T,typename Op,UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,int N>
struct Foreach<T,Op,tp_mat<UPLO,TRANS,DIAG,N> > {
  Foreach(const tp_mat<UPLO,TRANS,DIAG,N>& _ma,T* _a,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i]);
  }
};

  // MISSING: ForeachT

template <typename T,typename Op,UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,
          int NB>
struct Foreach2<T,Op,tp_mat<UPLO,TRANS,DIAG,NA>,tp_mat<UPLO,TRANS,DIAG,NB> > {
  Foreach2(const tp_mat<UPLO,TRANS,DIAG,NA>& _ma,T* _a,
           const tp_mat<UPLO,TRANS,DIAG,NB>& _mb,T* _b,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i],_b[i]);
  }
};

template <typename T,typename Op,UpperLowerFlag UPLO, TransposeFlag TRANS,DiagonalFlag DIAG,
          int NA,
          int NB,
          int NC>
struct Foreach3<T,Op,
                tp_mat<UPLO,TRANS,DIAG,NA>,
                tp_mat<UPLO,TRANS,DIAG,NB>,
                tp_mat<UPLO,TRANS,DIAG,NC> > {
  Foreach3(const tp_mat<UPLO,TRANS,DIAG,NA>& _ma,T* _a,
           const tp_mat<UPLO,TRANS,DIAG,NB>& _mb,T* _b,
           const tp_mat<UPLO,TRANS,DIAG,NC>& _mc,T* _c,Op _op) {
    _VC_DBG_BLAS_SCOPE();

    assert(_ma.n()==_mb.n());
    assert(_ma.n()==_mc.n());

    int n=_ma.n()*(_ma.n()+1)/2;
    for (int i=0;i<n;++i)
      _op(_a[i],_b[i],_c[i]);
  }
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------

# endif // DOXYGEN_SKIP

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
