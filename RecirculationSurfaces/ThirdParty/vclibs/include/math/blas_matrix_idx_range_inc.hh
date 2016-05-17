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

    \arg Provide ranges of indices with support for "end".

    Support `x(3:5)` as `x(_(2,4))`
    `A(2:3,4:5)` as `A(_(1,2),_(3,4))`,
    `A(:,(k+1):end-1)` as `A($,_(k,end-1))`,
    `A((3:4)+1,:)` as `A(_(2,3)+1,:)`,...

    This file is independent of most other stuff.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

# ifndef DOXYGEN_SKIP

template <bool END1,bool END2> struct idx_range {
private:
  // undefined ("A(end-2,1)")
};

template <> struct idx_range<false,false> {
  idx_range(int _i1=0,int _i2=0) : i1(_i1), i2(_i2) {}   // default: empty
  int i1,i2;
  int eval_i1(int /*_n*/) const { return i1; }
  int eval_i2(int /*_n*/) const { return i2; }
  idx_range<false,false> operator+(int _i) const {
    return idx_range<false,false>(i1+_i,i2+_i);
  }
  idx_range<false,false> operator-(int _i) const {
    return idx_range<false,false>(i1-_i,i2-_i);
  }
};
template <> struct idx_range<false,true> {
  idx_range(int _i1=0,int _i2=0) : i1(_i1), i2(_i2) {}  // default: all
  int i1,i2;
  int eval_i1(int /*_n*/) const { return i1; }
  int eval_i2(int _n)     const { return _n-i2-1; }
};
template <> struct idx_range<true,true> {
  idx_range(int _i1=0, int _i2=0) : i1(_i1), i2(_i2) {}  // default: empty
  int i1,i2;
  int eval_i1(int _n) const { return _n-i1-1; }
  int eval_i2(int _n) const { return _n-i2-1; }
};
// idx_range<true,false> ~ "A(end-1,2)" is undefined

struct idx_end { // global "end"
  idx_end() : offset(0) {}
  int offset;
};

inline idx_end operator-(const idx_end& _end,int _i) {
  idx_end e; e.offset=_end.offset+_i; return e;
}

struct make_idx_range { // global _
  idx_range<false,true> operator()(void) const { // all
    return idx_range<false,true>();
  }
  idx_range<false,false> operator()(int _i1,int _i2) const {
    return idx_range<false,false>(_i1,_i2);
  }
  idx_range<false,false> operator()(int _i1) const {
    return idx_range<false,false>(_i1,_i1);
  }
  idx_range<false,true> operator()(int _i1,const idx_end& _ie2) const {
    return idx_range<false,true>(_i1,_ie2.offset);
  }
  idx_range<true,true> operator()(const idx_end& _ie1,const idx_end& _ie2) const {
    return idx_range<true,true>(_ie1.offset,_ie2.offset);
  }
  idx_range<false,true> operator()(const idx_end& _ie1) const {
    return idx_range<false,true>(_ie1.offset,_ie1.offset);
  }

  make_idx_range() : dummy(0) {}
  int dummy;
};

template <bool END1,bool END2> // debug
inline std::ostream&
operator<<(std::ostream& _out,const idx_range<END1,END2>& r) {
  return _out << "idx_range<"
              << (END1?'T':'F') << ',' << (END2?'T':'F')
              << ">(" << r.i1 << ',' << r.i2 << ")";
}

# endif // DOXYGEN_SKIP

/** \defgroup vc_blas_index Matlab-like indexing of matrices and vectors.
    \ingroup vc_blas

    Indices are created from special symbols `"_"`, `"end"` (synonym
    `"END"`), `"$"`, e.g., `A(_(2,end-2),$)` refers to `A(3:end-2),:)`
    in Matlab notation. (We use 0-based indices.)
 */

/// create index range (usage: `_(i,j)`) \ingroup vc_blas_index
extern const make_idx_range _;
/// index range ":" ("all") \ingroup vc_blas


# ifdef __clang__
//#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wdollar-in-identifier-extension"
# endif

/** \var $
    create index range "all" (`x(:)` in Matlab) 
    \ingroup vc_blas_index
*/
extern const idx_range<false,true> $;

//# ifdef __clang__
//#  pragma GCC diagnostic pop
//# endif

/// "end" specifier for indexing (usage `end-int`) \ingroup vc_blas_index
extern const idx_end end;
/** alias for VC::math::blas::end to avoid need for `using` (disambiguation) 
    \ingroup vc_blas_index
*/
extern const idx_end& END;

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
