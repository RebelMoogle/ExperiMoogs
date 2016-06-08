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

    \internal

    \arg
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

#ifndef VC_MATH_BLAS_MATRIX_EXPR_HH
# error "don't include directly"
#endif

/* NOTE
 * ====
 *
 * Rules are defined for (combinations of) vector_reference_t and
 * vector_const_reference_t or matrix_reference_t and
 * matrix_const_reference_t, respectively.
 *
 * Rules could be formulated more elegantly relying on C++11 with
 * return types inferred from arguments. (This would save some
 * explicit enumeration of cases, especially for "a-b == a+(-b)".
 * (The idea is to have a single rule "expr-ANY => expr+(-ANY)".)
 *
 * Currently, many rules ("a-b") are implemented in a non-elegant but
 * more efficient (disregarding optimization) "copy&paste" fashion in
 * favor of relying on "a+(-b)" other rules.
 */

//
// sums of expressions
//

enum Expr_sum { EXPR_sum=100 };

template <typename T,typename A,typename B>
struct Expression<Type_vector,Expr_sum,T,A,B> {
 Expression(const A& _a,const B& _b)
 : a(_a), b(_b) { }

 template <typename LHS> void assign(LHS& _lhs) const {
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// a=" << a << std::endl;
   std::cerr << "// a.assign(_lhs);\\n";
# endif
   a.assign(_lhs);
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// b=" << b << std::endl;
   std::cerr << "// b.plus_assign(_lhs);\\n";
# endif

   b.plus_assign(_lhs);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// a=" << a << std::endl;
   std::cerr << "// _a.plus_assign(_lhs);\\n";
# endif
   a.plus_assign(_lhs);
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// b=" << b << std::endl;
   std::cerr << "// _b.plus_assign(_lhs);\\n";
# endif
   b.plus_assign(_lhs);
 }

 A a;
 B b;
};

template <typename T,typename EA,typename EB>
struct Expression<Type_matrix,Expr_sum,T,EA,EB> {
 Expression(const EA& _A,const EB& _B)
 : A(_A), B(_B) { }

 template <typename LHS> void assign(LHS& _lhs) const {
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// A=" << A << std::endl;
   std::cerr << "// A.assign(_lhs);\\n";
# endif
   A.assign(_lhs);
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// B=" << B << std::endl;
   std::cerr << "// B.plus_assign(_lhs);\\n";
# endif
   B.plus_assign(_lhs);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// A=" << A << std::endl;
   std::cerr << "// _A.plus_assign(_lhs);\\n";
# endif
   A.plus_assign(_lhs);
# ifdef _VC_BLAS_TRACE_EXPR
   std::cerr << "// B=" << B << std::endl;
   std::cerr << "// _B.plus_assign(_lhs);\\n";
# endif
   B.plus_assign(_lhs);
 }

 EA A;
 EB B;
};

//
// unary +sum
//

template <typename TYPE,typename T,typename A1,typename A2>
Expression<TYPE,Expr_sum,T,A1,A2> operator+(const Expression<TYPE,Expr_sum,T,A1,A2>& _sum) {
  return _sum; // Note: -sum yields an error
}

// NOTE: there is no unary "-sum", -(A+B) yields an error

//
// --- vector expressions ---
//

//
// vexpr+vexpr, vexpr-vexpr
//

// vexpr+vexpr
template <typename T,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,EXPRB,T,B1,B2> >
operator+(const Expression<Type_vector,EXPRA,T,A1,A2>& _a,
          const Expression<Type_vector,EXPRB,T,B1,B2>& _b) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,EXPRB,T,B1,B2> >(_a,_b);
}

// vexpr-vexpr
template <typename T,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,EXPRB,T,B1,B2> >
operator-(const Expression<Type_vector,EXPRA,T,A1,A2>& _a,
          const Expression<Type_vector,EXPRB,T,B1,B2>& _b) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,EXPRB,T,B1,B2> >(_a,-_b);
}

//
// vexpr+init, vexpr-init
//

// vexpr+init
template <typename T,typename S,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_initializer,EXPRB,S,B1,B2> >
operator+(const Expression<Type_vector,EXPRA,T,A1,A2>& _a,
          const Expression<Type_initializer,EXPRB,S,B1,B2>& _b) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_initializer,EXPRB,S,B1,B2> >(_a,_b);
}

// vexpr-init
template <typename T,typename S,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_initializer,EXPRB,S,B1,B2> >
operator-(const Expression<Type_vector,EXPRA,T,A1,A2>& _a,
          const Expression<Type_initializer,EXPRB,S,B1,B2>& _b) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_initializer,EXPRB,S,B1,B2> >(_a,-_b);
}


//
// vexpr+x, x+vexpr
//

// vexpr+x
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator+(const Expression<Type_vector,EXPRA,T,A1,A2>& _expr,
          const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(1),_y));
}
// vexpr+x
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator+(const Expression<Type_vector,EXPRA,T,A1,A2>& _expr,
          const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()));
}
// x+vexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_vector,EXPRA,T,A1,A2> >
operator+(const vector_const_reference_t<T,VY>& _y,
          const Expression<Type_vector,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_vector,EXPRA,T,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y),_expr);
}
// x+vexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_vector,EXPRA,T,A1,A2> >
operator+(const vector_reference_t<T,VY>& _y,
          const Expression<Type_vector,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_vector,EXPRA,T,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()),_expr);
}

//
// vexpr-x, x-vexpr
//

// vexpr-x
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator-(const Expression<Type_vector,EXPRA,T,A1,A2>& _expr,
          const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y));
}
// vexpr-x
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,EXPRA,T,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator-(const Expression<Type_vector,EXPRA,T,A1,A2>& _expr,
          const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,EXPRA,T,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y.const_ref()));
}
// x-vexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_vector,EXPRA,T,A1,A2> >
operator-(const vector_const_reference_t<T,VY>& _y,
          const Expression<Type_vector,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_vector,EXPRA,T,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y),-_expr);
}
// x-vexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_vector,EXPRA,T,A1,A2> >
operator-(const vector_reference_t<T,VY>& _y,
          const Expression<Type_vector,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_vector,EXPRA,T,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()),-_expr);
}


//
// x+y
//


// x+y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator+(const vector_const_reference_t<T,VX>& _x,const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x)+Expression<Type_vector,Expr_ax,T,VY>(T(1),_y);
}
// x+y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator+(const vector_reference_t<T,VX>& _x,const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x.const_ref())+Expression<Type_vector,Expr_ax,T,VY>(T(1),_y);
}
// x+y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator+(const vector_const_reference_t<T,VX>& _x,const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x)+Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref());
}
// x+y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator+(const vector_reference_t<T,VX>& _x,const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x.const_ref())+Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref());
}

//
// x-y
//

// x-y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator-(const vector_const_reference_t<T,VX>& _x,const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x)+Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y);
}
// x-y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator-(const vector_reference_t<T,VX>& _x,const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x.const_ref())+Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y);
}
// x-y
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator-(const vector_const_reference_t<T,VX>& _x,const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x)+Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y.const_ref());
}
template <typename T,typename VX,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VX>,Expression<Type_vector,Expr_ax,T,VY> >
operator-(const vector_reference_t<T,VX>& _x,const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_ax,T,VX>(T(1),_x.const_ref())+Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y.const_ref());
}

//
// init+x, x+init
//

// init+x
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator+(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(1),_y));
}
// init+x
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator+(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()));
}
// x+init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator+(const vector_const_reference_t<T,VY>& _y,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y),_expr);
}
// x+init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator+(const vector_reference_t<T,VY>& _y,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()),_expr);
}

//
// init-x, x-init (if "-init" is defined)
//

// init-x
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator-(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const vector_const_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y));
}
// init-x
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_vector,Expr_ax,T,VY> >
operator-(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const vector_reference_t<T,VY>& _y) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_vector,Expr_ax,T,VY> >
   (_expr,Expression<Type_vector,Expr_ax,T,VY>(T(-1),_y.const_ref()));
}
// x-init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator-(const vector_const_reference_t<T,VY>& _y,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y),-_expr);
}
// x-init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename VY>
Expression<Type_vector,Expr_sum,T,
           Expression<Type_vector,Expr_ax,T,VY>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator-(const vector_reference_t<T,VY>& _y,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_vector,Expr_sum,T,
                    Expression<Type_vector,Expr_ax,T,VY>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_vector,Expr_ax,T,VY>(T(1),_y.const_ref()),-_expr);
}


//
// --- matrix expressions ---
//

//
// mexpr+mexpr, mexpr-mexpr
//

// mexpr+mexpr
template <typename T,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,EXPRB,T,B1,B2> >
operator+(const Expression<Type_matrix,EXPRA,T,A1,A2>& _a,
          const Expression<Type_matrix,EXPRB,T,B1,B2>& _b) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,EXPRB,T,B1,B2> >(_a,_b);
}

// mexpr-mexpr
template <typename T,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,EXPRB,T,B1,B2> >
operator-(const Expression<Type_matrix,EXPRA,T,A1,A2>& _a,
          const Expression<Type_matrix,EXPRB,T,B1,B2>& _b) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,EXPRB,T,B1,B2> >(_a,-_b);
}

// mexpr+init
template <typename T,typename S,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_initializer,EXPRB,S,B1,B2> >
operator+(const Expression<Type_matrix,EXPRA,T,A1,A2>& _a,
          const Expression<Type_initializer,EXPRB,S,B1,B2>& _b) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_initializer,EXPRB,S,B1,B2> >(_a,_b);
}

// mexpr-init
template <typename T,typename S,typename EXPRA,typename EXPRB,
          typename A1,typename A2,typename B1,typename B2>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_initializer,EXPRB,S,B1,B2> >
operator-(const Expression<Type_matrix,EXPRA,T,A1,A2>& _a,
          const Expression<Type_initializer,EXPRB,S,B1,B2>& _b) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_initializer,EXPRB,S,B1,B2> >(_a,-_b);
}

//
// mexpr+A, A+mexpr
//

// mexpr+A
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr,
          const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B));
}
// mexpr+A
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr,
          const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref()));
}
// A+mexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_matrix,EXPRA,T,A1,A2> >
operator+(const matrix_const_reference_t<T,MB>& _B,
          const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_matrix,EXPRA,T,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B),_expr);
}
// A+mexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_matrix,EXPRA,T,A1,A2> >
operator+(const matrix_reference_t<T,MB>& _B,
          const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_matrix,EXPRA,T,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref()),_expr);
}

//
// mexpr-A, A-mexpr
//

// mexpr-A
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr,
          const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B));
}
// mexpr-A
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,EXPRA,T,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr,
          const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,EXPRA,T,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B.const_ref()));
}

// A-mexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_matrix,EXPRA,T,A1,A2> >
operator-(const matrix_const_reference_t<T,MB>& _B,
          const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr) {
  return _B+(-_expr);
}
// A-mexpr
template <typename T,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_matrix,EXPRA,T,A1,A2> >
operator-(const matrix_reference_t<T,MB>& _B,
          const Expression<Type_matrix,EXPRA,T,A1,A2>& _expr) {
  return _B.const_ref()+(-_expr);
}

//
// A+B
//

// A+B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const matrix_const_reference_t<T,MA>& _A,const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A)+Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B);
}
// A+B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const matrix_reference_t<T,MA>& _A,const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A.const_ref())+Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B);
}
// A+B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const matrix_const_reference_t<T,MA>& _A,const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A)+Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref());
}
// A+B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const matrix_reference_t<T,MA>& _A,const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A.const_ref())+Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref());
}

//
// A-B
//

// A-B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const matrix_const_reference_t<T,MA>& _A,const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A)+Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B);
}
// A-B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const matrix_reference_t<T,MA>& _A,const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A.const_ref())+Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B);
}
// A-B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const matrix_const_reference_t<T,MA>& _A,const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A)+Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B.const_ref());
}
// A-B
template <typename T,typename MA,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MA>,Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const matrix_reference_t<T,MA>& _A,const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_aA,T,MA>(T(1),_A.const_ref())+Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B.const_ref());
}

//
// init+A, A+init
//

// init+A
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B));
}
// init+A
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator+(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref()));
}
// A+init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator+(const matrix_const_reference_t<T,MB>& _B,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B),_expr);
}
// A+init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator+(const matrix_reference_t<T,MB>& _B,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref()),_expr);
}

//
// init-A, A-init (if "-init" is defined)
//

// init-A
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const matrix_const_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B));
}
// init-A
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_initializer,EXPRA,S,A1,A2>,
           Expression<Type_matrix,Expr_aA,T,MB> >
operator-(const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr,
          const matrix_reference_t<T,MB>& _B) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_initializer,EXPRA,S,A1,A2>,
                    Expression<Type_matrix,Expr_aA,T,MB> >
   (_expr,Expression<Type_matrix,Expr_aA,T,MB>(T(-1),_B.const_ref()));
}
// A-init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator-(const matrix_const_reference_t<T,MB>& _B,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B),-_expr);
}
// A-init
template <typename T,typename S,typename EXPRA,typename A1,typename A2,typename MB>
Expression<Type_matrix,Expr_sum,T,
           Expression<Type_matrix,Expr_aA,T,MB>,
           Expression<Type_initializer,EXPRA,S,A1,A2> >
operator-(const matrix_reference_t<T,MB>& _B,
          const Expression<Type_initializer,EXPRA,S,A1,A2>& _expr) {
  return Expression<Type_matrix,Expr_sum,T,
                    Expression<Type_matrix,Expr_aA,T,MB>,
                    Expression<Type_initializer,EXPRA,S,A1,A2> >
   (Expression<Type_matrix,Expr_aA,T,MB>(T(1),_B.const_ref()),-_expr);
}

//
// --- few more rules that are easier to specify here ---
//


//
// -zeros, +zeros
//

inline
Expression<Type_initializer,Expr_zeros>
operator-(const Expression<Type_initializer,Expr_zeros>& _zeros) { return _zeros; }

inline
Expression<Type_initializer,Expr_zeros>
operator+(const Expression<Type_initializer,Expr_zeros>& _zeros) { return _zeros; }

//
// -ones, +ones
//

inline
Expression<Type_initializer,Expr_all,int>
operator-(const Expression<Type_initializer,Expr_ones>& /*_ones*/) {
  return Expression<Type_initializer,Expr_all,int>(-1);
}
inline
Expression<Type_initializer,Expr_ones>
operator+(const Expression<Type_initializer,Expr_ones>& _ones) { return _ones; }


//
// --- undefined operations are explicitly forbidden ---
//

struct Expr_ERROR {};

// forbidden: vector*vector
template <typename T,typename EXPRA,typename EXPRB,typename A1,typename A2,typename B1,typename B2>
Expr_ERROR operator*(const Expression<Type_vector,EXPRA,T,A1,A2>&,const Expression<Type_vector,EXPRB,T,B1,B2>&) {
  static_assertion(!std::is_same<T,T>::value,"undefined operation vector*vector");
  return Expr_ERROR();
}
// forbidden: vector*matrix
template <typename T,typename EXPRA,typename EXPRB,typename A1,typename A2,typename B1,typename B2>
Expr_ERROR operator*(const Expression<Type_vector,EXPRA,T,A1,A2>&,const Expression<Type_matrix,EXPRB,T,B1,B2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation vector*matrix");
 return Expr_ERROR();
}
// foirbidden: ANY*sum
template <typename T,typename ANY,typename TYPE,typename A1,typename A2>
Expr_ERROR operator*(const ANY&,const Expression<TYPE,Expr_sum,T,A1,A2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation ?*sum");
 return Expr_ERROR(); // could be solved by decltype(): (a+b)*s -> (s*a)+(s*b) -- BUT inefficient
}
// forbidden sum*ANY
template <typename T,typename ANY,typename TYPE,typename A1,typename A2>
Expr_ERROR operator*(const Expression<TYPE,Expr_sum,T,A1,A2>&,const ANY&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation *");
 return Expr_ERROR(); // could be solved by decltype(): s*(a+b) -> (s*a)+(s*b) -- BUT inefficient
}
// forbidden: ANY/sum
template <typename T,typename ANY,typename TYPE,typename A1,typename A2>
Expr_ERROR operator/(const ANY&,const Expression<TYPE,Expr_sum,T,A1,A2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation ?/sum");
 return Expr_ERROR();
}
// forbidden -sum
template <typename T,typename TYPE,typename A1,typename A2>
Expr_ERROR operator-(const Expression<TYPE,Expr_sum,T,A1,A2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation -sum");
 return Expr_ERROR(); // could be solved by decltype(): -(a+x) -> (-a)+(-b) -- BUT inefficient
}
// forbidden sum/ANY
template <typename T,typename ANY,typename TYPE,typename A1,typename A2>
Expr_ERROR operator/(const Expression<TYPE,Expr_sum,T,A1,A2>&,const ANY&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation sum/?");
 return Expr_ERROR();
}
// forbidden ANY/vec
template <typename T,typename ANY,typename EXPR,typename A1,typename A2>
Expr_ERROR operator/(const ANY&,const Expression<Type_vector,EXPR,T,A1,A2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation ?/vector"); 
 return Expr_ERROR();
}
// forbidden ANY/mat
template <typename T,typename ANY,typename EXPR,typename A1,typename A2>
Expr_ERROR operator/(const ANY&,const Expression<Type_matrix,EXPR,T,A1,A2>&) {
static_assertion(!std::is_same<T,T>::value,"undefined operation ?/matrix");
 return Expr_ERROR();
}

//
// --- debug code used #ifdef _VC_BLAS_TRACE_EXPR ---
//

template <typename TYPE,typename EXPR,typename T,typename A1,typename A2>
std::ostream& operator<<(std::ostream& _out,const Expression<TYPE,EXPR,T,A1,A2>& _e) {
  return _out << "Expression<*>";
}
template <typename TYPE,typename T,typename A1,typename A2>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_vector,Expr_sum,T,A1,A2>& _e) {
  return _out << "Expression<Type_vector,Expr_sum,...>["
              << _e.a << "," << _e.b << "]";
}
template <typename TYPE,typename T,typename A1,typename A2>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,Expr_sum,T,A1,A2>& _e) {
  return _out << "Expression<Type_matrix,Expr_sum,...>["
              << _e.EA << "," << _e.EB << "]";
}
template <typename EXPR,typename T,typename A1,typename A2>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,EXPR,T,A1,A2>& _e) {
  return _out << "Expression<Type_matrix,...>";
}
template <typename EXPR,typename T,typename A1,typename A2>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_vector,EXPR,T,A1,A2>& _e) {
  return _out << "Expression<Type_vector,...>";
}
template <typename T,typename VX>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_vector,Expr_ax,T,VX>& _e) {
  return _out << "ax:(" << _e.a << "," << _e.x << ")";
}
template <typename T,typename MA,typename VX>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_vector,Expr_aAx,T,MA,VX>& _e) {
  return _out << "aAx:(" << _e.a << "," << _e.A << "," << _e.x << ")";
}
template <typename T,typename MA>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,Expr_aA,T,MA>& _e) {
  return _out << "aAx:(" << _e.a << "," << _e.A << ")";
}
template <typename T,typename MA,typename MB>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _e) {
  return _out << "aAB:(" << _e.a << "," << _e.A << "," << _e.B << ")";
}
template <typename T,typename VX>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,Expr_axxt,T,VX>& _e) {
  return _out << "axxt:(" << _e.a << "," << _e.x << ")";
}
template <typename T,typename MA>
std::ostream& operator<<(std::ostream& _out,const Expression<Type_matrix,Expr_aAtA,T,MA>& _e) {
  return _out << "aAtA:(" << _e.a << "," << _e.A << ")";
}
