
//
// THIS CODE WAS GENERATED AUTOMATICALLY. DO NOT EDIT !!!
//

# ifdef DOXYGEN_SKIP

/** Expression templates.
    ingroup vc_blas
    iternal
    Expression<TYPE,EXPR,T,ARG1,ARG2> provides an interface to
    evaluating terms which include matrix and vector expressions.
    This enables the use of certain arithmetic expressions and
    their mapping to BLAS calls.
    (The details remain undocumented.)
     a matrix_reference_t, vector_reference_t
*/
typename <typename TYPE,typename EXPR,typename T,typename ARG1,typename ARG2>
Expression {};

# endif

# ifndef DOXYGEN_SKIP

//
// expression types
//

enum Type_error { TYPE_error=0 };
enum Type_vector { TYPE_vector=1 };
enum Type_matrix { TYPE_matrix=2 };
enum Type_initializer { TYPE_initializer=3 };

//
// the expression (will be specialized)
//

template <typename TYPE,typename EXPR,
          typename ARG1=void,
          typename ARG2=void,
          typename ARG3=void>
struct Expression {
 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n static_assert(!std::is_same<LHS,LHS>::value,\"undefined evaluation method assign\")";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   static_assert(!std::is_same<LHS,LHS>::value,"undefined evaluation method assign");
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n static_assert(!std::is_same<LHS,LHS>::value,\"undefined evaluation method plus_assign\")";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   static_assert(!std::is_same<LHS,LHS>::value,"undefined evaluation method plus_assign");
 }
private:
  Expression();
};

//
// expressions (specializations)
//
                                             
enum Expr_ax { EXPR_ax=19 };

template <typename T,typename VX>
struct Expression<Type_vector,Expr_ax,T,VX> {
 typedef T value_t;
 Expression(const T& _a,
            const vector_const_reference_t<T,VX>& _x)
 : a(_a), x(_x) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n if (a==T(1)) _lhs.copy(x); else _lhs.copy_scal(a,x);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   if (a==T(1)) _lhs.copy(x); else _lhs.copy_scal(a,x);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.axpy(a,x);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.axpy(a,x);
 }

 T a;
 vector_const_reference_t<T,VX> x;
};
              
enum Expr_aAx { EXPR_aAx=20 };

template <typename T,typename MA,typename VX>
struct Expression<Type_vector,Expr_aAx,T,MA,VX> {
 typedef T value_t;
 Expression(const T& _a,
            const matrix_const_reference_t<T,MA>& _A,
            const vector_const_reference_t<T,VX>& _x)
 : a(_a), A(_A), x(_x) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.mv(a,A,x,T(0));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.mv(a,A,x,T(0));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.mv(a,A,x,T(1));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.mv(a,A,x,T(1));
 }

 T a;
 matrix_const_reference_t<T,MA> A;
 vector_const_reference_t<T,VX> x;
};
              
enum Expr_aA { EXPR_aA=21 };

template <typename T,typename MA>
struct Expression<Type_matrix,Expr_aA,T,MA> {
 typedef T value_t;
 Expression(const T& _a,
            const matrix_const_reference_t<T,MA>& _A)
 : a(_a), A(_A) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n cp(A,_lhs); if (a!=T(1)) _lhs.mscal(a);";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   cp(A,_lhs); if (a!=T(1)) _lhs.mscal(a);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n assert(_lhs.data()!=A.data()); _lhs.madd(a,A);";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   assert(_lhs.data()!=A.data()); _lhs.madd(a,A);
 }

 T a;
 matrix_const_reference_t<T,MA> A;
};
              
enum Expr_aAB { EXPR_aAB=22 };

template <typename T,typename MA,typename MB>
struct Expression<Type_matrix,Expr_aAB,T,MA,MB> {
 typedef T value_t;
 Expression(const T& _a,
            const matrix_const_reference_t<T,MA>& _A,
            const matrix_const_reference_t<T,MB>& _B)
 : a(_a), A(_A), B(_B) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.mm(a,A,B,T(0));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A << ',' << "B=" << B;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.mm(a,A,B,T(0));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.mm(a,A,B,T(1));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A << ',' << "B=" << B;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.mm(a,A,B,T(1));
 }

 T a;
 matrix_const_reference_t<T,MA> A;
 matrix_const_reference_t<T,MB> B;
};
              
enum Expr_aAtA { EXPR_aAtA=23 };

template <typename T,typename MA>
struct Expression<Type_matrix,Expr_aAtA,T,MA> {
 typedef T value_t;
 Expression(const T& _a,
            const matrix_const_reference_t<T,MA>& _A)
 : a(_a), A(_A) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.rk(a,trans(A),T(0));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.rk(a,trans(A),T(0));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.rk(a,trans(A),T(1));";
    std::cerr << " // "<< "a=" << a << ',' << "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.rk(a,trans(A),T(1));
 }

 T a;
 matrix_const_reference_t<T,MA> A;
};
              
enum Expr_axxt { EXPR_axxt=24 };

template <typename T,typename VX>
struct Expression<Type_matrix,Expr_axxt,T,VX> {
 typedef T value_t;
 Expression(const T& _a,
            const vector_const_reference_t<T,VX>& _x)
 : a(_a), x(_x) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_zero(); _lhs.r(a,x);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_zero(); _lhs.r(a,x);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.r(a,x);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.r(a,x);
 }

 T a;
 vector_const_reference_t<T,VX> x;
};
              
enum Expr_aABt { EXPR_aABt=25 };

template <typename T,typename VX,typename VY>
struct Expression<Type_matrix,Expr_aABt,T,VX,VY> {
 typedef T value_t;
 Expression(const T& _a,
            const vector_const_reference_t<T,VX>& _x,
            const vector_const_reference_t<T,VY>& _y)
 : a(_a), x(_x), y(_y) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_zero(); _lhs.r2(a,x,y);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x << ',' << "y=" << y;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_zero(); _lhs.r2(a,x,y);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.r2(a,x,y);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x << ',' << "y=" << y;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.r2(a,x,y);
 }

 T a;
 vector_const_reference_t<T,VX> x;
 vector_const_reference_t<T,VY> y;
};
              
enum Expr_zeros { EXPR_zeros=26 };

template <>
struct Expression<Type_initializer,Expr_zeros> {
 Expression()
  { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_zero();";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_zero();
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n ";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   ;
 }

};
              
enum Expr_ones { EXPR_ones=27 };

template <>
struct Expression<Type_initializer,Expr_ones> {
 Expression()
  { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_all(typename LHS::value_t(1));";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_all(typename LHS::value_t(1));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.adds(typename LHS::value_t(1))";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.adds(typename LHS::value_t(1));
 }

};
              
enum Expr_random_values { EXPR_random_values=28 };

template <>
struct Expression<Type_initializer,Expr_random_values> {
 Expression()
  { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_rand();";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_rand();
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n static_assert(!std::is_same<LHS,LHS>::value,\"undefined evaluation method plus_assign\")";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   static_assert(!std::is_same<LHS,LHS>::value,"undefined evaluation method plus_assign");
 }

};
              
enum Expr_all { EXPR_all=29 };

template <typename T>
struct Expression<Type_initializer,Expr_all,T> {
 typedef T value_t;
 Expression(const T& _a)
 : a(_a) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_all(typename LHS::value_t(a));";
    std::cerr << " // "<< "a=" << a;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_all(typename LHS::value_t(a));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.adds(typename LHS::value_t(a))";
    std::cerr << " // "<< "a=" << a;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.adds(typename LHS::value_t(a));
 }

 T a;
};
              
enum Expr_eye { EXPR_eye=30 };

template <>
struct Expression<Type_initializer,Expr_eye> {
 Expression()
  { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_eye();";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_eye();
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n diag(_lhs).adds(typename LHS::value_t(1))";
    std::cerr << " // ";
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   diag(_lhs).adds(typename LHS::value_t(1));
 }

};
              
enum Expr_seye { EXPR_seye=31 };

template <>
struct Expression<Type_initializer,Expr_seye> {
 Expression(const double& _d)
 : d(_d) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_eye(typename LHS::value_t(d));";
    std::cerr << " // "<< "d=" << d;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_eye(typename LHS::value_t(d));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n diag(_lhs).adds(typename LHS::value_t(d))";
    std::cerr << " // "<< "d=" << d;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   diag(_lhs).adds(typename LHS::value_t(d));
 }

 double d;
};
              
enum Expr_unit { EXPR_unit=32 };

template <>
struct Expression<Type_initializer,Expr_unit> {
 Expression(const int& _i,
            const double& _d)
 : i(_i), d(_d) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_unit(i,typename LHS::value_t(d))";
    std::cerr << " // "<< "i=" << i << ',' << "d=" << d;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_unit(i,typename LHS::value_t(d));
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs(i)+=typename LHS::value_t(d)";
    std::cerr << " // "<< "i=" << i << ',' << "d=" << d;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs(i)+=typename LHS::value_t(d);
 }

 int i;
 double d;
};
              
enum Expr_cp_column { EXPR_cp_column=33 };

template <typename T,typename MA>
struct Expression<Type_vector,Expr_cp_column,T,MA> {
 typedef T value_t;
 Expression(const matrix_const_reference_t<T,MA>& _A,
            const int& _i)
 : A(_A), i(_i) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n cp_column(A,_lhs,i);";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   cp_column(A,_lhs,i);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.axpy(T(1),column(A,i)); /* GE only */";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.axpy(T(1),column(A,i)); /* GE only */;
 }

 matrix_const_reference_t<T,MA> A;
 int i;
};
              
enum Expr_cp_row { EXPR_cp_row=34 };

template <typename T,typename MA>
struct Expression<Type_vector,Expr_cp_row,T,MA> {
 typedef T value_t;
 Expression(const matrix_const_reference_t<T,MA>& _A,
            const int& _i)
 : A(_A), i(_i) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n cp_row(A,_lhs,i);";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   cp_row(A,_lhs,i);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.axpy(T(1),row(A,i)); /* GE only */";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.axpy(T(1),row(A,i)); /* GE only */;
 }

 matrix_const_reference_t<T,MA> A;
 int i;
};
              
enum Expr_cp_diag { EXPR_cp_diag=35 };

template <typename T,typename MA>
struct Expression<Type_vector,Expr_cp_diag,T,MA> {
 typedef T value_t;
 Expression(const matrix_const_reference_t<T,MA>& _A)
 : A(_A) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n cp_diag(A,_lhs);";
    std::cerr << " // "<< "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   cp_diag(A,_lhs);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.axpy(T(1),diag(A)); /* GE only */";
    std::cerr << " // "<< "A=" << A;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.axpy(T(1),diag(A)); /* GE only */;
 }

 matrix_const_reference_t<T,MA> A;
};
              
enum Expr_cp_block { EXPR_cp_block=36 };

template <typename T,typename MA>
struct Expression<Type_matrix,Expr_cp_block,T,MA> {
 typedef T value_t;
 Expression(const matrix_const_reference_t<T,MA>& _A,
            const int& _i,
            const int& _j,
            const int& _i1,
            const int& _j1)
 : A(_A), i(_i), j(_j), i1(_i1), j1(_j1) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n cp_block(A,i,i1,j,j1,_lhs);";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i << ',' << "j=" << j << ',' << "i1=" << i1 << ',' << "j1=" << j1;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   cp_block(A,i,i1,j,j1,_lhs);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n assert(_lhs.data()!=A.data()); _lhs.madd(T(1),block(A,i,i1,j,j1)); /* GE only */";
    std::cerr << " // "<< "A=" << A << ',' << "i=" << i << ',' << "j=" << j << ',' << "i1=" << i1 << ',' << "j1=" << j1;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   assert(_lhs.data()!=A.data()); _lhs.madd(T(1),block(A,i,i1,j,j1)); /* GE only */;
 }

 matrix_const_reference_t<T,MA> A;
 int i;
 int j;
 int i1;
 int j1;
};
              
enum Expr_diag_mat { EXPR_diag_mat=37 };

template <typename T,typename VX>
struct Expression<Type_matrix,Expr_diag_mat,T,VX> {
 typedef T value_t;
 Expression(const T& _a,
            const vector_const_reference_t<T,VX>& _x)
 : a(_a), x(_x) { }

 template <typename LHS> void assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n _lhs.ld_diag(x,a);";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   _lhs.ld_diag(x,a);
 }
 template <typename LHS> void plus_assign(LHS& _lhs) const {
   
# ifdef _VC_BLAS_TRACE_EXPR
    std::cerr << ">> _VC_BLAS_TRACE_EXPR:\n diag(_lhs).axpy(a,x); /* GE only */";
    std::cerr << " // "<< "a=" << a << ',' << "x=" << x;
    std::cerr << "\n<< _VC_BLAS_TRACE_EXPR\n";
# endif // _VC_BLAS_TRACE_EXPR

   diag(_lhs).axpy(a,x); /* GE only */;
 }

 T a;
 vector_const_reference_t<T,VX> x;
};
       
//
// rules
//

/// operator+:x=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator+(const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_ax,T,VX>(1,_x);
}

/// operator+:xx=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator+(const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_ax,T,VX>(1,_xx.const_ref());
}

/// operator-:x=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator-(const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_ax,T,VX>(-1,_x);
}

/// operator-:xx=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator-(const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_ax,T,VX>(-1,_xx.const_ref());
}

/// operator+:ax=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator+(const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a,_ax.x);
}

/// operator-:ax=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator-(const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_ax,T,VX>(-_ax.a,_ax.x);
}

/// operator+:A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator+(const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(1,_A);
}

/// operator+:AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator+(const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(1,_AA);
}

/// operator-:A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator-(const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(-1,_A);
}

/// operator-:AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator-(const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(-1,_AA);
}

/// operator+:aAx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator+(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a,_aAx.A,_aAx.x);
}

/// operator-:aAx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator-(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(-_aAx.a,_aAx.A,_aAx.x);
}

/// operator+:aA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator+(const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a,_aA.A);
}

/// operator-:aA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator-(const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(-_aA.a,_aA.A);
}

/// operator+:aAB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator+(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a,_aAB.A,_aAB.B);
}

/// operator-:aAB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator-(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(-_aAB.a,_aAB.A,_aAB.B);
}

/// operator+:axxt=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator+(const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a,_axxt.x);
}

/// operator-:axxt=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator-(const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(-_axxt.a,_axxt.x);
}

/// operator+:aAtA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator+(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a,_aAtA.A);
}

/// operator-:aAtA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator-(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(-_aAtA.a,_aAtA.A);
}

/// operator*:f,x=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const float& _f,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_f),_x);
}

/// operator*:x,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_const_reference_t<T,VX>& _x,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_f),_x);
}

/// operator*:f,xx=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const float& _f,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_f),_xx.const_ref());
}

/// operator*:xx,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_reference_t<T,VX>& _xx,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_f),_xx.const_ref());
}

/// operator/:x,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_const_reference_t<T,VX>& _x,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_f),_x);
}

/// operator/:xx,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_reference_t<T,VX>& _xx,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_f),_xx.const_ref());
}

/// operator*:f,ax=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const float& _f,const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_f),_ax.x);
}

/// operator*:ax,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_f),_ax.x);
}

/// operator/:ax,f=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const float& _f) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a/T(_f),_ax.x);
}

/// operator*:f,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const float& _f,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_f),_A);
}

/// operator*:A,f=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_const_reference_t<T,MA>& _A,const float& _f) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_f),_A);
}

/// operator*:f,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const float& _f,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_f),_AA.const_ref());
}

/// operator*:AA,f=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_reference_t<T,MA>& _AA,const float& _f) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_f),_AA.const_ref());
}

/// operator/:f,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const float& _f,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_f),_A);
}

/// operator/:f,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const float& _f,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_f),_AA.const_ref());
}

/// operator/:AA,f=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const matrix_reference_t<T,MA>& _AA,const float& _f) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_f),_AA.const_ref());
}

/// operator*:f,aAB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const float& _f,const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_f),_aAB.A,_aAB.B);
}

/// operator*:aAB,f=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const float& _f) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_f),_aAB.A,_aAB.B);
}

/// operator/:aAB,f=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator/(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const float& _f) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a/T(_f),_aAB.A,_aAB.B);
}

/// operator*:f,aAx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const float& _f,const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_f),_aAx.A,_aAx.x);
}

/// operator*:aAx,f=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const float& _f) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_f),_aAx.A,_aAx.x);
}

/// operator/:aAx,f=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator/(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const float& _f) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a/T(_f),_aAx.A,_aAx.x);
}

/// operator*:f,aA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const float& _f,const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_f),_aA.A);
}

/// operator*:aA,f=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const float& _f) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_f),_aA.A);
}

/// operator/:aA,f=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const float& _f) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a/T(_f),_aA.A);
}

/// operator*:f,aAtA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const float& _f,const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_f),_aAtA.A);
}

/// operator*:aAtA,f=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const float& _f) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_f),_aAtA.A);
}

/// operator/:aAtA,f=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator/(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const float& _f) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a/T(_f),_aAtA.A);
}

/// operator*:f,axxt=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const float& _f,const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_f),_axxt.x);
}

/// operator*:axxt,f=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt,const float& _f) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_f),_axxt.x);
}

/// operator*:d,x=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const double& _d,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_d),_x);
}

/// operator*:x,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_const_reference_t<T,VX>& _x,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_d),_x);
}

/// operator*:d,xx=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const double& _d,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_d),_xx.const_ref());
}

/// operator*:xx,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_reference_t<T,VX>& _xx,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_d),_xx.const_ref());
}

/// operator/:x,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_const_reference_t<T,VX>& _x,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_d),_x);
}

/// operator/:xx,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_reference_t<T,VX>& _xx,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_d),_xx.const_ref());
}

/// operator*:d,ax=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const double& _d,const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_d),_ax.x);
}

/// operator*:ax,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_d),_ax.x);
}

/// operator/:ax,d=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const double& _d) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a/T(_d),_ax.x);
}

/// operator*:d,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const double& _d,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_d),_A);
}

/// operator*:A,d=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_const_reference_t<T,MA>& _A,const double& _d) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_d),_A);
}

/// operator*:d,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const double& _d,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_d),_AA.const_ref());
}

/// operator*:AA,d=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_reference_t<T,MA>& _AA,const double& _d) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_d),_AA.const_ref());
}

/// operator/:d,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const double& _d,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_d),_A);
}

/// operator/:d,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const double& _d,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_d),_AA.const_ref());
}

/// operator/:AA,d=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const matrix_reference_t<T,MA>& _AA,const double& _d) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_d),_AA.const_ref());
}

/// operator*:d,aAB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const double& _d,const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_d),_aAB.A,_aAB.B);
}

/// operator*:aAB,d=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const double& _d) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_d),_aAB.A,_aAB.B);
}

/// operator/:aAB,d=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator/(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const double& _d) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a/T(_d),_aAB.A,_aAB.B);
}

/// operator*:d,aAx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const double& _d,const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_d),_aAx.A,_aAx.x);
}

/// operator*:aAx,d=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const double& _d) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_d),_aAx.A,_aAx.x);
}

/// operator/:aAx,d=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator/(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const double& _d) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a/T(_d),_aAx.A,_aAx.x);
}

/// operator*:d,aA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const double& _d,const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_d),_aA.A);
}

/// operator*:aA,d=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const double& _d) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_d),_aA.A);
}

/// operator/:aA,d=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const double& _d) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a/T(_d),_aA.A);
}

/// operator*:d,aAtA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const double& _d,const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_d),_aAtA.A);
}

/// operator*:aAtA,d=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const double& _d) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_d),_aAtA.A);
}

/// operator/:aAtA,d=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator/(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const double& _d) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a/T(_d),_aAtA.A);
}

/// operator*:d,axxt=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const double& _d,const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_d),_axxt.x);
}

/// operator*:axxt,d=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt,const double& _d) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_d),_axxt.x);
}

/// operator*:i,x=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const int& _i,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_i),_x);
}

/// operator*:x,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_const_reference_t<T,VX>& _x,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_i),_x);
}

/// operator*:i,xx=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const int& _i,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_i),_xx.const_ref());
}

/// operator*:xx,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const vector_reference_t<T,VX>& _xx,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(_i),_xx.const_ref());
}

/// operator/:x,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_const_reference_t<T,VX>& _x,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_i),_x);
}

/// operator/:xx,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const vector_reference_t<T,VX>& _xx,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(T(1)/T(_i),_xx.const_ref());
}

/// operator*:i,ax=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const int& _i,const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_i),_ax.x);
}

/// operator*:ax,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator*(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a*T(_i),_ax.x);
}

/// operator/:ax,i=>ax
template <typename T,typename VX>Expression<Type_vector,Expr_ax,T,VX>
operator/(const Expression<Type_vector,Expr_ax,T,VX>& _ax,const int& _i) {
 return Expression<Type_vector,Expr_ax,T,VX>(_ax.a/T(_i),_ax.x);
}

/// operator*:i,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const int& _i,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_i),_A);
}

/// operator*:A,i=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_const_reference_t<T,MA>& _A,const int& _i) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_i),_A);
}

/// operator*:i,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const int& _i,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_i),_AA.const_ref());
}

/// operator*:AA,i=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const matrix_reference_t<T,MA>& _AA,const int& _i) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(_i),_AA.const_ref());
}

/// operator/:i,A=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const int& _i,const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_i),_A);
}

/// operator/:i,AA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const int& _i,const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_i),_AA.const_ref());
}

/// operator/:AA,i=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const matrix_reference_t<T,MA>& _AA,const int& _i) {
 return Expression<Type_matrix,Expr_aA,T,MA>(T(1)/T(_i),_AA.const_ref());
}

/// operator*:i,aAB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const int& _i,const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_i),_aAB.A,_aAB.B);
}

/// operator*:aAB,i=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const int& _i) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a*T(_i),_aAB.A,_aAB.B);
}

/// operator/:aAB,i=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator/(const Expression<Type_matrix,Expr_aAB,T,MA,MB>& _aAB,const int& _i) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aAB.a/T(_i),_aAB.A,_aAB.B);
}

/// operator*:i,aAx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const int& _i,const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_i),_aAx.A,_aAx.x);
}

/// operator*:aAx,i=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const int& _i) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a*T(_i),_aAx.A,_aAx.x);
}

/// operator/:aAx,i=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator/(const Expression<Type_vector,Expr_aAx,T,MA,VX>& _aAx,const int& _i) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aAx.a/T(_i),_aAx.A,_aAx.x);
}

/// operator*:i,aA=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const int& _i,const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_i),_aA.A);
}

/// operator*:aA,i=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const int& _i) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a*T(_i),_aA.A);
}

/// operator/:aA,i=>aA
template <typename T,typename MA>Expression<Type_matrix,Expr_aA,T,MA>
operator/(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const int& _i) {
 return Expression<Type_matrix,Expr_aA,T,MA>(_aA.a/T(_i),_aA.A);
}

/// operator*:i,aAtA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const int& _i,const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_i),_aAtA.A);
}

/// operator*:aAtA,i=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator*(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const int& _i) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a*T(_i),_aAtA.A);
}

/// operator/:aAtA,i=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
operator/(const Expression<Type_matrix,Expr_aAtA,T,MA>& _aAtA,const int& _i) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(_aAtA.a/T(_i),_aAtA.A);
}

/// operator*:i,axxt=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const int& _i,const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_i),_axxt.x);
}

/// operator*:axxt,i=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
operator*(const Expression<Type_matrix,Expr_axxt,T,VX>& _axxt,const int& _i) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(_axxt.a*T(_i),_axxt.x);
}

/// operator*:A,x=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const matrix_const_reference_t<T,MA>& _A,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(T(1),_A,_x);
}

/// operator*:A,xx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const matrix_const_reference_t<T,MA>& _A,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(T(1),_A,_xx.const_ref());
}

/// operator*:AA,x=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const matrix_reference_t<T,MA>& _AA,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(T(1),_AA.const_ref(),_x);
}

/// operator*:AA,xx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const matrix_reference_t<T,MA>& _AA,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(T(1),_AA.const_ref(),_xx.const_ref());
}

/// operator*:aA,x=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aA.a,_aA.A,_x);
}

/// operator*:aA,xx=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aA.a,_aA.A,_xx.const_ref());
}

/// operator*:aA,ax=>aAx
template <typename T,typename MA,typename VX>Expression<Type_vector,Expr_aAx,T,MA,VX>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const Expression<Type_vector,Expr_ax,T,VX>& _ax) {
 return Expression<Type_vector,Expr_aAx,T,MA,VX>(_aA.a*_ax.a,_aA.a,_ax.x);
}

/// operator*:A,B=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_const_reference_t<T,MA>& _A,const matrix_const_reference_t<T,MB>& _B) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(1,_A,_B);
}

/// operator*:AA,B=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_reference_t<T,MA>& _AA,const matrix_const_reference_t<T,MB>& _B) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(1,_AA.const_ref(),_B);
}

/// operator*:A,BB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_const_reference_t<T,MA>& _A,const matrix_reference_t<T,MB>& _BB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(1,_A,_BB.const_ref());
}

/// operator*:AA,BB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_reference_t<T,MA>& _AA,const matrix_reference_t<T,MB>& _BB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(1,_AA.const_ref(),_BB.const_ref());
}

/// operator*:aA,B=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const matrix_const_reference_t<T,MB>& _B) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aA.a,_aA.A,_B);
}

/// operator*:B,aA=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_const_reference_t<T,MB>& _B,const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aA.a,_aA.A,_B);
}

/// operator*:aA,BB=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const Expression<Type_matrix,Expr_aA,T,MA>& _aA,const matrix_reference_t<T,MB>& _BB) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aA.a,_aA.A,_BB.const_ref());
}

/// operator*:BB,aA=>aAB
template <typename T,typename MA,typename MB>Expression<Type_matrix,Expr_aAB,T,MA,MB>
operator*(const matrix_reference_t<T,MB>& _BB,const Expression<Type_matrix,Expr_aA,T,MA>& _aA) {
 return Expression<Type_matrix,Expr_aAB,T,MA,MB>(_aA.a,_aA.A,_BB.const_ref());
}

/// outer_prod:x=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
outer_prod(const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(T(1),_x);
}

/// outer_prod:xx=>axxt
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
outer_prod(const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_matrix,Expr_axxt,T,VX>(T(1),_xx.const_ref());
}

/// ata:A=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
ata(const matrix_const_reference_t<T,MA>& _A) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(T(1),_A);
}

/// ata:AA=>aAtA
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
ata(const matrix_reference_t<T,MA>& _AA) {
 return Expression<Type_matrix,Expr_aAtA,T,MA>(T(1),_AA.const_ref());
}

/// all:a=>all
template <typename T>Expression<Type_initializer,Expr_all,T>
all(const T& _a) {
 return Expression<Type_initializer,Expr_all,T>(_a);
}

/// operator-:all=>all
template <typename T>Expression<Type_initializer,Expr_all,T>
operator-(const Expression<Type_initializer,Expr_all,T>& _all) {
 return Expression<Type_initializer,Expr_all,T>(-_all.a);
}

/// unit:i=>unit
inline Expression<Type_initializer,Expr_unit>
unit(const int& _i) {
 return Expression<Type_initializer,Expr_unit>(_i,1.0);
}

/// unit:i,a=>unit
template <typename T>Expression<Type_initializer,Expr_unit>
unit(const int& _i,const T& _a) {
 return Expression<Type_initializer,Expr_unit>(_i,_a);
}

/// operator-:unit=>unit
inline Expression<Type_initializer,Expr_unit>
operator-(const Expression<Type_initializer,Expr_unit>& _unit) {
 return Expression<Type_initializer,Expr_unit>(_unit.i,-_unit.d);
}

/// operator*:unit,d=>unit
inline Expression<Type_initializer,Expr_unit>
operator*(const Expression<Type_initializer,Expr_unit>& _unit,const double& _d) {
 return Expression<Type_initializer,Expr_unit>(_unit.i,_unit.d*_d);
}

/// operator*:d,unit=>unit
inline Expression<Type_initializer,Expr_unit>
operator*(const double& _d,const Expression<Type_initializer,Expr_unit>& _unit) {
 return Expression<Type_initializer,Expr_unit>(_unit.i,_unit.d*_d);
}

/// operator/:unit,d=>unit
inline Expression<Type_initializer,Expr_unit>
operator/(const Expression<Type_initializer,Expr_unit>& _unit,const double& _d) {
 return Expression<Type_initializer,Expr_unit>(_unit.i,_unit.d/_d);
}

/// operator-:eye=>seye
inline Expression<Type_initializer,Expr_seye>
operator-(const Expression<Type_initializer,Expr_eye>& _eye) {
 return Expression<Type_initializer,Expr_seye>(-1.0);
}

/// operator*:eye,d=>seye
inline Expression<Type_initializer,Expr_seye>
operator*(const Expression<Type_initializer,Expr_eye>& _eye,const double& _d) {
 return Expression<Type_initializer,Expr_seye>(_d);
}

/// operator*:d,eye=>seye
inline Expression<Type_initializer,Expr_seye>
operator*(const double& _d,const Expression<Type_initializer,Expr_eye>& _eye) {
 return Expression<Type_initializer,Expr_seye>(_d);
}

/// operator/:eye,d=>seye
inline Expression<Type_initializer,Expr_seye>
operator/(const Expression<Type_initializer,Expr_eye>& _eye,const double& _d) {
 return Expression<Type_initializer,Expr_seye>(1.0/_d);
}

/// operator-:seye=>seye
inline Expression<Type_initializer,Expr_seye>
operator-(const Expression<Type_initializer,Expr_seye>& _seye) {
 return Expression<Type_initializer,Expr_seye>(-_seye.d);
}

/// operator*:seye,d=>seye
inline Expression<Type_initializer,Expr_seye>
operator*(const Expression<Type_initializer,Expr_seye>& _seye,const double& _d) {
 return Expression<Type_initializer,Expr_seye>(_seye.d*_d);
}

/// operator*:d,seye=>seye
inline Expression<Type_initializer,Expr_seye>
operator*(const double& _d,const Expression<Type_initializer,Expr_seye>& _seye) {
 return Expression<Type_initializer,Expr_seye>(_seye.d*_d);
}

/// operator/:seye,d=>seye
inline Expression<Type_initializer,Expr_seye>
operator/(const Expression<Type_initializer,Expr_seye>& _seye,const double& _d) {
 return Expression<Type_initializer,Expr_seye>(_seye.d/_d);
}

/// diag:x=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
diag(const vector_const_reference_t<T,VX>& _x) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(1),_x);
}

/// diag:xx=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
diag(const vector_reference_t<T,VX>& _xx) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(1),_xx);
}

/// operator-:diag_mat=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator-(const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(-_diag_mat.a,_diag_mat.x);
}

/// operator*:diag_mat,i=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator*(const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat,const int& _i) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(_i)*_diag_mat.a,_diag_mat.x);
}

/// operator*:i,diag_mat=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator*(const int& _i,const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(_i)*_diag_mat.a,_diag_mat.x);
}

/// operator/:diag_mat,i=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator/(const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat,const int& _i) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(_diag_mat.a/T(_i),_diag_mat.x);
}

/// operator/:i,diag_mat=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator/(const int& _i,const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(_diag_mat.a/T(_i),_diag_mat.x);
}

/// operator*:diag_mat,d=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator*(const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat,const double& _d) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(_d)*_diag_mat.a,_diag_mat.x);
}

/// operator*:d,diag_mat=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator*(const double& _d,const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(T(_d)*_diag_mat.a,_diag_mat.x);
}

/// operator/:diag_mat,d=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator/(const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat,const double& _d) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(_diag_mat.a/T(_d),_diag_mat.x);
}

/// operator/:d,diag_mat=>diag_mat
template <typename T,typename VX>Expression<Type_matrix,Expr_diag_mat,T,VX>
operator/(const double& _d,const Expression<Type_matrix,Expr_diag_mat,T,VX>& _diag_mat) {
 return Expression<Type_matrix,Expr_diag_mat,T,VX>(_diag_mat.a/T(_d),_diag_mat.x);
}


# define VC_MATH_BLAS_MATRIX_EXPR_HH
#  include "blas_matrix_expr_rules_inc.hh"
# undef VC_MATH_BLAS_MATRIX_EXPR_HH

//
// THIS CODE WAS GENERATED AUTOMATICALLY. DO NOT EDIT !!!
//

# endif // DOXYGEN_SKIP

