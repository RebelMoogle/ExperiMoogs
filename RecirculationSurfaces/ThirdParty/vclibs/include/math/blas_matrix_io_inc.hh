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

    \arg Provide output to ostream:
    operator<<(ostream&,matrix_const_reference_t), etc.

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================

#ifndef DOXYGEN_SKIP

template <typename T,typename A>
std::ostream& operator<<(std::ostream& _os,
                         const matrix_const_reference_t<T,A> _a) {
  return _os << Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >(_a);
}
template <typename T,typename A>
std::ostream& operator<<(std::ostream& _os,
                         const matrix_reference_t<T,A> _a) {
  return _os << Expression<Type_matrix,Expr_prettyprint,matrix_const_reference_t<T,A> >(_a);
}

//-----------------------------------------------------------------------------

template <typename A>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_matrix,Expr_print,A>& _p) {

  const A& a=_p.a;

  for (int i=0;i<a.m_rows();++i) {
    for (int j=0;j<a.n_cols();++j)
      _os << at(a,i,j) << (j!=a.n_cols()-1 ? ' ' : (i!=a.m_rows()-1 ? ';' : ' '));;
  }
  return _os;
}

//-----------------------------------------------------------------------------

template <typename A>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_matrix,Expr_prettyprint,A>& _p) {

  const A& a=_p.a;

  static const char* type[]={"GE","GB","SY","SB","SP","SP","TR","TB","TP"};

  _os << type[a.matrix_type()] << ':' << a.m() << 'x' << a.n()
      << ':' << char(a.trans())
      << "  [\n";

  for (int i=0;i<a.m_rows();++i) {
    _os << "   ";
    for (int j=0;j<a.n_cols();++j)
      _os << at(a,i,j) << ' ';
    if (i<a.m_rows()-1) _os << std::endl;
  }

  return _os << " ]";
}

//-----------------------------------------------------------------------------

template <typename A>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_matrix,Expr_pmatlab,A>& _p) {

  const A& a=_p.a;

  _os << '[';

  int i,j;
  for (i=0;i<a.m_rows();++i) {
    for (j=0;j<a.n_cols()-1;++j)
      _os << at(a,i,j) << ',';
    _os << at(a,i,j) << (i!=a.m_rows()-1 ? ';' : ']');
  }

  if (a.m()==0 && a.n()==0)
    _os << " ]";

  return _os;
}

//-----------------------------------------------------------------------------

#ifdef _MSC_VER
# pragma warning(push)
# pragma warning(disable:4512)
#endif

struct _PDataOp {
  _PDataOp(std::ostream& _out) : out(_out) {}
  template <typename T> void operator()(const T& _a) { out << _a << ' '; }
  std::ostream& out;
};

#ifdef _MSC_VER
# pragma warning(pop)
#endif

template <typename A>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_matrix,Expr_pdata,A>& _p) {
  Foreach<typename A::value_t,_PDataOp,typename A::matrix_t>
    (_p.a.matrix(),(typename A::value_t*) _p.a.data(),_PDataOp(_os));
  return _os;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template <typename T,typename V>
std::ostream& operator<<(std::ostream& _os,const vector_const_reference_t<T,V>& _x) {
  return _os << Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,V> >(_x);
}
template <typename T,typename V>
std::ostream& operator<<(std::ostream& _os,const vector_reference_t<T,V>& _x) {
  return _os << Expression<Type_vector,Expr_prettyprint,vector_const_reference_t<T,V> >(_x);
}

//-----------------------------------------------------------------------------

template <typename V>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_vector,Expr_print,V>& _p) {

  V x=_p.a;

  for (int i=0;i<x.n();++i)
    _os << _at(x,i) << ' ';

  return _os;
}

//-----------------------------------------------------------------------------

template <typename V>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_vector,Expr_prettyprint,V>& _p) {

  V x=_p.a;

  _os << "vector:" << x.n()
      << "[ ";

  for (int i=0;i<x.n();++i)
    _os << _at(x,i) << ' ';

  return _os << " ]";
}

//-----------------------------------------------------------------------------

template <typename V>
std::ostream& operator<<(std::ostream& _os,
                         const Expression<Type_vector,Expr_pmatlab,V>& _p) {

  V x=_p.a;

  _os << '[';

  int i;
  for (i=0;i<x.n()-1;++i)
    _os << _at(x,i) << ';';

  return _os << _at(x,i) << ']';
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/** Write matrix in MATLAB V4 format (Level 4 MAT).
    \ingroup vc_blas_io
    Simple output to MAT-file, supports a subset of the Level 4 file format.
    Multiple matrices can be output to the same file.
    \param _out output stream
    \param _a matrix
    \param _name matlab name of matrix
    \b Note: The corresponding \c read_mat4() is implemented as a method of
    GE_Matrix.
    \return _out
    \sa GE_Matrix<T>::read_mat4, [\ref vc_mat4]
    Note: this function calls mat4::write_matrix() and will raise an exception
    on IO error.
 */
template <typename T,typename A>
std::ostream& write_mat4(std::ostream& _out,
                         const matrix_const_reference_t<T,A>& _a,
                         const std::string& _name) {
  return mat4::write_matrix(_out,_name,_a.m(),_a.n(),_a.ld(),_a.data());
}

/// provided for convenience: short for write_mat4(trans(trans(_x)),_name) \ingroup vc_blas_io
template <typename T,typename V>
std::ostream& write_mat4(std::ostream& _out,
                         const vector_const_reference_t<T,V>& _x,
                         const std::string& _name) {
  return write_mat4(_out,trans(trans(_x)),_name);
}

# ifndef DOXYGEN_SKIP

template <typename T,typename A>
std::ostream& write_mat4(std::ostream& _out,
                         const matrix_reference_t<T,A>& _a,
                         const std::string& _name) {
  return write_mat4(_out,_a.const_ref(),_name);
}

template <typename T,typename A>
std::ostream& write_mat4(std::ostream& _out,
                         const vector_reference_t<T,A>& _a,
                         const std::string& _name) {
  return write_mat4(_out,_a.const_ref(),_name);
}

# endif

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC