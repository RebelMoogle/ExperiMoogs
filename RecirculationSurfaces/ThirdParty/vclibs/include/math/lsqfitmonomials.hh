//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_MATH_POLYNOMIALS_LSQFITMONOMIALS_HH
#define VC_MATH_POLYNOMIALS_LSQFITMONOMIALS_HH


//== INCLUDES =================================================================

#include <cassert>
#include <limits>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "lsqfit.hh"
#include "mvmonomials.hh"


//== CLASS DEFINITION =========================================================
namespace VC {
namespace math {
namespace polynomials {

/** \class LsqFitMonomials mvmonomials.hh
    \ingroup vc_math
    \brief (Weighted) least-squares fit of polynomial in monomial form.
    \tparam T real type for domain points ans results
    \tparam D dimension of domain
    \tparam Q maxmimum (total) polynomial degree
    \sa MVMonomialBasis, VC::math::LsqFitContext
*/
template <typename T,unsigned D,unsigned Q>
class LsqFitMonomials {
public:
  typedef T value_type; //!< scalar type

  enum { Dimension=D }; //!< dimension of domain
  enum { MaxDegree=Q }; //!< maxmimum degree

  typedef MVMonomialBasis<T,D,Q> basis_t; //!< monimal basis

  enum { NCoeff=basis_t::NCoeff }; //!< dimension of basis

  ///
  LsqFitMonomials(unsigned _data_dim,LsqFitContext<T>* _context)
    : context(_context) {
    m_data_dim=_data_dim;
  }

  LsqFitContext<T>* context;

  /// get dimension of data (=# of right hand sides)
  unsigned data_dim() const { return m_data_dim; }

  /// preallocate storage for handling _m_max samples
  void reserve(unsigned _m_max) {
    context->setup_gelsd(_m_max,NCoeff,m_data_dim);
  }

  /** Least-squares fit of _m samples (_points,_values) using SVD.

      All arrays are stored in C-style.
      \param _m number of samples
      \param _points position of samples
      \param _values data_dim() dimensional data on right-hand-side
      \param _weights weighted least-squares fit if `!=0`
      \param _coeff coefficients of monomial basis as \a output C-order
      \param _condition maxmium condition number of system matrix
      \return false on failure, i.e., if `context.info!=0` or
      `context->rank<min(NCoeff,m)` (see status information context)
   */
  bool
  wfit_svd(unsigned _m,const T* _points,const T* _values,const T* _weights,
           T* _coeff,T _condition) {
    context->setup_gelsd(_m,NCoeff,m_data_dim);

    unsigned ldb=std::max(_m,unsigned(NCoeff));

    /// copy rhs
    lapack::cp_c2f(_values,&(context->b[0]),_m,m_data_dim,ldb);
    /// evaluate basis
    for (unsigned j=0;j<_m;++j)
      basis_t::evaluate(_points+j*m_data_dim,&context->A[j],_m);

    if (_weights!=0)
      apply_weights(_m,_weights);

    if (context->solve_gelsd(_condition)) {
      /// copy solution
      lapack::cp_f2c(&(context->b[0]),_coeff,NCoeff,m_data_dim,ldb);
      return true;
    }
    else
      return false;
  }

  /// Least-squares fit: same wfit_svd() w/o _weights
  bool
  fit_svd(unsigned _m,const T* _points,const T* _values,
          T* _coeff,T _condition) {
    return wfit_svd(_m,_points,_values,(const T*) 0,_coeff,_condition);
  }

  /** Least-squares fit of _m samples (_points,_values) using QR-decomposition.

      All arrays are stored in C-style.
      \param _m number of samples
      \param _points position of samples
      \param _values data_dim() dimensional data on right-hand-side
      \param _weights weighted least-squares fit if `!=0`
      \param _coeff coefficients of monomial basis as \a output C-order
      \param _condition maxmium condition number of system matrix
      \return false on failure, i.e., if `context.info!=0` or
      `context->rank<min(NCoeff,m)` (see status information context)
   */
  bool
  wfit_qr(unsigned _m,const T* _points,const T* _values,const T* _weights,
          T* _coeff,T _condition) {
    context->setup_gelsy(_m,NCoeff,m_data_dim);

    unsigned ldb=std::max(_m,unsigned(NCoeff));

    /// copy rhs
    lapack::cp_c2f(_values,&(context->b[0]),_m,m_data_dim,ldb);
    /// evaluate basis
    for (unsigned j=0;j<_m;++j)
      basis_t::evaluate(_points+j*m_data_dim,&context->A[j],_m);

    if (_weights!=0)
      apply_weights(_m,_weights);

    if (context->solve_gelsy(_condition)) {
      /// copy solution
      lapack::cp_f2c(&(context->b[0]),_coeff,NCoeff,m_data_dim,ldb);
      return true;
    }
    else
      return false;
  }

  /// Least-squares fit: same wfit_qr() w/o _weights.
  bool
  fit_qr(unsigned _m,const T* _points,const T* _values,
         T* _coeff,T _condition) {
    return wfit_qr(_m,_points,_values,(const T*) 0,_coeff,_condition);
  }

protected:

  void apply_weights(unsigned _m,const T* _w) {
    assert(_w!=0);
    T* a=&context->A[0];
    for (unsigned i=0;i<NCoeff;++i) {
      const T* w=_w;
      for (unsigned j=0;j<_m;++j,++w,++a)
        *a*=*w;
    }
    T* b=&context->b[0];
    for (unsigned i=0;i<m_data_dim;++i) {
      const T* w=_w;
      for (unsigned j=0;j<_m;++j,++w,++b)
        *b*=*w;
    }
  }

  unsigned                      m_data_dim; //!< dimension of data
  std::vector<T>                m_A;        //!< evaluation of basis functions
  std::vector<T>                m_sigma;    //!< singular values
  std::vector<T>                m_b;        //!< right hand side
  std::vector<T>                m_work;     //!< for LAPACK
  std::vector<lapack::INTEGER>  m_iwork; //!< for LAPACK
};

//-----------------------------------------------------------------------------

/** \class LsqFitMonomialsDQ mvmonomials.hh
    \ingroup vc_math
    \brief Least-squares fit of polynomial in monomial form.

    Same as LsqFitMonomials but with variable dimension and degree
    \sa LsqFitMonomials, MVMonomialBasisDQ
*/
template <typename T>
class LsqFitMonomialsDQ {
public:
  typedef T value_type; //!< scalar type

  ///
  LsqFitMonomialsDQ(unsigned _dimension,unsigned _data_dim,
                    LsqFitContext<T>* _context)
    : curr_basis(_dimension,0), context(_context) {
    m_dim=_dimension;
    m_data_dim=_data_dim;
  }

  /// get dimension of domain
  unsigned dim() const { return m_dim; }
  /// get dimension of data (=# of right hand sides)
  unsigned data_dim() const { return m_data_dim; }

  /// preallocate storage for handling _m_max samples
  void reserve(unsigned _m_max,unsigned _max_degree=6) {
    context->setup_gelsd(_m_max,ncoeff(m_dim,_max_degree),m_data_dim);
  }

  MVMonomialBasisDQ<T> curr_basis; //!< current basis
  LsqFitContext<T>*    context;

  /** Least-squares fit of _m samples (_points,_values) using SVD.

      All arrays are stored in C-style.
      \param _degree polynomial degree
      \param _m number of samples
      \param _points position of samples
      \param _values data_dim() dimensional data on right-hand-side
      \param _weights weighted least-squares fit if `!=0`
      \param _coeff coefficients of monomial basis as \a output C-order
      \param _condition maxmium condition number of system matrix
      \return false on failure, i.e., if `context.info!=0` or
      `context->rank<min(NCoeff,m)` (see status information context)
   */
  bool
  wfit_svd(unsigned _degree,
           unsigned _m,const T* _points,const T* _values,const T* _weights,
           T* _coeff,T _condition) {
    if (curr_basis.degree()!=_degree)
      curr_basis=MVMonomialBasisDQ<T>(m_dim,_degree);

    unsigned ncoeff=curr_basis.ncoeff();
    context->setup_gelsd(_m,ncoeff,m_data_dim);

    unsigned ldb=std::max(_m,ncoeff);

    /// copy rhs
    lapack::cp_c2f(_values,&(context->b[0]),_m,m_data_dim,ldb);
    /// evaluate basis
    curr_basis.evaluate(_points,&context->A[0],_m,_m);

    if (_weights!=0)
      apply_weights(_m,ncoeff,_weights);

    if (context->solve_gelsd(_condition)) {
      /// copy solution
      lapack::cp_f2c(&(context->b[0]),_coeff,ncoeff,m_data_dim,ldb);
      return true;
    }
    else
      return false;
  }

  /// Least-squares fitL: same as wfit_svd() w/o _weights.
  bool
  fit_svd(unsigned _degree,unsigned _m,const T* _points,const T* _values,
          T* _coeff,T _condition) {
    return wfit_svd(_degree,_m,_points,_values,0,_coeff,_condition);
  }

  /** Least-squares fit of _m samples (_points,_values) using QR-decomposition.

      All arrays are stored in C-style.
      \param _degree polynomial degree
      \param _m number of samples
      \param _points position of samples
      \param _values data_dim() dimensional data on right-hand-side
      \param _weights weighted least-squares fit if `!=0`
      \param _coeff coefficients of monomial basis as \a output C-order
      \param _condition maxmium condition number of system matrix
      \return false on failure, i.e., if `context.info!=0` or
      `context->rank<min(NCoeff,m)` (see status information context)
   */
  bool
  wfit_qr(unsigned _degree,
          unsigned _m,const T* _points,const T* _values,const T* _weights,
          T* _coeff,T _condition) {
    if (curr_basis.degree()!=_degree)
      curr_basis=MVMonomialBasisDQ<T>(m_dim,_degree);

    unsigned ncoeff=curr_basis.ncoeff();
    context->setup_gelsy(_m,ncoeff,m_data_dim);

    unsigned ldb=std::max(_m,ncoeff);

    /// copy rhs
    lapack::cp_c2f(_values,&(context->b[0]),_m,m_data_dim,ldb);
    /// evaluate basis
    curr_basis.evaluate(_points,&context->A[0],_m,_m);

    if (_weights!=0)
      apply_weights(_m,ncoeff,_weights);

    if (context->solve_gelsy(_condition)) {
      /// copy solution
      lapack::cp_f2c(&(context->b[0]),_coeff,ncoeff,m_data_dim,ldb);
      return true;
    }
    else
      return false;
  }

  /// Least-squares fit: same as wfit_qr() w/o _weights
  bool
  fit_qr(unsigned _degree,unsigned _m,const T* _points,const T* _values,
          T* _coeff,T _condition) {
    return wfit_qr(_degree,_m,_points,_values,0,_coeff,_condition);
  }

protected:

  void apply_weights(unsigned _m,unsigned _ncoeff,const T* _w) {
    assert(_w!=0);
    T* a=&context->A[0];
    for (unsigned i=0;i<_ncoeff;++i) {
      const T* w=_w;
      for (unsigned j=0;j<_m;++j,++w,++a)
        *a*=*w;
    }
    T* b=&context->b[0];
    for (unsigned i=0;i<m_data_dim;++i) {
      const T* w=_w;
      for (unsigned j=0;j<_m;++j,++w,++b)
        *b*=*w;
    }
  }

  unsigned                      m_dim;      //!< dimension of domain
  unsigned                      m_data_dim; //!< dimension of data
  std::vector<T>                m_A;        //!< evaluation of basis functions
  std::vector<T>                m_sigma;    //!< singular values
  std::vector<T>                m_b;        //!< right hand side
  std::vector<T>                m_work;     //!< for LAPACK
  std::vector<lapack::INTEGER>  m_iwork; //!< for LAPACK
};

//=============================================================================
} // namespace polynomials
} // namespace math
} // namespace VC
//=============================================================================
#endif // VC_MATH_POLYNOMIALS_LSQFITMONOMIALS_HH defined
