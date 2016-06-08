//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_MATH_MVMONOMIAL_HH
#define VC_MATH_MVMONOMIAL_HH


//== INCLUDES =================================================================

#include <cassert>

#include "../base/memory.hh"


//== CLASS DEFINITION =========================================================
namespace VC {
namespace math {


/// Polynomials: evaluation of bases, fitting, evaluation.
namespace polynomials {

# ifndef DOXYGEN_SKIP

/** \class _MVMonomialBasis mvmonomials.hh
    \ingroup vc_math
    \sa MVMonomialBasis
*/
template <typename T,unsigned D,unsigned Q>
class _MVMonomialBasis {
public:
  typedef T value_type; //!< scalar type

  enum { Dimension=D }; //!< dimension of domain
  enum { MaxDegree=Q }; //!< maxmimum degree

  /// number of coefficients/basis functions = dimension of basis
  enum { NCoeff=_MVMonomialBasis<T,D,Q-1>::NCoeff+_MVMonomialBasis<T,D-1,Q>::NCoeff };
  /// number of basis functions with degree == MaxDegree
  enum { NCoeffMaxDegree=NCoeff-_MVMonomialBasis<T,D,Q-1>::NCoeff };

  typedef _MVMonomialBasis<T,D,Q-1> LowerDegree;

  template <unsigned Di,unsigned J>
  class _Recursion {
  public:
    typedef _Recursion<Di-1,J+_MVMonomialBasis<T,Di,Q>::NCoeffMaxDegree> _RecursionStep;

    static void recurse(const T* _x,T* _b,unsigned _ldb=1) {
      // multiply
      T* bo=_b+J*_ldb;
      T* bi=_b+(NCoeff-_MVMonomialBasis<T,Di,Q>::NCoeffMaxDegree)*_ldb;
      for (unsigned j=0;j<_MVMonomialBasis<T,Di,Q>::NCoeffMaxDegree;++j,bo+=_ldb,bi+=_ldb)
        *bo=(*bi)*(*_x);
        //_b[(J+j)*_ldb]=_b[(NCoeff-_MVMonomialBasis<T,Di,Q>::NCoeffMaxDegree+j)*_ldb]*(*_x);

      // recurse in dimensions (next _x - coordinate)
      _RecursionStep::recurse(_x+1,_b,_ldb);
    }
  };

  typedef _Recursion<D,NCoeff> _RecursionStart;

  template <unsigned J>
  class _Recursion<0,J> {
  public:
  static void recurse(const T*,T*,unsigned) {}
  };

  /** Evaluate all basis functions at _x and write result to _b.
      \param _x domain point in R(T)^Dimension
      \param _b array storing result (NCoeff values)
      \param _ldb leading dimension of _b (e.g., _b is matrix in FORTRAN order)
   */
  static void evaluate(const T* _x,T* _b,unsigned _ldb=1) {
    LowerDegree::evaluate(_x,_b,_ldb);
    LowerDegree::_RecursionStart::recurse(_x,_b,_ldb);
  }
};

template <typename T,unsigned Q>
    class _MVMonomialBasis<T,0,Q> {
public:
  typedef T value_type;
  enum { Dimension=0 };
  enum { MaxDegree=Q };
  enum { NCoeff=1 };
};

template <typename T,unsigned D>
class _MVMonomialBasis<T,D,0> {
public:
  typedef T value_type;
  enum { Dimension=D };
  enum { MaxDegree=0 };
  enum { NCoeff=1 };

  template <unsigned Di,unsigned J>
  class _Recursion {
  public:
    static void recurse(const T* _x,T* _b,unsigned _ldb) {
      for (unsigned j=0;j<Dimension;++j)
        _b[(1+j*_ldb)]=_x[j]; // linear terms
    }
  };
  typedef _Recursion<0,0> _RecursionStart;

  static void evaluate(const T* /*_x*/,T* _b,unsigned /*_ldb=1*/) { *_b=T(1); }
  static void evaluate(const T* /*_x*/,T* _b) { *_b=T(1); }
};

template <typename T>
class _MVMonomialBasis<T,0,0> {
public:
  typedef T value_type;
  enum { Dimension=0 };
  enum { MaxDegree=0 };
  enum { NCoeff=1 };
};

# endif // DOXYGEN_SKIP

//-----------------------------------------------------------------------------

/** \class MVMonomialBasis mvmonomials.hh
    \ingroup vc_math
    \brief Evaluation of multivariate polynomial basis in monomial form.
    \arg T real type for domain points ans results
    \arg D dimension of domain
    \arg Q maxmimum (total) polynomial degree
*/
template <typename T,unsigned D,unsigned Q>
class MVMonomialBasis {
public:
  typedef T value_type; //!< scalar type

  enum { Dimension=D }; //!< dimension of domain
  enum { MaxDegree=Q }; //!< maxmimum degree

  /// number of coefficients/basis functions = dimension of basis
  enum { NCoeff=_MVMonomialBasis<T,D,Q-1>::NCoeff+_MVMonomialBasis<T,D-1,Q>::NCoeff };
  /// number of basis functions with degree == MaxDegree
  enum { NCoeffMaxDegree=NCoeff-_MVMonomialBasis<T,D,Q-1>::NCoeff };

  /** Evaluate all basis functions at _x and write result to _b.
      \param _x domain point in R(T)^Dimension
      \param _b array storing result (NCoeff values)
      \param _ldb leading dimension of _b (e.g., _b is matrix in FORTRAN order)
   */
  static void evaluate(const T* _x,T* _b,unsigned _ldb=1) {
    _MVMonomialBasis<T,D,Q>::evaluate(_x,_b,_ldb);
  }
};

# ifndef DOXYGEN_SKIP
#  include "mvmonomials.inc" // specializations
# endif

//-----------------------------------------------------------------------------

/** Get number of coefficients for multivariate polynomial.
    \ingroup vc_math
    \param _dimension of domain
    \param _degree maxmimum (total) degree
    \return dimension of polynomial basis
 */
unsigned ncoeff(unsigned _dimension,unsigned _degree);

/// returned by d_mvbasis_evaluation, evaluate in FORTRAN order
typedef void (*d_mvbasis_evaluation_t)
             (const double* _x,double* _b,unsigned _n,unsigned _ldb);
/// returned by f_mvbasis_evaluation, evaluate in FORTRAN order
typedef void (*f_mvbasis_evaluation_t)
             (const float* _x,float* _b,unsigned _n,unsigned _ldb);

/// return fuction to evaluate polynomial at multiple points (or null)
d_mvbasis_evaluation_t d_mvbasis_evaluation(unsigned _dimension,unsigned _degree);
/// return fuction to evaluate polynomial at multiple points (or null)
f_mvbasis_evaluation_t f_mvbasis_evaluation(unsigned _dimension,unsigned _degree);


/** \class MVMonomialBasisDQ_base mvmonomials.hh
    \ingroup vc_math
    \brief Base class and interface to MVMonomialBasisDQ.
    \arg T real type for domain points ans results
    \sa MVMonomialBasisDQ
*/
template <typename T>
class MVMonomialBasisDQ_base {
public:
  unsigned ncoeff() const { return m_ncoeff; }       //!< get number of coefficients
  unsigned dimension() const { return m_dimension; } //!< get dimension
  unsigned degree() const { return m_degree; }       //!< get degree
  /// evaluate _n points _x to _b
  inline void evaluate(const T* _x,T* _b,unsigned _n=1,unsigned _ldb=1) {
    assert(m_eval!=0);
    m_eval(_x,_b,_n,_ldb);
  }
protected:
  MVMonomialBasisDQ_base(unsigned _dimension,unsigned _degree)
    : m_eval(0), m_dimension(_dimension), m_degree(_degree),
      m_ncoeff(::VC::math::polynomials::ncoeff(_dimension,_degree)) {}

  typedef void (*mvbasis_evaluation_t)(const T*,T*,unsigned,unsigned);

  mvbasis_evaluation_t m_eval;
  unsigned             m_dimension;
  unsigned             m_degree;
  unsigned             m_ncoeff;
};

/** \class MVMonomialBasisDQ mvmonomials.hh
    \ingroup vc_math
    \brief Like MVMonomialBasis but variable dimension and degree.
    \arg T real type for domain points ans results
    \sa MVMonomialBasisDQ_base


    \todo evaluate polynomial (simple)

    \todo enumerate get transformations (weights and indices) to get partials
    (especially Jacobian)
*/
template <typename T>
class MVMonomialBasisDQ : public MVMonomialBasisDQ_base<T> {
};

# ifndef DOXYGEN_SKIP

template <>
class MVMonomialBasisDQ<double> : public MVMonomialBasisDQ_base<double> {
public:
  typedef double value_t;
  typedef MVMonomialBasisDQ_base<value_t>::mvbasis_evaluation_t
  mvbasis_evaluation_t;

  MVMonomialBasisDQ(unsigned _dimension=1,unsigned _degree=1)
    : MVMonomialBasisDQ_base<value_t>(_dimension,_degree) {
    m_eval=d_mvbasis_evaluation(_dimension,_degree);
    assert(m_eval!=0 && "not implemented");
  }
};

template <>
class MVMonomialBasisDQ<float> : public MVMonomialBasisDQ_base<float> {
public:
  typedef float value_t;
  typedef MVMonomialBasisDQ_base<value_t>::mvbasis_evaluation_t
  mvbasis_evaluation_t;

  MVMonomialBasisDQ(unsigned _dimension=1,unsigned _degree=1)
    : MVMonomialBasisDQ_base<value_t>(_dimension,_degree) {
    m_eval=f_mvbasis_evaluation(_dimension,_degree);
    assert(m_eval!=0 && "not implemented");
  }
};

# endif // DOXYGEN_SKIP

//=============================================================================
} // namespace polynomials
} // namespace math
} // namespace VC
//=============================================================================
#endif // VC_MATH_MVMONOMIAL_HH defined