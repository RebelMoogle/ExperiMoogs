//=============================================================================
// $TEMPLATE_HEADLINE$
// ----------------------------------------------------------------------------
// $Id$
// $Revision$
//
//=============================================================================

#ifndef VC_MATH_LSQFIT_HH
#define VC_MATH_LSQFIT_HH


//== INCLUDES =================================================================

#include <cassert>
#include <vector>

#include "lapack.hh"


//== CLASS DEFINITION =========================================================
namespace VC {
namespace math {

/** \class LsqFitContext lsqfit.hh
    \ingroup vc_math
    \brief Store context information and provide least-squares fit.

    This context information stores the problem and its solution as
    well as all temporary memory required for solving.
    Instances of, e.g., polynomials::LsqFitMonomials may share conetxts (but
    not use them concurrently).

    \todo implement normal equations

    \sa polynomials::LsqFitMonomials
 */
template <typename T>
class LsqFitContext {
public:
  LsqFitContext()
  : curr_setup(UNDEF_SETUP), prev_setup(UNDEF_SETUP),
    m(0), n(0), nrhs(0) {}

  /// Allocate storage and get dimensions for subsequent call to gelsd().
  void setup_gelsd(unsigned _m,unsigned _n,unsigned _nrhs) {
    A.resize(_m*_n);
    b.resize(_m*_nrhs);
    sigma.resize(std::min(_m,_n));

    if (_m>m || _n>n || _nrhs>nrhs ||
        (prev_setup!=GELSD_SETUP && curr_setup!=GELSD_SETUP)) {
      // estimate workspace
      unsigned lwork=lapack::gelsd_lwork(_m,_n,_nrhs);

      if (work.size()<lwork)
        work.resize(lwork);

      unsigned liwork=lapack::gelsd_liwork(_m,_n,_nrhs);

      if (iwork.size()<liwork)
        iwork.resize(liwork);
    }
    m=_m; n=_n; nrhs=_nrhs; curr_setup=GELSD_SETUP;
  }
  /// Allocate storage and get dimensions for subsequent call to gelsy().
  void setup_gelsy(unsigned _m,unsigned _n,unsigned _nrhs) {
    A.resize(_m*_n);
    b.resize(_m*_nrhs);
    sigma.resize(0);

    if (_m>m || _n>n || _nrhs>nrhs ||
        (prev_setup!=GELSY_SETUP || curr_setup!=GELSY_SETUP)) {
      // estimate workspace

      unsigned lwork=lapack::gelsy_lwork(_m,_n,_nrhs);

      if (work.size()<lwork)
        work.resize(lwork);

      if (iwork.size()<_n) // permutation JPVT
        iwork.resize(_n);
    }
    m=_m; n=_n; nrhs=_nrhs; curr_setup=GELSY_SETUP;
  }

  /** Find least-squares solution of A*x=b.
      Calls VC::math::lapack::gelsd() (SVD).
  */
  bool solve_gelsd(double _condition) {
    assert(curr_setup==GELSD_SETUP);

    int irank;
    info=lapack::gelsd(m,n,nrhs, &A[0],m, &b[0],m, &sigma[0],
                       T(1)/_condition,irank,&work[0],int(work.size()),&iwork[0]);

    rank=irank;

    if (unsigned(work[0])>work.capacity())
      work.resize(std::max(work.size()*2,
                                    size_t(work[0])));

    unsigned exp_rank=std::min(m,n);

    prev_setup=curr_setup;
    curr_setup=UNDEF_SETUP;

    return info==0 && rank==exp_rank;
  }

  /** Find least-squares/least-norom solution of A*x=b.
      Calls VC::math::lapack::gelsy() (QR factorization).
  */
  bool solve_gelsy(double _condition) {
    assert(curr_setup==GELSY_SETUP);

    int irank;
    info=lapack::gelsy(m,n,nrhs, &A[0],m, &b[0],m, &iwork[0],
                       T(1)/_condition,irank,&work[0],int(work.size()));

    rank=irank;

    if (unsigned(work[0])>work.capacity())
      work.resize(std::max(work.size()*2,
                                    size_t(work[0])));

    unsigned exp_rank=std::min(m,n);

    prev_setup=curr_setup;
    curr_setup=UNDEF_SETUP;

    return info==0 && rank==exp_rank;
  }

  enum Setup {
    UNDEF_SETUP, GELSD_SETUP, GELSY_SETUP
  };

  Setup    curr_setup;    //!< set be setup_*(), reset by solve_*()
  Setup    prev_setup;    //!< previous setup
  unsigned m,n;           //!< dimensions of system
  unsigned nrhs;          //!< number of right-hand-sides
  std::vector<T> A;       //!< system matrix
  std::vector<T> b;       //!< right-hand-side (input) / solution (output)
  std::vector<T> sigma;   //!< singular values set by gelsd
  unsigned rank;          //!< effctive rank w.r.t. prescribed condition
  int      info;          //!< LAPACK return code

  /// get condition number of A from sigma
  T condition() const {
    assert(!sigma.empty() && "only solve_gelsd() yields singular values");
    return sigma.front()/sigma.back();
  }

  std::vector<T>                work;  //!< floating point storage
  std::vector<lapack::INTEGER>  iwork; //!< integer storage
};

//-----------------------------------------------------------------------------

//=============================================================================
} // namespace math
} // namespace VC
//=============================================================================
#endif // VC_MATH_LSQFIT_HH defined
