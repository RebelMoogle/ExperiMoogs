#ifndef VC_MATH_BLAS_PROTOTYPES_HH
#define VC_MATH_BLAS_PROTOTYPES_HH

/** \file blas_prototypes.hh Math: BLAS prototypes
 */

/** \defgroup vc_blas_core BLAS core interface
    \ingroup vc_blas
    C/C++ interface to <a href="http://www.netlib.org/blas/">BLAS</a> subroutines.

    An quick reference is found
    <a href="http://www.netlib.org/lapack/lug/node145.html">here</a>
    <a href="http://www.netlib.org/blas/blasqr.ps">[ps]</a>.

    Only a subset of prototypes is provided. Especially no functions dealing with
    complex values.

    <b>Naming conventions</b>

    Number types
    \arg D: double
    \arg S: float
    \arg (Z,C)

    Matrix type
    \arg GE: general
    \arg SY: symmetric
    \arg SB: symmetric banded
    \arg SP: symmetric packed
    \arg TR: triangular
    \arg TB: triangular banded
    \arg TP: triangular packed
    \arg (GB,HE,HB,HP)

    Operations
    \arg MV: matrix-vector multiplication
    \arg SV: solve linear system Ax=b (A is triangular)
    \arg MM: matrix-matrix multiplication
    \arg SM: solve linear system Ax=B (A is triangular)
    \arg ...

    Parameters to switches
    \arg _side: "L" (_A*_B, _A is left) or "R" (_B*_A)
    \arg _uplo: "U" (upper triangular) or "L" (lower triangular)
    \arg _trans: "N" (not transposed) or "T" (transposed) [or "C"]
    \arg _diag: "U" (unit diagonal) or "N" (no unit diagonal)
 */

# if defined(_USE_INTEL_FORTRAN_LIBS)
#  define _BLAS_FUNCTION(uname,lname) uname
# endif

# if !defined(_BLAS_FUNCTION)
#  define _BLAS_FUNCTION(uname,lname) lname ## _
# endif

/** \def _BLAS_FUNCTION
    \ingroup vc_blas
    Map names of BLAS function to system/library dependent symbols.
    Takes upper case and lower case names (w/o underscores) as arguements
    and returns symbol as used in BLAS library.
    For example
    \code
    // with
    #  define _BLAS_FUNCTION(uname,lname) lname ## _
    // and
    # define DSWAP _BLAS_FUNCTION(DSWAP,dswap)
    // "DSWAP" will be replaced by the symbol "dswap_".
    \endcode
    \sa _LAPACK_FUNCTION
*/
# ifdef DOXYGEN_SKIP
#  define _BLAS_FUNCTION "map linker symbols to upper case names"
#  error "doxygen only"
# endif

//
// double precision
//

//=============================================================================

//
// level 1
//

# define DSWAP _BLAS_FUNCTION(DSWAP,dswap)

/** BLAS swap _x and _y
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dswap.f">[reference]</a>
*/
extern "C"
void DSWAP(const ::VC::math::blas::INTEGER* _n,
           double* _x,const ::VC::math::blas::INTEGER* _incx,
           double* _y,const ::VC::math::blas::INTEGER* _incy);

# define DSCAL _BLAS_FUNCTION(DSCAL,dscal)

/** BLAS scal _x = _alpha*_x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dscal.f">[reference]</a>
*/
extern "C"
void DSCAL(const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           double* _x,const ::VC::math::blas::INTEGER* _incx);

# define DROT _BLAS_FUNCTION(DROT,drot)

/** BLAS rot
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/drot.f">[reference]</a>
 */
extern "C"
void DROT(const ::VC::math::blas::INTEGER* _n,
          double* _x,const ::VC::math::blas::INTEGER* _incx,
          double* _y,const ::VC::math::blas::INTEGER* _incy,
          const double* _c,const double* _s);

# define DROTG _BLAS_FUNCTION(DROTG,drotg)

/** BLAS rotg
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/drotg.f">[reference]</a>
 */
extern "C"
void DROTG(double* _a,double* _b,double* _c,double* _s);

# define DCOPY _BLAS_FUNCTION(DCOPY,dcopy)

/** BLAS copy _y = _x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dcopy.f">[reference]</a>
*/
extern "C"
void DCOPY(const ::VC::math::blas::INTEGER* _n,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           double* _y,const ::VC::math::blas::INTEGER* _incy);

# define DAXPY _BLAS_FUNCTION(DAXPY,daxpy)

/** BLAS AXPY _y = _alpha*_x+_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/daxpy.f">[reference]</a>
*/
extern "C"
void DAXPY(const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           double* _y, const ::VC::math::blas::INTEGER* _incy);

# define DDOT _BLAS_FUNCTION(DDOT,ddot)

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
extern "C"
double DDOT(const ::VC::math::blas::INTEGER* _n,
            const double* _x,const ::VC::math::blas::INTEGER* _incx,
            const double* _y,const ::VC::math::blas::INTEGER* _incy);

# define DSDOT _BLAS_FUNCTION(DSDOT,dsdot)

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
extern "C"
double DSDOT(const ::VC::math::blas::INTEGER* _n,
            const float* _x,const ::VC::math::blas::INTEGER* _incx,
            const float* _y,const ::VC::math::blas::INTEGER* _incy);

# define DNRM2 _BLAS_FUNCTION(DNRM2,dnrm2)

/** BLAS 2-norm ||_x||_2
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dnrm2.f">[reference]</a>
*/
extern "C"
double DNRM2(const ::VC::math::blas::INTEGER* _n,
             const double* _x,const ::VC::math::blas::INTEGER* _incx);

# define DASUM _BLAS_FUNCTION(DASUM,dasum)

/** BLAS 1-norm || _x || _1
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dasum.f">[reference]</a>
*/
extern "C"
double DASUM(const ::VC::math::blas::INTEGER* _n,
             const double* _x,const ::VC::math::blas::INTEGER* _incx);



//
// level 2
//

# define DGEMV _BLAS_FUNCTION(DGEMV,dgemv)

/** BLAS DGEMV _y=_alpha*op(_A)*_x+_beta*_y.
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/dgemv.f">[reference]</a>
*/
extern "C"
void DGEMV(const char* _trans,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _beta,
           double* _y, const ::VC::math::blas::INTEGER* _incy);

# define DSYMV _BLAS_FUNCTION(DSYMV,dsymv)

/** BLAS DSYMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsymv.f">[reference]</a>
*/
extern "C"
void DSYMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _beta,
           double* _y, const ::VC::math::blas::INTEGER* _incy);

# define DSBMV _BLAS_FUNCTION(DSBMV,dsbmv)

/** BLAS DSBMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsbmv.f">[reference]</a>
*/
extern "C"
void DSBMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const ::VC::math::blas::INTEGER* _k,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _beta,
           double* _y, const ::VC::math::blas::INTEGER* _incy);

# define DSPMV _BLAS_FUNCTION(DSPMV,dspmv)

/** BLAS DSPMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dspmv.f">[reference]</a>
*/
extern "C"
void DSPMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,const double* _Ap,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _beta,
           double* _y, const ::VC::math::blas::INTEGER* _incy);

# define DTRMV _BLAS_FUNCTION(DTRMV,dtrmv)

/** BLAS DTRMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtrmv.f">[reference]</a>
*/
extern "C"
void DTRMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _x,const ::VC::math::blas::INTEGER* _incx);

# define DTBMV _BLAS_FUNCTION(DTBMV,dtbmv)

/** BLAS DTBMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtbmv.f">[reference]</a>
*/
extern "C"
void DTBMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _x,const ::VC::math::blas::INTEGER* _incx);

# define DTPMV _BLAS_FUNCTION(DTPMV,dtpmv)

/** BLAS DTPMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtpmv.f">[reference]</a>
*/
extern "C"
void DTPMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const double* _Ap,
           double* _x,const ::VC::math::blas::INTEGER* _incx);




# define DTRSV _BLAS_FUNCTION(DTRSV,dtrsv)

/** BLAS DTRSV _x=op(_A)^(-1)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtrsv.f">[reference]</a>
*/
extern "C"
void DTRSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _x,const ::VC::math::blas::INTEGER* _incx);


# define DTBSV _BLAS_FUNCTION(DTBSV,dtbsv)

/** BLAS DTBSV _x=op(_A)^(-1)*_x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtbsv.f">[reference]</a>
*/
extern "C"
void DTBSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _x,const ::VC::math::blas::INTEGER* _incx);

# define DTPSV _BLAS_FUNCTION(DTPSV,dtpsv)

/** BLAS DTPSV _A=_x*_x^(T)+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtpsv.f">[reference]</a>
 */
extern "C"
void DTPSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const double* _Ap,
           double* _x,const ::VC::math::blas::INTEGER* _incx);



# define DSYR _BLAS_FUNCTION(DSYR,dsyr)

/** BLAS DSYR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsyr.f">[reference]</a>
*/
extern "C"
void DSYR(const char* _uplo,
          const ::VC::math::blas::INTEGER* _n,
          const double* _alpha,
          const double* _x,const ::VC::math::blas::INTEGER* _incx,
          double* _A,const ::VC::math::blas::INTEGER* _lda);

# define DSPR _BLAS_FUNCTION(DSPR,dspr)

/** BLAS DSPR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dspr.f">[reference]</a>
*/
extern "C"
void DSPR(const char* _uplo,
          const ::VC::math::blas::INTEGER* _n,
          const double* _alpha,
          const double* _x,const ::VC::math::blas::INTEGER* _incx,
          double* _Ap);

# define DSYR2 _BLAS_FUNCTION(DSYR2,dsyr2)

/** BLAS DSYR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsyr2.f">[reference]</a>
*/
extern "C"
void DSYR2(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _y,const ::VC::math::blas::INTEGER* _incy,
           double* _A,const ::VC::math::blas::INTEGER* _lda);

# define DSPR2 _BLAS_FUNCTION(DSPR2,dspr2)

/** BLAS DSPR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dspr2.f">[reference]</a>
*/
extern "C"
void DSPR2(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _x,const ::VC::math::blas::INTEGER* _incx,
           const double* _y,const ::VC::math::blas::INTEGER* _incy,
           double* _Ap);

//
// level 3
//

# define DGEMM _BLAS_FUNCTION(DGEMM,dgemm)

/** BLAS DGEMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dgemm.f">[reference]</a>
*/
extern "C"
void DGEMM(const char* _transA,const char* _transB,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _B,const ::VC::math::blas::INTEGER* _ldb,
           const double* _beta,
           double* _C,const ::VC::math::blas::INTEGER* _ldc);

# define DSYMM _BLAS_FUNCTION(DSYMM,dsymm)

/** BLAS DSYMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsymm.f">[reference]</a>
*/
extern "C"
void DSYMM(const char* _side,const char* _uplo,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _B,const ::VC::math::blas::INTEGER* _ldb,
           const double* _beta,
           double* _C,const ::VC::math::blas::INTEGER* _ldc);

# define DSYRK _BLAS_FUNCTION(DSYRK,dsyrk)

/** BLAS DSYRK _C=_alpha*_A*_A^T+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsyrk.f">[reference]</a>
*/
extern "C"
void DSYRK(const char* _uplo,const char* _trans,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           const double* _beta,
           double* _C,const ::VC::math::blas::INTEGER* _ldc);

# define DSYR2K _BLAS_FUNCTION(DSYR2K,dsyr2k)

/** BLAS DSYR2K _C=_alpha*_A*_B^T+_alpha*_B*_A^T+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dsyr2k.f">[reference]</a>
*/
extern "C"
void DSYR2K(const char* _uplo,const char* _trans,
            const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
            const double* _alpha,
            const double* _A,const ::VC::math::blas::INTEGER* _lda,
            const double* _B,const ::VC::math::blas::INTEGER* _ldb,
            const double* _beta,
            double* _C,const ::VC::math::blas::INTEGER* _ldc);


# define DTRMM _BLAS_FUNCTION(DTRMM,dtrmm)

/** BLAS DTRMM _B=_alpha*op(_A)*op(_B).
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtrmm.f">[reference]</a>
*/
extern "C"
void DTRMM(const char* _side,const char* _uplo,const char* _transA,const char* _diag,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _B,const ::VC::math::blas::INTEGER* _ldb);


# define DTRSM _BLAS_FUNCTION(DTRSM,dtrsm)

/** BLAS DTRSM _B=_alpha*op(_A^-1)*op(_B).
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dtrsm.f">[reference]</a>
*/
extern "C"
void DTRSM(const char* _side,const char* _uplo,const char* _transA,const char* _diag,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const double* _alpha,
           const double* _A,const ::VC::math::blas::INTEGER* _lda,
           double* _B,const ::VC::math::blas::INTEGER* _ldb);

//=============================================================================
//
// single precision
//
//=============================================================================

//
// level 1
//

# define SSWAP _BLAS_FUNCTION(SSWAP,sswap)

/** BLAS swap _x and _y
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dswap.f">[reference]</a>
*/
extern "C"
void SSWAP(const ::VC::math::blas::INTEGER* _n,
           float* _x,const ::VC::math::blas::INTEGER* _incx,
           float* _y,const ::VC::math::blas::INTEGER* _incy);

# define SSCAL _BLAS_FUNCTION(SSCAL,sscal)

/** BLAS scal _x = _alpha*_x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dscal.f">[reference]</a>
*/
extern "C"
void SSCAL(const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           float* _x,const ::VC::math::blas::INTEGER* _incx);

# define SROT _BLAS_FUNCTION(SROT,srot)

/** BLAS rot
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/srot.f">[reference]</a>
 */
extern "C"
void SROT(const ::VC::math::blas::INTEGER* _n,
          float* _x,const ::VC::math::blas::INTEGER* _incx,
          float* _y,const ::VC::math::blas::INTEGER* _incy,
          const float* _c,const float* _s);

# define SROTG _BLAS_FUNCTION(SROTG,srotg)

/** BLAS rotg
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/srotg.f">[reference]</a>
 */
extern "C"
void SROTG(float* _a,float* _b,float* _c,float* _s);

# define SCOPY _BLAS_FUNCTION(SCOPY,scopy)

/** BLAS copy _y = _x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dcopy.f">[reference]</a>
*/
extern "C"
void SCOPY(const ::VC::math::blas::INTEGER* _n,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           float* _y,const ::VC::math::blas::INTEGER* _incy);

# define SAXPY _BLAS_FUNCTION(SAXPY,saxpy)

/** BLAS AXPY _y = _alpha*_x+_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/daxpy.f">[reference]</a>
*/
extern "C"
void SAXPY(const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           float* _y, const ::VC::math::blas::INTEGER* _incy);

# define SDOT _BLAS_FUNCTION(SDOT,sdot)

/** BLAS dot product _x^T * _y
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ddot.f">[reference]</a>
*/
extern "C"
float SDOT(const ::VC::math::blas::INTEGER* _n,
            const float* _x,const ::VC::math::blas::INTEGER* _incx,
            const float* _y,const ::VC::math::blas::INTEGER* _incy);

# define SNRM2 _BLAS_FUNCTION(SNRM2,snrm2)

/** BLAS 2-norm ||_x||_2
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dnrm2.f">[reference]</a>
*/
extern "C"
float SNRM2(const ::VC::math::blas::INTEGER* _n,
             const float* _x,const ::VC::math::blas::INTEGER* _incx);

# define SASUM _BLAS_FUNCTION(SASUM,sasum)

/** BLAS 1-norm || _x || _1
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/dasum.f">[reference]</a>
*/
extern "C"
float SASUM(const ::VC::math::blas::INTEGER* _n,
             const float* _x,const ::VC::math::blas::INTEGER* _incx);


//
// level 2
//

# define SGEMV _BLAS_FUNCTION(SGEMV,sgemv)

/** BLAS DGEMV _y=_alpha*op(_A)*_x+_beta*_y.
    \ingroup vc_blas_core
     <a href="http://www.netlib.org/blas/sgemv.f">[reference]</a>
*/
extern "C"
void SGEMV(const char* _trans,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _beta,
           float* _y, const ::VC::math::blas::INTEGER* _incy);

# define SSYMV _BLAS_FUNCTION(SSYMV,ssymv)

/** BLAS DSYMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssymv.f">[reference]</a>
*/
extern "C"
void SSYMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _beta,
           float* _y, const ::VC::math::blas::INTEGER* _incy);

# define SSBMV _BLAS_FUNCTION(SSBMV,ssbmv)

/** BLAS DSBMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssbmv.f">[reference]</a>
*/
extern "C"
void SSBMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const ::VC::math::blas::INTEGER* _k,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _beta,
           float* _y, const ::VC::math::blas::INTEGER* _incy);

# define SSPMV _BLAS_FUNCTION(SSPMV,sspmv)

/** BLAS DSPMV _y=_alpha*A*_x+_beta*_y.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/sspmv.f">[reference]</a>
*/
extern "C"
void SSPMV(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,const float* _Ap,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _beta,
           float* _y, const ::VC::math::blas::INTEGER* _incy);

# define STRMV _BLAS_FUNCTION(STRMV,strmv)

/** BLAS DTRMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/strmv.f">[reference]</a>
*/
extern "C"
void STRMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _x,const ::VC::math::blas::INTEGER* _incx);

# define STBMV _BLAS_FUNCTION(STBMV,stbmv)

/** BLAS DTBMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/stbmv.f">[reference]</a>
*/
extern "C"
void STBMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _x,const ::VC::math::blas::INTEGER* _incx);

# define STPMV _BLAS_FUNCTION(STPMV,stpmv)

/** BLAS DTPMV _x=op(_A)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/stpmv.f">[reference]</a>
*/
extern "C"
void STPMV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const float* _Ap,
           float* _x,const ::VC::math::blas::INTEGER* _incx);




# define STRSV _BLAS_FUNCTION(STRSV,strsv)

/** BLAS DTRSV _x=op(_A)^(-1)*_x.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/strsv.f">[reference]</a>
*/
extern "C"
void STRSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _x,const ::VC::math::blas::INTEGER* _incx);


# define STBSV _BLAS_FUNCTION(STBSV,stbsv)

/** BLAS DTBSV _x=op(_A)^(-1)*_x
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/stbsv.f">[reference]</a>
*/
extern "C"
void STBSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _x,const ::VC::math::blas::INTEGER* _incx);

# define STPSV _BLAS_FUNCTION(STPSV,stpsv)

/** BLAS DTPSV _A=_x*_x^(T)+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/stpsv.f">[reference]</a>
 */
extern "C"
void STPSV(const char* _uplo,
           const char* _trans,
           const char* _diag,
           const ::VC::math::blas::INTEGER* _n,
           const float* _Ap,
           float* _x,const ::VC::math::blas::INTEGER* _incx);



# define SSYR _BLAS_FUNCTION(SSYR,ssyr)

/** BLAS DSYR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssyr.f">[reference]</a>
*/
extern "C"
void SSYR(const char* _uplo,
          const ::VC::math::blas::INTEGER* _n,
          const float* _alpha,
          const float* _x,const ::VC::math::blas::INTEGER* _incx,
          float* _A,const ::VC::math::blas::INTEGER* _lda);

# define SSPR _BLAS_FUNCTION(SSPR,sspr)

/** BLAS DSPR _A=_alpha*_x*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/sspr.f">[reference]</a>
*/
extern "C"
void SSPR(const char* _uplo,
          const ::VC::math::blas::INTEGER* _n,
          const float* _alpha,
          const float* _x,const ::VC::math::blas::INTEGER* _incx,
          float* _Ap);

# define SSYR2 _BLAS_FUNCTION(SSYR2,ssyr2)

/** BLAS DSYR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssyr2.f">[reference]</a>
*/
extern "C"
void SSYR2(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _y,const ::VC::math::blas::INTEGER* _incy,
           float* _A,const ::VC::math::blas::INTEGER* _lda);

# define SSPR2 _BLAS_FUNCTION(SSPR2,sspr2)

/** BLAS DSPR2 _A=_alpha*_x*_y^T+_alpha*_y*_x^T+_A.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/sspr2.f">[reference]</a>
*/
extern "C"
void SSPR2(const char* _uplo,
           const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _x,const ::VC::math::blas::INTEGER* _incx,
           const float* _y,const ::VC::math::blas::INTEGER* _incy,
           float* _Ap);

//
// level 3
//

# define SGEMM _BLAS_FUNCTION(SGEMM,sgemm)

/** BLAS DGEMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/sgemm.f">[reference]</a>
*/
extern "C"
void SGEMM(const char* _transA,const char* _transB,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _B,const ::VC::math::blas::INTEGER* _ldb,
           const float* _beta,
           float* _C,const ::VC::math::blas::INTEGER* _ldc);

# define SSYMM _BLAS_FUNCTION(SSYMM,ssymm)

/** BLAS DSYMM _C=_alpha*op(_A)*op(_B)+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssymm.f">[reference]</a>
*/
extern "C"
void SSYMM(const char* _side,const char* _uplo,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _B,const ::VC::math::blas::INTEGER* _ldb,
           const float* _beta,
           float* _C,const ::VC::math::blas::INTEGER* _ldc);

# define SSYRK _BLAS_FUNCTION(SSYRK,ssyrk)

/** BLAS DSYRK _C=_alpha*_A*_A^T+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssyrk.f">[reference]</a>
*/
extern "C"
void SSYRK(const char* _uplo,const char* _trans,
           const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           const float* _beta,
           float* _C,const ::VC::math::blas::INTEGER* _ldc);

# define SSYR2K _BLAS_FUNCTION(SSYR2K,ssyr2k)

/** BLAS DSYR2K _C=_alpha*_A*_B^T+_alpha*_B*_A^T+_beta*_C.
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/ssyr2k.f">[reference]</a>
*/
extern "C"
void SSYR2K(const char* _uplo,const char* _trans,
            const ::VC::math::blas::INTEGER* _n,const ::VC::math::blas::INTEGER* _k,
            const float* _alpha,
            const float* _A,const ::VC::math::blas::INTEGER* _lda,
            const float* _B,const ::VC::math::blas::INTEGER* _ldb,
            const float* _beta,
            float* _C,const ::VC::math::blas::INTEGER* _ldc);


# define STRMM _BLAS_FUNCTION(STRMM,strmm)

/** BLAS DTRMM _B=_alpha*op(_A)*op(_B).
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/strmm.f">[reference]</a>
*/
extern "C"
void STRMM(const char* _side,const char* _uplo,const char* _transA,const char* _diag,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _B,const ::VC::math::blas::INTEGER* _ldb);


# define STRSM _BLAS_FUNCTION(STRSM,strsm)

/** BLAS DTRSM _B=_alpha*op(_A^-1)*op(_B).
    \ingroup vc_blas_core
    <a href="http://www.netlib.org/blas/strsm.f">[reference]</a>
*/
extern "C"
void STRSM(const char* _side,const char* _uplo,const char* _transA,const char* _diag,
           const ::VC::math::blas::INTEGER* _m,const ::VC::math::blas::INTEGER* _n,
           const float* _alpha,
           const float* _A,const ::VC::math::blas::INTEGER* _lda,
           float* _B,const ::VC::math::blas::INTEGER* _ldb);

//=============================================================================

#endif // VC_MATH_BLAS_PROTOTYPES_HH