#ifndef VC_MATH_LAPACK_PROTOTYPES_HH
#define VC_MATH_LAPACK_PROTOTYPES_HH

/** \file lapack_prototypes.hh Math: LAPACK prototypes
 */

#include "lapack_types.hh"

/** \defgroup vc_lapack_core LAPACK core interface
    \ingroup vc_lapack
    C/C++ interface to <a href="http://www.netlib.org/lapack/">LAPACK</a> subroutines.

    \arg <a href="http://www.netlib.org/lapack/lug/index.html">LAPACK User Guide</a>

    \arg <a href="">Search</a> for LAPACK routines by problem domain.

    \arg Overview of double precision routines
    <a href="http://www.netlib.org/lapack/double/">http://www.netlib.org/lapack/double/</a>

    \arg Overview of single precision routines
    <a href="http://www.netlib.org/lapack/single/">http://www.netlib.org/lapack/single/</a>
 */

# if defined(_USE_INTEL_FORTRAN_LIBS)
#  define _LAPACK_FUNCTION(uname,lname) uname
# endif

# if !defined(_LAPACK_FUNCTION)
#  define _LAPACK_FUNCTION(uname,lname) lname ## _
# endif

/** \def _LAPACK_FUNCTION
    \ingroup vc_lapack
    Map names of LAPACK function to  system/library dependent symbols.
    Takes upper case and lower case names (w/o underscores) as arguements
    and returns symbol as used in LAPACK library.
    For example
    \code
    // with
    #  define _LAPACK_FUNCTION(uname,lname) lname ## _
    // and
    #  define XERBLA _LAPACK_FUNCTION(XERBLA,xerbla)
    // "XERBLA" will be replaced by the symbol "xerbla_".
    \endcode
    \sa _BLAS_FUNCTION
*/
# ifdef DOXYGEN_SKIP
#  define _LAPACK_FUNCTION "map linker symbols to upper case names"
#  error "doxygen only"
# endif

# define XERBLA _LAPACK_FUNCTION(XERBLA,xerbla)



#include "lapack_prototypes_inc.hh"

extern "C" {

# define DGEQR2 _LAPACK_FUNCTION(DGEQR2,dgeqr2)


/** LAPACK DGEQR2.
Computes a QR factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgeqr2.f">http://www.netlib.org/lapack/double/dgeqr2.f</a>]
\ingroup vc_lapack_core
*/
void DGEQR2(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGEQR2 _LAPACK_FUNCTION(SGEQR2,sgeqr2)


/** LAPACK SGEQR2.
Computes a QR factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgeqr2.f">http://www.netlib.org/lapack/single/sgeqr2.f</a>]
\ingroup vc_lapack_core
*/
void SGEQR2(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);



# define DLASR _LAPACK_FUNCTION(DLASR,dlasr)


/** LAPACK DLASR.
DLASR applies a sequence of plane rotations to a real matrix A, from either the left or the right.
[<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlasr.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlasr.f</a>]
\ingroup vc_lapack_core
*/
void DLASR(const char* _side,
           const char* _pivot,
           const char* _direct,
           ::VC::math::lapack::INTEGER* _m,
           ::VC::math::lapack::INTEGER* _n,
           const double* _c,
           const double* _s,
           double* _a,
           ::VC::math::lapack::INTEGER* _lda);

# define SLASR _LAPACK_FUNCTION(SLASR,slasr)


/** LAPACK SLASR.
SLASR applies a sequence of plane rotations to a real matrix A, from either the left or the right.
[<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/slasr.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/slasr.f</a>]
\ingroup vc_lapack_core
*/
void SLASR(const char* _side,
           const char* _pivot,
           const char* _direct,
           ::VC::math::lapack::INTEGER* _m,
           ::VC::math::lapack::INTEGER* _n,
           const float* _c,
           const float* _s,
           float* _a,
           ::VC::math::lapack::INTEGER* _lda);

# define DLARF _LAPACK_FUNCTION(DLARF,dlarf)


/** LAPACK DLARF.
    DLARF applies a real elementary reflector H to a real m by n matrix C,
    from either the left or the right. H is represented in the form
    \code
        H = I - tau * v * v'
    \endcode
    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.

    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarf.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarf.f</a>]
\ingroup vc_lapack_core
*/
void DLARF(const char* _side,
           ::VC::math::lapack::INTEGER* _m,
           ::VC::math::lapack::INTEGER* _n,
           const double* _v,
           ::VC::math::lapack::INTEGER* _incv,
           const double* _tau,
           double* _c,
           ::VC::math::lapack::INTEGER* _ldc,
           double* _work);


# define SLARF _LAPACK_FUNCTION(SLARF,slarf)


/** LAPACK SLARF.
    SLARF applies a real elementary reflector H to a real m by n matrix C,
    from either the left or the right. H is represented in the form
    \code
        H = I - tau * v * v'
    \endcode
    where tau is a real scalar and v is a real vector.

    If tau = 0, then H is taken to be the unit matrix.

    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarf.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarf.f</a>]
\ingroup vc_lapack_core
*/
void SLARF(const char* _side,
           ::VC::math::lapack::INTEGER* _m,
           ::VC::math::lapack::INTEGER* _n,
           const float* _v,
           ::VC::math::lapack::INTEGER* _incv,
           const float* _tau,
           float* _c,
           ::VC::math::lapack::INTEGER* _ldc,
           float* _work);


# define DLARFG _LAPACK_FUNCTION(DLARFG,dlarfg)


/** LAPACK DLARFG.
    DLARFG generates a real elementary reflector H of order n, such that
    \code
         H * ( alpha ) = ( beta ),   H' * H = I.
             (   x   )   (   0  )
    \endcode
    where alpha and beta are scalars, and x is an (n-1)-element real
    vector. H is represented in the form
    \code
        H = I - tau * ( 1 ) * ( 1 v' ) ,
                      ( v )
    \endcode
    where tau is a real scalar and v is a real (n-1)-element
    vector.

    If the elements of x are all zero, then tau = 0 and H is taken to be
    the unit matrix.

    Otherwise  1 <= tau <= 2.

    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfg.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfg.f</a>]
\ingroup vc_lapack_core
*/
void DLARFG(::VC::math::lapack::INTEGER* _n,
            double* _alpha,
            double* _x,
            ::VC::math::lapack::INTEGER* _incx,
            double* _tau);

# define SLARFG _LAPACK_FUNCTION(SLARFG,slarfg)


/** LAPACK SLARFG.
    SLARFG generates a real elementary reflector H of order n, such that
    \code
         H * ( alpha ) = ( beta ),   H' * H = I.
             (   x   )   (   0  )
    \endcode
    where alpha and beta are scalars, and x is an (n-1)-element real
    vector. H is represented in the form
    \code
        H = I - tau * ( 1 ) * ( 1 v' ) ,
                      ( v )
    \endcode
    where tau is a real scalar and v is a real (n-1)-element
    vector.

    If the elements of x are all zero, then tau = 0 and H is taken to be
    the unit matrix.

    Otherwise  1 <= tau <= 2.

    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfg.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfg.f</a>]
\ingroup vc_lapack_core
*/
void SLARFG(::VC::math::lapack::INTEGER* _n,
            float* _alpha,
            float* _x,
            ::VC::math::lapack::INTEGER* _incx,
            float* _tau);


# define DLARFT _LAPACK_FUNCTION(DLARFT,dlarft)

/** LAPACK DLARFT.
    DLARFT applies a real block reflector H or its transpose H' to a
    real m by n matrix C, from either the left or the right.
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarft.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarft.f</a>]
\ingroup vc_lapack_core
 */
void DLARFT(const char* _direct,
            const char* _strorev,
            const ::VC::math::lapack::INTEGER* _n,
            const ::VC::math::lapack::INTEGER* _k,
            double* _v,
            const ::VC::math::lapack::INTEGER* _ldv,
            const double* _tau,
            double* _t,
            const ::VC::math::lapack::INTEGER* _ldt);


# define SLARFT _LAPACK_FUNCTION(SLARFT,slarft)

/** LAPACK SLARFT.
    SLARFT applies a real block reflector H or its transpose H' to a
    real m by n matrix C, from either the left or the right.
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarft.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarft.f</a>]
\ingroup vc_lapack_core
 */
void SLARFT(const char* _direct,
            const char* _strorev,
            const ::VC::math::lapack::INTEGER* _n,
            const ::VC::math::lapack::INTEGER* _k,
            float* _v,
            const ::VC::math::lapack::INTEGER* _ldv,
            const float* _tau,
            float* _t,
            const ::VC::math::lapack::INTEGER* _ldt);


# define DLARFB _LAPACK_FUNCTION(DLARFB,dlarfb)

/** LAPACK DLARFB.
    DLARFB applies a real block reflector H or its transpose H' to a
    real m by n matrix C, from either the left or the right.
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfb.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlarfb.f</a>]
\ingroup vc_lapack_core
 */
void DLARFB(const char* _side,
            const char* _trans,
            const char* _direct,
            const char* _strorev,
            const ::VC::math::lapack::INTEGER* _m,
            const ::VC::math::lapack::INTEGER* _n,
            const ::VC::math::lapack::INTEGER* _k,
            const double* _v,
            const ::VC::math::lapack::INTEGER* _ldv,
            const double* _t,
            const ::VC::math::lapack::INTEGER* _ldt,
            double* _c,
            const ::VC::math::lapack::INTEGER* _ldc,
            double* _work,
            const ::VC::math::lapack::INTEGER* _ldwork);


# define SLARFB _LAPACK_FUNCTION(SLARFB,slarfb)

/** LAPACK SLARFB.
    SLARFB applies a real block reflector H or its transpose H' to a
    real m by n matrix C, from either the left or the right.
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarfb.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/slarfb.f</a>]
\ingroup vc_lapack_core
 */
void SLARFB(const char* _side,
            const char* _trans,
            const char* _direct,
            const char* _strorev,
            const ::VC::math::lapack::INTEGER* _m,
            const ::VC::math::lapack::INTEGER* _n,
            const ::VC::math::lapack::INTEGER* _k,
            const float* _v,
            const ::VC::math::lapack::INTEGER* _ldv,
            const float* _t,
            const ::VC::math::lapack::INTEGER* _ldt,
            float* _c,
            const ::VC::math::lapack::INTEGER* _ldc,
            float* _work,
            const ::VC::math::lapack::INTEGER* _ldwork);


# define DLANGE _LAPACK_FUNCTION(DLANGE,dlange)

/** LAPACK DLANGE.
    DLANGE - return the value of the one norm, or the Frobenius
    norm, or the infinity norm, or the element of largest absolute
    value of a real matrix A
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlange.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/dlange.f</a>]
\ingroup vc_lapack_core
 */
double DLANGE(const char* _norm,
              ::VC::math::lapack::INTEGER* _m,
              ::VC::math::lapack::INTEGER* _n,
              const double* _a,::VC::math::lapack::INTEGER* _lda,
              double* _work);


# define SLANGE _LAPACK_FUNCTION(SLANGE,slange)

/** LAPACK SLANGE.
    SLANGE - return the value of the one norm, or the Frobenius
    norm, or the infinity norm, or the element of largest absolute
    value of a real matrix A
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/slange.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/slange.f</a>]
\ingroup vc_lapack_core
 */
float SLANGE(const char* _norm,
              ::VC::math::lapack::INTEGER* _m,
              ::VC::math::lapack::INTEGER* _n,
              const float* _a,::VC::math::lapack::INTEGER* _lda,
             float* _work);

#define ILAENV _LAPACK_FUNCTION(ILAENV,ilaenv)

/** LAPACK ILAENV
    ILAENV is called from the LAPACK routines to choose problem-dependent
    parameters for the local environment.
    [<a href="http://www.netlib.org/lapack/lapack-3.1.1/SRC/ilaenv.f">http://www.netlib.org/lapack/lapack-3.1.1/SRC/ilaenv.f</a>]
\ingroup vc_lapack_core
 */
::VC::math::lapack::INTEGER
ILAENV(const ::VC::math::lapack::INTEGER* _spec,
       const char* _name,
       const char* _opts,
       const ::VC::math::lapack::INTEGER* _n1,
       const ::VC::math::lapack::INTEGER* _n2,
       const ::VC::math::lapack::INTEGER* _n3,
       const ::VC::math::lapack::INTEGER* _n4,
       ::VC::math::lapack::INTEGER _len_name, // THIS IS POTENTIALLY
       ::VC::math::lapack::INTEGER _len_opts  // ***NOT PORTABLE*** (works for gfortan)
       );

} // extern "C"

#endif // VC_MATH_LAPACK_PROTOTYPES_HH
