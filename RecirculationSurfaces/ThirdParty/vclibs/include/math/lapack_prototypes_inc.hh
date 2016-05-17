extern "C" {

# define DPTRFS _LAPACK_FUNCTION(DPTRFS,dptrfs)


/** LAPACK DPTRFS.
Improves the computed solution to a symmetric positive definite tridiagonal system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dptrfs.f">http://www.netlib.org/lapack/double/dptrfs.f</a>]
\ingroup vc_lapack_core
*/
void DPTRFS(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _d,
        const double* _e,
        const double* _df,
        const double* _ef,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SPTRFS _LAPACK_FUNCTION(SPTRFS,sptrfs)


/** LAPACK SPTRFS.
Improves the computed solution to a symmetric positive definite tridiagonal system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/sptrfs.f">http://www.netlib.org/lapack/single/sptrfs.f</a>]
\ingroup vc_lapack_core
*/
void SPTRFS(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _d,
        const float* _e,
        const float* _df,
        const float* _ef,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSYGVD _LAPACK_FUNCTION(DSYGVD,dsygvd)


/** LAPACK DSYGVD.
Computes all eigenvalues and the eigenvectors of  a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x. If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dsygvd.f">http://www.netlib.org/lapack/double/dsygvd.f</a>]
\ingroup vc_lapack_core
*/
void DSYGVD(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _w,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYGVD _LAPACK_FUNCTION(SSYGVD,ssygvd)


/** LAPACK SSYGVD.
Computes all eigenvalues and the eigenvectors of  a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x. If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/ssygvd.f">http://www.netlib.org/lapack/single/ssygvd.f</a>]
\ingroup vc_lapack_core
*/
void SSYGVD(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _w,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEBRD _LAPACK_FUNCTION(DGEBRD,dgebrd)


/** LAPACK DGEBRD.
Reduces a general rectangular matrix to real bidiagonal form by an orthogonal transformation.
[<a href="http://www.netlib.org/lapack/double/dgebrd.f">http://www.netlib.org/lapack/double/dgebrd.f</a>]
\ingroup vc_lapack_core
*/
void DGEBRD(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _d,
        double* _e,
        double* _tauq,
        double* _taup,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEBRD _LAPACK_FUNCTION(SGEBRD,sgebrd)


/** LAPACK SGEBRD.
Reduces a general rectangular matrix to real bidiagonal form by an orthogonal transformation.
[<a href="http://www.netlib.org/lapack/single/sgebrd.f">http://www.netlib.org/lapack/single/sgebrd.f</a>]
\ingroup vc_lapack_core
*/
void SGEBRD(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _d,
        float* _e,
        float* _tauq,
        float* _taup,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORGBR _LAPACK_FUNCTION(DORGBR,dorgbr)


/** LAPACK DORGBR.
Generates the orthogonal transformation matrices from a reduction to bidiagonal form determined by DGEBRD.
[<a href="http://www.netlib.org/lapack/double/dorgbr.f">http://www.netlib.org/lapack/double/dorgbr.f</a>]
\ingroup vc_lapack_core
*/
void DORGBR(
        const char* _vect,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGBR _LAPACK_FUNCTION(SORGBR,sorgbr)


/** LAPACK SORGBR.
Generates the orthogonal transformation matrices from a reduction to bidiagonal form determined by DGEBRD.
[<a href="http://www.netlib.org/lapack/single/sorgbr.f">http://www.netlib.org/lapack/single/sorgbr.f</a>]
\ingroup vc_lapack_core
*/
void SORGBR(
        const char* _vect,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEVD _LAPACK_FUNCTION(DSTEVD,dstevd)


/** LAPACK DSTEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dstevd.f">http://www.netlib.org/lapack/double/dstevd.f</a>]
\ingroup vc_lapack_core
*/
void DSTEVD(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEVD _LAPACK_FUNCTION(SSTEVD,sstevd)


/** LAPACK SSTEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/sstevd.f">http://www.netlib.org/lapack/single/sstevd.f</a>]
\ingroup vc_lapack_core
*/
void SSTEVD(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPPEQU _LAPACK_FUNCTION(DPPEQU,dppequ)


/** LAPACK DPPEQU.
Computes row and column scalings to equilibrate a symmetric positive definite matrix in packed storage and reduce its condition number.
[<a href="http://www.netlib.org/lapack/double/dppequ.f">http://www.netlib.org/lapack/double/dppequ.f</a>]
\ingroup vc_lapack_core
*/
void DPPEQU(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        double* _s,
        double* _scond,
        double* _amax,
        ::VC::math::lapack::INTEGER* _info);

# define SPPEQU _LAPACK_FUNCTION(SPPEQU,sppequ)


/** LAPACK SPPEQU.
Computes row and column scalings to equilibrate a symmetric positive definite matrix in packed storage and reduce its condition number.
[<a href="http://www.netlib.org/lapack/single/sppequ.f">http://www.netlib.org/lapack/single/sppequ.f</a>]
\ingroup vc_lapack_core
*/
void SPPEQU(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        float* _s,
        float* _scond,
        float* _amax,
        ::VC::math::lapack::INTEGER* _info);


# define DPPCON _LAPACK_FUNCTION(DPPCON,dppcon)


/** LAPACK DPPCON.
Estimates the reciprocal of the condition number of a symmetric positive definite matrix in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/double/dppcon.f">http://www.netlib.org/lapack/double/dppcon.f</a>]
\ingroup vc_lapack_core
*/
void DPPCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPPCON _LAPACK_FUNCTION(SPPCON,sppcon)


/** LAPACK SPPCON.
Estimates the reciprocal of the condition number of a symmetric positive definite matrix in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/single/sppcon.f">http://www.netlib.org/lapack/single/sppcon.f</a>]
\ingroup vc_lapack_core
*/
void SPPCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPOTRF _LAPACK_FUNCTION(DPOTRF,dpotrf)


/** LAPACK DPOTRF.
Computes the Cholesky factorization of a symmetric positive definite matrix.
[<a href="http://www.netlib.org/lapack/double/dpotrf.f">http://www.netlib.org/lapack/double/dpotrf.f</a>]
\ingroup vc_lapack_core
*/
void DPOTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);

# define SPOTRF _LAPACK_FUNCTION(SPOTRF,spotrf)


/** LAPACK SPOTRF.
Computes the Cholesky factorization of a symmetric positive definite matrix.
[<a href="http://www.netlib.org/lapack/single/spotrf.f">http://www.netlib.org/lapack/single/spotrf.f</a>]
\ingroup vc_lapack_core
*/
void SPOTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);


# define DSYTRS _LAPACK_FUNCTION(DSYTRS,dsytrs)


/** LAPACK DSYTRS.
Solves a real symmetric indefinite system of linear equations AX=B, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/double/dsytrs.f">http://www.netlib.org/lapack/double/dsytrs.f</a>]
\ingroup vc_lapack_core
*/
void DSYTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SSYTRS _LAPACK_FUNCTION(SSYTRS,ssytrs)


/** LAPACK SSYTRS.
Solves a real symmetric indefinite system of linear equations AX=B, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/single/ssytrs.f">http://www.netlib.org/lapack/single/ssytrs.f</a>]
\ingroup vc_lapack_core
*/
void SSYTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGEQP3 dgeqp3_

/** LAPACK DGEQP3.
Computes a QR factorization with column pivoting of a general rectangular matrix using Level 3 BLAS.
[<a href="http://www.netlib.org/lapack/double/dgeqp3.f">http://www.netlib.org/lapack/double/dgeqp3.f</a>]
\ingroup vc_lapack_core
*/
void DGEQP3(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _jpvt,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEQP3 sgeqp3_

/** LAPACK SGEQP3.
Computes a QR factorization with column pivoting of a general rectangular matrix using Level 3 BLAS.
[<a href="http://www.netlib.org/lapack/single/sgeqp3.f">http://www.netlib.org/lapack/single/sgeqp3.f</a>]
\ingroup vc_lapack_core
*/
void SGEQP3(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _jpvt,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGSVP _LAPACK_FUNCTION(DGGSVP,dggsvp)


/** LAPACK DGGSVP.

[<a href="http://www.netlib.org/lapack/double/dggsvp.f">http://www.netlib.org/lapack/double/dggsvp.f</a>]
\ingroup vc_lapack_core
*/
void DGGSVP(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _tola,
        const double* _tolb,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        ::VC::math::lapack::INTEGER* _iwork,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGGSVP _LAPACK_FUNCTION(SGGSVP,sggsvp)


/** LAPACK SGGSVP.

[<a href="http://www.netlib.org/lapack/single/sggsvp.f">http://www.netlib.org/lapack/single/sggsvp.f</a>]
\ingroup vc_lapack_core
*/
void SGGSVP(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _tola,
        const float* _tolb,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        ::VC::math::lapack::INTEGER* _iwork,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGELS _LAPACK_FUNCTION(DGELS,dgels)


/** LAPACK DGELS.
Computes the least squares solution to an over-determined system of linear equations, A X=B or A**H X=B,  or the minimum norm solution of an under-determined system, where A is a general rectangular matrix of full rank,  using a QR or LQ factorization of A.
[<a href="http://www.netlib.org/lapack/double/dgels.f">http://www.netlib.org/lapack/double/dgels.f</a>]
\ingroup vc_lapack_core
*/
void DGELS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGELS _LAPACK_FUNCTION(SGELS,sgels)


/** LAPACK SGELS.
Computes the least squares solution to an over-determined system of linear equations, A X=B or A**H X=B,  or the minimum norm solution of an under-determined system, where A is a general rectangular matrix of full rank,  using a QR or LQ factorization of A.
[<a href="http://www.netlib.org/lapack/single/sgels.f">http://www.netlib.org/lapack/single/sgels.f</a>]
\ingroup vc_lapack_core
*/
void SGELS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSYSVX _LAPACK_FUNCTION(DSYSVX,dsysvx)


/** LAPACK DSYSVX.
Solves a real symmetric indefinite system  of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dsysvx.f">http://www.netlib.org/lapack/double/dsysvx.f</a>]
\ingroup vc_lapack_core
*/
void DSYSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYSVX _LAPACK_FUNCTION(SSYSVX,ssysvx)


/** LAPACK SSYSVX.
Solves a real symmetric indefinite system  of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/ssysvx.f">http://www.netlib.org/lapack/single/ssysvx.f</a>]
\ingroup vc_lapack_core
*/
void SSYSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPGVD _LAPACK_FUNCTION(DSPGVD,dspgvd)


/** LAPACK DSPGVD.
Computes all eigenvalues and eigenvectors of  a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage. If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dspgvd.f">http://www.netlib.org/lapack/double/dspgvd.f</a>]
\ingroup vc_lapack_core
*/
void DSPGVD(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _bp,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSPGVD _LAPACK_FUNCTION(SSPGVD,sspgvd)


/** LAPACK SSPGVD.
Computes all eigenvalues and eigenvectors of  a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage. If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/sspgvd.f">http://www.netlib.org/lapack/single/sspgvd.f</a>]
\ingroup vc_lapack_core
*/
void SSPGVD(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _bp,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPPTRF _LAPACK_FUNCTION(DPPTRF,dpptrf)


/** LAPACK DPPTRF.
Computes the Cholesky factorization of a symmetric positive definite matrix in packed storage.
[<a href="http://www.netlib.org/lapack/double/dpptrf.f">http://www.netlib.org/lapack/double/dpptrf.f</a>]
\ingroup vc_lapack_core
*/
void DPPTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        ::VC::math::lapack::INTEGER* _info);

# define SPPTRF _LAPACK_FUNCTION(SPPTRF,spptrf)


/** LAPACK SPPTRF.
Computes the Cholesky factorization of a symmetric positive definite matrix in packed storage.
[<a href="http://www.netlib.org/lapack/single/spptrf.f">http://www.netlib.org/lapack/single/spptrf.f</a>]
\ingroup vc_lapack_core
*/
void SPPTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        ::VC::math::lapack::INTEGER* _info);


# define DPORFS _LAPACK_FUNCTION(DPORFS,dporfs)


/** LAPACK DPORFS.
Improves the computed solution to a symmetric positive definite system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dporfs.f">http://www.netlib.org/lapack/double/dporfs.f</a>]
\ingroup vc_lapack_core
*/
void DPORFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPORFS _LAPACK_FUNCTION(SPORFS,sporfs)


/** LAPACK SPORFS.
Improves the computed solution to a symmetric positive definite system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/sporfs.f">http://www.netlib.org/lapack/single/sporfs.f</a>]
\ingroup vc_lapack_core
*/
void SPORFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORMRZ _LAPACK_FUNCTION(DORMRZ,dormrz)


/** LAPACK DORMRZ.

[<a href="http://www.netlib.org/lapack/double/dormrz.f">http://www.netlib.org/lapack/double/dormrz.f</a>]
\ingroup vc_lapack_core
*/
void DORMRZ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMRZ _LAPACK_FUNCTION(SORMRZ,sormrz)


/** LAPACK SORMRZ.

[<a href="http://www.netlib.org/lapack/single/sormrz.f">http://www.netlib.org/lapack/single/sormrz.f</a>]
\ingroup vc_lapack_core
*/
void SORMRZ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEQLF _LAPACK_FUNCTION(DGEQLF,dgeqlf)


/** LAPACK DGEQLF.
Computes a QL factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgeqlf.f">http://www.netlib.org/lapack/double/dgeqlf.f</a>]
\ingroup vc_lapack_core
*/
void DGEQLF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEQLF _LAPACK_FUNCTION(SGEQLF,sgeqlf)


/** LAPACK SGEQLF.
Computes a QL factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgeqlf.f">http://www.netlib.org/lapack/single/sgeqlf.f</a>]
\ingroup vc_lapack_core
*/
void SGEQLF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGTSVX _LAPACK_FUNCTION(DGTSVX,dgtsvx)


/** LAPACK DGTSVX.
Solves a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number  and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dgtsvx.f">http://www.netlib.org/lapack/double/dgtsvx.f</a>]
\ingroup vc_lapack_core
*/
void DGTSVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _dl,
        const double* _d,
        const double* _du,
        double* _dlf,
        double* _df,
        double* _duf,
        double* _du2,
        ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGTSVX _LAPACK_FUNCTION(SGTSVX,sgtsvx)


/** LAPACK SGTSVX.
Solves a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number  and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sgtsvx.f">http://www.netlib.org/lapack/single/sgtsvx.f</a>]
\ingroup vc_lapack_core
*/
void SGTSVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _dl,
        const float* _d,
        const float* _du,
        float* _dlf,
        float* _df,
        float* _duf,
        float* _du2,
        ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGTTRS _LAPACK_FUNCTION(DGTTRS,dgttrs)


/** LAPACK DGTTRS.
Solves a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGTTRF.
[<a href="http://www.netlib.org/lapack/double/dgttrs.f">http://www.netlib.org/lapack/double/dgttrs.f</a>]
\ingroup vc_lapack_core
*/
void DGTTRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _dl,
        const double* _d,
        const double* _du,
        const double* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGTTRS _LAPACK_FUNCTION(SGTTRS,sgttrs)


/** LAPACK SGTTRS.
Solves a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGTTRF.
[<a href="http://www.netlib.org/lapack/single/sgttrs.f">http://www.netlib.org/lapack/single/sgttrs.f</a>]
\ingroup vc_lapack_core
*/
void SGTTRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _dl,
        const float* _d,
        const float* _du,
        const float* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DORGQL _LAPACK_FUNCTION(DORGQL,dorgql)


/** LAPACK DORGQL.
Generates all or part of the orthogonal matrix Q from a QL factorization determined by DGEQLF.
[<a href="http://www.netlib.org/lapack/double/dorgql.f">http://www.netlib.org/lapack/double/dorgql.f</a>]
\ingroup vc_lapack_core
*/
void DORGQL(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGQL _LAPACK_FUNCTION(SORGQL,sorgql)


/** LAPACK SORGQL.
Generates all or part of the orthogonal matrix Q from a QL factorization determined by DGEQLF.
[<a href="http://www.netlib.org/lapack/single/sorgql.f">http://www.netlib.org/lapack/single/sorgql.f</a>]
\ingroup vc_lapack_core
*/
void SORGQL(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPTRS _LAPACK_FUNCTION(DSPTRS,dsptrs)


/** LAPACK DSPTRS.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/double/dsptrs.f">http://www.netlib.org/lapack/double/dsptrs.f</a>]
\ingroup vc_lapack_core
*/
void DSPTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SSPTRS _LAPACK_FUNCTION(SSPTRS,ssptrs)


/** LAPACK SSPTRS.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/single/ssptrs.f">http://www.netlib.org/lapack/single/ssptrs.f</a>]
\ingroup vc_lapack_core
*/
void SSPTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DPOTRI _LAPACK_FUNCTION(DPOTRI,dpotri)


/** LAPACK DPOTRI.
Computes the inverse of a symmetric positive definite matrix, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/double/dpotri.f">http://www.netlib.org/lapack/double/dpotri.f</a>]
\ingroup vc_lapack_core
*/
void DPOTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);

# define SPOTRI _LAPACK_FUNCTION(SPOTRI,spotri)


/** LAPACK SPOTRI.
Computes the inverse of a symmetric positive definite matrix, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/single/spotri.f">http://www.netlib.org/lapack/single/spotri.f</a>]
\ingroup vc_lapack_core
*/
void SPOTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);


# define DORMR3 dormr3_

/** LAPACK DORMR3.

[<a href="http://www.netlib.org/lapack/double/dormr3.f">http://www.netlib.org/lapack/double/dormr3.f</a>]
\ingroup vc_lapack_core
*/
void DORMR3(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SORMR3 sormr3_

/** LAPACK SORMR3.

[<a href="http://www.netlib.org/lapack/single/sormr3.f">http://www.netlib.org/lapack/single/sormr3.f</a>]
\ingroup vc_lapack_core
*/
void SORMR3(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGELQF _LAPACK_FUNCTION(DGELQF,dgelqf)


/** LAPACK DGELQF.
Computes an LQ factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgelqf.f">http://www.netlib.org/lapack/double/dgelqf.f</a>]
\ingroup vc_lapack_core
*/
void DGELQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGELQF _LAPACK_FUNCTION(SGELQF,sgelqf)


/** LAPACK SGELQF.
Computes an LQ factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgelqf.f">http://www.netlib.org/lapack/single/sgelqf.f</a>]
\ingroup vc_lapack_core
*/
void SGELQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBEVD _LAPACK_FUNCTION(DSBEVD,dsbevd)


/** LAPACK DSBEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric band matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dsbevd.f">http://www.netlib.org/lapack/double/dsbevd.f</a>]
\ingroup vc_lapack_core
*/
void DSBEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSBEVD _LAPACK_FUNCTION(SSBEVD,ssbevd)


/** LAPACK SSBEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric band matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/ssbevd.f">http://www.netlib.org/lapack/single/ssbevd.f</a>]
\ingroup vc_lapack_core
*/
void SSBEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORGHR _LAPACK_FUNCTION(DORGHR,dorghr)


/** LAPACK DORGHR.
Generates the orthogonal transformation matrix from a reduction to Hessenberg form determined by DGEHRD.
[<a href="http://www.netlib.org/lapack/double/dorghr.f">http://www.netlib.org/lapack/double/dorghr.f</a>]
\ingroup vc_lapack_core
*/
void DORGHR(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGHR _LAPACK_FUNCTION(SORGHR,sorghr)


/** LAPACK SORGHR.
Generates the orthogonal transformation matrix from a reduction to Hessenberg form determined by DGEHRD.
[<a href="http://www.netlib.org/lapack/single/sorghr.f">http://www.netlib.org/lapack/single/sorghr.f</a>]
\ingroup vc_lapack_core
*/
void SORGHR(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPSVX _LAPACK_FUNCTION(DSPSVX,dspsvx)


/** LAPACK DSPSVX.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dspsvx.f">http://www.netlib.org/lapack/double/dspsvx.f</a>]
\ingroup vc_lapack_core
*/
void DSPSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        double* _afp,
        ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSPSVX _LAPACK_FUNCTION(SSPSVX,sspsvx)


/** LAPACK SSPSVX.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sspsvx.f">http://www.netlib.org/lapack/single/sspsvx.f</a>]
\ingroup vc_lapack_core
*/
void SSPSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        float* _afp,
        ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPGV _LAPACK_FUNCTION(DSPGV,dspgv)


/** LAPACK DSPGV.
Computes all eigenvalues and eigenvectors of  a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage.
[<a href="http://www.netlib.org/lapack/double/dspgv.f">http://www.netlib.org/lapack/double/dspgv.f</a>]
\ingroup vc_lapack_core
*/
void DSPGV(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _bp,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSPGV _LAPACK_FUNCTION(SSPGV,sspgv)


/** LAPACK SSPGV.
Computes all eigenvalues and eigenvectors of  a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage.
[<a href="http://www.netlib.org/lapack/single/sspgv.f">http://www.netlib.org/lapack/single/sspgv.f</a>]
\ingroup vc_lapack_core
*/
void SSPGV(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _bp,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEGR _LAPACK_FUNCTION(DSTEGR,dstegr)


/** LAPACK DSTEGR.

[<a href="http://www.netlib.org/lapack/double/dstegr.f">http://www.netlib.org/lapack/double/dstegr.f</a>]
\ingroup vc_lapack_core
*/
void DSTEGR(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEGR _LAPACK_FUNCTION(SSTEGR,sstegr)


/** LAPACK SSTEGR.

[<a href="http://www.netlib.org/lapack/single/sstegr.f">http://www.netlib.org/lapack/single/sstegr.f</a>]
\ingroup vc_lapack_core
*/
void SSTEGR(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPPRFS _LAPACK_FUNCTION(DPPRFS,dpprfs)


/** LAPACK DPPRFS.
Improves the computed solution to a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dpprfs.f">http://www.netlib.org/lapack/double/dpprfs.f</a>]
\ingroup vc_lapack_core
*/
void DPPRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        const double* _afp,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPPRFS _LAPACK_FUNCTION(SPPRFS,spprfs)


/** LAPACK SPPRFS.
Improves the computed solution to a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/spprfs.f">http://www.netlib.org/lapack/single/spprfs.f</a>]
\ingroup vc_lapack_core
*/
void SPPRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        const float* _afp,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPBEQU _LAPACK_FUNCTION(DPBEQU,dpbequ)


/** LAPACK DPBEQU.
Computes row and column scalings to equilibrate a symmetric positive definite band matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/double/dpbequ.f">http://www.netlib.org/lapack/double/dpbequ.f</a>]
\ingroup vc_lapack_core
*/
void DPBEQU(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _s,
        double* _scond,
        double* _amax,
        ::VC::math::lapack::INTEGER* _info);

# define SPBEQU _LAPACK_FUNCTION(SPBEQU,spbequ)


/** LAPACK SPBEQU.
Computes row and column scalings to equilibrate a symmetric positive definite band matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/single/spbequ.f">http://www.netlib.org/lapack/single/spbequ.f</a>]
\ingroup vc_lapack_core
*/
void SPBEQU(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _s,
        float* _scond,
        float* _amax,
        ::VC::math::lapack::INTEGER* _info);


# define DPPTRI _LAPACK_FUNCTION(DPPTRI,dpptri)


/** LAPACK DPPTRI.
Computes the inverse of a symmetric positive definite matrix in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/double/dpptri.f">http://www.netlib.org/lapack/double/dpptri.f</a>]
\ingroup vc_lapack_core
*/
void DPPTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        ::VC::math::lapack::INTEGER* _info);

# define SPPTRI _LAPACK_FUNCTION(SPPTRI,spptri)


/** LAPACK SPPTRI.
Computes the inverse of a symmetric positive definite matrix in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/single/spptri.f">http://www.netlib.org/lapack/single/spptri.f</a>]
\ingroup vc_lapack_core
*/
void SPPTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        ::VC::math::lapack::INTEGER* _info);


# define DPBCON _LAPACK_FUNCTION(DPBCON,dpbcon)


/** LAPACK DPBCON.
Estimates the reciprocal of the condition number of a symmetric positive definite band matrix, using the Cholesky factorization computed by DPBTRF.
[<a href="http://www.netlib.org/lapack/double/dpbcon.f">http://www.netlib.org/lapack/double/dpbcon.f</a>]
\ingroup vc_lapack_core
*/
void DPBCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPBCON _LAPACK_FUNCTION(SPBCON,spbcon)


/** LAPACK SPBCON.
Estimates the reciprocal of the condition number of a symmetric positive definite band matrix, using the Cholesky factorization computed by DPBTRF.
[<a href="http://www.netlib.org/lapack/single/spbcon.f">http://www.netlib.org/lapack/single/spbcon.f</a>]
\ingroup vc_lapack_core
*/
void SPBCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGBSVX _LAPACK_FUNCTION(DGBSVX,dgbsvx)


/** LAPACK DGBSVX.
Solves a general banded system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dgbsvx.f">http://www.netlib.org/lapack/double/dgbsvx.f</a>]
\ingroup vc_lapack_core
*/
void DGBSVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        ::VC::math::lapack::INTEGER* _ipiv,
        char* _equed,
        double* _r,
        double* _c,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGBSVX _LAPACK_FUNCTION(SGBSVX,sgbsvx)


/** LAPACK SGBSVX.
Solves a general banded system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sgbsvx.f">http://www.netlib.org/lapack/single/sgbsvx.f</a>]
\ingroup vc_lapack_core
*/
void SGBSVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        ::VC::math::lapack::INTEGER* _ipiv,
        char* _equed,
        float* _r,
        float* _c,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGESVX _LAPACK_FUNCTION(DGESVX,dgesvx)


/** LAPACK DGESVX.
Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dgesvx.f">http://www.netlib.org/lapack/double/dgesvx.f</a>]
\ingroup vc_lapack_core
*/
void DGESVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        ::VC::math::lapack::INTEGER* _ipiv,
        char* _equed,
        double* _r,
        double* _c,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGESVX _LAPACK_FUNCTION(SGESVX,sgesvx)


/** LAPACK SGESVX.
Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sgesvx.f">http://www.netlib.org/lapack/single/sgesvx.f</a>]
\ingroup vc_lapack_core
*/
void SGESVX(
        const char* _fact,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        ::VC::math::lapack::INTEGER* _ipiv,
        char* _equed,
        float* _r,
        float* _c,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGBTRS _LAPACK_FUNCTION(DGBTRS,dgbtrs)


/** LAPACK DGBTRS.
Solves a general banded system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGBTRF.
[<a href="http://www.netlib.org/lapack/double/dgbtrs.f">http://www.netlib.org/lapack/double/dgbtrs.f</a>]
\ingroup vc_lapack_core
*/
void DGBTRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGBTRS _LAPACK_FUNCTION(SGBTRS,sgbtrs)


/** LAPACK SGBTRS.
Solves a general banded system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGBTRF.
[<a href="http://www.netlib.org/lapack/single/sgbtrs.f">http://www.netlib.org/lapack/single/sgbtrs.f</a>]
\ingroup vc_lapack_core
*/
void SGBTRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGETRS _LAPACK_FUNCTION(DGETRS,dgetrs)


/** LAPACK DGETRS.
Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/double/dgetrs.f">http://www.netlib.org/lapack/double/dgetrs.f</a>]
\ingroup vc_lapack_core
*/
void DGETRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGETRS _LAPACK_FUNCTION(SGETRS,sgetrs)


/** LAPACK SGETRS.
Solves a general system of linear equations AX=B, A**T X=B or A**H X=B, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/single/sgetrs.f">http://www.netlib.org/lapack/single/sgetrs.f</a>]
\ingroup vc_lapack_core
*/
void SGETRS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DPBTRF _LAPACK_FUNCTION(DPBTRF,dpbtrf)


/** LAPACK DPBTRF.
Computes the Cholesky factorization of a symmetric positive definite band matrix.
[<a href="http://www.netlib.org/lapack/double/dpbtrf.f">http://www.netlib.org/lapack/double/dpbtrf.f</a>]
\ingroup vc_lapack_core
*/
void DPBTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _info);

# define SPBTRF _LAPACK_FUNCTION(SPBTRF,spbtrf)


/** LAPACK SPBTRF.
Computes the Cholesky factorization of a symmetric positive definite band matrix.
[<a href="http://www.netlib.org/lapack/single/spbtrf.f">http://www.netlib.org/lapack/single/spbtrf.f</a>]
\ingroup vc_lapack_core
*/
void SPBTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _info);


# define DGBBRD _LAPACK_FUNCTION(DGBBRD,dgbbrd)


/** LAPACK DGBBRD.

[<a href="http://www.netlib.org/lapack/double/dgbbrd.f">http://www.netlib.org/lapack/double/dgbbrd.f</a>]
\ingroup vc_lapack_core
*/
void DGBBRD(
        const char* _vect,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ncc,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _d,
        double* _e,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _pt,
        ::VC::math::lapack::INTEGER* _ldpt,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGBBRD _LAPACK_FUNCTION(SGBBRD,sgbbrd)


/** LAPACK SGBBRD.

[<a href="http://www.netlib.org/lapack/single/sgbbrd.f">http://www.netlib.org/lapack/single/sgbbrd.f</a>]
\ingroup vc_lapack_core
*/
void SGBBRD(
        const char* _vect,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ncc,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _d,
        float* _e,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _pt,
        ::VC::math::lapack::INTEGER* _ldpt,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DHSEQR _LAPACK_FUNCTION(DHSEQR,dhseqr)


/** LAPACK DHSEQR.
Computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the multishift QR algorithm.
[<a href="http://www.netlib.org/lapack/double/dhseqr.f">http://www.netlib.org/lapack/double/dhseqr.f</a>]
\ingroup vc_lapack_core
*/
void DHSEQR(
        const char* _job,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        double* _wr,
        double* _wi,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SHSEQR _LAPACK_FUNCTION(SHSEQR,shseqr)


/** LAPACK SHSEQR.
Computes the eigenvalues and Schur factorization of an upper Hessenberg matrix, using the multishift QR algorithm.
[<a href="http://www.netlib.org/lapack/single/shseqr.f">http://www.netlib.org/lapack/single/shseqr.f</a>]
\ingroup vc_lapack_core
*/
void SHSEQR(
        const char* _job,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        float* _wr,
        float* _wi,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBGVD _LAPACK_FUNCTION(DSBGVD,dsbgvd)


/** LAPACK DSBGVD.

[<a href="http://www.netlib.org/lapack/double/dsbgvd.f">http://www.netlib.org/lapack/double/dsbgvd.f</a>]
\ingroup vc_lapack_core
*/
void DSBGVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSBGVD _LAPACK_FUNCTION(SSBGVD,ssbgvd)


/** LAPACK SSBGVD.

[<a href="http://www.netlib.org/lapack/single/ssbgvd.f">http://www.netlib.org/lapack/single/ssbgvd.f</a>]
\ingroup vc_lapack_core
*/
void SSBGVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEQRF _LAPACK_FUNCTION(DGEQRF,dgeqrf)


/** LAPACK DGEQRF.
Computes a QR factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgeqrf.f">http://www.netlib.org/lapack/double/dgeqrf.f</a>]
\ingroup vc_lapack_core
*/
void DGEQRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEQRF _LAPACK_FUNCTION(SGEQRF,sgeqrf)


/** LAPACK SGEQRF.
Computes a QR factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgeqrf.f">http://www.netlib.org/lapack/single/sgeqrf.f</a>]
\ingroup vc_lapack_core
*/
void SGEQRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGELSS _LAPACK_FUNCTION(DGELSS,dgelss)


/** LAPACK DGELSS.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B,  using the singular value decomposition of A.
[<a href="http://www.netlib.org/lapack/double/dgelss.f">http://www.netlib.org/lapack/double/dgelss.f</a>]
\ingroup vc_lapack_core
*/
void DGELSS(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _s,
        const double* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGELSS _LAPACK_FUNCTION(SGELSS,sgelss)


/** LAPACK SGELSS.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B,  using the singular value decomposition of A.
[<a href="http://www.netlib.org/lapack/single/sgelss.f">http://www.netlib.org/lapack/single/sgelss.f</a>]
\ingroup vc_lapack_core
*/
void SGELSS(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _s,
        const float* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSYSV _LAPACK_FUNCTION(DSYSV,dsysv)


/** LAPACK DSYSV.
Solves a real symmetric indefinite system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dsysv.f">http://www.netlib.org/lapack/double/dsysv.f</a>]
\ingroup vc_lapack_core
*/
void DSYSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYSV _LAPACK_FUNCTION(SSYSV,ssysv)


/** LAPACK SSYSV.
Solves a real symmetric indefinite system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/ssysv.f">http://www.netlib.org/lapack/single/ssysv.f</a>]
\ingroup vc_lapack_core
*/
void SSYSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTREXC _LAPACK_FUNCTION(DTREXC,dtrexc)


/** LAPACK DTREXC.
Reorders the Schur factorization of a matrix by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/double/dtrexc.f">http://www.netlib.org/lapack/double/dtrexc.f</a>]
\ingroup vc_lapack_core
*/
void DTREXC(
        const char* _compq,
        ::VC::math::lapack::INTEGER* _n,
        double* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        ::VC::math::lapack::INTEGER* _ifst,
        ::VC::math::lapack::INTEGER* _ilst,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define STREXC _LAPACK_FUNCTION(STREXC,strexc)


/** LAPACK STREXC.
Reorders the Schur factorization of a matrix by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/single/strexc.f">http://www.netlib.org/lapack/single/strexc.f</a>]
\ingroup vc_lapack_core
*/
void STREXC(
        const char* _compq,
        ::VC::math::lapack::INTEGER* _n,
        float* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        ::VC::math::lapack::INTEGER* _ifst,
        ::VC::math::lapack::INTEGER* _ilst,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGEBAK _LAPACK_FUNCTION(DGEBAK,dgebak)


/** LAPACK DGEBAK.
Transforms eigenvectors of a balanced matrix to those of the original matrix supplied to DGEBAL.
[<a href="http://www.netlib.org/lapack/double/dgebak.f">http://www.netlib.org/lapack/double/dgebak.f</a>]
\ingroup vc_lapack_core
*/
void DGEBAK(
        const char* _job,
        const char* _side,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const double* _scale,
        ::VC::math::lapack::INTEGER* _m,
        double* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        ::VC::math::lapack::INTEGER* _info);

# define SGEBAK _LAPACK_FUNCTION(SGEBAK,sgebak)


/** LAPACK SGEBAK.
Transforms eigenvectors of a balanced matrix to those of the original matrix supplied to DGEBAL.
[<a href="http://www.netlib.org/lapack/single/sgebak.f">http://www.netlib.org/lapack/single/sgebak.f</a>]
\ingroup vc_lapack_core
*/
void SGEBAK(
        const char* _job,
        const char* _side,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const float* _scale,
        ::VC::math::lapack::INTEGER* _m,
        float* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        ::VC::math::lapack::INTEGER* _info);


# define DGEHRD _LAPACK_FUNCTION(DGEHRD,dgehrd)


/** LAPACK DGEHRD.
Reduces a general matrix to upper Hessenberg form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/double/dgehrd.f">http://www.netlib.org/lapack/double/dgehrd.f</a>]
\ingroup vc_lapack_core
*/
void DGEHRD(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEHRD _LAPACK_FUNCTION(SGEHRD,sgehrd)


/** LAPACK SGEHRD.
Reduces a general matrix to upper Hessenberg form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/single/sgehrd.f">http://www.netlib.org/lapack/single/sgehrd.f</a>]
\ingroup vc_lapack_core
*/
void SGEHRD(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORMBR _LAPACK_FUNCTION(DORMBR,dormbr)


/** LAPACK DORMBR.
Multiplies a general matrix by one of the orthogonal transformation  matrices from a reduction to bidiagonal form determined by DGEBRD.
[<a href="http://www.netlib.org/lapack/double/dormbr.f">http://www.netlib.org/lapack/double/dormbr.f</a>]
\ingroup vc_lapack_core
*/
void DORMBR(
        const char* _vect,
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMBR _LAPACK_FUNCTION(SORMBR,sormbr)


/** LAPACK SORMBR.
Multiplies a general matrix by one of the orthogonal transformation  matrices from a reduction to bidiagonal form determined by DGEBRD.
[<a href="http://www.netlib.org/lapack/single/sormbr.f">http://www.netlib.org/lapack/single/sormbr.f</a>]
\ingroup vc_lapack_core
*/
void SORMBR(
        const char* _vect,
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTRSNA _LAPACK_FUNCTION(DTRSNA,dtrsna)


/** LAPACK DTRSNA.
Estimates the reciprocal condition numbers (sensitivities) of selected eigenvalues and eigenvectors of an upper quasi-triangular matrix.
[<a href="http://www.netlib.org/lapack/double/dtrsna.f">http://www.netlib.org/lapack/double/dtrsna.f</a>]
\ingroup vc_lapack_core
*/
void DTRSNA(
        const char* _job,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const double* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        const double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        const double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        double* _s,
        double* _sep,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        double* _work,
        ::VC::math::lapack::INTEGER* _ldwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STRSNA _LAPACK_FUNCTION(STRSNA,strsna)


/** LAPACK STRSNA.
Estimates the reciprocal condition numbers (sensitivities) of selected eigenvalues and eigenvectors of an upper quasi-triangular matrix.
[<a href="http://www.netlib.org/lapack/single/strsna.f">http://www.netlib.org/lapack/single/strsna.f</a>]
\ingroup vc_lapack_core
*/
void STRSNA(
        const char* _job,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const float* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        const float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        const float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        float* _s,
        float* _sep,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        float* _work,
        ::VC::math::lapack::INTEGER* _ldwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPBRFS _LAPACK_FUNCTION(DPBRFS,dpbrfs)


/** LAPACK DPBRFS.
Improves the computed solution to a symmetric positive definite banded system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dpbrfs.f">http://www.netlib.org/lapack/double/dpbrfs.f</a>]
\ingroup vc_lapack_core
*/
void DPBRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const double* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPBRFS _LAPACK_FUNCTION(SPBRFS,spbrfs)


/** LAPACK SPBRFS.
Improves the computed solution to a symmetric positive definite banded system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/spbrfs.f">http://www.netlib.org/lapack/single/spbrfs.f</a>]
\ingroup vc_lapack_core
*/
void SPBRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const float* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEBAL _LAPACK_FUNCTION(DGEBAL,dgebal)


/** LAPACK DGEBAL.
Balances a general matrix in order to improve the accuracy of computed eigenvalues.
[<a href="http://www.netlib.org/lapack/double/dgebal.f">http://www.netlib.org/lapack/double/dgebal.f</a>]
\ingroup vc_lapack_core
*/
void DGEBAL(
        const char* _job,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _scale,
        ::VC::math::lapack::INTEGER* _info);

# define SGEBAL _LAPACK_FUNCTION(SGEBAL,sgebal)


/** LAPACK SGEBAL.
Balances a general matrix in order to improve the accuracy of computed eigenvalues.
[<a href="http://www.netlib.org/lapack/single/sgebal.f">http://www.netlib.org/lapack/single/sgebal.f</a>]
\ingroup vc_lapack_core
*/
void SGEBAL(
        const char* _job,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _scale,
        ::VC::math::lapack::INTEGER* _info);


# define DHSEIN _LAPACK_FUNCTION(DHSEIN,dhsein)


/** LAPACK DHSEIN.
Computes specified right and/or left eigenvectors of an upper Hessenberg matrix by inverse iteration.
[<a href="http://www.netlib.org/lapack/double/dhsein.f">http://www.netlib.org/lapack/double/dhsein.f</a>]
\ingroup vc_lapack_core
*/
void DHSEIN(
        const char* _side,
        const char* _eigsrc,
        const char* _initv,
        ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const double* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        double* _wr,
        const double* _wi,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        double* _work,
        ::VC::math::lapack::INTEGER* _ifaill,
        ::VC::math::lapack::INTEGER* _ifailr,
        ::VC::math::lapack::INTEGER* _info);

# define SHSEIN _LAPACK_FUNCTION(SHSEIN,shsein)


/** LAPACK SHSEIN.
Computes specified right and/or left eigenvectors of an upper Hessenberg matrix by inverse iteration.
[<a href="http://www.netlib.org/lapack/single/shsein.f">http://www.netlib.org/lapack/single/shsein.f</a>]
\ingroup vc_lapack_core
*/
void SHSEIN(
        const char* _side,
        const char* _eigsrc,
        const char* _initv,
        ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const float* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        float* _wr,
        const float* _wi,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        float* _work,
        ::VC::math::lapack::INTEGER* _ifaill,
        ::VC::math::lapack::INTEGER* _ifailr,
        ::VC::math::lapack::INTEGER* _info);


# define DSYEV _LAPACK_FUNCTION(DSYEV,dsyev)


/** LAPACK DSYEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix.
[<a href="http://www.netlib.org/lapack/double/dsyev.f">http://www.netlib.org/lapack/double/dsyev.f</a>]
\ingroup vc_lapack_core
*/
void DSYEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _w,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYEV _LAPACK_FUNCTION(SSYEV,ssyev)


/** LAPACK SSYEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix.
[<a href="http://www.netlib.org/lapack/single/ssyev.f">http://www.netlib.org/lapack/single/ssyev.f</a>]
\ingroup vc_lapack_core
*/
void SSYEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _w,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORGQR _LAPACK_FUNCTION(DORGQR,dorgqr)


/** LAPACK DORGQR.
Generates all or part of the orthogonal matrix Q from a QR factorization determined by DGEQRF.
[<a href="http://www.netlib.org/lapack/double/dorgqr.f">http://www.netlib.org/lapack/double/dorgqr.f</a>]
\ingroup vc_lapack_core
*/
void DORGQR(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGQR _LAPACK_FUNCTION(SORGQR,sorgqr)


/** LAPACK SORGQR.
Generates all or part of the orthogonal matrix Q from a QR factorization determined by DGEQRF.
[<a href="http://www.netlib.org/lapack/single/sorgqr.f">http://www.netlib.org/lapack/single/sorgqr.f</a>]
\ingroup vc_lapack_core
*/
void SORGQR(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTPRFS _LAPACK_FUNCTION(DTPRFS,dtprfs)


/** LAPACK DTPRFS.
Provides forward and backward error bounds for the solution of a triangular system of linear equations AX=B, A**T X=B or A**H X=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/double/dtprfs.f">http://www.netlib.org/lapack/double/dtprfs.f</a>]
\ingroup vc_lapack_core
*/
void DTPRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STPRFS _LAPACK_FUNCTION(STPRFS,stprfs)


/** LAPACK STPRFS.
Provides forward and backward error bounds for the solution of a triangular system of linear equations AX=B, A**T X=B or A**H X=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/single/stprfs.f">http://www.netlib.org/lapack/single/stprfs.f</a>]
\ingroup vc_lapack_core
*/
void STPRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPTTRS _LAPACK_FUNCTION(DPTTRS,dpttrs)


/** LAPACK DPTTRS.
Solves a symmetric positive definite tridiagonal system of linear equations, using the LDL**H factorization computed by DPTTRF.
[<a href="http://www.netlib.org/lapack/double/dpttrs.f">http://www.netlib.org/lapack/double/dpttrs.f</a>]
\ingroup vc_lapack_core
*/
void DPTTRS(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _d,
        const double* _e,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPTTRS _LAPACK_FUNCTION(SPTTRS,spttrs)


/** LAPACK SPTTRS.
Solves a symmetric positive definite tridiagonal system of linear equations, using the LDL**H factorization computed by DPTTRF.
[<a href="http://www.netlib.org/lapack/single/spttrs.f">http://www.netlib.org/lapack/single/spttrs.f</a>]
\ingroup vc_lapack_core
*/
void SPTTRS(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _d,
        const float* _e,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DORGTR _LAPACK_FUNCTION(DORGTR,dorgtr)


/** LAPACK DORGTR.
Generates the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSYTRD.
[<a href="http://www.netlib.org/lapack/double/dorgtr.f">http://www.netlib.org/lapack/double/dorgtr.f</a>]
\ingroup vc_lapack_core
*/
void DORGTR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGTR _LAPACK_FUNCTION(SORGTR,sorgtr)


/** LAPACK SORGTR.
Generates the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSYTRD.
[<a href="http://www.netlib.org/lapack/single/sorgtr.f">http://www.netlib.org/lapack/single/sorgtr.f</a>]
\ingroup vc_lapack_core
*/
void SORGTR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGESDD _LAPACK_FUNCTION(DGESDD,dgesdd)


/** LAPACK DGESDD.
Computes the singular value decomposition (SVD) of a general rectangular matrix using divide-and-conquer.
[<a href="http://www.netlib.org/lapack/double/dgesdd.f">http://www.netlib.org/lapack/double/dgesdd.f</a>]
\ingroup vc_lapack_core
*/
void DGESDD(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _s,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGESDD _LAPACK_FUNCTION(SGESDD,sgesdd)


/** LAPACK SGESDD.
Computes the singular value decomposition (SVD) of a general rectangular matrix using divide-and-conquer.
[<a href="http://www.netlib.org/lapack/single/sgesdd.f">http://www.netlib.org/lapack/single/sgesdd.f</a>]
\ingroup vc_lapack_core
*/
void SGESDD(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _s,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPTSVX _LAPACK_FUNCTION(DPTSVX,dptsvx)


/** LAPACK DPTSVX.
Solves a symmetric positive definite tridiagonal system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dptsvx.f">http://www.netlib.org/lapack/double/dptsvx.f</a>]
\ingroup vc_lapack_core
*/
void DPTSVX(
        const char* _fact,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _d,
        const double* _e,
        double* _df,
        double* _ef,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SPTSVX _LAPACK_FUNCTION(SPTSVX,sptsvx)


/** LAPACK SPTSVX.
Solves a symmetric positive definite tridiagonal system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sptsvx.f">http://www.netlib.org/lapack/single/sptsvx.f</a>]
\ingroup vc_lapack_core
*/
void SPTSVX(
        const char* _fact,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _d,
        const float* _e,
        float* _df,
        float* _ef,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSYEVR _LAPACK_FUNCTION(DSYEVR,dsyevr)


/** LAPACK DSYEVR.
Computes selected eigenvalues, and optionally, eigenvectors of a real symmetric matrix.  Eigenvalues are computed by the dqds algorithm, and eigenvectors are computed from various "good" LDL^T representations (also known as Relatively Robust Representations).
[<a href="http://www.netlib.org/lapack/double/dsyevr.f">http://www.netlib.org/lapack/double/dsyevr.f</a>]
\ingroup vc_lapack_core
*/
void DSYEVR(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYEVR _LAPACK_FUNCTION(SSYEVR,ssyevr)


/** LAPACK SSYEVR.
Computes selected eigenvalues, and optionally, eigenvectors of a real symmetric matrix.  Eigenvalues are computed by the dqds algorithm, and eigenvectors are computed from various "good" LDL^T representations (also known as Relatively Robust Representations).
[<a href="http://www.netlib.org/lapack/single/ssyevr.f">http://www.netlib.org/lapack/single/ssyevr.f</a>]
\ingroup vc_lapack_core
*/
void SSYEVR(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGQRF _LAPACK_FUNCTION(DGGQRF,dggqrf)


/** LAPACK DGGQRF.

[<a href="http://www.netlib.org/lapack/double/dggqrf.f">http://www.netlib.org/lapack/double/dggqrf.f</a>]
\ingroup vc_lapack_core
*/
void DGGQRF(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _taua,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _taub,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGQRF _LAPACK_FUNCTION(SGGQRF,sggqrf)


/** LAPACK SGGQRF.

[<a href="http://www.netlib.org/lapack/single/sggqrf.f">http://www.netlib.org/lapack/single/sggqrf.f</a>]
\ingroup vc_lapack_core
*/
void SGGQRF(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _taua,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _taub,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEGS _LAPACK_FUNCTION(DGEGS,dgegs)


/** LAPACK DGEGS.

[<a href="http://www.netlib.org/lapack/double/dgegs.f">http://www.netlib.org/lapack/double/dgegs.f</a>]
\ingroup vc_lapack_core
*/
void DGEGS(
        const char* _jobvsl,
        const char* _jobvsr,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        double* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEGS _LAPACK_FUNCTION(SGEGS,sgegs)


/** LAPACK SGEGS.

[<a href="http://www.netlib.org/lapack/single/sgegs.f">http://www.netlib.org/lapack/single/sgegs.f</a>]
\ingroup vc_lapack_core
*/
void SGEGS(
        const char* _jobvsl,
        const char* _jobvsr,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        float* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTRCON _LAPACK_FUNCTION(DTRCON,dtrcon)


/** LAPACK DTRCON.
Estimates the reciprocal of the condition number of a triangular matrix, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/double/dtrcon.f">http://www.netlib.org/lapack/double/dtrcon.f</a>]
\ingroup vc_lapack_core
*/
void DTRCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STRCON _LAPACK_FUNCTION(STRCON,strcon)


/** LAPACK STRCON.
Estimates the reciprocal of the condition number of a triangular matrix, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/single/strcon.f">http://www.netlib.org/lapack/single/strcon.f</a>]
\ingroup vc_lapack_core
*/
void STRCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTPTRI _LAPACK_FUNCTION(DTPTRI,dtptri)


/** LAPACK DTPTRI.
Computes the inverse of a triangular matrix in packed storage.
[<a href="http://www.netlib.org/lapack/double/dtptri.f">http://www.netlib.org/lapack/double/dtptri.f</a>]
\ingroup vc_lapack_core
*/
void DTPTRI(
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        ::VC::math::lapack::INTEGER* _info);

# define STPTRI _LAPACK_FUNCTION(STPTRI,stptri)


/** LAPACK STPTRI.
Computes the inverse of a triangular matrix in packed storage.
[<a href="http://www.netlib.org/lapack/single/stptri.f">http://www.netlib.org/lapack/single/stptri.f</a>]
\ingroup vc_lapack_core
*/
void STPTRI(
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        ::VC::math::lapack::INTEGER* _info);


# define DPPSV _LAPACK_FUNCTION(DPPSV,dppsv)


/** LAPACK DPPSV.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/double/dppsv.f">http://www.netlib.org/lapack/double/dppsv.f</a>]
\ingroup vc_lapack_core
*/
void DPPSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ap,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPPSV _LAPACK_FUNCTION(SPPSV,sppsv)


/** LAPACK SPPSV.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/single/sppsv.f">http://www.netlib.org/lapack/single/sppsv.f</a>]
\ingroup vc_lapack_core
*/
void SPPSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ap,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSPSV _LAPACK_FUNCTION(DSPSV,dspsv)


/** LAPACK DSPSV.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/double/dspsv.f">http://www.netlib.org/lapack/double/dspsv.f</a>]
\ingroup vc_lapack_core
*/
void DSPSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ap,
        ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SSPSV _LAPACK_FUNCTION(SSPSV,sspsv)


/** LAPACK SSPSV.
Solves a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/single/sspsv.f">http://www.netlib.org/lapack/single/sspsv.f</a>]
\ingroup vc_lapack_core
*/
void SSPSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ap,
        ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGELSX _LAPACK_FUNCTION(DGELSX,dgelsx)


/** LAPACK DGELSX.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B, using a complete orthogonal factorization of A.
[<a href="http://www.netlib.org/lapack/double/dgelsx.f">http://www.netlib.org/lapack/double/dgelsx.f</a>]
\ingroup vc_lapack_core
*/
void DGELSX(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _jpvt,
        const double* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGELSX _LAPACK_FUNCTION(SGELSX,sgelsx)


/** LAPACK SGELSX.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B, using a complete orthogonal factorization of A.
[<a href="http://www.netlib.org/lapack/single/sgelsx.f">http://www.netlib.org/lapack/single/sgelsx.f</a>]
\ingroup vc_lapack_core
*/
void SGELSX(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _jpvt,
        const float* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DBDSQR _LAPACK_FUNCTION(DBDSQR,dbdsqr)


/** LAPACK DBDSQR.
Computes the singular value decomposition (SVD) of a real bidiagonal matrix, using the bidiagonal QR algorithm.
[<a href="http://www.netlib.org/lapack/double/dbdsqr.f">http://www.netlib.org/lapack/double/dbdsqr.f</a>]
\ingroup vc_lapack_core
*/
void DBDSQR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ncvt,
        ::VC::math::lapack::INTEGER* _nru,
        ::VC::math::lapack::INTEGER* _ncc,
        double* _d,
        double* _e,
        double* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SBDSQR _LAPACK_FUNCTION(SBDSQR,sbdsqr)


/** LAPACK SBDSQR.
Computes the singular value decomposition (SVD) of a real bidiagonal matrix, using the bidiagonal QR algorithm.
[<a href="http://www.netlib.org/lapack/single/sbdsqr.f">http://www.netlib.org/lapack/single/sbdsqr.f</a>]
\ingroup vc_lapack_core
*/
void SBDSQR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ncvt,
        ::VC::math::lapack::INTEGER* _nru,
        ::VC::math::lapack::INTEGER* _ncc,
        float* _d,
        float* _e,
        float* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DTGEXC _LAPACK_FUNCTION(DTGEXC,dtgexc)


/** LAPACK DTGEXC.

[<a href="http://www.netlib.org/lapack/double/dtgexc.f">http://www.netlib.org/lapack/double/dtgexc.f</a>]
\ingroup vc_lapack_core
*/
void DTGEXC(
        const ::VC::math::lapack::LOGICAL* _wantq,
        const ::VC::math::lapack::LOGICAL* _wantz,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _ifst,
        ::VC::math::lapack::INTEGER* _ilst,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define STGEXC _LAPACK_FUNCTION(STGEXC,stgexc)


/** LAPACK STGEXC.

[<a href="http://www.netlib.org/lapack/single/stgexc.f">http://www.netlib.org/lapack/single/stgexc.f</a>]
\ingroup vc_lapack_core
*/
void STGEXC(
        const ::VC::math::lapack::LOGICAL* _wantq,
        const ::VC::math::lapack::LOGICAL* _wantz,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _ifst,
        ::VC::math::lapack::INTEGER* _ilst,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPBSTF _LAPACK_FUNCTION(DPBSTF,dpbstf)


/** LAPACK DPBSTF.

[<a href="http://www.netlib.org/lapack/double/dpbstf.f">http://www.netlib.org/lapack/double/dpbstf.f</a>]
\ingroup vc_lapack_core
*/
void DPBSTF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _info);

# define SPBSTF _LAPACK_FUNCTION(SPBSTF,spbstf)


/** LAPACK SPBSTF.

[<a href="http://www.netlib.org/lapack/single/spbstf.f">http://www.netlib.org/lapack/single/spbstf.f</a>]
\ingroup vc_lapack_core
*/
void SPBSTF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _info);


# define DGELSY _LAPACK_FUNCTION(DGELSY,dgelsy)


/** LAPACK DGELSY.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B, using a complete orthogonal factorization of A.
[<a href="http://www.netlib.org/lapack/double/dgelsy.f">http://www.netlib.org/lapack/double/dgelsy.f</a>]
\ingroup vc_lapack_core
*/
void DGELSY(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _jpvt,
        const double* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGELSY _LAPACK_FUNCTION(SGELSY,sgelsy)


/** LAPACK SGELSY.
Computes the minimum norm least squares solution to an over- or under-determined system of linear equations A X=B, using a complete orthogonal factorization of A.
[<a href="http://www.netlib.org/lapack/single/sgelsy.f">http://www.netlib.org/lapack/single/sgelsy.f</a>]
\ingroup vc_lapack_core
*/
void SGELSY(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _jpvt,
        const float* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DBDSDC _LAPACK_FUNCTION(DBDSDC,dbdsdc)


/** LAPACK DBDSDC.
Computes the singular value decomposition (SVD) of a real bidiagonal matrix, using a divide and conquer method.
[<a href="http://www.netlib.org/lapack/double/dbdsdc.f">http://www.netlib.org/lapack/double/dbdsdc.f</a>]
\ingroup vc_lapack_core
*/
void DBDSDC(
        const char* _uplo,
        const char* _compq,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        double* _q,
        ::VC::math::lapack::INTEGER* _iq,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SBDSDC _LAPACK_FUNCTION(SBDSDC,sbdsdc)


/** LAPACK SBDSDC.
Computes the singular value decomposition (SVD) of a real bidiagonal matrix, using a divide and conquer method.
[<a href="http://www.netlib.org/lapack/single/sbdsdc.f">http://www.netlib.org/lapack/single/sbdsdc.f</a>]
\ingroup vc_lapack_core
*/
void SBDSDC(
        const char* _uplo,
        const char* _compq,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        float* _q,
        ::VC::math::lapack::INTEGER* _iq,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORMHR _LAPACK_FUNCTION(DORMHR,dormhr)


/** LAPACK DORMHR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to Hessenberg form determined by DGEHRD.
[<a href="http://www.netlib.org/lapack/double/dormhr.f">http://www.netlib.org/lapack/double/dormhr.f</a>]
\ingroup vc_lapack_core
*/
void DORMHR(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMHR _LAPACK_FUNCTION(SORMHR,sormhr)


/** LAPACK SORMHR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to Hessenberg form determined by DGEHRD.
[<a href="http://www.netlib.org/lapack/single/sormhr.f">http://www.netlib.org/lapack/single/sormhr.f</a>]
\ingroup vc_lapack_core
*/
void SORMHR(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGBAK _LAPACK_FUNCTION(DGGBAK,dggbak)


/** LAPACK DGGBAK.

[<a href="http://www.netlib.org/lapack/double/dggbak.f">http://www.netlib.org/lapack/double/dggbak.f</a>]
\ingroup vc_lapack_core
*/
void DGGBAK(
        const char* _job,
        const char* _side,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const double* _lscale,
        const double* _rscale,
        ::VC::math::lapack::INTEGER* _m,
        double* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        ::VC::math::lapack::INTEGER* _info);

# define SGGBAK _LAPACK_FUNCTION(SGGBAK,sggbak)


/** LAPACK SGGBAK.

[<a href="http://www.netlib.org/lapack/single/sggbak.f">http://www.netlib.org/lapack/single/sggbak.f</a>]
\ingroup vc_lapack_core
*/
void SGGBAK(
        const char* _job,
        const char* _side,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        const float* _lscale,
        const float* _rscale,
        ::VC::math::lapack::INTEGER* _m,
        float* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        ::VC::math::lapack::INTEGER* _info);


# define DORMQL _LAPACK_FUNCTION(DORMQL,dormql)


/** LAPACK DORMQL.
Multiplies a general matrix by the orthogonal matrix from a QL factorization determined by DGEQLF.
[<a href="http://www.netlib.org/lapack/double/dormql.f">http://www.netlib.org/lapack/double/dormql.f</a>]
\ingroup vc_lapack_core
*/
void DORMQL(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMQL _LAPACK_FUNCTION(SORMQL,sormql)


/** LAPACK SORMQL.
Multiplies a general matrix by the orthogonal matrix from a QL factorization determined by DGEQLF.
[<a href="http://www.netlib.org/lapack/single/sormql.f">http://www.netlib.org/lapack/single/sormql.f</a>]
\ingroup vc_lapack_core
*/
void SORMQL(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGHRD _LAPACK_FUNCTION(DGGHRD,dgghrd)


/** LAPACK DGGHRD.

[<a href="http://www.netlib.org/lapack/double/dgghrd.f">http://www.netlib.org/lapack/double/dgghrd.f</a>]
\ingroup vc_lapack_core
*/
void DGGHRD(
        const char* _compq,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _info);

# define SGGHRD _LAPACK_FUNCTION(SGGHRD,sgghrd)


/** LAPACK SGGHRD.

[<a href="http://www.netlib.org/lapack/single/sgghrd.f">http://www.netlib.org/lapack/single/sgghrd.f</a>]
\ingroup vc_lapack_core
*/
void SGGHRD(
        const char* _compq,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _info);


# define DTGSNA _LAPACK_FUNCTION(DTGSNA,dtgsna)


/** LAPACK DTGSNA.

[<a href="http://www.netlib.org/lapack/double/dtgsna.f">http://www.netlib.org/lapack/double/dtgsna.f</a>]
\ingroup vc_lapack_core
*/
void DTGSNA(
        const char* _job,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        const double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        double* _s,
        double* _dif,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STGSNA _LAPACK_FUNCTION(STGSNA,stgsna)


/** LAPACK STGSNA.

[<a href="http://www.netlib.org/lapack/single/stgsna.f">http://www.netlib.org/lapack/single/stgsna.f</a>]
\ingroup vc_lapack_core
*/
void STGSNA(
        const char* _job,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        const float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        float* _s,
        float* _dif,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSYGST _LAPACK_FUNCTION(DSYGST,dsygst)


/** LAPACK DSYGST.
Reduces a symmetric-definite generalized eigenproblem Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard form, where B has been factorized by DPOTRF.
[<a href="http://www.netlib.org/lapack/double/dsygst.f">http://www.netlib.org/lapack/double/dsygst.f</a>]
\ingroup vc_lapack_core
*/
void DSYGST(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SSYGST _LAPACK_FUNCTION(SSYGST,ssygst)


/** LAPACK SSYGST.
Reduces a symmetric-definite generalized eigenproblem Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard form, where B has been factorized by DPOTRF.
[<a href="http://www.netlib.org/lapack/single/ssygst.f">http://www.netlib.org/lapack/single/ssygst.f</a>]
\ingroup vc_lapack_core
*/
void SSYGST(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DPOTRS _LAPACK_FUNCTION(DPOTRS,dpotrs)


/** LAPACK DPOTRS.
Solves a symmetric positive definite system of linear equations AX=B, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/double/dpotrs.f">http://www.netlib.org/lapack/double/dpotrs.f</a>]
\ingroup vc_lapack_core
*/
void DPOTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPOTRS _LAPACK_FUNCTION(SPOTRS,spotrs)


/** LAPACK SPOTRS.
Solves a symmetric positive definite system of linear equations AX=B, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/single/spotrs.f">http://www.netlib.org/lapack/single/spotrs.f</a>]
\ingroup vc_lapack_core
*/
void SPOTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEVR _LAPACK_FUNCTION(DSTEVR,dstevr)


/** LAPACK DSTEVR.
Computes selected eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.  Eigenvalues are computed by the dqds algorithm, and eigenvectors are computed from various "good" LDL^T representations (also known as Relatively Robust Representations).
[<a href="http://www.netlib.org/lapack/double/dstevr.f">http://www.netlib.org/lapack/double/dstevr.f</a>]
\ingroup vc_lapack_core
*/
void DSTEVR(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEVR _LAPACK_FUNCTION(SSTEVR,sstevr)


/** LAPACK SSTEVR.
Computes selected eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.  Eigenvalues are computed by the dqds algorithm, and eigenvectors are computed from various "good" LDL^T representations (also known as Relatively Robust Representations).
[<a href="http://www.netlib.org/lapack/single/sstevr.f">http://www.netlib.org/lapack/single/sstevr.f</a>]
\ingroup vc_lapack_core
*/
void SSTEVR(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _isuppz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGBAL _LAPACK_FUNCTION(DGGBAL,dggbal)


/** LAPACK DGGBAL.

[<a href="http://www.netlib.org/lapack/double/dggbal.f">http://www.netlib.org/lapack/double/dggbal.f</a>]
\ingroup vc_lapack_core
*/
void DGGBAL(
        const char* _job,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _lscale,
        double* _rscale,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGGBAL _LAPACK_FUNCTION(SGGBAL,sggbal)


/** LAPACK SGGBAL.

[<a href="http://www.netlib.org/lapack/single/sggbal.f">http://www.netlib.org/lapack/single/sggbal.f</a>]
\ingroup vc_lapack_core
*/
void SGGBAL(
        const char* _job,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _lscale,
        float* _rscale,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DPOSVX _LAPACK_FUNCTION(DPOSVX,dposvx)


/** LAPACK DPOSVX.
Solves a symmetric positive definite system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dposvx.f">http://www.netlib.org/lapack/double/dposvx.f</a>]
\ingroup vc_lapack_core
*/
void DPOSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        char* _equed,
        double* _s,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPOSVX _LAPACK_FUNCTION(SPOSVX,sposvx)


/** LAPACK SPOSVX.
Solves a symmetric positive definite system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sposvx.f">http://www.netlib.org/lapack/single/sposvx.f</a>]
\ingroup vc_lapack_core
*/
void SPOSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        char* _equed,
        float* _s,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEGV _LAPACK_FUNCTION(DGEGV,dgegv)


/** LAPACK DGEGV.

[<a href="http://www.netlib.org/lapack/double/dgegv.f">http://www.netlib.org/lapack/double/dgegv.f</a>]
\ingroup vc_lapack_core
*/
void DGEGV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEGV _LAPACK_FUNCTION(SGEGV,sgegv)


/** LAPACK SGEGV.

[<a href="http://www.netlib.org/lapack/single/sgegv.f">http://www.netlib.org/lapack/single/sgegv.f</a>]
\ingroup vc_lapack_core
*/
void SGEGV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGERQF _LAPACK_FUNCTION(DGERQF,dgerqf)


/** LAPACK DGERQF.
Computes an RQ factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgerqf.f">http://www.netlib.org/lapack/double/dgerqf.f</a>]
\ingroup vc_lapack_core
*/
void DGERQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGERQF _LAPACK_FUNCTION(SGERQF,sgerqf)


/** LAPACK SGERQF.
Computes an RQ factorization of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgerqf.f">http://www.netlib.org/lapack/single/sgerqf.f</a>]
\ingroup vc_lapack_core
*/
void SGERQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTRRFS _LAPACK_FUNCTION(DTRRFS,dtrrfs)


/** LAPACK DTRRFS.
Provides forward and backward error bounds for the solution of a triangular system of linear equations A X=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/double/dtrrfs.f">http://www.netlib.org/lapack/double/dtrrfs.f</a>]
\ingroup vc_lapack_core
*/
void DTRRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STRRFS _LAPACK_FUNCTION(STRRFS,strrfs)


/** LAPACK STRRFS.
Provides forward and backward error bounds for the solution of a triangular system of linear equations A X=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/single/strrfs.f">http://www.netlib.org/lapack/single/strrfs.f</a>]
\ingroup vc_lapack_core
*/
void STRRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTPCON _LAPACK_FUNCTION(DTPCON,dtpcon)


/** LAPACK DTPCON.
Estimates the reciprocal of the condition number of a triangular matrix in packed storage, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/double/dtpcon.f">http://www.netlib.org/lapack/double/dtpcon.f</a>]
\ingroup vc_lapack_core
*/
void DTPCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STPCON _LAPACK_FUNCTION(STPCON,stpcon)


/** LAPACK STPCON.
Estimates the reciprocal of the condition number of a triangular matrix in packed storage, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/single/stpcon.f">http://www.netlib.org/lapack/single/stpcon.f</a>]
\ingroup vc_lapack_core
*/
void STPCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DOPGTR _LAPACK_FUNCTION(DOPGTR,dopgtr)


/** LAPACK DOPGTR.
Generates the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSPTRD.
[<a href="http://www.netlib.org/lapack/double/dopgtr.f">http://www.netlib.org/lapack/double/dopgtr.f</a>]
\ingroup vc_lapack_core
*/
void DOPGTR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        const double* _tau,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SOPGTR _LAPACK_FUNCTION(SOPGTR,sopgtr)


/** LAPACK SOPGTR.
Generates the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSPTRD.
[<a href="http://www.netlib.org/lapack/single/sopgtr.f">http://www.netlib.org/lapack/single/sopgtr.f</a>]
\ingroup vc_lapack_core
*/
void SOPGTR(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        const float* _tau,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSPEV _LAPACK_FUNCTION(DSPEV,dspev)


/** LAPACK DSPEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix in packed storage.
[<a href="http://www.netlib.org/lapack/double/dspev.f">http://www.netlib.org/lapack/double/dspev.f</a>]
\ingroup vc_lapack_core
*/
void DSPEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSPEV _LAPACK_FUNCTION(SSPEV,sspev)


/** LAPACK SSPEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix in packed storage.
[<a href="http://www.netlib.org/lapack/single/sspev.f">http://www.netlib.org/lapack/single/sspev.f</a>]
\ingroup vc_lapack_core
*/
void SSPEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DTRTRI _LAPACK_FUNCTION(DTRTRI,dtrtri)


/** LAPACK DTRTRI.
Computes the inverse of a triangular matrix.
[<a href="http://www.netlib.org/lapack/double/dtrtri.f">http://www.netlib.org/lapack/double/dtrtri.f</a>]
\ingroup vc_lapack_core
*/
void DTRTRI(
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);

# define STRTRI _LAPACK_FUNCTION(STRTRI,strtri)


/** LAPACK STRTRI.
Computes the inverse of a triangular matrix.
[<a href="http://www.netlib.org/lapack/single/strtri.f">http://www.netlib.org/lapack/single/strtri.f</a>]
\ingroup vc_lapack_core
*/
void STRTRI(
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _info);


# define DTRSYL _LAPACK_FUNCTION(DTRSYL,dtrsyl)


/** LAPACK DTRSYL.
Solves the Sylvester matrix equation A X +/- X B=C where A and B are upper quasi-triangular, and may be transposed.
[<a href="http://www.netlib.org/lapack/double/dtrsyl.f">http://www.netlib.org/lapack/double/dtrsyl.f</a>]
\ingroup vc_lapack_core
*/
void DTRSYL(
        const char* _trana,
        const char* _tranb,
        ::VC::math::lapack::INTEGER* _isgn,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _scale,
        ::VC::math::lapack::INTEGER* _info);

# define STRSYL _LAPACK_FUNCTION(STRSYL,strsyl)


/** LAPACK STRSYL.
Solves the Sylvester matrix equation A X +/- X B=C where A and B are upper quasi-triangular, and may be transposed.
[<a href="http://www.netlib.org/lapack/single/strsyl.f">http://www.netlib.org/lapack/single/strsyl.f</a>]
\ingroup vc_lapack_core
*/
void STRSYL(
        const char* _trana,
        const char* _tranb,
        ::VC::math::lapack::INTEGER* _isgn,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _scale,
        ::VC::math::lapack::INTEGER* _info);


# define DPPTRS _LAPACK_FUNCTION(DPPTRS,dpptrs)


/** LAPACK DPPTRS.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/double/dpptrs.f">http://www.netlib.org/lapack/double/dpptrs.f</a>]
\ingroup vc_lapack_core
*/
void DPPTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPPTRS _LAPACK_FUNCTION(SPPTRS,spptrs)


/** LAPACK SPPTRS.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, using the Cholesky factorization computed by DPPTRF.
[<a href="http://www.netlib.org/lapack/single/spptrs.f">http://www.netlib.org/lapack/single/spptrs.f</a>]
\ingroup vc_lapack_core
*/
void SPPTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DPPSVX _LAPACK_FUNCTION(DPPSVX,dppsvx)


/** LAPACK DPPSVX.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dppsvx.f">http://www.netlib.org/lapack/double/dppsvx.f</a>]
\ingroup vc_lapack_core
*/
void DPPSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ap,
        double* _afp,
        char* _equed,
        double* _s,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPPSVX _LAPACK_FUNCTION(SPPSVX,sppsvx)


/** LAPACK SPPSVX.
Solves a symmetric positive definite system of linear equations AX=B, where A is held in packed storage, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/sppsvx.f">http://www.netlib.org/lapack/single/sppsvx.f</a>]
\ingroup vc_lapack_core
*/
void SPPSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ap,
        float* _afp,
        char* _equed,
        float* _s,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTRSEN _LAPACK_FUNCTION(DTRSEN,dtrsen)


/** LAPACK DTRSEN.
Reorders the Schur factorization of a matrix in order to find an orthonormal basis of a right invariant subspace corresponding to selected eigenvalues, and returns reciprocal condition numbers (sensitivities) of the average of the cluster of eigenvalues and of the invariant subspace.
[<a href="http://www.netlib.org/lapack/double/dtrsen.f">http://www.netlib.org/lapack/double/dtrsen.f</a>]
\ingroup vc_lapack_core
*/
void DTRSEN(
        const char* _job,
        const char* _compq,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        double* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _wr,
        double* _wi,
        ::VC::math::lapack::INTEGER* _m,
        double* _s,
        double* _sep,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define STRSEN _LAPACK_FUNCTION(STRSEN,strsen)


/** LAPACK STRSEN.
Reorders the Schur factorization of a matrix in order to find an orthonormal basis of a right invariant subspace corresponding to selected eigenvalues, and returns reciprocal condition numbers (sensitivities) of the average of the cluster of eigenvalues and of the invariant subspace.
[<a href="http://www.netlib.org/lapack/single/strsen.f">http://www.netlib.org/lapack/single/strsen.f</a>]
\ingroup vc_lapack_core
*/
void STRSEN(
        const char* _job,
        const char* _compq,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        float* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _wr,
        float* _wi,
        ::VC::math::lapack::INTEGER* _m,
        float* _s,
        float* _sep,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPGST _LAPACK_FUNCTION(DSPGST,dspgst)


/** LAPACK DSPGST.
Reduces a symmetric-definite generalized eigenproblem Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard form,  where A and B are held in packed storage, and B has been factorized by DPPTRF.
[<a href="http://www.netlib.org/lapack/double/dspgst.f">http://www.netlib.org/lapack/double/dspgst.f</a>]
\ingroup vc_lapack_core
*/
void DSPGST(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        const double* _bp,
        ::VC::math::lapack::INTEGER* _info);

# define SSPGST _LAPACK_FUNCTION(SSPGST,sspgst)


/** LAPACK SSPGST.
Reduces a symmetric-definite generalized eigenproblem Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard form,  where A and B are held in packed storage, and B has been factorized by DPPTRF.
[<a href="http://www.netlib.org/lapack/single/sspgst.f">http://www.netlib.org/lapack/single/sspgst.f</a>]
\ingroup vc_lapack_core
*/
void SSPGST(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        const float* _bp,
        ::VC::math::lapack::INTEGER* _info);


# define DSYEVX _LAPACK_FUNCTION(DSYEVX,dsyevx)


/** LAPACK DSYEVX.
Computes selected eigenvalues and eigenvectors of a symmetric matrix.
[<a href="http://www.netlib.org/lapack/double/dsyevx.f">http://www.netlib.org/lapack/double/dsyevx.f</a>]
\ingroup vc_lapack_core
*/
void DSYEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSYEVX _LAPACK_FUNCTION(SSYEVX,ssyevx)


/** LAPACK SSYEVX.
Computes selected eigenvalues and eigenvectors of a symmetric matrix.
[<a href="http://www.netlib.org/lapack/single/ssyevx.f">http://www.netlib.org/lapack/single/ssyevx.f</a>]
\ingroup vc_lapack_core
*/
void SSYEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DSGESV _LAPACK_FUNCTION(DSGESV,dsgesv)


/** LAPACK DSGESV.

[<a href="http://www.netlib.org/lapack/double/dsgesv.f">http://www.netlib.org/lapack/double/dsgesv.f</a>]
\ingroup vc_lapack_core
*/
void DSGESV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _work,
        double* _swork,
        ::VC::math::lapack::INTEGER* _iter,
        ::VC::math::lapack::INTEGER* _info);

# define SSGESV _LAPACK_FUNCTION(SSGESV,ssgesv)


/** LAPACK SSGESV.

[<a href="http://www.netlib.org/lapack/single/ssgesv.f">http://www.netlib.org/lapack/single/ssgesv.f</a>]
\ingroup vc_lapack_core
*/
void SSGESV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _work,
        float* _swork,
        ::VC::math::lapack::INTEGER* _iter,
        ::VC::math::lapack::INTEGER* _info);


# define DHGEQZ _LAPACK_FUNCTION(DHGEQZ,dhgeqz)


/** LAPACK DHGEQZ.

[<a href="http://www.netlib.org/lapack/double/dhgeqz.f">http://www.netlib.org/lapack/double/dhgeqz.f</a>]
\ingroup vc_lapack_core
*/
void DHGEQZ(
        const char* _job,
        const char* _compq,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        double* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SHGEQZ _LAPACK_FUNCTION(SHGEQZ,shgeqz)


/** LAPACK SHGEQZ.

[<a href="http://www.netlib.org/lapack/single/shgeqz.f">http://www.netlib.org/lapack/single/shgeqz.f</a>]
\ingroup vc_lapack_core
*/
void SHGEQZ(
        const char* _job,
        const char* _compq,
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _h,
        ::VC::math::lapack::INTEGER* _ldh,
        float* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGRQF _LAPACK_FUNCTION(DGGRQF,dggrqf)


/** LAPACK DGGRQF.

[<a href="http://www.netlib.org/lapack/double/dggrqf.f">http://www.netlib.org/lapack/double/dggrqf.f</a>]
\ingroup vc_lapack_core
*/
void DGGRQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _taua,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _taub,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGRQF _LAPACK_FUNCTION(SGGRQF,sggrqf)


/** LAPACK SGGRQF.

[<a href="http://www.netlib.org/lapack/single/sggrqf.f">http://www.netlib.org/lapack/single/sggrqf.f</a>]
\ingroup vc_lapack_core
*/
void SGGRQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _taua,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _taub,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEQPF _LAPACK_FUNCTION(DGEQPF,dgeqpf)


/** LAPACK DGEQPF.
Computes a QR factorization with column pivoting of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgeqpf.f">http://www.netlib.org/lapack/double/dgeqpf.f</a>]
\ingroup vc_lapack_core
*/
void DGEQPF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _jpvt,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SGEQPF _LAPACK_FUNCTION(SGEQPF,sgeqpf)


/** LAPACK SGEQPF.
Computes a QR factorization with column pivoting of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgeqpf.f">http://www.netlib.org/lapack/single/sgeqpf.f</a>]
\ingroup vc_lapack_core
*/
void SGEQPF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _jpvt,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DORGLQ _LAPACK_FUNCTION(DORGLQ,dorglq)


/** LAPACK DORGLQ.
Generates all or part of the orthogonal matrix Q from an LQ factorization determined by DGELQF.
[<a href="http://www.netlib.org/lapack/double/dorglq.f">http://www.netlib.org/lapack/double/dorglq.f</a>]
\ingroup vc_lapack_core
*/
void DORGLQ(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGLQ _LAPACK_FUNCTION(SORGLQ,sorglq)


/** LAPACK SORGLQ.
Generates all or part of the orthogonal matrix Q from an LQ factorization determined by DGELQF.
[<a href="http://www.netlib.org/lapack/single/sorglq.f">http://www.netlib.org/lapack/single/sorglq.f</a>]
\ingroup vc_lapack_core
*/
void SORGLQ(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPTSV _LAPACK_FUNCTION(DPTSV,dptsv)


/** LAPACK DPTSV.
Solves a symmetric positive definite tridiagonal system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dptsv.f">http://www.netlib.org/lapack/double/dptsv.f</a>]
\ingroup vc_lapack_core
*/
void DPTSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _d,
        double* _e,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPTSV _LAPACK_FUNCTION(SPTSV,sptsv)


/** LAPACK SPTSV.
Solves a symmetric positive definite tridiagonal system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/sptsv.f">http://www.netlib.org/lapack/single/sptsv.f</a>]
\ingroup vc_lapack_core
*/
void SPTSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _d,
        float* _e,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DTBCON _LAPACK_FUNCTION(DTBCON,dtbcon)


/** LAPACK DTBCON.
Estimates the reciprocal of the condition number of a triangular band matrix, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/double/dtbcon.f">http://www.netlib.org/lapack/double/dtbcon.f</a>]
\ingroup vc_lapack_core
*/
void DTBCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STBCON _LAPACK_FUNCTION(STBCON,stbcon)


/** LAPACK STBCON.
Estimates the reciprocal of the condition number of a triangular band matrix, in either the 1-norm or the infinity-norm.
[<a href="http://www.netlib.org/lapack/single/stbcon.f">http://www.netlib.org/lapack/single/stbcon.f</a>]
\ingroup vc_lapack_core
*/
void STBCON(
        const char* _norm,
        const char* _uplo,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORMQR _LAPACK_FUNCTION(DORMQR,dormqr)


/** LAPACK DORMQR.
Multiplies a general matrix by the orthogonal matrix from a QR factorization determined by DGEQRF.
[<a href="http://www.netlib.org/lapack/double/dormqr.f">http://www.netlib.org/lapack/double/dormqr.f</a>]
\ingroup vc_lapack_core
*/
void DORMQR(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMQR _LAPACK_FUNCTION(SORMQR,sormqr)


/** LAPACK SORMQR.
Multiplies a general matrix by the orthogonal matrix from a QR factorization determined by DGEQRF.
[<a href="http://www.netlib.org/lapack/single/sormqr.f">http://www.netlib.org/lapack/single/sormqr.f</a>]
\ingroup vc_lapack_core
*/
void SORMQR(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTREVC _LAPACK_FUNCTION(DTREVC,dtrevc)


/** LAPACK DTREVC.
Computes some or all of the right and/or left eigenvectors of an upper quasi-triangular matrix.
[<a href="http://www.netlib.org/lapack/double/dtrevc.f">http://www.netlib.org/lapack/double/dtrevc.f</a>]
\ingroup vc_lapack_core
*/
void DTREVC(
        const char* _side,
        const char* _howmny,
        ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const double* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define STREVC _LAPACK_FUNCTION(STREVC,strevc)


/** LAPACK STREVC.
Computes some or all of the right and/or left eigenvectors of an upper quasi-triangular matrix.
[<a href="http://www.netlib.org/lapack/single/strevc.f">http://www.netlib.org/lapack/single/strevc.f</a>]
\ingroup vc_lapack_core
*/
void STREVC(
        const char* _side,
        const char* _howmny,
        ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const float* _t,
        ::VC::math::lapack::INTEGER* _ldt,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DTGSYL _LAPACK_FUNCTION(DTGSYL,dtgsyl)


/** LAPACK DTGSYL.

[<a href="http://www.netlib.org/lapack/double/dtgsyl.f">http://www.netlib.org/lapack/double/dtgsyl.f</a>]
\ingroup vc_lapack_core
*/
void DTGSYL(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _ijob,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        const double* _d,
        ::VC::math::lapack::INTEGER* _ldd,
        const double* _e,
        ::VC::math::lapack::INTEGER* _lde,
        double* _f,
        ::VC::math::lapack::INTEGER* _ldf,
        double* _scale,
        double* _dif,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STGSYL _LAPACK_FUNCTION(STGSYL,stgsyl)


/** LAPACK STGSYL.

[<a href="http://www.netlib.org/lapack/single/stgsyl.f">http://www.netlib.org/lapack/single/stgsyl.f</a>]
\ingroup vc_lapack_core
*/
void STGSYL(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _ijob,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        const float* _d,
        ::VC::math::lapack::INTEGER* _ldd,
        const float* _e,
        ::VC::math::lapack::INTEGER* _lde,
        float* _f,
        ::VC::math::lapack::INTEGER* _ldf,
        float* _scale,
        float* _dif,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPBTRS _LAPACK_FUNCTION(DPBTRS,dpbtrs)


/** LAPACK DPBTRS.
Solves a symmetric positive definite banded system of linear equations AX=B, using the Cholesky factorization computed by DPBTRF.
[<a href="http://www.netlib.org/lapack/double/dpbtrs.f">http://www.netlib.org/lapack/double/dpbtrs.f</a>]
\ingroup vc_lapack_core
*/
void DPBTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPBTRS _LAPACK_FUNCTION(SPBTRS,spbtrs)


/** LAPACK SPBTRS.
Solves a symmetric positive definite banded system of linear equations AX=B, using the Cholesky factorization computed by DPBTRF.
[<a href="http://www.netlib.org/lapack/single/spbtrs.f">http://www.netlib.org/lapack/single/spbtrs.f</a>]
\ingroup vc_lapack_core
*/
void SPBTRS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DPBSVX _LAPACK_FUNCTION(DPBSVX,dpbsvx)


/** LAPACK DPBSVX.
Solves a symmetric positive definite banded system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/double/dpbsvx.f">http://www.netlib.org/lapack/double/dpbsvx.f</a>]
\ingroup vc_lapack_core
*/
void DPBSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        char* _equed,
        double* _s,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _rcond,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPBSVX _LAPACK_FUNCTION(SPBSVX,spbsvx)


/** LAPACK SPBSVX.
Solves a symmetric positive definite banded system of linear equations AX=B, and provides an estimate of the condition number and error bounds on the solution.
[<a href="http://www.netlib.org/lapack/single/spbsvx.f">http://www.netlib.org/lapack/single/spbsvx.f</a>]
\ingroup vc_lapack_core
*/
void SPBSVX(
        const char* _fact,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        char* _equed,
        float* _s,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _rcond,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSYGVX _LAPACK_FUNCTION(DSYGVX,dsygvx)


/** LAPACK DSYGVX.
Computes selected eigenvalues, and optionally, the eigenvectors of a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
[<a href="http://www.netlib.org/lapack/double/dsygvx.f">http://www.netlib.org/lapack/double/dsygvx.f</a>]
\ingroup vc_lapack_core
*/
void DSYGVX(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSYGVX _LAPACK_FUNCTION(SSYGVX,ssygvx)


/** LAPACK SSYGVX.
Computes selected eigenvalues, and optionally, the eigenvectors of a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
[<a href="http://www.netlib.org/lapack/single/ssygvx.f">http://www.netlib.org/lapack/single/ssygvx.f</a>]
\ingroup vc_lapack_core
*/
void SSYGVX(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DSBGV _LAPACK_FUNCTION(DSBGV,dsbgv)


/** LAPACK DSBGV.

[<a href="http://www.netlib.org/lapack/double/dsbgv.f">http://www.netlib.org/lapack/double/dsbgv.f</a>]
\ingroup vc_lapack_core
*/
void DSBGV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSBGV _LAPACK_FUNCTION(SSBGV,ssbgv)


/** LAPACK SSBGV.

[<a href="http://www.netlib.org/lapack/single/ssbgv.f">http://www.netlib.org/lapack/single/ssbgv.f</a>]
\ingroup vc_lapack_core
*/
void SSBGV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSPEVX _LAPACK_FUNCTION(DSPEVX,dspevx)


/** LAPACK DSPEVX.
Computes selected eigenvalues and eigenvectors of a symmetric matrix in packed storage.
[<a href="http://www.netlib.org/lapack/double/dspevx.f">http://www.netlib.org/lapack/double/dspevx.f</a>]
\ingroup vc_lapack_core
*/
void DSPEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSPEVX _LAPACK_FUNCTION(SSPEVX,sspevx)


/** LAPACK SSPEVX.
Computes selected eigenvalues and eigenvectors of a symmetric matrix in packed storage.
[<a href="http://www.netlib.org/lapack/single/sspevx.f">http://www.netlib.org/lapack/single/sspevx.f</a>]
\ingroup vc_lapack_core
*/
void SSPEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DGEESX _LAPACK_FUNCTION(DGEESX,dgeesx)


/** LAPACK DGEESX.
Computes the eigenvalues and Schur factorization of a general matrix, orders the factorization so that selected eigenvalues are at the top left of the Schur form, and computes reciprocal condition numbers for the average of the selected eigenvalues, and for the associated right invariant subspace.
[<a href="http://www.netlib.org/lapack/double/dgeesx.f">http://www.netlib.org/lapack/double/dgeesx.f</a>]
\ingroup vc_lapack_core
*/
void DGEESX(
        const char* _jobvs,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _select,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _sdim,
        double* _wr,
        double* _wi,
        double* _vs,
        ::VC::math::lapack::INTEGER* _ldvs,
        double* _rconde,
        double* _rcondv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEESX _LAPACK_FUNCTION(SGEESX,sgeesx)


/** LAPACK SGEESX.
Computes the eigenvalues and Schur factorization of a general matrix, orders the factorization so that selected eigenvalues are at the top left of the Schur form, and computes reciprocal condition numbers for the average of the selected eigenvalues, and for the associated right invariant subspace.
[<a href="http://www.netlib.org/lapack/single/sgeesx.f">http://www.netlib.org/lapack/single/sgeesx.f</a>]
\ingroup vc_lapack_core
*/
void SGEESX(
        const char* _jobvs,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _select,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _sdim,
        float* _wr,
        float* _wi,
        float* _vs,
        ::VC::math::lapack::INTEGER* _ldvs,
        float* _rconde,
        float* _rcondv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBGST _LAPACK_FUNCTION(DSBGST,dsbgst)


/** LAPACK DSBGST.

[<a href="http://www.netlib.org/lapack/double/dsbgst.f">http://www.netlib.org/lapack/double/dsbgst.f</a>]
\ingroup vc_lapack_core
*/
void DSBGST(
        const char* _vect,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const double* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSBGST _LAPACK_FUNCTION(SSBGST,ssbgst)


/** LAPACK SSBGST.

[<a href="http://www.netlib.org/lapack/single/ssbgst.f">http://www.netlib.org/lapack/single/ssbgst.f</a>]
\ingroup vc_lapack_core
*/
void SSBGST(
        const char* _vect,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const float* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DORMTR _LAPACK_FUNCTION(DORMTR,dormtr)


/** LAPACK DORMTR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSYTRD.
[<a href="http://www.netlib.org/lapack/double/dormtr.f">http://www.netlib.org/lapack/double/dormtr.f</a>]
\ingroup vc_lapack_core
*/
void DORMTR(
        const char* _side,
        const char* _uplo,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMTR _LAPACK_FUNCTION(SORMTR,sormtr)


/** LAPACK SORMTR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSYTRD.
[<a href="http://www.netlib.org/lapack/single/sormtr.f">http://www.netlib.org/lapack/single/sormtr.f</a>]
\ingroup vc_lapack_core
*/
void SORMTR(
        const char* _side,
        const char* _uplo,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DORGRQ _LAPACK_FUNCTION(DORGRQ,dorgrq)


/** LAPACK DORGRQ.
Generates all or part of the orthogonal matrix Q from an RQ factorization determined by DGERQF.
[<a href="http://www.netlib.org/lapack/double/dorgrq.f">http://www.netlib.org/lapack/double/dorgrq.f</a>]
\ingroup vc_lapack_core
*/
void DORGRQ(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORGRQ _LAPACK_FUNCTION(SORGRQ,sorgrq)


/** LAPACK SORGRQ.
Generates all or part of the orthogonal matrix Q from an RQ factorization determined by DGERQF.
[<a href="http://www.netlib.org/lapack/single/sorgrq.f">http://www.netlib.org/lapack/single/sorgrq.f</a>]
\ingroup vc_lapack_core
*/
void SORGRQ(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEV _LAPACK_FUNCTION(DSTEV,dstev)


/** LAPACK DSTEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/double/dstev.f">http://www.netlib.org/lapack/double/dstev.f</a>]
\ingroup vc_lapack_core
*/
void DSTEV(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEV _LAPACK_FUNCTION(SSTEV,sstev)


/** LAPACK SSTEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/single/sstev.f">http://www.netlib.org/lapack/single/sstev.f</a>]
\ingroup vc_lapack_core
*/
void SSTEV(
        const char* _jobz,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEVX _LAPACK_FUNCTION(DSTEVX,dstevx)


/** LAPACK DSTEVX.
Computes selected eigenvalues and eigenvectors of a real symmetric tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/double/dstevx.f">http://www.netlib.org/lapack/double/dstevx.f</a>]
\ingroup vc_lapack_core
*/
void DSTEVX(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEVX _LAPACK_FUNCTION(SSTEVX,sstevx)


/** LAPACK SSTEVX.
Computes selected eigenvalues and eigenvectors of a real symmetric tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/single/sstevx.f">http://www.netlib.org/lapack/single/sstevx.f</a>]
\ingroup vc_lapack_core
*/
void SSTEVX(
        const char* _jobz,
        const char* _range,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DTPTRS _LAPACK_FUNCTION(DTPTRS,dtptrs)


/** LAPACK DTPTRS.
Solves a triangular system of linear equations AX=B, A**T X=B or A**H X=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/double/dtptrs.f">http://www.netlib.org/lapack/double/dtptrs.f</a>]
\ingroup vc_lapack_core
*/
void DTPTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define STPTRS _LAPACK_FUNCTION(STPTRS,stptrs)


/** LAPACK STPTRS.
Solves a triangular system of linear equations AX=B, A**T X=B or A**H X=B, where A is held in packed storage.
[<a href="http://www.netlib.org/lapack/single/stptrs.f">http://www.netlib.org/lapack/single/stptrs.f</a>]
\ingroup vc_lapack_core
*/
void STPTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DTGSEN _LAPACK_FUNCTION(DTGSEN,dtgsen)


/** LAPACK DTGSEN.

[<a href="http://www.netlib.org/lapack/double/dtgsen.f">http://www.netlib.org/lapack/double/dtgsen.f</a>]
\ingroup vc_lapack_core
*/
void DTGSEN(
        ::VC::math::lapack::INTEGER* _ijob,
        const ::VC::math::lapack::LOGICAL* _wantq,
        const ::VC::math::lapack::LOGICAL* _wantz,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _m,
        double* _pl,
        double* _pr,
        double* _dif,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define STGSEN _LAPACK_FUNCTION(STGSEN,stgsen)


/** LAPACK STGSEN.

[<a href="http://www.netlib.org/lapack/single/stgsen.f">http://www.netlib.org/lapack/single/stgsen.f</a>]
\ingroup vc_lapack_core
*/
void STGSEN(
        ::VC::math::lapack::INTEGER* _ijob,
        const ::VC::math::lapack::LOGICAL* _wantq,
        const ::VC::math::lapack::LOGICAL* _wantz,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        ::VC::math::lapack::INTEGER* _m,
        float* _pl,
        float* _pr,
        float* _dif,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEBZ _LAPACK_FUNCTION(DSTEBZ,dstebz)


/** LAPACK DSTEBZ.
Computes selected eigenvalues of a real symmetric tridiagonal matrix by bisection.
[<a href="http://www.netlib.org/lapack/double/dstebz.f">http://www.netlib.org/lapack/double/dstebz.f</a>]
\ingroup vc_lapack_core
*/
void DSTEBZ(
        const char* _range,
        const char* _order,
        ::VC::math::lapack::INTEGER* _n,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        const double* _d,
        const double* _e,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _nsplit,
        double* _w,
        ::VC::math::lapack::INTEGER* _iblock,
        ::VC::math::lapack::INTEGER* _isplit,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEBZ _LAPACK_FUNCTION(SSTEBZ,sstebz)


/** LAPACK SSTEBZ.
Computes selected eigenvalues of a real symmetric tridiagonal matrix by bisection.
[<a href="http://www.netlib.org/lapack/single/sstebz.f">http://www.netlib.org/lapack/single/sstebz.f</a>]
\ingroup vc_lapack_core
*/
void SSTEBZ(
        const char* _range,
        const char* _order,
        ::VC::math::lapack::INTEGER* _n,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        const float* _d,
        const float* _e,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _nsplit,
        float* _w,
        ::VC::math::lapack::INTEGER* _iblock,
        ::VC::math::lapack::INTEGER* _isplit,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGESV _LAPACK_FUNCTION(DGESV,dgesv)


/** LAPACK DGESV.
Solves a general system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dgesv.f">http://www.netlib.org/lapack/double/dgesv.f</a>]
\ingroup vc_lapack_core
*/
void DGESV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGESV _LAPACK_FUNCTION(SGESV,sgesv)


/** LAPACK SGESV.
Solves a general system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/sgesv.f">http://www.netlib.org/lapack/single/sgesv.f</a>]
\ingroup vc_lapack_core
*/
void SGESV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGEEVX _LAPACK_FUNCTION(DGEEVX,dgeevx)


/** LAPACK DGEEVX.
Computes the eigenvalues and left and right eigenvectors of a general matrix,  with preliminary balancing of the matrix, and computes reciprocal condition numbers for the eigenvalues and right eigenvectors.
[<a href="http://www.netlib.org/lapack/double/dgeevx.f">http://www.netlib.org/lapack/double/dgeevx.f</a>]
\ingroup vc_lapack_core
*/
void DGEEVX(
        const char* _balanc,
        const char* _jobvl,
        const char* _jobvr,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _wr,
        double* _wi,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _scale,
        double* _abnrm,
        double* _rconde,
        double* _rcondv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEEVX _LAPACK_FUNCTION(SGEEVX,sgeevx)


/** LAPACK SGEEVX.
Computes the eigenvalues and left and right eigenvectors of a general matrix,  with preliminary balancing of the matrix, and computes reciprocal condition numbers for the eigenvalues and right eigenvectors.
[<a href="http://www.netlib.org/lapack/single/sgeevx.f">http://www.netlib.org/lapack/single/sgeevx.f</a>]
\ingroup vc_lapack_core
*/
void SGEEVX(
        const char* _balanc,
        const char* _jobvl,
        const char* _jobvr,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _wr,
        float* _wi,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _scale,
        float* _abnrm,
        float* _rconde,
        float* _rcondv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGESVD _LAPACK_FUNCTION(DGESVD,dgesvd)


/** LAPACK DGESVD.
Computes the singular value decomposition (SVD) of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/double/dgesvd.f">http://www.netlib.org/lapack/double/dgesvd.f</a>]
\ingroup vc_lapack_core
*/
void DGESVD(
        const char* _jobu,
        const char* _jobvt,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _s,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGESVD _LAPACK_FUNCTION(SGESVD,sgesvd)


/** LAPACK SGESVD.
Computes the singular value decomposition (SVD) of a general rectangular matrix.
[<a href="http://www.netlib.org/lapack/single/sgesvd.f">http://www.netlib.org/lapack/single/sgesvd.f</a>]
\ingroup vc_lapack_core
*/
void SGESVD(
        const char* _jobu,
        const char* _jobvt,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _s,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _vt,
        ::VC::math::lapack::INTEGER* _ldvt,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEES _LAPACK_FUNCTION(DGEES,dgees)


/** LAPACK DGEES.
Computes the eigenvalues and Schur factorization of a general matrix, and orders the factorization so that selected eigenvalues are at the top left of the Schur form.
[<a href="http://www.netlib.org/lapack/double/dgees.f">http://www.netlib.org/lapack/double/dgees.f</a>]
\ingroup vc_lapack_core
*/
void DGEES(
        const char* _jobvs,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _select,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _sdim,
        double* _wr,
        double* _wi,
        double* _vs,
        ::VC::math::lapack::INTEGER* _ldvs,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEES _LAPACK_FUNCTION(SGEES,sgees)


/** LAPACK SGEES.
Computes the eigenvalues and Schur factorization of a general matrix, and orders the factorization so that selected eigenvalues are at the top left of the Schur form.
[<a href="http://www.netlib.org/lapack/single/sgees.f">http://www.netlib.org/lapack/single/sgees.f</a>]
\ingroup vc_lapack_core
*/
void SGEES(
        const char* _jobvs,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _select,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _sdim,
        float* _wr,
        float* _wi,
        float* _vs,
        ::VC::math::lapack::INTEGER* _ldvs,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTBRFS _LAPACK_FUNCTION(DTBRFS,dtbrfs)


/** LAPACK DTBRFS.
Provides forward and backward error bounds for the solution of a triangular banded system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/double/dtbrfs.f">http://www.netlib.org/lapack/double/dtbrfs.f</a>]
\ingroup vc_lapack_core
*/
void DTBRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define STBRFS _LAPACK_FUNCTION(STBRFS,stbrfs)


/** LAPACK STBRFS.
Provides forward and backward error bounds for the solution of a triangular banded system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/single/stbrfs.f">http://www.netlib.org/lapack/single/stbrfs.f</a>]
\ingroup vc_lapack_core
*/
void STBRFS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSYTRD _LAPACK_FUNCTION(DSYTRD,dsytrd)


/** LAPACK DSYTRD.
Reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/double/dsytrd.f">http://www.netlib.org/lapack/double/dsytrd.f</a>]
\ingroup vc_lapack_core
*/
void DSYTRD(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _d,
        double* _e,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYTRD _LAPACK_FUNCTION(SSYTRD,ssytrd)


/** LAPACK SSYTRD.
Reduces a symmetric matrix to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/single/ssytrd.f">http://www.netlib.org/lapack/single/ssytrd.f</a>]
\ingroup vc_lapack_core
*/
void SSYTRD(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _d,
        float* _e,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPGVX _LAPACK_FUNCTION(DSPGVX,dspgvx)


/** LAPACK DSPGVX.
Computes selected eigenvalues, and optionally, eigenvectors of a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage.
[<a href="http://www.netlib.org/lapack/double/dspgvx.f">http://www.netlib.org/lapack/double/dspgvx.f</a>]
\ingroup vc_lapack_core
*/
void DSPGVX(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _bp,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSPGVX _LAPACK_FUNCTION(SSPGVX,sspgvx)


/** LAPACK SSPGVX.
Computes selected eigenvalues, and optionally, eigenvectors of a generalized symmetric-definite generalized eigenproblem,  Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed storage.
[<a href="http://www.netlib.org/lapack/single/sspgvx.f">http://www.netlib.org/lapack/single/sspgvx.f</a>]
\ingroup vc_lapack_core
*/
void SSPGVX(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _bp,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DTGEVC _LAPACK_FUNCTION(DTGEVC,dtgevc)


/** LAPACK DTGEVC.
Computes some or all of the right and/or left generalized eigenvectors of a pair of upper triangular matrices.
[<a href="http://www.netlib.org/lapack/double/dtgevc.f">http://www.netlib.org/lapack/double/dtgevc.f</a>]
\ingroup vc_lapack_core
*/
void DTGEVC(
        const char* _side,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const double* _s,
        ::VC::math::lapack::INTEGER* _lds,
        const double* _p,
        ::VC::math::lapack::INTEGER* _ldp,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define STGEVC _LAPACK_FUNCTION(STGEVC,stgevc)


/** LAPACK STGEVC.
Computes some or all of the right and/or left generalized eigenvectors of a pair of upper triangular matrices.
[<a href="http://www.netlib.org/lapack/single/stgevc.f">http://www.netlib.org/lapack/single/stgevc.f</a>]
\ingroup vc_lapack_core
*/
void STGEVC(
        const char* _side,
        const char* _howmny,
        const ::VC::math::lapack::LOGICAL* _select,
        ::VC::math::lapack::INTEGER* _n,
        const float* _s,
        ::VC::math::lapack::INTEGER* _lds,
        const float* _p,
        ::VC::math::lapack::INTEGER* _ldp,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _mm,
        ::VC::math::lapack::INTEGER* _m,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSYCON _LAPACK_FUNCTION(DSYCON,dsycon)


/** LAPACK DSYCON.
Estimates the reciprocal of the condition number of a real symmetric indefinite matrix, using the factorization computed by DSYTRF.
[<a href="http://www.netlib.org/lapack/double/dsycon.f">http://www.netlib.org/lapack/double/dsycon.f</a>]
\ingroup vc_lapack_core
*/
void DSYCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYCON _LAPACK_FUNCTION(SSYCON,ssycon)


/** LAPACK SSYCON.
Estimates the reciprocal of the condition number of a real symmetric indefinite matrix, using the factorization computed by DSYTRF.
[<a href="http://www.netlib.org/lapack/single/ssycon.f">http://www.netlib.org/lapack/single/ssycon.f</a>]
\ingroup vc_lapack_core
*/
void SSYCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSTERF _LAPACK_FUNCTION(DSTERF,dsterf)


/** LAPACK DSTERF.
Computes all eigenvalues of a real symmetric tridiagonal matrix, using a root-free variant of the QL or QR algorithm.
[<a href="http://www.netlib.org/lapack/double/dsterf.f">http://www.netlib.org/lapack/double/dsterf.f</a>]
\ingroup vc_lapack_core
*/
void DSTERF(
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        ::VC::math::lapack::INTEGER* _info);

# define SSTERF _LAPACK_FUNCTION(SSTERF,ssterf)


/** LAPACK SSTERF.
Computes all eigenvalues of a real symmetric tridiagonal matrix, using a root-free variant of the QL or QR algorithm.
[<a href="http://www.netlib.org/lapack/single/ssterf.f">http://www.netlib.org/lapack/single/ssterf.f</a>]
\ingroup vc_lapack_core
*/
void SSTERF(
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEQR _LAPACK_FUNCTION(DSTEQR,dsteqr)


/** LAPACK DSTEQR.
Computes all eigenvalues and eigenvectors of a real symmetric tridiagonal matrix, using the implicit QL or QR algorithm.
[<a href="http://www.netlib.org/lapack/double/dsteqr.f">http://www.netlib.org/lapack/double/dsteqr.f</a>]
\ingroup vc_lapack_core
*/
void DSTEQR(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEQR _LAPACK_FUNCTION(SSTEQR,ssteqr)


/** LAPACK SSTEQR.
Computes all eigenvalues and eigenvectors of a real symmetric tridiagonal matrix, using the implicit QL or QR algorithm.
[<a href="http://www.netlib.org/lapack/single/ssteqr.f">http://www.netlib.org/lapack/single/ssteqr.f</a>]
\ingroup vc_lapack_core
*/
void SSTEQR(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEDC _LAPACK_FUNCTION(DSTEDC,dstedc)


/** LAPACK DSTEDC.

[<a href="http://www.netlib.org/lapack/double/dstedc.f">http://www.netlib.org/lapack/double/dstedc.f</a>]
\ingroup vc_lapack_core
*/
void DSTEDC(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEDC _LAPACK_FUNCTION(SSTEDC,sstedc)


/** LAPACK SSTEDC.

[<a href="http://www.netlib.org/lapack/single/sstedc.f">http://www.netlib.org/lapack/single/sstedc.f</a>]
\ingroup vc_lapack_core
*/
void SSTEDC(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DOPMTR _LAPACK_FUNCTION(DOPMTR,dopmtr)


/** LAPACK DOPMTR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSPTRD.
[<a href="http://www.netlib.org/lapack/double/dopmtr.f">http://www.netlib.org/lapack/double/dopmtr.f</a>]
\ingroup vc_lapack_core
*/
void DOPMTR(
        const char* _side,
        const char* _uplo,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SOPMTR _LAPACK_FUNCTION(SOPMTR,sopmtr)


/** LAPACK SOPMTR.
Multiplies a general matrix by the orthogonal transformation matrix from a reduction to tridiagonal form determined by DSPTRD.
[<a href="http://www.netlib.org/lapack/single/sopmtr.f">http://www.netlib.org/lapack/single/sopmtr.f</a>]
\ingroup vc_lapack_core
*/
void SOPMTR(
        const char* _side,
        const char* _uplo,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGGESX _LAPACK_FUNCTION(DGGESX,dggesx)


/** LAPACK DGGESX.

[<a href="http://www.netlib.org/lapack/double/dggesx.f">http://www.netlib.org/lapack/double/dggesx.f</a>]
\ingroup vc_lapack_core
*/
void DGGESX(
        const char* _jobvsl,
        const char* _jobvsr,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _selctg,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _sdim,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        double* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        double* _rconde,
        double* _rcondv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGESX _LAPACK_FUNCTION(SGGESX,sggesx)


/** LAPACK SGGESX.

[<a href="http://www.netlib.org/lapack/single/sggesx.f">http://www.netlib.org/lapack/single/sggesx.f</a>]
\ingroup vc_lapack_core
*/
void SGGESX(
        const char* _jobvsl,
        const char* _jobvsr,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _selctg,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _sdim,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        float* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        float* _rconde,
        float* _rcondv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBEVX _LAPACK_FUNCTION(DSBEVX,dsbevx)


/** LAPACK DSBEVX.
Computes selected eigenvalues and eigenvectors of a symmetric band matrix.
[<a href="http://www.netlib.org/lapack/double/dsbevx.f">http://www.netlib.org/lapack/double/dsbevx.f</a>]
\ingroup vc_lapack_core
*/
void DSBEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSBEVX _LAPACK_FUNCTION(SSBEVX,ssbevx)


/** LAPACK SSBEVX.
Computes selected eigenvalues and eigenvectors of a symmetric band matrix.
[<a href="http://www.netlib.org/lapack/single/ssbevx.f">http://www.netlib.org/lapack/single/ssbevx.f</a>]
\ingroup vc_lapack_core
*/
void SSBEVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DSYTRF _LAPACK_FUNCTION(DSYTRF,dsytrf)


/** LAPACK DSYTRF.
Computes the factorization of a real symmetric-indefinite matrix, using the diagonal pivoting method.
[<a href="http://www.netlib.org/lapack/double/dsytrf.f">http://www.netlib.org/lapack/double/dsytrf.f</a>]
\ingroup vc_lapack_core
*/
void DSYTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYTRF _LAPACK_FUNCTION(SSYTRF,ssytrf)


/** LAPACK SSYTRF.
Computes the factorization of a real symmetric-indefinite matrix, using the diagonal pivoting method.
[<a href="http://www.netlib.org/lapack/single/ssytrf.f">http://www.netlib.org/lapack/single/ssytrf.f</a>]
\ingroup vc_lapack_core
*/
void SSYTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGTCON _LAPACK_FUNCTION(DGTCON,dgtcon)


/** LAPACK DGTCON.
Estimates the reciprocal of the condition number of a general tridiagonal matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGTTRF.
[<a href="http://www.netlib.org/lapack/double/dgtcon.f">http://www.netlib.org/lapack/double/dgtcon.f</a>]
\ingroup vc_lapack_core
*/
void DGTCON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        const double* _dl,
        const double* _d,
        const double* _du,
        const double* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGTCON _LAPACK_FUNCTION(SGTCON,sgtcon)


/** LAPACK SGTCON.
Estimates the reciprocal of the condition number of a general tridiagonal matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGTTRF.
[<a href="http://www.netlib.org/lapack/single/sgtcon.f">http://www.netlib.org/lapack/single/sgtcon.f</a>]
\ingroup vc_lapack_core
*/
void SGTCON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        const float* _dl,
        const float* _d,
        const float* _du,
        const float* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEEV _LAPACK_FUNCTION(DGEEV,dgeev)


/** LAPACK DGEEV.
Computes the eigenvalues and left and right eigenvectors of a general matrix.
[<a href="http://www.netlib.org/lapack/double/dgeev.f">http://www.netlib.org/lapack/double/dgeev.f</a>]
\ingroup vc_lapack_core
*/
void DGEEV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _wr,
        double* _wi,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGEEV _LAPACK_FUNCTION(SGEEV,sgeev)


/** LAPACK SGEEV.
Computes the eigenvalues and left and right eigenvectors of a general matrix.
[<a href="http://www.netlib.org/lapack/single/sgeev.f">http://www.netlib.org/lapack/single/sgeev.f</a>]
\ingroup vc_lapack_core
*/
void SGEEV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _wr,
        float* _wi,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGES _LAPACK_FUNCTION(DGGES,dgges)


/** LAPACK DGGES.

[<a href="http://www.netlib.org/lapack/double/dgges.f">http://www.netlib.org/lapack/double/dgges.f</a>]
\ingroup vc_lapack_core
*/
void DGGES(
        const char* _jobvsl,
        const char* _jobvsr,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _selctg,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _sdim,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        double* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGES _LAPACK_FUNCTION(SGGES,sgges)


/** LAPACK SGGES.

[<a href="http://www.netlib.org/lapack/single/sgges.f">http://www.netlib.org/lapack/single/sgges.f</a>]
\ingroup vc_lapack_core
*/
void SGGES(
        const char* _jobvsl,
        const char* _jobvsr,
        const char* _sort,
        ::VC::math::lapack::EXTERNAL_PROCEDURE _selctg,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _sdim,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vsl,
        ::VC::math::lapack::INTEGER* _ldvsl,
        float* _vsr,
        ::VC::math::lapack::INTEGER* _ldvsr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTRTRS _LAPACK_FUNCTION(DTRTRS,dtrtrs)


/** LAPACK DTRTRS.
Solves a triangular system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/double/dtrtrs.f">http://www.netlib.org/lapack/double/dtrtrs.f</a>]
\ingroup vc_lapack_core
*/
void DTRTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define STRTRS _LAPACK_FUNCTION(STRTRS,strtrs)


/** LAPACK STRTRS.
Solves a triangular system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/single/strtrs.f">http://www.netlib.org/lapack/single/strtrs.f</a>]
\ingroup vc_lapack_core
*/
void STRTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSPTRD _LAPACK_FUNCTION(DSPTRD,dsptrd)


/** LAPACK DSPTRD.
Reduces a symmetric matrix in packed storage to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/double/dsptrd.f">http://www.netlib.org/lapack/double/dsptrd.f</a>]
\ingroup vc_lapack_core
*/
void DSPTRD(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _d,
        double* _e,
        double* _tau,
        ::VC::math::lapack::INTEGER* _info);

# define SSPTRD _LAPACK_FUNCTION(SSPTRD,ssptrd)


/** LAPACK SSPTRD.
Reduces a symmetric matrix in packed storage to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/single/ssptrd.f">http://www.netlib.org/lapack/single/ssptrd.f</a>]
\ingroup vc_lapack_core
*/
void SSPTRD(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _d,
        float* _e,
        float* _tau,
        ::VC::math::lapack::INTEGER* _info);


# define DGGSVD _LAPACK_FUNCTION(DGGSVD,dggsvd)


/** LAPACK DGGSVD.

[<a href="http://www.netlib.org/lapack/double/dggsvd.f">http://www.netlib.org/lapack/double/dggsvd.f</a>]
\ingroup vc_lapack_core
*/
void DGGSVD(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alpha,
        double* _beta,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGSVD _LAPACK_FUNCTION(SGGSVD,sggsvd)


/** LAPACK SGGSVD.

[<a href="http://www.netlib.org/lapack/single/sggsvd.f">http://www.netlib.org/lapack/single/sggsvd.f</a>]
\ingroup vc_lapack_core
*/
void SGGSVD(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alpha,
        float* _beta,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGTSV _LAPACK_FUNCTION(DGTSV,dgtsv)


/** LAPACK DGTSV.
Solves a general tridiagonal system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dgtsv.f">http://www.netlib.org/lapack/double/dgtsv.f</a>]
\ingroup vc_lapack_core
*/
void DGTSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _dl,
        double* _d,
        double* _du,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGTSV _LAPACK_FUNCTION(SGTSV,sgtsv)


/** LAPACK SGTSV.
Solves a general tridiagonal system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/sgtsv.f">http://www.netlib.org/lapack/single/sgtsv.f</a>]
\ingroup vc_lapack_core
*/
void SGTSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _dl,
        float* _d,
        float* _du,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGTTRF _LAPACK_FUNCTION(DGTTRF,dgttrf)


/** LAPACK DGTTRF.
Computes an LU factorization of a general tridiagonal matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/double/dgttrf.f">http://www.netlib.org/lapack/double/dgttrf.f</a>]
\ingroup vc_lapack_core
*/
void DGTTRF(
        ::VC::math::lapack::INTEGER* _n,
        double* _dl,
        double* _d,
        double* _du,
        double* _du2,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);

# define SGTTRF _LAPACK_FUNCTION(SGTTRF,sgttrf)


/** LAPACK SGTTRF.
Computes an LU factorization of a general tridiagonal matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/single/sgttrf.f">http://www.netlib.org/lapack/single/sgttrf.f</a>]
\ingroup vc_lapack_core
*/
void SGTTRF(
        ::VC::math::lapack::INTEGER* _n,
        float* _dl,
        float* _d,
        float* _du,
        float* _du2,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);


# define DGGEVX _LAPACK_FUNCTION(DGGEVX,dggevx)


/** LAPACK DGGEVX.

[<a href="http://www.netlib.org/lapack/double/dggevx.f">http://www.netlib.org/lapack/double/dggevx.f</a>]
\ingroup vc_lapack_core
*/
void DGGEVX(
        const char* _balanc,
        const char* _jobvl,
        const char* _jobvr,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        double* _lscale,
        double* _rscale,
        double* _abnrm,
        double* _bbnrm,
        double* _rconde,
        double* _rcondv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGEVX _LAPACK_FUNCTION(SGGEVX,sggevx)


/** LAPACK SGGEVX.

[<a href="http://www.netlib.org/lapack/single/sggevx.f">http://www.netlib.org/lapack/single/sggevx.f</a>]
\ingroup vc_lapack_core
*/
void SGGEVX(
        const char* _balanc,
        const char* _jobvl,
        const char* _jobvr,
        const char* _sense,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        ::VC::math::lapack::INTEGER* _ilo,
        ::VC::math::lapack::INTEGER* _ihi,
        float* _lscale,
        float* _rscale,
        float* _abnrm,
        float* _bbnrm,
        float* _rconde,
        float* _rcondv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::LOGICAL* _bwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSTEIN _LAPACK_FUNCTION(DSTEIN,dstein)


/** LAPACK DSTEIN.
Computes selected eigenvectors of a real symmetric tridiagonal matrix by inverse iteration.
[<a href="http://www.netlib.org/lapack/double/dstein.f">http://www.netlib.org/lapack/double/dstein.f</a>]
\ingroup vc_lapack_core
*/
void DSTEIN(
        ::VC::math::lapack::INTEGER* _n,
        const double* _d,
        const double* _e,
        ::VC::math::lapack::INTEGER* _m,
        const double* _w,
        const ::VC::math::lapack::INTEGER* _iblock,
        const ::VC::math::lapack::INTEGER* _isplit,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSTEIN _LAPACK_FUNCTION(SSTEIN,sstein)


/** LAPACK SSTEIN.
Computes selected eigenvectors of a real symmetric tridiagonal matrix by inverse iteration.
[<a href="http://www.netlib.org/lapack/single/sstein.f">http://www.netlib.org/lapack/single/sstein.f</a>]
\ingroup vc_lapack_core
*/
void SSTEIN(
        ::VC::math::lapack::INTEGER* _n,
        const float* _d,
        const float* _e,
        ::VC::math::lapack::INTEGER* _m,
        const float* _w,
        const ::VC::math::lapack::INTEGER* _iblock,
        const ::VC::math::lapack::INTEGER* _isplit,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DSPCON _LAPACK_FUNCTION(DSPCON,dspcon)


/** LAPACK DSPCON.
Estimates the reciprocal of the condition number of a real symmetric indefinite matrix in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/double/dspcon.f">http://www.netlib.org/lapack/double/dspcon.f</a>]
\ingroup vc_lapack_core
*/
void DSPCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSPCON _LAPACK_FUNCTION(SSPCON,sspcon)


/** LAPACK SSPCON.
Estimates the reciprocal of the condition number of a real symmetric indefinite matrix in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/single/sspcon.f">http://www.netlib.org/lapack/single/sspcon.f</a>]
\ingroup vc_lapack_core
*/
void SSPCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBGVX _LAPACK_FUNCTION(DSBGVX,dsbgvx)


/** LAPACK DSBGVX.

[<a href="http://www.netlib.org/lapack/double/dsbgvx.f">http://www.netlib.org/lapack/double/dsbgvx.f</a>]
\ingroup vc_lapack_core
*/
void DSBGVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        const double* _vl,
        const double* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const double* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);

# define SSBGVX _LAPACK_FUNCTION(SSBGVX,ssbgvx)


/** LAPACK SSBGVX.

[<a href="http://www.netlib.org/lapack/single/ssbgvx.f">http://www.netlib.org/lapack/single/ssbgvx.f</a>]
\ingroup vc_lapack_core
*/
void SSBGVX(
        const char* _jobz,
        const char* _range,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _ka,
        ::VC::math::lapack::INTEGER* _kb,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _bb,
        ::VC::math::lapack::INTEGER* _ldbb,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        const float* _vl,
        const float* _vu,
        ::VC::math::lapack::INTEGER* _il,
        ::VC::math::lapack::INTEGER* _iu,
        const float* _abstol,
        ::VC::math::lapack::INTEGER* _m,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _ifail,
        ::VC::math::lapack::INTEGER* _info);


# define DORMLQ _LAPACK_FUNCTION(DORMLQ,dormlq)


/** LAPACK DORMLQ.
Multiplies a general matrix by the orthogonal matrix from an LQ factorization determined by DGELQF.
[<a href="http://www.netlib.org/lapack/double/dormlq.f">http://www.netlib.org/lapack/double/dormlq.f</a>]
\ingroup vc_lapack_core
*/
void DORMLQ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMLQ _LAPACK_FUNCTION(SORMLQ,sormlq)


/** LAPACK SORMLQ.
Multiplies a general matrix by the orthogonal matrix from an LQ factorization determined by DGELQF.
[<a href="http://www.netlib.org/lapack/single/sormlq.f">http://www.netlib.org/lapack/single/sormlq.f</a>]
\ingroup vc_lapack_core
*/
void SORMLQ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPBSV _LAPACK_FUNCTION(DPBSV,dpbsv)


/** LAPACK DPBSV.
Solves a symmetric positive definite banded system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dpbsv.f">http://www.netlib.org/lapack/double/dpbsv.f</a>]
\ingroup vc_lapack_core
*/
void DPBSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPBSV _LAPACK_FUNCTION(SPBSV,spbsv)


/** LAPACK SPBSV.
Solves a symmetric positive definite banded system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/spbsv.f">http://www.netlib.org/lapack/single/spbsv.f</a>]
\ingroup vc_lapack_core
*/
void SPBSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSYRFS _LAPACK_FUNCTION(DSYRFS,dsyrfs)


/** LAPACK DSYRFS.
Improves the computed solution to a real symmetric indefinite system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dsyrfs.f">http://www.netlib.org/lapack/double/dsyrfs.f</a>]
\ingroup vc_lapack_core
*/
void DSYRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYRFS _LAPACK_FUNCTION(SSYRFS,ssyrfs)


/** LAPACK SSYRFS.
Improves the computed solution to a real symmetric indefinite system of linear equations AX=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/ssyrfs.f">http://www.netlib.org/lapack/single/ssyrfs.f</a>]
\ingroup vc_lapack_core
*/
void SSYRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPTRF _LAPACK_FUNCTION(DSPTRF,dsptrf)


/** LAPACK DSPTRF.
Computes the factorization of a real symmetric-indefinite matrix in packed storage, using the diagonal pivoting method.
[<a href="http://www.netlib.org/lapack/double/dsptrf.f">http://www.netlib.org/lapack/double/dsptrf.f</a>]
\ingroup vc_lapack_core
*/
void DSPTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);

# define SSPTRF _LAPACK_FUNCTION(SSPTRF,ssptrf)


/** LAPACK SSPTRF.
Computes the factorization of a real symmetric-indefinite matrix in packed storage, using the diagonal pivoting method.
[<a href="http://www.netlib.org/lapack/single/ssptrf.f">http://www.netlib.org/lapack/single/ssptrf.f</a>]
\ingroup vc_lapack_core
*/
void SSPTRF(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);


# define DGELSD _LAPACK_FUNCTION(DGELSD,dgelsd)


/** LAPACK DGELSD.
Computes the least squares solution to an over-determined system of linear equations, A X=B or A**H X=B,  or the minimum norm solution of an under-determined system, using a divide and conquer method, where A is a general rectangular matrix of full rank, using a QR or LQ factorization of A.
[<a href="http://www.netlib.org/lapack/double/dgelsd.f">http://www.netlib.org/lapack/double/dgelsd.f</a>]
\ingroup vc_lapack_core
*/
void DGELSD(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _s,
        const double* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGELSD _LAPACK_FUNCTION(SGELSD,sgelsd)


/** LAPACK SGELSD.
Computes the least squares solution to an over-determined system of linear equations, A X=B or A**H X=B,  or the minimum norm solution of an under-determined system, using a divide and conquer method, where A is a general rectangular matrix of full rank, using a QR or LQ factorization of A.
[<a href="http://www.netlib.org/lapack/single/sgelsd.f">http://www.netlib.org/lapack/single/sgelsd.f</a>]
\ingroup vc_lapack_core
*/
void SGELSD(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _s,
        const float* _rcond,
        ::VC::math::lapack::INTEGER* _rank,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGEV _LAPACK_FUNCTION(DGGEV,dggev)


/** LAPACK DGGEV.

[<a href="http://www.netlib.org/lapack/double/dggev.f">http://www.netlib.org/lapack/double/dggev.f</a>]
\ingroup vc_lapack_core
*/
void DGGEV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _alphar,
        double* _alphai,
        double* _beta,
        double* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        double* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGEV _LAPACK_FUNCTION(SGGEV,sggev)


/** LAPACK SGGEV.

[<a href="http://www.netlib.org/lapack/single/sggev.f">http://www.netlib.org/lapack/single/sggev.f</a>]
\ingroup vc_lapack_core
*/
void SGGEV(
        const char* _jobvl,
        const char* _jobvr,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _alphar,
        float* _alphai,
        float* _beta,
        float* _vl,
        ::VC::math::lapack::INTEGER* _ldvl,
        float* _vr,
        ::VC::math::lapack::INTEGER* _ldvr,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGEEQU _LAPACK_FUNCTION(DGEEQU,dgeequ)


/** LAPACK DGEEQU.
Computes row and column scalings to equilibrate a general rectangular matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/double/dgeequ.f">http://www.netlib.org/lapack/double/dgeequ.f</a>]
\ingroup vc_lapack_core
*/
void DGEEQU(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _r,
        double* _c,
        double* _rowcnd,
        double* _colcnd,
        double* _amax,
        ::VC::math::lapack::INTEGER* _info);

# define SGEEQU _LAPACK_FUNCTION(SGEEQU,sgeequ)


/** LAPACK SGEEQU.
Computes row and column scalings to equilibrate a general rectangular matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/single/sgeequ.f">http://www.netlib.org/lapack/single/sgeequ.f</a>]
\ingroup vc_lapack_core
*/
void SGEEQU(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _r,
        float* _c,
        float* _rowcnd,
        float* _colcnd,
        float* _amax,
        ::VC::math::lapack::INTEGER* _info);


# define DSYTRI _LAPACK_FUNCTION(DSYTRI,dsytri)


/** LAPACK DSYTRI.
Computes the inverse of a real symmetric indefinite matrix, using the factorization computed by DSYTRF.
[<a href="http://www.netlib.org/lapack/double/dsytri.f">http://www.netlib.org/lapack/double/dsytri.f</a>]
\ingroup vc_lapack_core
*/
void DSYTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSYTRI _LAPACK_FUNCTION(SSYTRI,ssytri)


/** LAPACK SSYTRI.
Computes the inverse of a real symmetric indefinite matrix, using the factorization computed by DSYTRF.
[<a href="http://www.netlib.org/lapack/single/ssytri.f">http://www.netlib.org/lapack/single/ssytri.f</a>]
\ingroup vc_lapack_core
*/
void SSYTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DPTEQR _LAPACK_FUNCTION(DPTEQR,dpteqr)


/** LAPACK DPTEQR.
Computes all eigenvalues and eigenvectors of a real symmetric positive definite tridiagonal matrix, by computing the SVD of its bidiagonal Cholesky factor.
[<a href="http://www.netlib.org/lapack/double/dpteqr.f">http://www.netlib.org/lapack/double/dpteqr.f</a>]
\ingroup vc_lapack_core
*/
void DPTEQR(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SPTEQR _LAPACK_FUNCTION(SPTEQR,spteqr)


/** LAPACK SPTEQR.
Computes all eigenvalues and eigenvectors of a real symmetric positive definite tridiagonal matrix, by computing the SVD of its bidiagonal Cholesky factor.
[<a href="http://www.netlib.org/lapack/single/spteqr.f">http://www.netlib.org/lapack/single/spteqr.f</a>]
\ingroup vc_lapack_core
*/
void SPTEQR(
        const char* _compz,
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DPOSV _LAPACK_FUNCTION(DPOSV,dposv)


/** LAPACK DPOSV.
Solves a symmetric positive definite system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dposv.f">http://www.netlib.org/lapack/double/dposv.f</a>]
\ingroup vc_lapack_core
*/
void DPOSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SPOSV _LAPACK_FUNCTION(SPOSV,sposv)


/** LAPACK SPOSV.
Solves a symmetric positive definite system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/sposv.f">http://www.netlib.org/lapack/single/sposv.f</a>]
\ingroup vc_lapack_core
*/
void SPOSV(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DGECON _LAPACK_FUNCTION(DGECON,dgecon)


/** LAPACK DGECON.
Estimates the reciprocal of the condition number of a general matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/double/dgecon.f">http://www.netlib.org/lapack/double/dgecon.f</a>]
\ingroup vc_lapack_core
*/
void DGECON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGECON _LAPACK_FUNCTION(SGECON,sgecon)


/** LAPACK SGECON.
Estimates the reciprocal of the condition number of a general matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/single/sgecon.f">http://www.netlib.org/lapack/single/sgecon.f</a>]
\ingroup vc_lapack_core
*/
void SGECON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGTRFS _LAPACK_FUNCTION(DGTRFS,dgtrfs)


/** LAPACK DGTRFS.
Improves the computed solution to a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dgtrfs.f">http://www.netlib.org/lapack/double/dgtrfs.f</a>]
\ingroup vc_lapack_core
*/
void DGTRFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _dl,
        const double* _d,
        const double* _du,
        const double* _dlf,
        const double* _df,
        const double* _duf,
        const double* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGTRFS _LAPACK_FUNCTION(SGTRFS,sgtrfs)


/** LAPACK SGTRFS.
Improves the computed solution to a general tridiagonal system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/sgtrfs.f">http://www.netlib.org/lapack/single/sgtrfs.f</a>]
\ingroup vc_lapack_core
*/
void SGTRFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _dl,
        const float* _d,
        const float* _du,
        const float* _dlf,
        const float* _df,
        const float* _duf,
        const float* _du2,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTZRQF _LAPACK_FUNCTION(DTZRQF,dtzrqf)


/** LAPACK DTZRQF.
Computes an RQ factorization of an upper trapezoidal matrix.
[<a href="http://www.netlib.org/lapack/double/dtzrqf.f">http://www.netlib.org/lapack/double/dtzrqf.f</a>]
\ingroup vc_lapack_core
*/
void DTZRQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        ::VC::math::lapack::INTEGER* _info);

# define STZRQF _LAPACK_FUNCTION(STZRQF,stzrqf)


/** LAPACK STZRQF.
Computes an RQ factorization of an upper trapezoidal matrix.
[<a href="http://www.netlib.org/lapack/single/stzrqf.f">http://www.netlib.org/lapack/single/stzrqf.f</a>]
\ingroup vc_lapack_core
*/
void STZRQF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        ::VC::math::lapack::INTEGER* _info);


# define DSPRFS _LAPACK_FUNCTION(DSPRFS,dsprfs)


/** LAPACK DSPRFS.
Improves the computed solution to a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dsprfs.f">http://www.netlib.org/lapack/double/dsprfs.f</a>]
\ingroup vc_lapack_core
*/
void DSPRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ap,
        const double* _afp,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSPRFS _LAPACK_FUNCTION(SSPRFS,ssprfs)


/** LAPACK SSPRFS.
Improves the computed solution to a real symmetric indefinite system of linear equations AX=B, where A is held in packed storage, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/ssprfs.f">http://www.netlib.org/lapack/single/ssprfs.f</a>]
\ingroup vc_lapack_core
*/
void SSPRFS(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ap,
        const float* _afp,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBTRD _LAPACK_FUNCTION(DSBTRD,dsbtrd)


/** LAPACK DSBTRD.
Reduces a symmetric band matrix to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/double/dsbtrd.f">http://www.netlib.org/lapack/double/dsbtrd.f</a>]
\ingroup vc_lapack_core
*/
void DSBTRD(
        const char* _vect,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _d,
        double* _e,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSBTRD _LAPACK_FUNCTION(SSBTRD,ssbtrd)


/** LAPACK SSBTRD.
Reduces a symmetric band matrix to real symmetric tridiagonal form by an orthogonal similarity transformation.
[<a href="http://www.netlib.org/lapack/single/ssbtrd.f">http://www.netlib.org/lapack/single/ssbtrd.f</a>]
\ingroup vc_lapack_core
*/
void SSBTRD(
        const char* _vect,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _d,
        float* _e,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGBTRF _LAPACK_FUNCTION(DGBTRF,dgbtrf)


/** LAPACK DGBTRF.
Computes an LU factorization of a general band matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/double/dgbtrf.f">http://www.netlib.org/lapack/double/dgbtrf.f</a>]
\ingroup vc_lapack_core
*/
void DGBTRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);

# define SGBTRF _LAPACK_FUNCTION(SGBTRF,sgbtrf)


/** LAPACK SGBTRF.
Computes an LU factorization of a general band matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/single/sgbtrf.f">http://www.netlib.org/lapack/single/sgbtrf.f</a>]
\ingroup vc_lapack_core
*/
void SGBTRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);


# define DGETRF _LAPACK_FUNCTION(DGETRF,dgetrf)


/** LAPACK DGETRF.
Computes an LU factorization of a general matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/double/dgetrf.f">http://www.netlib.org/lapack/double/dgetrf.f</a>]
\ingroup vc_lapack_core
*/
void DGETRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);

# define SGETRF _LAPACK_FUNCTION(SGETRF,sgetrf)


/** LAPACK SGETRF.
Computes an LU factorization of a general matrix, using partial pivoting with row interchanges.
[<a href="http://www.netlib.org/lapack/single/sgetrf.f">http://www.netlib.org/lapack/single/sgetrf.f</a>]
\ingroup vc_lapack_core
*/
void SGETRF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        ::VC::math::lapack::INTEGER* _ipiv,
        ::VC::math::lapack::INTEGER* _info);


# define DSPTRI _LAPACK_FUNCTION(DSPTRI,dsptri)


/** LAPACK DSPTRI.
Computes the inverse of a real symmetric indefinite matrix in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/double/dsptri.f">http://www.netlib.org/lapack/double/dsptri.f</a>]
\ingroup vc_lapack_core
*/
void DSPTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSPTRI _LAPACK_FUNCTION(SSPTRI,ssptri)


/** LAPACK SSPTRI.
Computes the inverse of a real symmetric indefinite matrix in packed storage, using the factorization computed by DSPTRF.
[<a href="http://www.netlib.org/lapack/single/ssptri.f">http://www.netlib.org/lapack/single/ssptri.f</a>]
\ingroup vc_lapack_core
*/
void SSPTRI(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DORMRQ _LAPACK_FUNCTION(DORMRQ,dormrq)


/** LAPACK DORMRQ.
Multiplies a general matrix by the orthogonal matrix from an RQ factorization determined by DGERQF.
[<a href="http://www.netlib.org/lapack/double/dormrq.f">http://www.netlib.org/lapack/double/dormrq.f</a>]
\ingroup vc_lapack_core
*/
void DORMRQ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _tau,
        double* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SORMRQ _LAPACK_FUNCTION(SORMRQ,sormrq)


/** LAPACK SORMRQ.
Multiplies a general matrix by the orthogonal matrix from an RQ factorization determined by DGERQF.
[<a href="http://www.netlib.org/lapack/single/sormrq.f">http://www.netlib.org/lapack/single/sormrq.f</a>]
\ingroup vc_lapack_core
*/
void SORMRQ(
        const char* _side,
        const char* _trans,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _tau,
        float* _c,
        ::VC::math::lapack::INTEGER* _ldc,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSBEV _LAPACK_FUNCTION(DSBEV,dsbev)


/** LAPACK DSBEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric band matrix.
[<a href="http://www.netlib.org/lapack/double/dsbev.f">http://www.netlib.org/lapack/double/dsbev.f</a>]
\ingroup vc_lapack_core
*/
void DSBEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SSBEV _LAPACK_FUNCTION(SSBEV,ssbev)


/** LAPACK SSBEV.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric band matrix.
[<a href="http://www.netlib.org/lapack/single/ssbev.f">http://www.netlib.org/lapack/single/ssbev.f</a>]
\ingroup vc_lapack_core
*/
void SSBEV(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGERFS _LAPACK_FUNCTION(DGERFS,dgerfs)


/** LAPACK DGERFS.
Improves the computed solution to a general system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dgerfs.f">http://www.netlib.org/lapack/double/dgerfs.f</a>]
\ingroup vc_lapack_core
*/
void DGERFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGERFS _LAPACK_FUNCTION(SGERFS,sgerfs)


/** LAPACK SGERFS.
Improves the computed solution to a general system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/sgerfs.f">http://www.netlib.org/lapack/single/sgerfs.f</a>]
\ingroup vc_lapack_core
*/
void SGERFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _af,
        ::VC::math::lapack::INTEGER* _ldaf,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPTCON _LAPACK_FUNCTION(DPTCON,dptcon)


/** LAPACK DPTCON.
Computes the reciprocal of the condition number of a symmetric positive definite tridiagonal matrix, using the LDL**H factorization computed by DPTTRF.
[<a href="http://www.netlib.org/lapack/double/dptcon.f">http://www.netlib.org/lapack/double/dptcon.f</a>]
\ingroup vc_lapack_core
*/
void DPTCON(
        ::VC::math::lapack::INTEGER* _n,
        const double* _d,
        const double* _e,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _info);

# define SPTCON _LAPACK_FUNCTION(SPTCON,sptcon)


/** LAPACK SPTCON.
Computes the reciprocal of the condition number of a symmetric positive definite tridiagonal matrix, using the LDL**H factorization computed by DPTTRF.
[<a href="http://www.netlib.org/lapack/single/sptcon.f">http://www.netlib.org/lapack/single/sptcon.f</a>]
\ingroup vc_lapack_core
*/
void SPTCON(
        ::VC::math::lapack::INTEGER* _n,
        const float* _d,
        const float* _e,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _info);


# define DGETRI _LAPACK_FUNCTION(DGETRI,dgetri)


/** LAPACK DGETRI.
Computes the inverse of a general matrix, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/double/dgetri.f">http://www.netlib.org/lapack/double/dgetri.f</a>]
\ingroup vc_lapack_core
*/
void DGETRI(
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGETRI _LAPACK_FUNCTION(SGETRI,sgetri)


/** LAPACK SGETRI.
Computes the inverse of a general matrix, using the LU factorization computed by DGETRF.
[<a href="http://www.netlib.org/lapack/single/sgetri.f">http://www.netlib.org/lapack/single/sgetri.f</a>]
\ingroup vc_lapack_core
*/
void SGETRI(
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const ::VC::math::lapack::INTEGER* _ipiv,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGBRFS _LAPACK_FUNCTION(DGBRFS,dgbrfs)


/** LAPACK DGBRFS.
Improves the computed solution to a general banded system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/double/dgbrfs.f">http://www.netlib.org/lapack/double/dgbrfs.f</a>]
\ingroup vc_lapack_core
*/
void DGBRFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const double* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        double* _ferr,
        double* _berr,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGBRFS _LAPACK_FUNCTION(SGBRFS,sgbrfs)


/** LAPACK SGBRFS.
Improves the computed solution to a general banded system of linear equations AX=B, A**T X=B or A**H X=B, and provides forward and backward error bounds for the solution.
[<a href="http://www.netlib.org/lapack/single/sgbrfs.f">http://www.netlib.org/lapack/single/sgbrfs.f</a>]
\ingroup vc_lapack_core
*/
void SGBRFS(
        const char* _trans,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const float* _afb,
        ::VC::math::lapack::INTEGER* _ldafb,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _x,
        ::VC::math::lapack::INTEGER* _ldx,
        float* _ferr,
        float* _berr,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGLSE _LAPACK_FUNCTION(DGGLSE,dgglse)


/** LAPACK DGGLSE.
Solves the LSE (Constrained Linear Least Squares Problem) using the GRQ (Generalized RQ) factorization
[<a href="http://www.netlib.org/lapack/double/dgglse.f">http://www.netlib.org/lapack/double/dgglse.f</a>]
\ingroup vc_lapack_core
*/
void DGGLSE(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _p,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _c,
        double* _d,
        double* _x,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGLSE _LAPACK_FUNCTION(SGGLSE,sgglse)


/** LAPACK SGGLSE.
Solves the LSE (Constrained Linear Least Squares Problem) using the GRQ (Generalized RQ) factorization
[<a href="http://www.netlib.org/lapack/single/sgglse.f">http://www.netlib.org/lapack/single/sgglse.f</a>]
\ingroup vc_lapack_core
*/
void SGGLSE(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _p,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _c,
        float* _d,
        float* _x,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DDISNA _LAPACK_FUNCTION(DDISNA,ddisna)


/** LAPACK DDISNA.

[<a href="http://www.netlib.org/lapack/double/ddisna.f">http://www.netlib.org/lapack/double/ddisna.f</a>]
\ingroup vc_lapack_core
*/
void DDISNA(
        const char* _job,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const double* _d,
        double* _sep,
        ::VC::math::lapack::INTEGER* _info);

# define SDISNA _LAPACK_FUNCTION(SDISNA,sdisna)


/** LAPACK SDISNA.

[<a href="http://www.netlib.org/lapack/single/sdisna.f">http://www.netlib.org/lapack/single/sdisna.f</a>]
\ingroup vc_lapack_core
*/
void SDISNA(
        const char* _job,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        const float* _d,
        float* _sep,
        ::VC::math::lapack::INTEGER* _info);


# define DTBTRS _LAPACK_FUNCTION(DTBTRS,dtbtrs)


/** LAPACK DTBTRS.
Solves a triangular banded system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/double/dtbtrs.f">http://www.netlib.org/lapack/double/dtbtrs.f</a>]
\ingroup vc_lapack_core
*/
void DTBTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define STBTRS _LAPACK_FUNCTION(STBTRS,stbtrs)


/** LAPACK STBTRS.
Solves a triangular banded system of linear equations AX=B, A**T X=B or A**H X=B.
[<a href="http://www.netlib.org/lapack/single/stbtrs.f">http://www.netlib.org/lapack/single/stbtrs.f</a>]
\ingroup vc_lapack_core
*/
void STBTRS(
        const char* _uplo,
        const char* _trans,
        const char* _diag,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kd,
        ::VC::math::lapack::INTEGER* _nrhs,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSYEVD _LAPACK_FUNCTION(DSYEVD,dsyevd)


/** LAPACK DSYEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dsyevd.f">http://www.netlib.org/lapack/double/dsyevd.f</a>]
\ingroup vc_lapack_core
*/
void DSYEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _w,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYEVD _LAPACK_FUNCTION(SSYEVD,ssyevd)


/** LAPACK SSYEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/ssyevd.f">http://www.netlib.org/lapack/single/ssyevd.f</a>]
\ingroup vc_lapack_core
*/
void SSYEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _w,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DPTTRF _LAPACK_FUNCTION(DPTTRF,dpttrf)


/** LAPACK DPTTRF.
Computes the LDL**H factorization of a symmetric positive definite tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/double/dpttrf.f">http://www.netlib.org/lapack/double/dpttrf.f</a>]
\ingroup vc_lapack_core
*/
void DPTTRF(
        ::VC::math::lapack::INTEGER* _n,
        double* _d,
        double* _e,
        ::VC::math::lapack::INTEGER* _info);

# define SPTTRF _LAPACK_FUNCTION(SPTTRF,spttrf)


/** LAPACK SPTTRF.
Computes the LDL**H factorization of a symmetric positive definite tridiagonal matrix.
[<a href="http://www.netlib.org/lapack/single/spttrf.f">http://www.netlib.org/lapack/single/spttrf.f</a>]
\ingroup vc_lapack_core
*/
void SPTTRF(
        ::VC::math::lapack::INTEGER* _n,
        float* _d,
        float* _e,
        ::VC::math::lapack::INTEGER* _info);


# define DGBEQU _LAPACK_FUNCTION(DGBEQU,dgbequ)


/** LAPACK DGBEQU.
Computes row and column scalings to equilibrate a general band matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/double/dgbequ.f">http://www.netlib.org/lapack/double/dgbequ.f</a>]
\ingroup vc_lapack_core
*/
void DGBEQU(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        double* _r,
        double* _c,
        double* _rowcnd,
        double* _colcnd,
        double* _amax,
        ::VC::math::lapack::INTEGER* _info);

# define SGBEQU _LAPACK_FUNCTION(SGBEQU,sgbequ)


/** LAPACK SGBEQU.
Computes row and column scalings to equilibrate a general band matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/single/sgbequ.f">http://www.netlib.org/lapack/single/sgbequ.f</a>]
\ingroup vc_lapack_core
*/
void SGBEQU(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        float* _r,
        float* _c,
        float* _rowcnd,
        float* _colcnd,
        float* _amax,
        ::VC::math::lapack::INTEGER* _info);


# define DTZRZF _LAPACK_FUNCTION(DTZRZF,dtzrzf)


/** LAPACK DTZRZF.
Computes an RZ factorization of an upper trapezoidal matrix (blocked version of DTZRQF).
[<a href="http://www.netlib.org/lapack/double/dtzrzf.f">http://www.netlib.org/lapack/double/dtzrzf.f</a>]
\ingroup vc_lapack_core
*/
void DTZRZF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _tau,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define STZRZF _LAPACK_FUNCTION(STZRZF,stzrzf)


/** LAPACK STZRZF.
Computes an RZ factorization of an upper trapezoidal matrix (blocked version of DTZRQF).
[<a href="http://www.netlib.org/lapack/single/stzrzf.f">http://www.netlib.org/lapack/single/stzrzf.f</a>]
\ingroup vc_lapack_core
*/
void STZRZF(
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _tau,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DTGSJA _LAPACK_FUNCTION(DTGSJA,dtgsja)


/** LAPACK DTGSJA.
Computes the generalized singular value decomposition of two real upper triangular (or trapezoidal) matrices as output by DGGSVP.
[<a href="http://www.netlib.org/lapack/double/dtgsja.f">http://www.netlib.org/lapack/double/dtgsja.f</a>]
\ingroup vc_lapack_core
*/
void DTGSJA(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const double* _tola,
        const double* _tolb,
        double* _alpha,
        double* _beta,
        double* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        double* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        double* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        double* _work,
        ::VC::math::lapack::INTEGER* _ncycle,
        ::VC::math::lapack::INTEGER* _info);

# define STGSJA _LAPACK_FUNCTION(STGSJA,stgsja)


/** LAPACK STGSJA.
Computes the generalized singular value decomposition of two real upper triangular (or trapezoidal) matrices as output by DGGSVP.
[<a href="http://www.netlib.org/lapack/single/stgsja.f">http://www.netlib.org/lapack/single/stgsja.f</a>]
\ingroup vc_lapack_core
*/
void STGSJA(
        const char* _jobu,
        const char* _jobv,
        const char* _jobq,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _k,
        ::VC::math::lapack::INTEGER* _l,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        const float* _tola,
        const float* _tolb,
        float* _alpha,
        float* _beta,
        float* _u,
        ::VC::math::lapack::INTEGER* _ldu,
        float* _v,
        ::VC::math::lapack::INTEGER* _ldv,
        float* _q,
        ::VC::math::lapack::INTEGER* _ldq,
        float* _work,
        ::VC::math::lapack::INTEGER* _ncycle,
        ::VC::math::lapack::INTEGER* _info);


# define DPOEQU _LAPACK_FUNCTION(DPOEQU,dpoequ)


/** LAPACK DPOEQU.
Computes row and column scalings to equilibrate a symmetric positive definite matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/double/dpoequ.f">http://www.netlib.org/lapack/double/dpoequ.f</a>]
\ingroup vc_lapack_core
*/
void DPOEQU(
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _s,
        double* _scond,
        double* _amax,
        ::VC::math::lapack::INTEGER* _info);

# define SPOEQU _LAPACK_FUNCTION(SPOEQU,spoequ)


/** LAPACK SPOEQU.
Computes row and column scalings to equilibrate a symmetric positive definite matrix and reduce its condition number.
[<a href="http://www.netlib.org/lapack/single/spoequ.f">http://www.netlib.org/lapack/single/spoequ.f</a>]
\ingroup vc_lapack_core
*/
void SPOEQU(
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _s,
        float* _scond,
        float* _amax,
        ::VC::math::lapack::INTEGER* _info);


# define DPOCON _LAPACK_FUNCTION(DPOCON,dpocon)


/** LAPACK DPOCON.
Estimates the reciprocal of the condition number of a symmetric positive definite matrix, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/double/dpocon.f">http://www.netlib.org/lapack/double/dpocon.f</a>]
\ingroup vc_lapack_core
*/
void DPOCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SPOCON _LAPACK_FUNCTION(SPOCON,spocon)


/** LAPACK SPOCON.
Estimates the reciprocal of the condition number of a symmetric positive definite matrix, using the Cholesky factorization computed by DPOTRF.
[<a href="http://www.netlib.org/lapack/single/spocon.f">http://www.netlib.org/lapack/single/spocon.f</a>]
\ingroup vc_lapack_core
*/
void SPOCON(
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        const float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGBSV _LAPACK_FUNCTION(DGBSV,dgbsv)


/** LAPACK DGBSV.
Solves a general banded system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/double/dgbsv.f">http://www.netlib.org/lapack/double/dgbsv.f</a>]
\ingroup vc_lapack_core
*/
void DGBSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _ipiv,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);

# define SGBSV _LAPACK_FUNCTION(SGBSV,sgbsv)


/** LAPACK SGBSV.
Solves a general banded system of linear equations AX=B.
[<a href="http://www.netlib.org/lapack/single/sgbsv.f">http://www.netlib.org/lapack/single/sgbsv.f</a>]
\ingroup vc_lapack_core
*/
void SGBSV(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        ::VC::math::lapack::INTEGER* _nrhs,
        float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        ::VC::math::lapack::INTEGER* _ipiv,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        ::VC::math::lapack::INTEGER* _info);


# define DSYGV _LAPACK_FUNCTION(DSYGV,dsygv)


/** LAPACK DSYGV.
Computes all eigenvalues and the eigenvectors of  a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
[<a href="http://www.netlib.org/lapack/double/dsygv.f">http://www.netlib.org/lapack/double/dsygv.f</a>]
\ingroup vc_lapack_core
*/
void DSYGV(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _w,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSYGV _LAPACK_FUNCTION(SSYGV,ssygv)


/** LAPACK SSYGV.
Computes all eigenvalues and the eigenvectors of  a generalized symmetric-definite generalized eigenproblem, Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
[<a href="http://www.netlib.org/lapack/single/ssygv.f">http://www.netlib.org/lapack/single/ssygv.f</a>]
\ingroup vc_lapack_core
*/
void SSYGV(
        ::VC::math::lapack::INTEGER* _itype,
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _w,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DSPEVD _LAPACK_FUNCTION(DSPEVD,dspevd)


/** LAPACK DSPEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix in packed storage.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/double/dspevd.f">http://www.netlib.org/lapack/double/dspevd.f</a>]
\ingroup vc_lapack_core
*/
void DSPEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        double* _ap,
        double* _w,
        double* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);

# define SSPEVD _LAPACK_FUNCTION(SSPEVD,sspevd)


/** LAPACK SSPEVD.
Computes all eigenvalues, and optionally, eigenvectors of a real symmetric matrix in packed storage.  If eigenvectors are desired, it uses a divide and conquer algorithm.
[<a href="http://www.netlib.org/lapack/single/sspevd.f">http://www.netlib.org/lapack/single/sspevd.f</a>]
\ingroup vc_lapack_core
*/
void SSPEVD(
        const char* _jobz,
        const char* _uplo,
        ::VC::math::lapack::INTEGER* _n,
        float* _ap,
        float* _w,
        float* _z,
        ::VC::math::lapack::INTEGER* _ldz,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _liwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGGGLM _LAPACK_FUNCTION(DGGGLM,dggglm)


/** LAPACK DGGGLM.

[<a href="http://www.netlib.org/lapack/double/dggglm.f">http://www.netlib.org/lapack/double/dggglm.f</a>]
\ingroup vc_lapack_core
*/
void DGGGLM(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        double* _a,
        ::VC::math::lapack::INTEGER* _lda,
        double* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        double* _d,
        double* _x,
        double* _y,
        double* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGGGLM _LAPACK_FUNCTION(SGGGLM,sggglm)


/** LAPACK SGGGLM.

[<a href="http://www.netlib.org/lapack/single/sggglm.f">http://www.netlib.org/lapack/single/sggglm.f</a>]
\ingroup vc_lapack_core
*/
void SGGGLM(
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _m,
        ::VC::math::lapack::INTEGER* _p,
        float* _a,
        ::VC::math::lapack::INTEGER* _lda,
        float* _b,
        ::VC::math::lapack::INTEGER* _ldb,
        float* _d,
        float* _x,
        float* _y,
        float* _work,
        ::VC::math::lapack::INTEGER* _lwork,
        ::VC::math::lapack::INTEGER* _info);


# define DGBCON _LAPACK_FUNCTION(DGBCON,dgbcon)


/** LAPACK DGBCON.
Estimates the reciprocal of the condition number of a general band matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGBTRF.
[<a href="http://www.netlib.org/lapack/double/dgbcon.f">http://www.netlib.org/lapack/double/dgbcon.f</a>]
\ingroup vc_lapack_core
*/
void DGBCON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        const double* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const double* _anorm,
        double* _rcond,
        double* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);

# define SGBCON _LAPACK_FUNCTION(SGBCON,sgbcon)


/** LAPACK SGBCON.
Estimates the reciprocal of the condition number of a general band matrix, in either the 1-norm or the infinity-norm, using the LU factorization computed by DGBTRF.
[<a href="http://www.netlib.org/lapack/single/sgbcon.f">http://www.netlib.org/lapack/single/sgbcon.f</a>]
\ingroup vc_lapack_core
*/
void SGBCON(
        const char* _norm,
        ::VC::math::lapack::INTEGER* _n,
        ::VC::math::lapack::INTEGER* _kl,
        ::VC::math::lapack::INTEGER* _ku,
        const float* _ab,
        ::VC::math::lapack::INTEGER* _ldab,
        const ::VC::math::lapack::INTEGER* _ipiv,
        const float* _anorm,
        float* _rcond,
        float* _work,
        ::VC::math::lapack::INTEGER* _iwork,
        ::VC::math::lapack::INTEGER* _info);


} // extern "C"