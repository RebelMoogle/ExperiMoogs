# ifdef DOXYGEN_SKIP

/** \file

    This file is included by blas_matrix.hh.

    DO NOT INCLUDE IT DIRECTLY!

    \arg This file contains only documentation!

    \internal
 */

#ifndef VC_MATH_BLAS_MATRIX_HH
# error "don't include directly"
#endif

namespace VC {
namespace math {
namespace blas {
//=============================================================================


/** \defgroup vc_blas_matrix_operations Overview of interface to BLAS vector and matrix operations
    \ingroup vc_blas

    Several "levels" of interfaces are provided.

    \section sec_blas_core BLAS core interface

    The namespace VC::math::blas defines C++ prototypes for BLAS functions
    and wrappers for C++-style calling with per value arguments if
    possible (e.g., for flags, dimensions, or scalar values). These
    wrappers are polymorphic and accept either \c float or \c double
    arguments.

    For example
    \arg DAXPY() and SAXPY() are mapped to
    \arg axpy(int,double,const double*,int,double*,int) and
    axpy(int,float,const float*,int,float*,int)

    The wrappers are defined for
    \arg [\ref vc_blas_1],
    \arg [\ref vc_blas_2], and
    \arg [\ref vc_blas_3].

    (Note that BLAS call may be replaced for optimization/inlining, see [\ref vc_vcblas].)

    Flags are
    \arg VC::math::blas::UpperLowerFlag,
    \arg VC::math::blas::TransposeFlag,
    \arg VC::math::blas::DiagonalFlag, and
    \arg VC::math::blas::SideFlag.


    \section sec_blas_mat BLAS matrix and vector structures

    BLAS deals with different kind of matrices, e.g., \a general
    matrices tagged VC::math::blas::GE.

    A summary on BAS matrix storage is provided here:
    http://www.netlib.org/lapack/lug/node121.html

    The parameters of the matrices like their dimension, flags (see
    above), etc. are collected in structures. For example, ge_mat
    captures the parameters to a general matrix.

    Parameters can be fixed at compile time (a certain number
    N!=VarInt) or free as run-time variables (template parameter
    VC::math::blas::VarInt).

    Note that for this reason parameters are \a checked at \a
    run-time, i.e., there is (generally) no compiler error for, e.g.,
    mismatch in dimensions! The operation will instead fail with an
    assertion at run-time.

    The next level of \a wrappers act on, e.g., a ge_mat and a pointer
    to the data. For these wrappers, \b specializations may be made
    available which, e.g., call inline code for fixed dimensions.


    Similarly, vector parameters (length and increment) are wrapped in
    VC::math::blas::vec.

    \sa [\ref vc_blas_mat] for a summary of matrix structures and BLAS
    function wrappers

    \section sec_blas_extensions Extensions to BLAS

    The following functions are defined as extensions to BLAS. They
    are equally defined for matrix_reference_t arguments and as
    methods of matrix_reference_t (see below).

    | method/function | expression                | comment                                                                  |
    |-----------------|---------------------------|--------------------------------------------------------------------------|
    | ld_zero()       | `A=VC::math::blas::zeros` | set all elements to  zero (all matrix types)                             |
    | ld_all()        | `A=VC::math::blas::all`   | set all element to constant (GE,SY,SP)                                   |
    | ld_eye()        | `A=VC::math::blas::eye`   | load identity (all matrix types)                                         |
    | ld_diag()       | `A=diag(x)`               | load diagonal matrix (GE, SY, no UnitDiag)                               |
    | mscal()         | `A*=x`                    | multiply by scalar                                                       |
    | adds()          | `A+=x`                    | elementwise add scalar (GE,SY,SP)                                        |
    | emul()          | `A*=B`                    | element-wise multiplication (vectors and GE only)                        |
    | ediv()          | `A/=B`                    | element-wise division (vectors and GE only)                              |
    | madd()          | `A+=B`                    | addition (same types, no VC::math::blas::UnitDiag)                       |
    | mscal_cols()    | `A*=diag(x)`              | scale columns = right multiplication by diagonal matrix (GE,GB,TR,TB,TP) |
    |                 | `trans(A)*=diag(x)`       | scale rows = left multiplication by diagonal matrix (GE,GB,TR,TB,TP)     |
    | cp()            | `A=B`e                    | copy matrix (any to GE, respect TransposeFlag; `SY<->SP`)                |
    | swap()          | `A=T; A=B; B=T;           | swap matrix **contents**                                                 |
    
    Element-wise maps can be defined by

    | function  | element-wise expression |
    |-----------|-------------------------|
    | map_f()   | a=op(a)                 |
    | map_f2()  | a=op(b)                 |
    | map_f12() | a=op(a,b)               |
    | map_f23() | a=op(b,c)               |

    
    for same matrix types and matching dimensions (for ge_mat, first
    argument to op may be transposed), processes only "stored"
    elements (i.e., disregards zeros (GB,SB,TR,TB,TP) or ones
    (TR,TB,TP)); may change unused elements for GB,SB. 

    \section sec_blas_indexing Indexing
    \anchor label_blas_indexing
    \arg The special symbols
    \code
    _(i1,i2)  ==> MATLAB "range" 'i1:i2'
    $         ==> MATLAB "range" ':', equivalent to empty "_()"
    end       ==> MATLAB specifier 'end' (supports index arithmetic)
    \endcode
    Are used for indexing with \c operator()() or \c block().
    (Note that "_", "$" and "end" live in the VC::math::blas namespace!)
    \arg operator()(...) for \a copying a submatrix (cp_block()).
    \code
    A=B(i1,i2,j1,j2);                   // A=B(i1:i2,j1:j2)
    A=B(_(i1,i2),_(j1,j2))              // A=B(i1:i2,j1:j2)
    A=B($,_(j1,j2)) == B(_(),_(j1,j2)); // A=B(:,j1:j2)
    A=B($,_(j1,end));                   // A=B(:,j1:end)
    A=B($,_(j1,end-k));                 // A=B(:,j1:end-k)

    x=B(i,_(j1,j2));                    // x=B(i,j1:j2)
    x=B(_(i1,end),end-2);               // x=B(i1:end,end-2)  %% vector

    a=A(i,j)                            // a=A(i,j)           %% scalar
    a=A(end-i,j)                        // a=A(end-i,j)
    a=A(i,end-j)                        // a=A(i,end-j)
    a=A(end-i,end-j)                    // a=A(end-i,end-j)
    \endcode
    \arg block() extracts a submatrix of a GE matrix as matrix_reference_t
    (or a subvector of a vector_reference_t).
    \arg see also column() row(), diag()
    \code
    matrix_reference_t<T,ge_mat<NoT,3,4> > A;
    block(A,i1,i2,j1,j2)=eye;
    // ...
    block(A,($,_(j1,end)) // ...
    \endcode
    \arg GE_Matrix A implements block(A,...) as A::operator()(...):
    \code
    GE_Matrix<T> A(3,4);
    A($,_(1,end))=eye;
    \endcode
    For vector_reference_t this is standard for operator()(...):
    \code
    vector_reference_t<T,V> x;
    x(i1,i2)=zeros;
    x(_(i1,i2))=zeros;
    x(_(i,end-1))=zeros;
    \endcode

    \section sec_blas_iterators Iterators and initialization
    \anchor label_blas_iterators

    Alternatively, iterators are provided for matrix access and \a
    initialization, e.g., VC::math::blas::GE_row_iterator.

    \code
    matrix_reference_t<double,ge_mat<...> > A=...;

    GE_row_iterator<double> ii=iterate_rows(A);
    ii << 1,2,3,i_fill_row(4),i_fill(5);

    GE_Matrix<double> B;
    GE_Matrix<double>::row_iterator_t jj=B.row_iterator();

    jj << i_fill(1.0,3),2,3,i_skip(4),5,i_skip_row,i_fill(6);

    \endcode

    \sa [\ref vc_blas_iter]

    \section sec_blas_matref BLAS matrix and vector references
    \anchor label_blas_matref

    The previously defined structures like, e.g., ge_mat, capture only
    information about "matrix shape". The do not reference any data.

    This is done by matrix_reference_t (matrix_const_reference_t) and
    vector_reference_t (vector_const_reference_t), which provide

    \arg the combination of matrix shape (e.g., ge_mat) with a pointer
    to data, and
    \arg operations on these data.

    You can use vec and mat to construct references or specify
    reference types conveniently. Here is an example:
    \code
    VecN<double,3> v;
    vec<3>::ref_t<double>::type vr=vec<3>::ref(v.data());  // fixed dimension
    vec<>::ref_t<double>::type  vs=vec<>::ref(v.data(),3); // variable dimension

    double A[4*3];
    mat<GE>::ref_t<double>::type     Ar=mat<GE>::ref(A,4,3);
    mat<GE,4,3>::ref_t<double>::type As=mat<GE>::ref(A);
    \endcode

    Note that

    \arg operations \a never allocate any (temporary) storage!

    \arg all operations modify at most the data referenced by the
    calling vector_reference_t or matrix_reference_t object!

    \arg not all operations may be defined for \a all matrix types, e.g.,
    a reference to a triangular matrix tr_mat cannot call
    \code
    matrix_reference_t<T,A>::mm(const T&,const matrix_const_reference_t<T,A>& _a,
    const matrix_const_reference_t<T,B>& _b,const T& _beta)
    \endcode
    because there is no way, trmm() could be called (via two levels of wrapping).
    The result is an error at \a compile \a time! -- Hence, it makes sense to
    know about available BLAS operations.

    There are various global functions on references, e.g., providing
    access to data (such as at()). Most of them are also available as
    methods of matrix_reference_t. Of particular interest is trans()
    which provides a \a reference to the transposed matrix (by
    changing the VC::math::blas::TransposeFlag, the data is not
    modified).

    \section sec_blas_matrices Matrix data structures

    matrix_reference_t (and matrix_const_reference_t) implement
    operations on matrices, however, no storage is
    allocated. Generally, you are advised to use these references,
    which provide linear algebra functionality.

    In addition, there are few classes provided for convenience which
    combine references with storage allocation, e.g., GE_Matrix. Note
    that this applies only to \a variable size matrices! (With the
    present implementation there would be a small but unacceptable
    storage overhead for fixed size matrices.)

    Allocators conforming to std::allocator are provided as template
    argument (can be replaced by, e.g., VC::base::small_vector_allocator).

    Vec_N, GE_Mat_MxN, SP_Mat_NxN, TP_Mat_NxN provide storage (and
    types) for small, fixed size matrices. In contrast to
    GE_Matrix,..., they are \a not derived from matrix_reference_t and
    provide no functionality themselves! Use GE_Mat_MxN::operator() to
    obtain a matrix_reference_t. (Remark: inheritance would "cost" an
    additional pointer for storing the reference.)

    See [\ref vc_blas_matrices].

    \section sec_blas_matrix_algorithms Matrix algorithms

    Matrix algorithms are mostly calls to respective LAPACK functions
    set up from matrix_reference_t. Some calls may require additional
    storage for (temporary) data.

    See [\ref vc_blas_algorithms], [\ref vc_lapack]

    \example vclibs/math/example_blas_matrix.cc
    This example demonstrates various features of
    [\ref vc_blas_matrix_operations].
 */

/** \defgroup vc_blas_mat BLAS matrix types
    \ingroup vc_blas

    BLAS matrix types template classes parameterized by the respective
    properties such as, e.g., matrix dimensions.

    \arg The parameters can be either fixed at compile time or
    run-time variables (specify VarInt) as template parameter.

    \arg Parameters are access via the \c mat_prop classes which
    define the respective attribute (e.g., mat_prop_m defines
    mat_prop_m::m() and hence, e.g., ge_mat::m()).

    For example ge_mat specifies a BLAS GEneral matrix. It is derived
    from various \c mat_prop classes defining methods ge_mat.m(),
    ge_mat.n(), etc.

    Available attributes are uplo(), trans(), diag(), m(), n(), k(),
    kl(), ku(), ld(), depending on the matrix type. Attributes might be
    defined as methods (run-time variable, template parameter is
    VC::math::blas::VarInt) or static class functions (attribute fixed
    at compile time). uplo(), trans(), and diag() are \a always fixed
    at compile time! \sa mat_prop

    Only the \a structure/shape of the matrix is specified -- not the
    data. matrix_reference_t combines this information with a pointer
    to data.

    Notes

    \arg Similarly, vectors are defined in vec (with attributes n()
    and inc()).

    \arg The TransposeFlag is relevant only to algorithms taking the
    matrix as input. It does neither transpose the matrix entries (in
    memory) nor does it swap dimensions m() and n(). \sa trans()

    \arg A summary on BLAS matrix storage is provided here:
    http://www.netlib.org/lapack/lug/node121.html

    \sa [\ref vc_blas_matrix_operations], vec, matrix_reference_t, vector_reference_t
 */

/** \defgroup vc_blas_matrices Matrix data structures
    \ingroup vc_blas

    Combines matrix_reference_t and matrix_const_reference_t with
    allocation of storage. Matrices implement all operators and
    methods of matrix_reference_t and matrix_const_reference_t by
    delegation via Matrix_delegate.

    For all matrices, resize() (e.g., GE_Matrix::resize()) \a
    invalidates any references to this matrix!

    \b Note:
    \arg There is \a no copy constructor!
    \arg Any call to the \a assignment operator is \a delegated to
    assignment operators on matrix_reference_t!
    \arg References to tr_mat, sy_mat are available as
    GE_Matrix::tr_lower(),... (similar for other matrices) as they don't
    appear on their own.
*/

/** \defgroup vc_blas_algorithms Matrix algorithms
    \ingroup vc_blas

    Algorithms on matrices represented as matrix_reference_t, mostly calls
    to LAPACK (see http://www.netlib.org/lapack/lug/).

    \arg LU decomposition (GE,GB): lu(), lu_subs(), lu_solve(), lu_invert()
    \arg Cholesky decomposition (SY,SP,SB): chol(), chol_subs(), chol_solve(),
    chol_invert()
    \arg Linear least squares problems: lsq_solve(), svd_solve()
    \arg Eigenvalues and eigenvectors: eig()
    \arg Singular value decomposition: svd()

    The following functions are provided for convenience (the argument is
    copied to temporary storage and is not destroyed).

    \arg rank() determines effective rank w.r.t. tolerance from svd()
    \arg norm2() compute L_2 norm from svd()
    \arg normf() compute Frobenius norm (GE)
    \arg cond2() compute condition numner w.r.t. L_2 norm using svd()
    \arg pinv() compute pseudo inverse from svd()

    Algorithms are implemented as global functions rather than
    methods of matrices. Most of them overwrite/destroy their
    operand(s).

    Temporary memory is either allocated from the stack via
    std::alloca() (no more than O(max(m,n)), or from
    lapack::Workspace(), or can be provided directly by
    the user (skips query for storage size lwork/iwork which is
    otherwise performed automatically).

    \todo lu_solve tridiagonal KL=KU=1 --> copy diagonals, gtrsv, (difference: _a unmodified)
    \todo lu/chol generalization / switch to tuned algorithms for small matrices

    \sa [\ref vc_lapack]
*/

/** \defgroup vc_blas_io Vector and matrix input/output
    \ingroup vc_blas

    Matrices and vectors can be output by std::ostream::operator<<.

    \section sec_blas_io_ascii ASCII output

    There are several output formats available
    \arg print(): prints all elements column-wise separated by space, rows
    are separated by semicolons
    \arg pp(): pretty-prints type, dimensions and elements in a
      multi-line representation \b (default)
    \arg pmatlab(): prints a single-line Matlab representation
    \arg pdata(): prints all stored data values separated by space,
     e.g., only triangular part w/o zeros for tr_mat. The specific order
     is defined implicitly by the matrix type. This is available only for
     matrices.

    Examples
    \code
    VC::math::blas::matrix_const_reference_t<T,A> A;

    cout << A;           // same as cout << print(A)
    cout << pp(A);
    cout << pmatlab(A);
    cout << pdata(A);
    \endcode

    \b Note: only vectors/matrices can be output, you \a cannot output
    vector/matrix expressions! You need to assign the expression to a
    (temporary) object first:

    \code
    VC::math::blas::vector_const_reference_t<T,v> x,y;

    cout << x+y; // fails (will output "Expression...")

    z=x+y;
    cout << z;
    \endcode

    \section sec_blas_io_binary Binary input/output

    \arg write_mat4() writes matrices (or vectors) in MATLAB V4 format
    (Level 4 MAT-file).
    \arg GE_Matrix<T>::read_mat4() reads matrices in MATLAB V4 format.

    Multiple matrices can be read or written sequentially.

    \code
    // write matrices
    VC::math::blas::matrix_const_reference_t<T,A> Aout;
    VC::math::blas::vector_const_reference_t<T,V> xout;
    ofstream out(...);

    write_mat4(Aout,"A").write_mat4(xout,"x");

    // read matrices
    ifstream in;
    VC::math::blas::GE_Matrix Ain, xin;             // note: no Vector input

    string nameA=Ain.read_mat4(in);
    string namex=xin.read_mat4(in); // access vector as column(x,1)

    assert(nameA=="A" && namex=="x");

    \endcode

    There is no Vector input, use GE_Matrix, check dimensions and
    access column().

    \sa [\ref vc_mat4] (provides more functionality)
 */

/** \defgroup vc_blas_details Details
    \ingroup vc_blas

    \todo diag() for SB,SY; simliar column(), row() (truncated
    vectors, varying dimensions)

    \todo max, min, sum, ...

    \todo specializations...

    \todo simplify code using Foreach, map (elementwise operations)

    \todo ForeachT for TR, TB, TP

    \todo permutations

    \todo matrix exponential, see Moler, van Loan; matlab expm1demo;
    also http://www.maths.uq.edu.au/expokit/
 */


/** Outer product `x*x'` in expression (for rank-1 update).
    \ingroup vc_blas
    Allows to write expressions
    
    C++                      | MATLAB
    -------------------------|-----------------
    S =alpha*outer_prod(x);  | S=alpha*x'*x;
    S+=alpha*outer_prod(x);  | S=S+alpha*x'*x;
    S-=alpha*outer_prod(x);  | S=S-alpha*x'*x;

    where `x` is a vector and `S` is a *symmetric* matrix.

    \sa VC::math::blas::syr
*/
template <typename T,typename VX>Expression<Type_matrix,Expr_axxt,T,VX>
outer_prod(const vector_const_reference_t<T,VX>& _x);


/** Expression `A'*A` for rank-k update.
    \ingroup vc_blas
    Allows to write expressions
    
    C++               | MATLAB
    ------------------|-----------------
    S =alpha*ata(A);  | S=alpha*A'*A;
    S+=alpha*ata(A);  | S=S+alpha*A'*A;
    S-=alpha*ata(A);  | S=S-alpha*A'*A;

    where `A` is a `GE` matrix and `S` is a *symmetric* matrix.        

    \sa VC::math::blas::syrk
*/
template <typename T,typename MA>Expression<Type_matrix,Expr_aAtA,T,MA>
ata(const matrix_const_reference_t<T,MA>& _A);



//
// SCRATCH
//

/*


** Load diagonal matrix (ld_diag)) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    VC::math::blas::vector_const_reference_t<T,V> d;
    A=diag(d);
    \endcode
*
template <typename T,typename V>
Closure_diag<T,V> diag(const vector_const_reference_t<T,V>& _x) {
  return Closure_diag<T,V>(_x);
}
** Load diagonal matrix (ld_diag)) \ingroup vc_blas
    Use as
    \code
    VC::math::blas::matrix_reference_t<T,A> A;
    VC::math::blas::vector_const_reference_t<T,V> d;
    A=diag(d);
    \endcode
*
template <typename T,typename V>
Closure_diag<T,V> diag(const vector_reference_t<T,V>& _x) {
  return Closure_diag<T,V>(_x.const_ref());
}

*/





//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC

# endif // DOXYGEN_SKIP
