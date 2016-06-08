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

    \arg Implementation of methods of VC::math::blas::vector_reference_t,
    VC::math::blas::matrix_reference_t,...

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

//
// Implementation of class methods vector_reference_t
//

// level 1

template <typename T,typename V>
template <typename W>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::swap(const vector_reference_t<T,W>& _x) const {
  ::VC::math::blas::swap(_x,*this);
  return *this;
}

template <typename T,typename V>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::scal(const T& _alpha) const {
  ::VC::math::blas::scal(_alpha,*this);
  return *this;
}

template <typename T,typename V>
template <typename W>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::copy(const vector_const_reference_t<T,W>& _x) const {
  ::VC::math::blas::copy(_x,*this);
  return *this;
}

template <typename T,typename V>
template <typename W>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::copy_scal(const T& _alpha,
                                   const vector_const_reference_t<T,W>& _x) const {
  ::VC::math::blas::copy_scal(_alpha,_x,*this);
  return *this;
}


template <typename T,typename V>
template <typename W>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::axpy(const T& _alpha,
                              const vector_const_reference_t<T,W>& _x) const {
  ::VC::math::blas::axpy(_alpha,_x,*this);
  return *this;
}

template <typename REFBASE>
template <typename W>
typename vector_const_reference_mixin<REFBASE>::value_t
vector_const_reference_mixin<REFBASE>::dot
(const vector_const_reference_t<typename vector_const_reference_mixin<REFBASE>::value_t,W>& _x) const {
  return ::VC::math::blas::dot(_x,this->const_ref());
}

template <typename REFBASE>
template <typename W>
double
vector_const_reference_mixin<REFBASE>::ddot
(const vector_const_reference_t<typename vector_const_reference_mixin<REFBASE>::value_t,W>& _x) const {
  return ::VC::math::blas::ddot(_x,this->const_ref());
}

template <typename REFBASE>
typename vector_const_reference_mixin<REFBASE>::value_t
vector_const_reference_mixin<REFBASE>::nrm2() const {
  return ::VC::math::blas::nrm2(this->const_ref());
}

template <typename REFBASE>
typename vector_const_reference_mixin<REFBASE>::value_t
vector_const_reference_mixin<REFBASE>::asum() const {
  return ::VC::math::blas::asum(this->const_ref());
}

// level 2

template <typename T,typename V>
template <typename A,typename W>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::mv
(const T& _alpha,
 const matrix_const_reference_t<T,A>& _A,
 const vector_const_reference_t<T,W>& _x,const T& _beta) const {
  static_assert(!std::is_same<typename A::matrix_id_t,TR_Mat>::value &&
                !std::is_same<typename A::matrix_id_t,TB_Mat>::value &&
                !std::is_same<typename A::matrix_id_t,TP_Mat>::value,
                "this matrix-vector product is not available for TR,TB,TP matrices");
  ::VC::math::blas::mv(_alpha,_A,_x,_beta,*this); // GE,GB,SY,SB,SP only
  return *this;
}

template <typename T,typename V>
template <typename A>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::mv
(const matrix_const_reference_t<T,A>& _A) const {
  static_assert(std::is_same<typename A::matrix_id_t,TR_Mat>::value ||
                std::is_same<typename A::matrix_id_t,TB_Mat>::value ||
                std::is_same<typename A::matrix_id_t,TP_Mat>::value,
                "this matrix-vector product is only available for TR,TB,TP matrix");
  ::VC::math::blas::mv(_A,*this); // TR, TB, TP only
  return *this;
}

template <typename T,typename V>
template <typename A>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::sv
(const matrix_const_reference_t<T,A>& _A) const {
  static_assert(std::is_same<typename A::matrix_id_t,TR_Mat>::value ||
                std::is_same<typename A::matrix_id_t,TB_Mat>::value ||
                std::is_same<typename A::matrix_id_t,TP_Mat>::value,
                "this matrix-vector product is only available for TR,TB,TP matrix");
  ::VC::math::blas::sv(_A,*this); // TR, TB, TP only
  return *this;
}

template <typename T,typename V>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::ld_all(const T& _a) const {
  ::VC::math::blas::ld_all(_a,*this);
  return *this;
}

template <typename T,typename V>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::ld_unit(int _i,const T& _a) const {
  ::VC::math::blas::ld_unit(_i,_a,*this);
  return *this;
}

template <typename T,typename V>
const vector_reference_t<T,V>&
vector_reference_t<T,V>::ld_rand() const {
  ::VC::math::blas::ld_rand(*this);
  return *this;
}

template <typename T,typename M>
const vector_reference_t<T,M>&
vector_reference_t<T,M>::adds(const T& _alpha) const {
  ::VC::math::blas::adds(_alpha,*this);
  return *this;
}


//-----------------------------------------------------------------------------

//
// Implementation of class methods matrix_reference_t
//

template <typename T,typename M>
template <typename A,typename B>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::mm
(const T& _alpha,
 const matrix_const_reference_t<T,A>& _a,
 const matrix_const_reference_t<T,B>& _b,
 const T& _beta) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,GE_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value,
                "matrix-matrix product requires GE or SY matrix");
  VC::math::blas::mm(_alpha,_a,_b,_beta,*this);
  return *this;
}

template <typename T,typename M>
template <typename A,typename B>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::mm
(SideFlag _side,
 const T& _alpha,
 const matrix_const_reference_t<T,A>& _a,
 const matrix_const_reference_t<T,B>& _b,
 const T& _beta) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,GE_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value,
                "this matrix-matrix product requires GE or SY matrix");
  VC::math::blas::mm(_side,_alpha,_a,_b,_beta,*this);
  return *this;
}

template <typename T,typename M>
template <typename A>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::mm
(const T& _alpha,const matrix_const_reference_t<T,A>& _a) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,TR_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,TB_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,TP_Mat>::value,
                "this matrix-matrix product only available for TR,TB,TP matrix");
  VC::math::blas::mm(_alpha,_a,*this);
  return *this;
}

template <typename T,typename M>
template <typename V>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::r
(const T& _alpha,const vector_const_reference_t<T,V>& _x) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,SP_Mat>::value,
                "rank-1 update requires SY or SP matrix");
  VC::math::blas::r(_alpha,_x,*this);
  return *this;
}

template <typename T,typename M>
template <typename V,typename W>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::r2
(const T& _alpha,
 const vector_const_reference_t<T,V>& _x,
 const vector_const_reference_t<T,W>& _y) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,SP_Mat>::value,
                "rank-2 update requires SY or SP matrix");
  VC::math::blas::r2(_alpha,_x,_y,*this);
  return *this;
}

template <typename T,typename M>
template <typename A>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::rk
(const T& _alpha,const matrix_const_reference_t<T,A>& _a,const T& _beta) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value,
                "rank-k update requires SY matrix");
  VC::math::blas::rk(_alpha,_a,_beta,*this);
  return *this;
}

template <typename T,typename M>
template <typename A,typename B>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::r2k
(const T& _alpha,
 const matrix_const_reference_t<T,A>& _a,
 const matrix_const_reference_t<T,B>& _b,const T& _beta) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,SY_Mat>::value,
                "rank-2 update requires SY matrix");
  VC::math::blas::r2k(_alpha,_a,_b,_beta,*this);
  return *this;
}

template <typename T,typename M>
template <typename A>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::sm
(const T& _alpha,const matrix_const_reference_t<T,A>& _a) const {
  static_assert(std::is_same<typename matrix_t::matrix_id_t,TR_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,TB_Mat>::value ||
                std::is_same<typename matrix_t::matrix_id_t,TP_Mat>::value,
                "solve (triangular) system only available for TR,TB,TP matrix");
  VC::math::blas::sm(_alpha,_a,*this);
  return *this;
}

//-----------------------------------------------------------------------------

template <typename T,typename M>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::ld_zero() const {
  ::VC::math::blas::ld_zero(*this);
  return *this;
}

template <typename T,typename M>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::ld_all(const T& _a) const {
  ::VC::math::blas::ld_all(_a,*this);
  return *this;
}

template <typename T,typename M>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::ld_eye(const T& _a) const {
  ::VC::math::blas::ld_eye(*this,_a);
  return *this;
}

template <typename T,typename M>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::ld_rand() const {
  ::VC::math::blas::ld_rand(*this);
  return *this;
}

template <typename T,typename M>
template <typename V>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::ld_diag(const vector_const_reference_t<T,V>& _d,
                                 const T& _alpha) const {
  ::VC::math::blas::ld_diag(_alpha,_d,*this);
  return *this;
}

//-----------------------------------------------------------------------------

template <typename T,typename M>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::mscal(const T& _alpha) const {
  ::VC::math::blas::mscal(_alpha,*this);
  return *this;
}

template <typename T,typename M>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::adds(const T& _alpha) const {
  ::VC::math::blas::adds(_alpha,*this);
  return *this;
}

template <typename T,typename M>
template <typename B>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::emul(const matrix_const_reference_t<T,B>& _b,
                              const T& _alpha) const {
  ::VC::math::blas::emul(_alpha,_b,*this);
  return *this;
}

template <typename T,typename M>
template <typename B>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::ediv(const matrix_const_reference_t<T,B>& _b,
                              const T& _alpha) const {
  ::VC::math::blas::ediv(_alpha,_b,*this);
  return *this;
}

template <typename T,typename M>
template <typename B>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::madd(const T& _alpha,
                              const matrix_const_reference_t<T,B>& _b) const {
  ::VC::math::blas::madd(_alpha,_b,*this);
  return *this;
}

template <typename T,typename M>
template <typename V>
const matrix_reference_t<T,M>&
matrix_reference_t<T,M>::scal_cols(const vector_const_reference_t<T,V>& _d,
                                    const T& _alpha) const {
  ::VC::math::blas::mscal_cols(_alpha,_d,*this);
  return *this;
}

template <typename T,typename M>
const matrix_reference_t<T,M>& matrix_reference_t<T,M>::add_trans(const T& _alpha) const {
  ::VC::math::blas::add_trans(_alpha,*this);
  return *this;
}

#endif // DOXYGEN_SKIP

//=============================================================================
} // namespace blas
} // namespace math
} // namespace VC
