/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef MATRIX_INVERSION_HPP
#define MATRIX_INVERSION_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>


#ifdef USE_LAPACK
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/blas/blas3.hpp>
namespace bnb = boost::numeric::bindings;
#endif

#ifdef USE_ATLAS
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
namespace bnb = boost::numeric::bindings;
#endif


namespace ub = boost::numeric::ublas;

namespace kv {

template <class T>
bool invert(const ub::matrix<T>& a, ub::matrix<T>& b) {
	ub::matrix<T> tmp(a);
	ub::permutation_matrix<> pm(tmp.size1());

	if (ub::lu_factorize(tmp, pm) != 0) return false;

	b = ub::identity_matrix<T>(tmp.size1());

	ub::lu_substitute(tmp, pm, b);

	return true;
}

// special version for double
#if defined(USE_LAPACK) || defined(USE_ATLAS)
template <>
bool invert(const ub::matrix<double>& a, ub::matrix<double>& b) {
	ub::matrix<double, ub::column_major> tmp(a);
	ub::permutation_matrix<int> pm(tmp.size1());

	#ifdef USE_LAPACK
	if (bnb::lapack::getrf(tmp, pm) != 0) return false;
	if (bnb::lapack::getri(tmp, pm) != 0) return false;
	#endif
	#ifdef USE_ATLAS
	if (bnb::atlas::getrf(tmp, pm) != 0) return false;
	if (bnb::atlas::getri(tmp, pm) != 0) return false;
	#endif

	b = tmp;

	return true;
}
#endif // defined(USE_LAPACK) || defined(USE_ATLAS)


template <class T>
bool linear_equation(const ub::matrix<T>& a, const ub::vector<T>& b, ub::vector<T>& x) {
	ub::matrix<T> tmp(a);
	ub::permutation_matrix<> pm(tmp.size1());

	if (ub::lu_factorize(tmp, pm) != 0) return false;

	x = b;

	ub::lu_substitute(tmp, pm, x);

	return true;
}

// special version for double
#if defined(USE_LAPACK) || defined(USE_ATLAS)
template <>
bool linear_equation(const ub::matrix<double>& a, const ub::vector<double>& b, ub::vector<double>& x) {
	ub::matrix<double, ub::column_major> tmp(a);
	// int i;
	// int size = tmp.size1();
	ub::matrix<double, ub::column_major> tmp2(tmp.size1(), 1);

	// for (i=0; i<size; i++) tmp2(i, 0) = b(i);
	// ub::column(tmp2, 0).assign(b);
	ub::column(tmp2, 0) = b;

	#ifdef USE_LAPACK
	if (bnb::lapack::gesv(tmp, tmp2) != 0) return false;
	#endif
	#ifdef USE_ATLAS
	if (bnb::atlas::gesv(tmp, tmp2) != 0) return false;
	#endif

	// for (i=0; i<size; i++) x(i) = tmp2(i, 0);
	x = ub::column(tmp2, 0);

	return true;
}
#endif // defined(USE_LAPACK) || defined(USE_ATLAS)

template <class T>
void mm_mult(const ub::matrix<T>& a, const ub::matrix<T>& b, ub::matrix<T>& c) {
	c = ub::prod(a, b);
}

// special version for double
#if defined(USE_LAPACK) || defined(USE_ATLAS)
template <>
void mm_mult(const ub::matrix<double>& a, const ub::matrix<double>& b, ub::matrix<double>& c) {
	ub::matrix<double, ub::column_major> ca(a);
	ub::matrix<double, ub::column_major> cb(a);
	ub::matrix<double, ub::column_major> cc(c);

	#ifdef USE_LAPACK
		bnb::blas::gemm(ca, cb, cc);
	#endif
	#ifdef USE_ATLAS
		bnb::atlas::gemm(ca, cb, cc);
	#endif

	c = cc;
}
#endif // defined(USE_LAPACK) || defined(USE_ATLAS)

} // namespace kv

#endif // MATRIX_INVERSION_HPP
