/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef KRAW_APPROX_HPP
#define KRAW_APPROX_HPP

// Krawczyk method using approximate solution
// Newton iteration can be applied in advance.

#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/autodif.hpp>
#include <kv/matrix-inversion.hpp>
#include <kv/make-candidate.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template <class T, class F>
bool
krawczyk_approx(F f, const ub::vector<T>& c, ub::vector< interval<T> >& result, int newton_max = 2, int verbose = 1)
{
	int s = c.size();

	ub::vector< interval<T> > I, fc, fi, Rfc, C, K;
	ub::matrix< interval<T> > fdc, fdi, M;
	ub::vector<T> c2, minus;
	ub::matrix<T> R;
	int i, j;
	bool r;
	ub::vector<T> newton_step;
	T tmp, tmp2;

	c2 = c;

	// Newton iteration
	// use interval<T> for argument of f
	// preparing for the case that f can not accept T.

	for (i=0; i<newton_max; i++) {
		C = c2;
		try {
			autodif< interval<T> >::split(f(autodif< interval<T> >::init(C)), fc, fdc);
		}
		catch (std::domain_error& e) {
			return false;
		}
		r = invert(mid(fdc), R);
		if (!r) return false;

		minus = prod(R, mid(fc));

		tmp = 1.;
		tmp2 = 0.;
		for (j=0; j<s; j++) {
			using std::abs;
			tmp = std::max(tmp, abs(c2(j)));
			tmp2 = std::max(tmp2, abs(minus(j)));
		}

		c2 = c2 - minus;
		if (verbose >= 1) {
			std::cout << "newton" << i << ": " << c2 << "\n";
		}
		if (tmp2 <= tmp * std::numeric_limits<T>::epsilon()) break;
	}

	C = c2;
	try {
		autodif< interval<T> >::split(f(autodif< interval<T> >::init(C)), fc, fdc);
	}
	catch (std::domain_error& e) {
		return false;
	}
	r = invert(mid(fdc), R);
	if (!r) return false;
	Rfc = prod(R, fc);

	newton_step.resize(s);
	for (i=0; i<s; i++) {
		newton_step(i) = norm(Rfc(i));
	}

	make_candidate(newton_step);

	I = C;
	for (i=0; i<s; i++) {
		tmp = std::numeric_limits<T>::epsilon() * norm(I(i)) * (s+1) * 2;
		tmp2 = std::numeric_limits<T>::min() * (s+1) * 2;
		if (newton_step(i) < tmp) newton_step(i) = tmp;
		if (newton_step(i) < tmp2) newton_step(i) = tmp2;
		I(i) += newton_step(i) * interval<T>(-1., 1.);
	}

	if (verbose >= 1) {
		std::cout << "I: " << I << "\n";
	}

	try {
		autodif< interval<T> >::split(f(autodif< interval<T> >::init(I)), fi, fdi);
	}
	catch (std::domain_error& e) {
		return false;
	}

	// M = ub::identity_matrix< interval<T> >(s) - prod(R, fdi);
	M = ub::identity_matrix< interval<T> >(s);
	M -= prod(R, fdi);

	K = C - Rfc +  prod(M, I - C);

	if (verbose >= 1) {
		std::cout << "K: " << K << "\n";
	}

	if (proper_subset(K, I)) {
		result = K;
		return true;
	} else {
		return false;
	}
}


namespace krawczyk_approx_sub {

// generate 1-d vector function from scalar function
template <class F>
struct MakeVec {
	F f;
	MakeVec(F f): f(f) {}

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		ub::vector<T> r(1);
		r(0) = f(x(0));
		return r;
	}
};

} // namespace krawczyk_approx_sub;


// 1 dimensional version

template <class T, class F>
bool
krawczyk_approx(F f, const T& c, interval<T>& result, int newton_max = 2, int verbose = 1)
{
	ub::vector<T> in(1);
	ub::vector< interval<T> > out;
	krawczyk_approx_sub::MakeVec<F> g(f);
	bool r;

	in(0) = c;
	r = krawczyk_approx(g, in, out, newton_max, verbose);
	if (!r) return r;
	result = out(0);
	return r;
}

} // namespace kv

#endif // KRAW_APPROX_HPP
