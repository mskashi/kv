/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/random.hpp>
#include <ctime>
#include <kv/autodif.hpp>

namespace ub = boost::numeric::ublas;


namespace kv {

template <class T, class F>
bool
newton(ub::vector<T>& x, F f, T tor = 1e2, int maxloop = 100)
{
	int s = x.size();
	int i, r;
	ub::vector<T> fx, nx;
	ub::matrix<T> fdx;
	T norm1, norm2;

	for (i=0; i<maxloop; i++) {
		ub::permutation_matrix<> pm(s);
		autodif<T>::split(f(autodif<T>::init(x)), fx, fdx);
		r = ub::lu_factorize(fdx, pm);
		if (r != 0) return false;
		ub::lu_substitute(fdx, pm, fx);

		norm1 = norm_inf(x);
		norm2 = norm_inf(fx);
		x = x - fx;
		if (norm2 <= std::numeric_limits<T>::min() * tor) return true;
		if (norm2 <= norm1*(std::numeric_limits<T>::epsilon()*tor)) return true;
	}
	return false;
}

template <class T, class F>
bool
rand_newton(ub::vector<T>& x, F f, T lower = -1., T upper = 1, T tol = 1e2, int maxloop = 100)
{
	int s = x.size();
	int i;

	using namespace boost;
	variate_generator< mt19937, uniform_real<T> > rand (mt19937(time(0)), uniform_real<T>(lower, upper));

	for (i=0; i<s; i++) x(i) = rand();
	return newton(x, f, tol, maxloop);
}

} // namespace kv

#endif // NEWTON_HPP
