/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef NEWTON_HPP
#define NEWTON_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <kv/autodif.hpp>
#include <limits>
#include <boost/random.hpp>
#include <ctime>


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T, class F>
bool
newton(F f, ub::vector<T>& x, T epsilon = std::numeric_limits<T>::epsilon(), int maxloop = 100)
{
	int s = x.size();
	int i, j, r;
	ub::vector<T> fx;
	ub::matrix<T> fdx;
	T norm1, norm2, tmp;

	for (i=0; i<maxloop; i++) {
		try {
			autodif<T>::split(f(autodif<T>::init(x)), fx, fdx);
			ub::permutation_matrix<> pm(s);
			r = ub::lu_factorize(fdx, pm);
			if (r != 0) return false;
			ub::lu_substitute(fdx, pm, fx);
		}
		catch (...) {
			return false;
		}

		norm1 = 1.;
		norm2 = 0.;
		for (j=0; j<s; j++) {
			tmp = (x(j) >= 0.) ? x(j) : -x(j);
			norm1 = std::max(norm1, tmp);
			tmp = (fx(j) >= 0.) ? fx(j) : -fx(j);
			norm2 = std::max(norm2, tmp);
		}

		x = x - fx;
		if (norm2 <= norm1 * epsilon) return true;
	}
	return false;
}

template <class T, class F>
bool
newton_random(F f, ub::vector<T>& x, T epsilon = std::numeric_limits<T>::epsilon(), int maxloop = 100)
{
	int s = x.size();
	int i;

	using namespace boost;
	variate_generator< mt19937, normal_distribution<> > rand (mt19937(time(0)), normal_distribution<>(0., 10.));

	for (i=0; i<s; i++) x(i) = rand();
	return newton(f, x, epsilon, maxloop);
}

} // namespace kv

#endif // NEWTON_HPP
