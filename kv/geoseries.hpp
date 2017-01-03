/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef GEOSERIES_HPP
#define GEOSERIES_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

/*
 * calculate infinite sum of interval geometric series
 *  1 + r_1 + r_2r_1 + r_3r_2r_1 + ...
 *  (r_i \in I)
 */

namespace kv {

template <class T> interval<T> geoseries(const interval<T>& I)
{
	interval<T> l, u, tmp;

	if (I >= 0.) {
		return 1. / (1. - I);
	}

	l = I.lower();
	u = I.upper();

	if (I <= 0.) {
		tmp = 1 - l * u;
		return interval<T>(((1+l)/tmp).lower(), ((1+u)/tmp).upper());
	}

	return interval<T>((1 + l/(1-u)).lower(), (1/(1-u)).upper());
}

} // namespace kv

#endif // GEOSERIES_HPP
