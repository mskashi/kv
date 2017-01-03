/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef BETA_HPP
#define BETA_HPP

#include <cmath>
#include <kv/defint-singular.hpp>

namespace kv {

template <class TT> struct Beta {
	TT x, y;
	Beta(TT x, TT y) : x(x), y(y) {}

	template <class T> T operator() (const T& t) {
		return pow(t, T(x)-1) * pow(1-t, T(y)-1);
	}
};

template <class TT> struct Beta_s1 {
	TT y;
	Beta_s1(TT y) : y(y) {}

	template <class T> T operator() (const T& t) {
		return pow(1-t, T(y)-1);
	}
};

template <class TT> struct Beta_s2 {
	TT x;
	Beta_s2(TT x) : x(x) {}

	template <class T> T operator() (const T& t) {
		return pow(t, T(x)-1);
	}
};

#if !defined(BETA_ORDER)
#define BETA_ORDER 18
#endif

template <class T>
interval<T> beta_r(const interval<T>& x, const interval<T>& y) {
	return	defint_power_autostep(Beta< interval<T> >(x,y), Beta_s1< interval<T> >(y), interval<T>(0.), interval<T>(0.5), BETA_ORDER, interval<T>(x-1))
		+ defint_power_autostep_r(Beta< interval<T> >(x,y), Beta_s2< interval<T> >(x), interval<T>(0.5), interval<T>(1.), BETA_ORDER, interval<T>(y-1));
}

template <class T>
interval<T> beta(const interval<T>& x, const interval<T>& y) {
	T n_x, n_y;
	int i;
	interval<T> r;

	using std::floor;
	n_x = floor(x.lower()) - 1.;
	n_y = floor(y.lower()) - 1.;

	r = beta_r(x - n_x, y - n_y);

	if (n_x >= 1.) {
		for (i=1; i<=n_x; i++) r *= x - i;
	} else {
		for (i=n_x+1.; i<=0.; i++) r /= x - i;
	}

	if (n_x + n_y >= 1.) {
		for (i=1; i<=n_x+n_y; i++) r /= x + y - i;
	} else {
		for (i=n_x+n_y+1.; i<=0.; i++) r *= x + y - i;
	}

	if (n_y >= 1.) {
		for (i=1; i<=n_y; i++) r *= y - i;
	} else {
		for (i=n_y+1.; i<=0.; i++) r /= y - i;
	}

	return r;
}

} // namespace kv

#endif // BETA_HPP
