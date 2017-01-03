/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef LOBACHEVSKY_HPP
#define LOBACHEVSKY_HPP

#include <kv/defint.hpp>
#include <kv/defint-singular.hpp>

// Lobachevsky function

namespace kv {

struct Lobachevsky {
	template <class T> T operator() (T x) {
		return log(2. * sin(x));
	}
};

struct Lobachevsky_nolog {
	public:
	template <class T> T operator() (T x) {
		return 2. * sin(x);
	}
};

// #define LOBACHEVSKY_TH 0.25

#ifndef LOBACHEVSKY_ORDER
#define LOBACHEVSKY_ORDER 14
#endif

// 0 \le x \le \pi/2
template <class T>
interval<T> lobachevsky_origin(const interval<T>& x) {
	if (x == 0.) return interval<T>(0.);
	#if 0
	if (x.lower() <= LOBA_TH) {
		return -defint_log2(Lobachevsky_nolog(), interval<T>(0.), x, LOBACHEVSKY_ORDER);
	} else {
		return -defint_log2(Lobachevsky_nolog(), interval<T>(0.), interval<T>(LOBACHEVSKY_TH), LOBACHEVSKY_ORDER) - defint_autostep(Lobachevsky(), interval<T>(LOBA_TH), x, LOBA_ORDER);
	}
	#endif
	return -defint_log2_autostep(Lobachevsky(), Lobachevsky_nolog(), interval<T>(0.), x, LOBACHEVSKY_ORDER);
}

template <class T>
interval<T> lobachevsky_point(const interval<T>& x) {
	static const interval<T> pi = constants< interval<T> >::pi();

	if (x.lower() >= (pi * 0.5).upper()) {
		return lobachevsky_point(x - pi);
	}

	if (zero_in(x)) {
		return interval<T>::hull(-lobachevsky_origin(interval<T>(0., -x.lower())), lobachevsky_origin(interval<T>(0.,x.upper())));
	}
	if (x.lower() > 0.) {
		return lobachevsky_origin(x);
	} else {
		return -lobachevsky_origin(-x);
	}
}

template <class T>
interval<T> lobachevsky(const interval<T>& x) {
	static const interval<T> pi = constants< interval<T> >::pi();
	static const interval<T> pi12 = pi * 0.5;
	static const interval<T> pi16 = pi / 6.;
	static const interval<T> pi56 = pi16 * 5.;
	static const interval<T> pi76 = pi16 * 7.;

	T n;
	interval<T> x2, r, m;
	m = lobachevsky_point(pi16);

	x2 = x;
	while (x2.lower() <= -pi12.upper() || x2.lower() >= pi12.upper()) {
		n = floor((x2.lower() / pi.lower()) + 0.5);
		x2 -= n * pi;
	}

	if (interval<T>(x2.upper()) - x2.lower() >= pi.upper()) {
		return interval<T>(-m.upper(), m.upper());
	}

	r = interval<T>::hull(lobachevsky_point(interval<T>(x2.lower())), lobachevsky_point(interval<T>(x2.upper())));

	if (subset(-pi16, x2)) {
		r = interval<T>::hull(r, -m);
	}
	if (subset(pi16, x2)) {
		r = interval<T>::hull(r, m);
	}
	if (subset(pi56, x2)) {
		r = interval<T>::hull(r, -m);
	}
	if (subset(pi76, x2)) {
		r = interval<T>::hull(r, m);
	}

	return r;
}

} // namespace kv

#endif // LOBACHEVSKY_HPP
