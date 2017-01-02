/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef LOBACHEVSKY_HPP
#define LOBACHEVSKY_HPP

#include <kv/defint.hpp>

// Lobachevsky function

namespace kv {

class Loba {
	public:
	template <class T> T operator() (T x) {
		return log(2. * sin(x));
	}
};

class Loba_nolog {
	public:
	template <class T> T operator() (T x) {
		return 2. * sin(x);
	}
};

template <class T, class F>
interval<T>
defint_log(F f, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
	int i;
	interval<T> step, c, z, result;
	psa< interval<T> > x, y;
	bool save_mode, save_uh, save_rh;

	step = end - start;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::mode() = 2;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;


	psa< interval<T> >::domain() = interval<T>(0, step.upper());


	c = start;
	x.v.resize(2);
	x.v(0) = c;
	x.v(1) = 1;
	x = setorder(x, order + singularity_order);

	y = f(x);
	for (i=0; i<singularity_order; i++) {
		if (! zero_in(x.v(i))) {
			std::cout << "defint_log: something wrong\n";
		}
	}

	for (i=singularity_order; i<x.v.size(); i++) {
		y.v(i-singularity_order) = y.v(i);
	}
	y.v.resize(y.v.size() - singularity_order);


	y = integrate(log(y));
	result = eval(y, step) + (T)singularity_order * step * (log(step) - 1.);

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

#define LOBA_TH 0.25
#define LOBA_ORDER 14

// 0 \le x \le \pi/2
template <class T>
interval<T> loba_origin(const interval<T>& x) {
	if (x == 0.) return interval<T>(0.);
	if (x.lower() <= LOBA_TH) {
		return -defint_log(Loba_nolog(), interval<T>(0.), x, LOBA_ORDER);
	} else {
		return -defint_log(Loba_nolog(), interval<T>(0.), interval<T>(LOBA_TH), LOBA_ORDER) - defint_autostep(Loba(), interval<T>(LOBA_TH), x, LOBA_ORDER);
	}
}

template <class T>
interval<T> loba_point(const interval<T>& x) {
	static const interval<T> pi = constants< interval<T> >::pi();

	if (x.lower() >= (pi * 0.5).upper()) {
		return loba_point(x - pi);
	}

	if (zero_in(x)) {
		return interval<T>::hull(-loba_origin(interval<T>(0., -x.lower())), loba_origin(interval<T>(0.,x.upper())));
	}
	if (x.lower() > 0.) {
		return loba_origin(x);
	} else {
		return -loba_origin(-x);
	}
}

template <class T>
interval<T> loba(const interval<T>& x) {
	static const interval<T> pi = constants< interval<T> >::pi();
	static const interval<T> pi12 = pi * 0.5;
	static const interval<T> pi16 = pi / 6.;
	static const interval<T> pi56 = pi16 * 5.;
	static const interval<T> pi76 = pi16 * 7.;

	T n;
	interval<T> x2, r, m;
	m = loba_point(pi16);

	x2 = x;
	while (x2.lower() <= -pi12.upper() || x2.lower() >= pi12.upper()) {
		n = floor((x2.lower() / pi.lower()) + 0.5);
		x2 -= n * pi;
	}

	if (interval<T>(x2.upper()) - x2.lower() >= pi.upper()) {
		return interval<T>(-m.upper(), m.upper());
	}

	r = interval<T>::hull(loba_point(interval<T>(x2.lower())), loba_point(interval<T>(x2.upper())));

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
