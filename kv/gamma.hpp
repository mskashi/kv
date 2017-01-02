/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef GAMMA_HPP
#define GAMMA_HPP

#include <kv/defint.hpp>

// Lobachevsky function

namespace kv {

template <class TT> class Gamma {
	public:
	TT x;
	
	Gamma(TT x) : x(x) {}

	template <class T> T operator() (T t) {
		return pow(t, (T)x - 1) * exp(-t);
	}
};

class Gamma_nopower {
	public:
	template <class T> T operator() (T t) {
		return exp(-t);
	}
};

// calculate \int_start^end x^power * f(x) dx

template <class T, class F>
interval<T>
defint_power(F f, interval<T> start, interval<T> end, int order, interval<T> power) {
	int i;
	interval<T> step, c, result, tx, tp;
	psa< interval<T> > x, y;
	bool save_mode, save_uh, save_rh;

	step = end - start;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::mode() = 2;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	psa< interval<T> >::domain() = interval<T>::hull(0, step);

	c = start;
	x.v.resize(2);
	x.v(0) = c;
	x.v(1) = 1;
	x = setorder(x, order);

	y = f(x);

	result = 0.;
	tp = power + 1.;
	tx = pow(step, tp);
	for (i=0; i<y.v.size(); i++) {
		result += y.v(i) * tx / tp;
		tp += 1.;
		tx *= step;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

#define GAMMA_TH1 0.25
#define GAMMA_TH2 40.
#define GAMMA_ORDER 14

// 1 \le x \le 2
template <class T>
interval<T> gamma_r(const interval<T>& x) {
	interval<T> result;

	result = defint_power(Gamma_nopower(), interval<T>(0.), interval<T>(GAMMA_TH1), GAMMA_ORDER, x - 1);

	result += defint_autostep(Gamma< interval<T> >(x), interval<T>(GAMMA_TH1), interval<T>(GAMMA_TH2), GAMMA_ORDER);

	result += interval<T>(0., (exp(-interval<T>(GAMMA_TH2)) * (interval<T>(GAMMA_TH2) + 1.)).upper());

	return result;
}

template <class T>
interval<T> gamma_pluspoint(const interval<T>& x) {
	T y;
	int i;
	interval<T> r;

	y = floor(x.lower());
	if (y == 0.) {
		return gamma_r(x + 1.) / x;
	} else if (y == 1.) {
		return gamma_r(x);
	} else {
		r = gamma_r(x - (y - 1.));
		for (i=y-1; i>=1; i--) {
			r *= x - i;
		}
		return r;
	}
}

template <class T>
interval<T> gamma_plus(const interval<T>& x) {
	static const interval<T> m = constants< interval<T> >::str("1.46163214496836234126265954232572132846819620400644635129598840859878644035380181024307499273372559");
	interval<T> tmp;

	tmp = interval<T>::hull(gamma_pluspoint(interval<T>(x.lower())), gamma_pluspoint(interval<T>(x.upper())));

	if (overlap(x, m)) {
		return interval<T>::hull(tmp, gamma_pluspoint(m));
	} else {
		return tmp;
	} 
}

template <class T>
interval<T> gamma(const interval<T>& x) {
	T tmp;
	interval<T> y;
	static const interval<T> pi = constants< interval<T> >::pi();

	if (in(0., x)) {
		return interval<T>::whole();
	}

	if (x.upper() < 0.) {
		tmp = -floor(-x.upper());
		y = x - tmp;
		if (in(0., y) || in(-1., y)) {
			return interval<T>::whole();
		}
		return - pi / (sin(pi * x) * gamma_plus(-x) * x);
	} else {
		return gamma_plus(x);
	}
}

} // namespace kv

#endif // GAMMA_HPP
