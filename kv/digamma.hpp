/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DIGAMMA_HPP
#define DIGAMMA_HPP

#include <kv/interval.hpp>
#include <kv/psa.hpp>
#include <kv/defint.hpp>

namespace kv {

template <class TT> struct Digamma_0 {
	TT x;
	Digamma_0(TT x) : x(x) {}
	template <class T> T operator() (const T& t) {
		return div_tn(exp(-t) - exp(-x * t) / div_tn(1-exp(-t), 1), 1);
	}
};

template <class TT> struct Digamma {
	TT x;
	Digamma(TT x) : x(x) {}
	template <class T> T operator() (const T& t) {
		return (exp(-t) - exp(-x * t) / ((1-exp(-t)) / t)) / t;
	}
};

#define DIGAMMA_TH1 0.125
#define DIGAMMA_TH2 35
#define DIGAMMA_ORDER 14

// work for x > 0
// return high precision output if x is in [1,2]
template <class T> interval<T> digamma_plus(const interval<T>& x) {
	interval<T> result, tmp;
	
	result = defint_singular(Digamma_0< interval<T> >(x), (interval<T>)0., (interval<T>)DIGAMMA_TH1, DIGAMMA_ORDER);

	result += defint_autostep(Digamma< interval<T> >(x), (interval<T>)DIGAMMA_TH1, (interval<T>)DIGAMMA_TH2, DIGAMMA_ORDER);

	tmp = exp(- interval<T>(DIGAMMA_TH2));
	result += interval<T>::hull(- 1. / (1. - tmp) * x * exp(- x * DIGAMMA_TH2), tmp);

	return result;
}


/*
 * digamma for small interval
 */

template <class T> interval<T> digamma_point(const interval<T>& x) {
	T tmp;
	int i;
	interval<T> r;

	tmp = floor(x.lower());

	if (tmp == 1.) {
		return digamma_plus(x);
	} 

	r = digamma_plus(x - tmp + 1.);
	if (tmp > 1.) {
		for (i = tmp-1; i>=1; i--) {
			r += 1. / (x - i);
		}
	} else {
		for (i = tmp-1; i<=-1; i++) {
			r -= 1. / (x - i - 1.);
		}
	}

	return r;
}

template <class T> interval<T> digamma(const interval<T>& x) {
	if (x.lower() < 0. && floor(x.lower()) != floor(x.upper())) {
		return interval<T>::whole();
	}
	return interval<T>::hull(digamma_point(interval<T>(x.lower())), digamma_point(interval<T>(x.upper())));
}

/*
 *  return x(t) / t^n assuming that xi = 0 (0 <= i < n) 
 */

template <class T> psa<T> div_tn(const psa<T>& x, int n) {
	int i;
	int s = x.v.size();
	psa<T> y;

	y.v.resize(s);

	for (i=0; i<s-n; i++) {
		y.v(i) = x.v(i+n);
	}
	for (i=s-n; i<s; i++) {
		y.v(i) = 0.;
	}

	return y;
}

/*
 *  return x(t) / y(t) assuming that xi = yi = 0 (0 <= i < n) 
 */

template <class T> psa<T> div_reduce(const psa<T>& x, const psa<T>& y, int n) {
	int i;
	int sx = x.v.size();
	int sy = y.v.size();
	psa<T> nx, ny;

	nx.v.resize(sx);
	for (i=0; i<sx-n; i++) {
		nx.v(i) = x.v(i+n);
	}
	for (i=sx-n; i<sx; i++) {
		nx.v(i) = 0.;
	}

	ny.v.resize(sy);
	for (i=0; i<sy-n; i++) {
		ny.v(i) = y.v(i+n);
	}
	for (i=sy-n; i<sy; i++) {
		ny.v(i) = 0.;
	}

	return nx / ny;
}


/*
 * defint for functions which has (removable) sigular point at "start".
 * almost same as defint except that expantion point is "start"
 *  (not center of interval).
 */

template <class T, class F>
interval<T>
defint_singular(F f, interval<T> start, interval<T> end, int order) {
	interval<T> step, result;
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

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1;
	x = setorder(x, order);
	y = integrate(f(x));

	result = eval(y, step);

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

} // namespace kv

#endif // DIGAMMA_HPP
