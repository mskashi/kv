/*
 * Copyright (c) 2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DEFINT_SINGULAR_HPP
#define DEFINT_SINGULAR_HPP

#include <kv/defint.hpp>

/*
 * helper functions for definite integral with singular point
 */

namespace kv {

// calculate \int_start^end (x-start)^power * f(x) dx
// (start < end) is needed.

template <class T, class F>
interval<T>
defint_power(F f, interval<T> start, interval<T> end, int order, interval<T> power) {
	int i;
	interval<T> step, result, tx, tp;
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

// calculate \int_start_end log(f(x)) dx
// (start < end) needed
// n = singularity order
// lim_{x->start} f(x)/(x-start)^i = 0 (0 \le i < n)
// lim_{x->start} f(x)/(x-start)^n != 0

template <class T, class F>
interval<T>
defint_log_singular(F f, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
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


	psa< interval<T> >::domain() = interval<T>::hull(0, step);


	c = start;
	x.v.resize(2);
	x.v(0) = c;
	x.v(1) = 1;
	x = setorder(x, order + singularity_order);

	y = f(x);
	for (i=0; i<singularity_order; i++) {
		if (! zero_in(y.v(i))) {
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


// calculate \int_start_end log(x)f(x) dx
// (start < end) needed

template <class T, class F>
interval<T>
defint_log(F f, interval<T> start, interval<T> end, int order) {
	int i;
	interval<T> step, c, result;
	interval<T> xn, logx;
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

	result = 0;
	logx = log(step);
	xn = step;
	for (i=0; i<y.v.size(); i++) {
		result += y.v(i) * xn * ((i+1) * logx - 1) / ((i+1) * (i+1));
		xn *= step;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

} // namespace kv

#endif // DEFINT_SINGULAR_HPP
