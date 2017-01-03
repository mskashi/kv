/*
 * Copyright (c) 2014-2015 Masahide Kashiwagi (kashi@waseda.jp)
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
	// set order preparing for constant function
	y = setorder(y, order);

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
	x = setorder(x, order + singularity_order);

	y = f(x);
	// set order preparing for constant function
	y = setorder(y, order + singularity_order);

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
	interval<T> step, result;
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

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1;
	x = setorder(x, order);

	y = f(x);
	// set order preparing for constant function
	y = setorder(y, order);

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


template <class T, class F1, class F2>
interval<T>
defint_singular_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t1, step, result;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	int n_rad;
	bool resized;
	int restart;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1.;
	x = setorder(x, order-1);

	psa< interval<T> >::mode() = 1;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	y = integrate(f2(x));

	// set order preparing for constant function
	y = setorder(y, order);

	radius = 0.;
	n_rad = 0;
	for (i=order; i>=1; i--) {
		m = norm(y.v(i));
		if (m == 0.) continue;
		radius_tmp = std::pow((double)m, 1./i);
		if (radius_tmp > radius) radius = radius_tmp;
		n_rad++;
		if (n_rad == 2) break;
	}
	radius = std::pow((double)epsilon, 1./order) / radius;

	psa< interval<T> >::mode() = 2;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = true;
	// psa< interval<T> >::record_history() = false;
	#endif

	resized = false;
	restart = 0;
	while (true) {
		step = end - start;

		if (subset(step, interval<T>(-radius, radius))) {
			flag = true;
			t1 = end;
		} else {
			flag = false;
			if (step.lower() > 0.) {
				t1 = mid(start + radius);
			} else {
				t1 = mid(start - radius);
			}
			step = t1 - start;
		}

		// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
		psa< interval<T> >::domain() = interval<T>::hull(0., step);

		try {
			y = integrate(f2(x));
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_singular_autostep: evaluation error");
			}
		}

		// result = eval(y, step) - eval(y, 0.);
		// eval(y, 0.) should be 0
		result = eval(y, step);

		if (resized == true) break;

		resized = true;
		m = rad(result) / epsilon;
		if (restart > 0) {
			radius /= std::max(1., std::pow((double)m, 1. / order));
		} else {
			radius /= std::pow((double)m, 1. / order);
		}
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	if (!flag) result += defint_autostep(f1, t1, end, order, epsilon);
	return result;
}


template <class T, class F1, class F2>
interval<T>
defint_power_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, interval<T> power, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t1, step, result, tx, tp;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	int n_rad;
	bool resized;
	int restart;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1.;
	x = setorder(x, order);

	psa< interval<T> >::mode() = 1;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	y = f2(x);

	// set order preparing for constant function
	y = setorder(y, order);

	radius = 0.;
	n_rad = 0;
	for (i=order; i>=1; i--) {
		m = norm(y.v(i)) / ((double)mid(power)+i+1);
		if (m == 0.) continue;
		radius_tmp = std::pow((double)m, 1./((double)mid(power)+i+1));
		if (radius_tmp > radius) radius = radius_tmp;
		n_rad++;
		if (n_rad == 2) break;
	}
	radius = std::pow((double)epsilon, 1./((double)mid(power)+order+1)) / radius;

	psa< interval<T> >::mode() = 2;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = true;
	// psa< interval<T> >::record_history() = false;
	#endif

	resized = false;
	restart = 0;
	while (true) {
		step = end - start;

		if (subset(step, interval<T>(-radius, radius))) {
			flag = true;
			t1 = end;
		} else {
			flag = false;
			if (step.lower() > 0.) {
				t1 = mid(start + radius);
			} else {
				t1 = mid(start - radius);
			}
			step = t1 - start;
		}

		// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
		psa< interval<T> >::domain() = interval<T>::hull(0., step);

		try {
			y = f2(x);
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_singular_autostep: evaluation error");
			}
		}

		result = 0.;
		tp = power + 1.;
		tx = pow(step, tp);
		for (i=0; i<y.v.size(); i++) {
			result += y.v(i) * tx / tp;
			tp += 1.;
			tx *= step;
		}

		if (resized == true) break;

		resized = true;
		m = rad(result) / epsilon;
		if (restart > 0) {
			radius /= std::max(1., std::pow((double)m, 1. /  ((double)mid(power)+order+1)));
		} else {
			radius /= std::pow((double)m, 1. / ((double)mid(power)+order+1));
		}
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	if (!flag) result += defint_autostep(f1, t1, end, order, epsilon);
	return result;
}


template <class T, class F1, class F2>
interval<T>
defint_log_singular_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t1, step, result;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	int n_rad;
	bool resized;
	int restart;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1.;
	x = setorder(x, order + singularity_order);

	psa< interval<T> >::mode() = 1;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif

	y = f2(x);
	// set order preparing for constant function
	y = setorder(y, order + singularity_order);

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

	radius = 0.;
	n_rad = 0;
	for (i=order+1; i>=1; i--) {
		m = norm(y.v(i));
		if (m == 0.) continue;
		radius_tmp = std::pow((double)m, 1./i);
		if (radius_tmp > radius) radius = radius_tmp;
		n_rad++;
		if (n_rad == 2) break;
	}
	radius = std::pow((double)epsilon, 1./order) / radius;

	psa< interval<T> >::mode() = 2;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = true;
	// psa< interval<T> >::record_history() = false;
	#endif

	resized = false;
	restart = 0;
	while (true) {
		step = end - start;

		if (subset(step, interval<T>(-radius, radius))) {
			flag = true;
			t1 = end;
		} else {
			flag = false;
			if (step.lower() > 0.) {
				t1 = mid(start + radius);
			} else {
				t1 = mid(start - radius);
			}
			step = t1 - start;
		}

		// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
		psa< interval<T> >::domain() = interval<T>::hull(0., step);

		try {
			y = f2(x);
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_singular_autostep: evaluation error");
			}
		}

		// set order preparing for constant function
		y = setorder(y, order + singularity_order);

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

		if (resized == true) break;

		resized = true;
		m = rad(result) / epsilon;
		if (restart > 0) {
			radius /= std::max(1., std::pow((double)m, 1. / (order+1)));
		} else {
			radius /= std::pow((double)m, 1. / (order+1));
		}
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	if (!flag) result += defint_autostep(f1, t1, end, order, epsilon);
	return result;
}


template <class T, class F1, class F2>
interval<T>
defint_log_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t1, step, result;
	interval<T> xn, logx;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	int n_rad;
	bool resized;
	int restart;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	x.v.resize(2);
	x.v(0) = start;
	x.v(1) = 1.;
	x = setorder(x, order);

	psa< interval<T> >::mode() = 1;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	y = f2(x);

	// set order preparing for constant function
	y = setorder(y, order);

	radius = 0.;
	n_rad = 0;
	for (i=order; i>=1; i--) {
		m = norm(y.v(i)) * std::fabs((i+1)*std::log(0.1)-1.)/(i+1.)/(i+1.);
		if (m == 0.) continue;
		radius_tmp = std::pow((double)m, 1./(i+1));
		if (radius_tmp > radius) radius = radius_tmp;
		n_rad++;
		if (n_rad == 2) break;
	}
	radius = std::pow((double)epsilon, 1./(order+1)) / radius;

	psa< interval<T> >::mode() = 2;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = true;
	// psa< interval<T> >::record_history() = false;
	#endif

	resized = false;
	restart = 0;
	while (true) {
		step = end - start;

		if (subset(step, interval<T>(-radius, radius))) {
			flag = true;
			t1 = end;
		} else {
			flag = false;
			if (step.lower() > 0.) {
				t1 = mid(start + radius);
			} else {
				t1 = mid(start - radius);
			}
			step = t1 - start;
		}

		// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
		psa< interval<T> >::domain() = interval<T>::hull(0., step);

		try {
			y = f2(x);
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_singular_autostep: evaluation error");
			}
		}

		result = 0;
		logx = log(step);
		xn = step;
		for (i=0; i<y.v.size(); i++) {
			result += y.v(i) * xn * ((i+1) * logx - 1) / ((i+1) * (i+1));
			xn *= step;
		}

		if (resized == true) break;

		resized = true;
		m = rad(result) / epsilon;
		if (restart > 0) {
			radius /= std::max(1., std::pow((double)m, 1. / (order+1)));
		} else {
			radius /= std::pow((double)m, 1. / (order+1));
		}
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	if (!flag) result += defint_autostep(f1, t1, end, order, epsilon);
	return result;
}

} // namespace kv

#endif // DEFINT_SINGULAR_HPP
