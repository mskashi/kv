/*
 * Copyright (c) 2014-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DEFINT_SINGULAR_HPP
#define DEFINT_SINGULAR_HPP

#include <kv/defint.hpp>

namespace kv {

/*
 * helper functions for definite integral with singular point
 */

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

// The function f(x) = x
struct Defint_x {
	template <class T> T operator()(const T& x) {
		return x;
	}
};

#if 0
// Not Used
// The function f(x) = -x
struct Defint_mx {
	template <class T> T operator()(const T& x) {
		return -x;
	}
};
#endif

// The function f(x) = 1
struct Defint_1 {
	template <class T> T operator()(const T& x) {
		return (T)1.;
	}
};

// Reverse Time : convert f(x) to f(-x)
template <class F> class Defint_Reverse {
	F f;
	public:
	Defint_Reverse(F f) : f(f) {}
	template <class T> T operator()(const T& x) {
		return f(-x);
	}
};

#if 0
// Not Used
template <class F1, class F2, class TT> class Defint_power_func {
	F1 f;
	F2 g;
	TT p;
	public:
	Defint_power_func(F1 f, F2 g, TT p) : f(f), g(g), p(p) {}
	template <class T> T operator() (const T& x) {
		return pow(f(x), (T)p) * g(x);
	}
};
#endif

#if 0
// Not Used
template <class F1, class F2> class Defint_log_func {
	F1 f;
	F2 g;
	public:
	Defint_log_func(F1 f, F2 g) : f(f), g(g) {}
	template <class T> T operator() (const T& x) {
		return log(f(x)) * g(x);
	}
};
#endif


/*
 * verified integrator with fixed step
 */

// calculate \int_start^end f(x)^power * g(x) dx
// (start < end) is needed.
// n = singularity order
// lim_{x->start} f(x)/(x-start)^i = 0 (0 \le i < n)
// lim_{x->start} f(x)/(x-start)^n != 0

template <class T, class F1, class F2>
interval<T>
defint_power3(F1 f, F2 g, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1) {
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

	y = pow(div_tn(f(x), singularity_order), power) * g(x);
	// set order preparing for constant function
	y = setorder(y, order);

	result = 0.;
	tp = power * singularity_order + 1.;
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

// calculate \int_start^end (x-start)^power * f(x) dx
// (start < end) is needed.

template <class T, class F>
interval<T>
defint_power(F f, interval<T> start, interval<T> end, int order, interval<T> power) {
	return defint_power3(Defint_x(), f, start, end, order, power, 1);
}

// calculate \int_start^end f(x)^power dx
// (start < end) is needed.
// n = singularity order
// lim_{x->start} f(x)/(x-start)^i = 0 (0 \le i < n)
// lim_{x->start} f(x)/(x-start)^n != 0

template <class T, class F>
interval<T>
defint_power2(F f, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1) {
	return defint_power3(f, Defint_1(), start, end, order, power, singularity_order);
}


// calculate \int_start^end log(f(x))g(x) dx
// (start < end) needed
// n = singularity order
// lim_{x->start} f(x)/(x-start)^i = 0 (0 \le i < n)
// lim_{x->start} f(x)/(x-start)^n != 0

template <class T, class F1, class F2>
interval<T>
defint_log3(F1 f, F2 g, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
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
	x = setorder(x, order + singularity_order);

	y = f(x);
	// set order preparing for constant function
	y = setorder(y, order + singularity_order);

	result = eval(integrate(log(div_tn(y, singularity_order)) * g(x)), step);

	x = setorder(x, order);
	y = singularity_order * g(x);
	// set order preparing for constant function
	y = setorder(y, order);

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


// calculate \int_start^end log(x-start)f(x) dx
// (start < end) needed

template <class T, class F>
interval<T>
defint_log(F f, interval<T> start, interval<T> end, int order) {
	return defint_log3(Defint_x(), f, start, end, order, 1);
}


// calculate \int_start^end log(f(x)) dx
// (start < end) needed
// n = singularity order
// lim_{x->start} f(x)/(x-start)^i = 0 (0 \le i < n)
// lim_{x->start} f(x)/(x-start)^n != 0

template <class T, class F>
interval<T>
defint_log2(F f, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
	return defint_log3(f, Defint_1(), start, end, order, singularity_order);
}


/*
 * defint for functions which has (removable) singular point at "start".
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

/*
 * verified integrator with automatic step size control
 */

template <class T, class F1, class F2, class F3>
interval<T>
defint_power3_autostep(F1 h, F2 f, F3 g, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
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

	y = pow(div_tn(f(x), singularity_order), power) * g(x);
	// set order preparing for constant function
	y = setorder(y, order);

	radius = 0.;
	n_rad = 0;
	for (i=order; i>=1; i--) {
		#ifdef DEFINT_STEPSIZE_MAG
		m = mag(y.v(i) / ((double)mid(power)+i+1));
		#else
		m = mig(y.v(i) / ((double)mid(power)+i+1));
		#endif
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
			y = pow(div_tn(f(x), singularity_order), power) * g(x);
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_power3_autostep: evaluation error");
			}
		}

		result = 0.;
		tp = power * singularity_order + 1.;
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

	#ifdef DEFINT_SHOW_STEPSIZE
	std::cout << "stepsize: " << step << "\n";
	#endif

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	// if (!flag) result += defint_autostep(Defint_power_func<F2, F3, interval<T> >(f, g, power), t1, end, order, epsilon);
	if (!flag) result += defint_autostep(h, t1, end, order, epsilon);
	return result;
}

template <class T, class F1, class F2>
interval<T>
defint_power_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, interval<T> power, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_power3_autostep(f1, Defint_x(), f2, start, end, order, power, 1, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_power2_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_power3_autostep(f1, f2, Defint_1(), start, end, order, power, singularity_order, epsilon);
}



template <class T, class F1, class F2, class F3>
interval<T>
defint_log3_autostep(F1 h, F2 f, F3 g, interval<T> start, interval<T> end, int order, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t1, step, result;
	interval<T> xn, logx;
	psa< interval<T> > x, y, y2;
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

	psa< interval<T> >::mode() = 1;
	#if DEFINT_FAST == 1
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif

	x = setorder(x, order + singularity_order);
	y = f(x);
	// set order preparing for constant function
	y = setorder(y, order + singularity_order);
	y = integrate(log(div_tn(y, singularity_order)) * g(x));

	x = setorder(x, order);
	y2 = singularity_order * g(x);
	// set order preparing for constant function
	y2 = setorder(y2, order);

	radius = 0.;
	n_rad = 0;
	for (i=order; i>=1; i--) {
		#ifdef DEFINT_STEPSIZE_MAG
		m = mag(y2.v(i) * ((i+1)*std::log(0.1)-1.)/(i+1.)/(i+1.));
		m = std::max(m, mag(y.v(i+1)));
		#else
		m = mig(y2.v(i) * ((i+1)*std::log(0.1)-1.)/(i+1.)/(i+1.));
		m = std::max(m, mig(y.v(i+1)));
		#endif
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
			x = setorder(x, order + singularity_order);
			y = f(x);
			result = eval(integrate(log(div_tn(y, singularity_order)) * g(x)), step);
			x = setorder(x, order);
			y = singularity_order * g(x);
		}
		catch (std::domain_error& e) {
			if (restart < RESTART_MAX) {
				psa< interval<T> >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
			} else {
				throw std::domain_error("defint_log3_autostep: evaluation error");
			}
		}

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

	#ifdef DEFINT_SHOW_STEPSIZE
	std::cout << "stepsize: " << step << "\n";
	#endif

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	// if (!flag) result += defint_autostep(Defint_log_func<F2, F3>(f, g), t1, end, order, epsilon);
	if (!flag) result += defint_autostep(h, t1, end, order, epsilon);
	return result;
}

template <class T, class F1, class F2>
interval<T>
defint_log_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_log3_autostep(f1, Defint_x(), f2, start, end, order, 1, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_log2_autostep(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_log3_autostep(f1, f2, Defint_1(), start, end, order, singularity_order, epsilon);
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
		#ifdef DEFINT_STEPSIZE_MAG
		m = mag(y.v(i));
		#else
		m = mig(y.v(i));
		#endif
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

	#ifdef DEFINT_SHOW_STEPSIZE
	std::cout << "stepsize: " << step << "\n";
	#endif

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	if (!flag) result += defint_autostep(f1, t1, end, order, epsilon);
	return result;
}


/*
 * integrators with right hand side singularity
 */

template <class T, class F1, class F2>
interval<T>
defint_power3_r(F1 f, F2 g, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1) {
	return defint_power3(Defint_Reverse<F1>(f), Defint_Reverse<F2>(g), -end, -start, order, power, singularity_order);
}

template <class T, class F>
interval<T>
defint_power_r(F f, interval<T> start, interval<T> end, int order, interval<T> power) {
	return defint_power3(Defint_x(), Defint_Reverse<F>(f), -end, -start, order, power, 1);
}

template <class T, class F>
interval<T>
defint_power2_r(F f, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1) {
	return defint_power3(Defint_Reverse<F>(f), Defint_1(), -end, -start, order, power, singularity_order);
}

template <class T, class F1, class F2>
interval<T>
defint_log3_r(F1 f, F2 g, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
	return defint_log3(Defint_Reverse<F1>(f), Defint_Reverse<F2>(g), -end, -start, order, singularity_order);
}

template <class T, class F>
interval<T>
defint_log_r(F f, interval<T> start, interval<T> end, int order) {
	return defint_log3(Defint_x(), Defint_Reverse<F>(f), -end, -start, order, 1);
}

template <class T, class F>
interval<T>
defint_log2_r(F f, interval<T> start, interval<T> end, int order, int singularity_order = 1) {
	return defint_log3(Defint_Reverse<F>(f), Defint_1(), -end, -start, order, singularity_order);
}

template <class T, class F>
interval<T>
defint_singular_r(F f, interval<T> start, interval<T> end, int order) {
	return defint_singular(Defint_Reverse<F>(f), -end, -start, order);
}

template <class T, class F1, class F2, class F3>
interval<T>
defint_power3_autostep_r(F1 h, F2 f, F3 g, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_power3_autostep(Defint_Reverse<F1>(h), Defint_Reverse<F2>(f), Defint_Reverse<F3>(g), -end, -start, order, power, singularity_order, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_power_autostep_r(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, interval<T> power, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_power3_autostep(Defint_Reverse<F1>(f1), Defint_x(), Defint_Reverse<F2>(f2), -end, -start, order, power, 1, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_power2_autostep_r(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, interval<T> power, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_power3_autostep(Defint_Reverse<F1>(f1), Defint_Reverse<F2>(f2), Defint_1(), -end, -start, order, power, singularity_order, epsilon);
}

template <class T, class F1, class F2, class F3>
interval<T>
defint_log3_autostep_r(F1 h, F2 f, F3 g, interval<T> start, interval<T> end, int order, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_log3_autostep(Defint_Reverse<F1>(h), Defint_Reverse<F2>(f), Defint_Reverse<F3>(g), -end, -start, order, singularity_order, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_log_autostep_r(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_log3_autostep(Defint_Reverse<F1>(f1), Defint_x(), Defint_Reverse<F2>(f2), -end, -start, order, 1, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_log2_autostep_r(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, int singularity_order = 1, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_log3_autostep(Defint_Reverse<F1>(f1), Defint_Reverse<F2>(f2), Defint_1(), -end, -start, order, singularity_order, epsilon);
}

template <class T, class F1, class F2>
interval<T>
defint_singular_autostep_r(F1 f1, F2 f2, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	return defint_singular_autostep(Defint_Reverse<F1>(f1), Defint_Reverse<F2>(f2), -end, -start, order, epsilon);
}


} // namespace kv

#endif // DEFINT_SINGULAR_HPP
