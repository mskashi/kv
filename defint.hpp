#ifndef DEFINT_HPP
#define DEFINT_HPP

// Definit Integration

#include <cmath>
#include <limits>
#include "interval.hpp"
#include "rdouble.hpp"
#include "interval-vector.hpp"
#include "psa.hpp"

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
interval<T>
defint(F f, interval<T> start, interval<T> end, int order, int n) {
	interval<T> r, step, c, z, result;
	psa< interval<T> > x, y;
	int i, j, s;
	bool save_mode, save_uh, save_rh;

	step = (end - start) / (T)n;
	r = step * 0.5;

	x.v.resize(2);
	// x.v(0) = 0.;
	x.v(1) = 1.;
	x = setorder(x, order);

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::mode() = 2;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	// psa< interval<T> >::domain() = interval<T>(-r.upper(), r.upper());
	psa< interval<T> >::domain() = interval<T>::hull(-r, r);

	result = 0.;

	for (i=0; i<n; i++) {
		c = start + (T)i * step + r;
		x.v(0) = c;
		y = integrate(f(x));

		z = eval(y, r) - eval(y, -r);

		#if 0
		s = y.v.size() - 1;
		z = (s % 2 == 0) ? y.v(s) - y.v(s) : y.v(s) * 2.;
		for (j=s-1; j>=0; j--) {
			if (j % 2 == 0) {
				z *= r;
			} else {
				z = z * r + y.v(j) * 2.;
			}
		}
		#endif
		
		result += z;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

template <class T, class F>
interval<T>
defint_autostep(F f, interval<T> start, interval<T> end, int order) {
	interval<T> t, t1, z, step, result;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	int n_rad;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	x.v.resize(2);
	// x.v(0) = 0.;
	x.v(1) = 1.;
	x = setorder(x, order-1);

	t = start;
	result = 0.;

	while (1) {
		flag = false;

		x.v(0) = t;
		psa< interval<T> >::mode() = 1;
		#ifdef DEFINT_FAST
		psa< interval<T> >::use_history() = false;
		psa< interval<T> >::record_history() = true;
		#endif
		y = integrate(f(x));

		// 念のため
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
#ifdef FORCE_STEP
		radius = FORCE_STEP;
#else
#ifdef EXP_STEP
		radius = std::exp(-EXP_STEP) / radius;
#else
		radius = std::pow((double)std::numeric_limits<T>::epsilon()*order/(order-1), 1./order) / radius;
#endif
#endif

		step = end - t;

		if (subset(step, interval<T>(-radius, radius))) {
			t1 = end;
			flag = true;
		} else {
			if (step.lower() > 0.) {
				t1 = mid(t + radius);
			} else {
				t1 = mid(t - radius);
			}
			step = t1 - t;
		}

		// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
		psa< interval<T> >::domain() = interval<T>::hull(0., step);

		x.v(0) = t;
		psa< interval<T> >::mode() = 2;
		#ifdef DEFINT_FAST
		psa< interval<T> >::use_history() = true;
		psa< interval<T> >::record_history() = false;
		#endif
		y = integrate(f(x));
		// z = eval(y, step) - eval(y, 0.);
		// eval(y, 0.) == 0のはず
		z = eval(y, step);
		result += z;
		if (flag) break;
		t = t1;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return result;
}

} // namespace kv

#endif // DEFINT_HPP
