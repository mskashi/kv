/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DEFINT_HPP
#define DEFINT_HPP

// Definit Integration

#include <cmath>
#include <limits>
#include <algorithm>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>


#ifndef DEFINT_FAST
#define DEFINT_FAST 1
#endif

#ifndef RESTART_MAX
#define RESTART_MAX 10
#endif


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
defint_autostep(F f, interval<T> start, interval<T> end, int order, T epsilon = std::numeric_limits<T>::epsilon()) {
	interval<T> t, t1, z, step, result;
	psa< interval<T> > x, y;
	int i;
	bool flag;
	bool save_mode, save_uh, save_rh;

	T radius, radius_tmp, m;
	T tolerance;
	int n_rad;
	bool resized;
	int restart;

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
		tolerance = std::max((T)1., norm(result)) * epsilon;

		x.v(0) = t;
		psa< interval<T> >::mode() = 1;
		#if DEFINT_FAST == 1
		psa< interval<T> >::use_history() = false;
		psa< interval<T> >::record_history() = true;
		psa< interval<T> >::history().clear();
		#endif
		y = integrate(f(x));

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
		radius = std::pow((double)tolerance, 1./order) / radius;

		// x.v(0) = t;
		psa< interval<T> >::mode() = 2;
		#if DEFINT_FAST == 1
		psa< interval<T> >::use_history() = true;
		// psa< interval<T> >::record_history() = false;
		#endif

		resized = false;
		restart = 0;
		while (true) {
			step = end - t;

			if (subset(step, interval<T>(-radius, radius))) {
				flag = true;
				t1 = end;
			} else {
				flag = false;
				if (step.lower() > 0.) {
					t1 = mid(t + radius);
				} else {
					t1 = mid(t - radius);
				}
				step = t1 - t;
			}

			// psa< interval<T> >::domain() = interval<T,P>(0., step.upper());
			psa< interval<T> >::domain() = interval<T>::hull(0., step);

			try {
				y = integrate(f(x));
			}
			catch (std::domain_error& e) {
				if (restart < RESTART_MAX) {
					psa< interval<T> >::use_history() = false;
					radius *= 0.5;
					restart++;
					continue;
				} else {
					throw std::domain_error("defint_autostep: evaluation error");
				}
			}

			// z = eval(y, step) - eval(y, 0.);
			// eval(y, 0.) should be 0
			z = eval(y, step);

			if (resized == true) break;

			resized = true;
			m = rad(z) / tolerance;
			if (restart > 0) {
				radius /= std::max(1., std::pow((double)m, 1. / order));
			} else {
				radius /= std::pow((double)m, 1. / order);
			}
		}
		#ifdef DEFINT_SHOW_STEPSIZE
		std::cout << "stepsize: " << step << "\n";
		#endif

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
