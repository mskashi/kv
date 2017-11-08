/*
 * Copyright (c) 2013-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_HPP
#define ODE_HPP

// ODE

#include <iostream>
#include <cmath>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/make-candidate.hpp>
#include <kv/psa.hpp>
#include <kv/ode-param.hpp>

#ifndef ODE_FAST
#define ODE_FAST 1
#endif

#ifndef ODE_STEP_COMPONENT
#define ODE_STEP_COMPONENT 0
#endif

#ifndef ODE_RESTART_RATIO
#define ODE_RESTART_RATIO 1
#endif

#ifndef ODE_COEF_MID
#define ODE_CORF_MID 0
#endif


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T, class F>
int
ode(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, const ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL) {
	int n = init.size();
	int i, j;

	ub::vector< psa< interval<T> > > x, y;
	psa< interval<T> > torg;
	psa< interval<T> > t;

	ub::vector< psa< interval<T> > > z, w;

	psa< interval<T> > temp;
	T m;
	ub::vector<T> newton_step;

	bool flag, resized;

	interval<T> deltat;
	ub::vector< interval<T> > result;

	T radius, radius_tmp;

	#if ODE_STEP_COMPONENT == 1
	ub::vector<T> tolerance(n);
	#else
	T tolerance;
	#endif

	int n_rad;

	#if ODE_RESTART_RATIO == 1
	T max_ratio;
	#endif

	int ret_val;
	interval<T> end2;
	int restart;

	bool save_mode, save_uh, save_rh;

	#if ODE_STEP_COMPONENT == 1
	for (i=0; i<n; i++) {
		tolerance(i) = std::max(T(1.), norm(init(i))) * p.epsilon;
	}
	#else
	m = 1.;
	for (i=0; i<n; i++) {
		m = std::max(m, norm(init(i)));
	}
	tolerance = m * p.epsilon;
	#endif

	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::mode() = 1;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;
	#if ODE_FAST == 1
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	for (j=0; j<p.order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< interval<T> >::use_history() = true;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) {
			y(i) = integrate(y(i));
			// set order preparing for constant function
			y(i) = setorder(y(i), j+1);
		}
		x = init + y;
	}

	if (p.autostep) {
		// use two non-zero coefficients of higher order term
		#if ODE_STEP_COMPONENT == 1
		radius = std::numeric_limits<T>::infinity();
		for (i=0; i<n; i++) {
			radius_tmp = 0.;
			n_rad = 0;
			for (j = p.order; j>=1; j--) {
				#if ODE_COEF_MID == 1
				using std::abs;
				m = abs(mid(x(i).v(j)));
				#else
				m = norm(x(i).v(j));
				#endif
				if (m == 0.) continue;
				radius_tmp = std::max(radius_tmp, (T)std::pow((double)m, 1./j));
				n_rad++;
				if (n_rad == 2) break;
			}
			radius = std::min(radius, std::pow((double)(tolerance(i)), 1./p.order) / radius_tmp);
		}
		#else // ODE_STEP_COMPONENT
		radius = 0.;
		n_rad = 0;
		for (j = p.order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				#if ODE_COEF_MID == 1
				using std::abs;
				m = std::max(m, abs(mid(x(i).v(j))));
				#else
				m = std::max(m, norm(x(i).v(j)));
				#endif
			}
			if (m == 0.) continue;
			radius = std::max(radius, (T)std::pow((double)m, 1./j));
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
		#endif // ODE_STEP_COMPONENT
	}

	psa< interval<T> >::mode() = 2;

	restart = 0;
	resized = false;

	while (true) {
		if (p.autostep) {
			end2 = mid(start + radius);
			if (end2 >= end.lower()) {
				end2 = end;
				ret_val = 2;
			} else {
				ret_val = 1;
			}
		} else {
			end2 = end;
			ret_val = 2;
		}
		deltat = end2 - start;

		psa< interval<T> >::domain() = interval<T>(0., deltat.upper());

		z = x;
		t = setorder(torg, p.order);

		try {
			w = f(z, t);
		}
		catch (std::domain_error& e) {
			if (restart < p.restart_max) {
				psa< interval<T> >::use_history() = false;
				if (p.verbose == 1) {
					std::cout << "ode: radius changed: " << radius;
				}
				radius *= 0.5;
				if (p.verbose == 1) {
					std::cout << " -> " << radius << "\n";
				}
				restart++;
				continue;
			} else {
				throw std::domain_error("ode: evaluation error");
			}
		}

		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, p.order);
		}
		w = init + w;

		newton_step.resize(n);
		for (i=0; i<n; i++) {
			newton_step(i) = norm(w(i).v(p.order) - z(i).v(p.order));
		}
		make_candidate(newton_step);
		for (i=0; i<n; i++) {
			z(i).v(p.order) += newton_step(i) * interval<T>(-1., 1.);
		}

		if (p.autostep && resized == false) {
			resized = true;
			m = 0.;
			#if ODE_STEP_COMPONENT == 1
			for (i=0; i<n; i++) {
				m = std::max(m, (rad(eval(z(i), deltat)) - rad(init(i))) / tolerance(i));
			}
			#else
			for (i=0; i<n; i++) {
				m = std::max(m, rad(eval(z(i), deltat)) - rad(init(i)));
			}
			m = m / tolerance;
			#endif
			if (restart > 0) {
				radius /= std::max(1., std::pow((double)m, 1. / p.order));
			} else {
				radius /= std::pow((double)m, 1. / p.order);
			}
			continue;
		}

		w = f(z, t);
		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, p.order);
		}
		w = init + w;

		flag = true;
		#if ODE_RESTART_RATIO == 1
		max_ratio = 0.;
		#endif
		for (i=0; i<n; i++) {
			#if ODE_RESTART_RATIO == 1
			max_ratio = std::max(max_ratio, width(w(i).v(p.order)) / width(z(i).v(p.order)));
			#endif
			flag = flag && subset(w(i).v(p.order), z(i).v(p.order));
		}
		if (flag) break;

		if (!p.autostep || restart >= p.restart_max) {
			ret_val = 0;
			break;
		}
		if (p.verbose == 1) {
			std::cout << "ode: radius changed: " << radius;
		}
		#if ODE_RESTART_RATIO == 1
		radius *= std::max(std::min((T)0.5, (T)0.5 / max_ratio), (T)0.125);
		#else
		radius *= 0.5;
		#endif
		if (p.verbose == 1) {
			std::cout << " -> " << radius << "\n";
		}
		restart++;
	}

	if (ret_val != 0) {
		for (j=0; j<p.iteration; j++) {
			z = w;
			w = f(z, t);
			for (i=0; i<n; i++) {
				temp = integrate(w(i));
				w(i) = setorder(temp, p.order);
			}
			w = init + w;
			for (i=0; i<n; i++) {
				w(i).v(p.order) = intersect(w(i).v(p.order), z(i).v(p.order));
			}
		}

		result.resize(n);
		for (i=0; i<n; i++) {
			result(i) = eval(w(i), deltat);
		}

		init = result;
		if (ret_val == 1) end = end2;
		if (result_psa != NULL) *result_psa = w;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return ret_val;
}

template <class T, class F>
int
odelong(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>()) {

	ub::vector< interval<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;

		r = ode(f, x, t, t1, p);
		if (r == 0) {
			if (ret_val == 1) {
				init = x;
				end = t;
			}
			return ret_val;
		}
		ret_val = 1;
		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << x << "\n";
		}
		if (r == 2) {
			init = x;
			return 2;
		}
		t = t1;
	}
}

} // namespace kv

#endif // ODE_HPP
