/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_LOHNER_HPP
#define ODE_LOHNER_HPP

// ODE by Lohner's method

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
#include <kv/autodif.hpp>
#include <kv/ode-param.hpp>


#ifndef ODE_FAST
#define ODE_FAST 1
#endif


namespace kv {

namespace ub = boost::numeric::ublas;

template <class T, class F>
int
ode_lohner(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL) {
	int n = init.size();
	int i, j;

	ub::vector< psa< interval<T> > > x, y;
	psa< interval<T> > torg;
	psa< interval<T> > t;

	ub::vector< psa< interval<T> > > z, w;

	psa< interval<T> > temp;
	T m;
	ub::vector<T> newton_step;

	bool flag;

	interval<T> deltat;
	ub::vector< interval<T> > result;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;

	int ret_val;
	interval<T> end2;
	int restart;

	bool save_mode, save_uh, save_rh;

	ub::vector< interval<T> > V, V2;
	interval<T> tste;

	m = 1.;
	for (i=0; i<n; i++) {
		m = std::max(m, norm(init(i)));
	}
	tolerance = m * p.epsilon;

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
	for (j=0; j<p.order-1; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< interval<T> >::use_history() = true;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) {
			y(i) = integrate(y(i));
			y(i) = setorder(y(i), j+1);
		}
		x = init + y;
	}

	if (p.autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = p.order-1; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m = std::max(m, norm(x(i).v(j)));
				using std::abs;
				m = std::max(m, abs(mid(x(i).v(j))));
			}
			if (m == 0.) continue;
			radius_tmp = std::pow((double)m, 1./j);
			// std::cout << j << " " << m << " " << radius_tmp << "\n";
			if (radius_tmp > radius) radius = radius_tmp;
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
	}
	
	restart = 0;

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

		tste = interval<T>(start.lower(), end2.upper());

		V = f(init, tste) * interval<T>(0., deltat.upper());
		newton_step = mag(V);
		make_candidate(newton_step);

		V = init;
		for (i=0; i<n; i++) {
			V(i) += newton_step(i) * interval<T>(-1., 1.);
		}
		V2 = init + f(V, tste) * interval<T>(0., deltat.upper());

		flag = true;
		for (i=0; i<n; i++) {
			flag = flag && subset(V2, V);
		}
		if (flag) break;

		if (!p.autostep || restart >= p.restart_max + 10) {
			ret_val = 0;
			break;
		}
		radius *= 0.5;
		restart++;
	}

	if (ret_val != 0) {
		for (j=0; j<p.iteration; j++) {
			V = V2;
			V2 = init + f(V, tste) * interval<T>(0., deltat.upper());
			for (i=0; i<n; i++) {
				V2 = intersect(V2, V);
			}
		}

		z = V2;
		torg.v.resize(2);
		torg.v(0) = tste; torg.v(1) = 1.;

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
			y = f(z, t);
			for (i=0; i<n; i++) {
				y(i) = integrate(y(i));
				y(i) = setorder(y(i), j+1);
			}
			z = init + y;
		}

		for (i=0; i<n; i++) {
			for (j=0; j<p.order; j++) {
				z(i).v(j) = x(i).v(j);
			}
		}

		result.resize(n);
		for (i=0; i<n; i++) {
			result(i) = eval(z(i), deltat);
		}

		init = result;
		if (ret_val == 1) end = end2;
		if (result_psa != NULL) *result_psa = z;
	}
	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return ret_val;
}


template <class T, class F>
int
ode_lohner(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< autodif< interval<T> > > >* result_psa = NULL) {
	int n = init.size();
	int i, j, k, km;

	ub::vector< psa< autodif< interval<T> > > > x, y;
	psa< autodif< interval<T> > > torg;
	psa< autodif< interval<T> > > t;

	ub::vector< psa< autodif< interval<T> > > > z, w;

	psa< autodif< interval<T> > > temp;
	T m;
	ub::vector<T> newton_step;

	bool flag;

	interval<T> deltat;
	ub::vector< autodif< interval<T> > > result, new_init;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;

	int ret_val;
	interval<T> end2;
	int restart;

	ub::matrix< interval<T> > save;

	bool save_mode, save_uh, save_rh;

	ub::vector< autodif< interval<T> > > V, V2;
	autodif< interval<T> > tste;


	new_init = autodif< interval<T> >::compress(init, save);

	m = 1.;
	for (i=0; i<n; i++) {
		m = std::max(m, norm(new_init(i).v));
		new_init(i).d.resize(n);
		for (j=0; j<n; j++) {
			m = std::max(m, norm(new_init(i).d(j)));
		}
	}
	tolerance = m * p.epsilon;

	x = new_init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< autodif< interval<T> > >::mode();
	save_uh = psa< autodif< interval<T> > >::use_history();
	save_rh = psa< autodif< interval<T> > >::record_history();
	psa< autodif< interval<T> > >::mode() = 1;
	psa< autodif< interval<T> > >::use_history() = false;
	psa< autodif< interval<T> > >::record_history() = false;
	#if ODE_FAST == 1
	psa< autodif< interval<T> > >::record_history() = true;
	psa< autodif< interval<T> > >::history().clear();
	#endif
	for (j=0; j<p.order-1; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< autodif< interval<T> > >::use_history() = true;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) {
			y(i) = integrate(y(i));
			y(i) = setorder(y(i), j+1);
		}
		x = new_init + y;
	}

	if (p.autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = p.order-1; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m = std::max(m, norm(x(i).v(j)));
				using std::abs;
				m = std::max(m, abs(mid(x(i).v(j).v)));
				#ifdef IGNORE_DIF_PART
				#else
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					using std::abs;
					m = std::max(m, abs(mid(x(i).v(j).d(k))));
				}
				#endif
			}
			if (m == 0.) continue;
			radius_tmp = std::pow((double)m, 1./j);
			// std::cout << j << " " << m << " " << radius_tmp << "\n";
			if (radius_tmp > radius) radius = radius_tmp;
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
	}
	

	restart = 0;

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

		tste = interval<T>(start.lower(), end2.upper());

		// V = f(new_init, tste) * interval<T>(0., deltat.upper());
		V = f(new_init, tste) * autodif< interval<T> >(interval<T>(0., deltat.upper())); // assist for VC++
		newton_step.resize(n + n * n);
		k = 0;
		for (i=0; i<n; i++) {
			newton_step(k++) = norm(V(i).v);
			km = V(i).d.size();
			for (j=0; j<km; j++) {
				newton_step(k++) = norm(V(i).d(j));
			}
			for (j=km; j<n; j++) newton_step(k++) = 0.;
		}
		make_candidate(newton_step);
		V = new_init;
		k = 0;
		for (i=0; i<n; i++) {
			V(i).v += newton_step(k++) * interval<T>(-1., 1.);
			V(i).d.resize(n);
			for (j=0; j<n; j++) {
				V(i).d(j) += newton_step(k++) * interval<T>(-1., 1.);
			}
		}
		// V2 = new_init + f(V, tste) * interval<T>(0., deltat.upper());
		V2 = new_init + f(V, tste) * autodif< interval<T> >(interval<T>(0., deltat.upper())); // assist for VC++

		flag = true;
		for (i=0; i<n; i++) {
			flag = flag && subset(V2(i).v, V(i).v);
			V2(i).d.resize(n);
			for (j=0; j<n; j++) {
				flag = flag && subset(V2(i).d(j), V(i).d(j));
			}
		}
		if (flag) break;

		if (!p.autostep || restart >= p.restart_max + 10) {
			ret_val = 0;
			break;
		}
		radius *= 0.5;
		restart++;
	}

	if (ret_val != 0) {
		for (j=0; j<p.iteration; j++) {
			V = V2;
			// V2 = new_init + f(V, tste) * interval<T>(0., deltat.upper());
			V2 = new_init + f(V, tste) * autodif< interval<T> >(interval<T>(0., deltat.upper())); // assist for VC++
			for (i=0; i<n; i++) {
				V2(i).v = intersect(V2(i).v, V(i).v);
				V2(i).d.resize(n);
				for (j=0; j<n; j++) {
					V2(i).d(j) = intersect(V2(i).d(j), V(i).d(j));
				}
			}
		}

		z = V2;
		torg.v.resize(2);
		torg.v(0) = tste; torg.v(1) = 1.;

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
			y = f(z, t);
			for (i=0; i<n; i++) {
				y(i) = integrate(y(i));
				y(i) = setorder(y(i), j+1);
			}
			z = new_init + y;
		}

		for (i=0; i<n; i++) {
			for (j=0; j<p.order; j++) {
				z(i).v(j) = x(i).v(j);
			}
		}

		for (i=0; i<n; i++) {
			for (j=0; j<=p.order; j++) {
				z(i).v(j).d.resize(n);
				z(i).v(j) = autodif< interval<T> >::expand(z(i).v(j), save);
			}
		}

		result.resize(n);
		for (i=0; i<n; i++) {
			result(i) = eval(z(i), (autodif< interval<T> >)deltat);
		}

		init = result;
		if (ret_val == 1) end = end2;
		if (result_psa != NULL) *result_psa = z;
	}

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;

	return ret_val;
}


template <class T, class F>
int
odelong_lohner(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>()) {

	ub::vector< interval<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;

		r = ode_lohner(f, x, t, t1, p);
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


template <class T, class F>
int
odelong_lohner(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>()) {

	ub::vector< autodif< interval<T> > > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;

		r = ode_lohner(f, x, t, t1, p);
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

#endif // ODE_LOHNER_HPP
