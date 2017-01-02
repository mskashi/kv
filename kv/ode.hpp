/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_HPP
#define ODE_HPP

// ODE

#include <iostream>
#include <cmath>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/make-candidate.hpp>
#include <kv/psa.hpp>

#ifndef ODE_FAST
#define ODE_FAST 1
#endif

// autostep==trueで失敗した場合のrestartの回数の上限
#ifndef RESTART_MAX
#define RESTART_MAX 1
#endif

#ifndef TOL1
#define TOL1 0.1
#endif

#if 0
#ifndef TOL2
#define TOL2 100.
#endif

#ifndef TOL3
#define TOL3 0.01
#endif
#endif


namespace ub = boost::numeric::ublas;


namespace kv {

template <class T, class F>
int
ode(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, int order, bool autostep = true, int iter_max = 2, ub::vector< psa< interval<T> > >* result_psa = NULL) {
	int n = init.size();
	int i, j;

	ub::vector< psa< interval<T> > > x, y;
	psa< interval<T> > torg;
	psa< interval<T> > t;

	ub::vector< psa< interval<T> > > z, w;

	psa< interval<T> > temp;
	T m, m_tmp;
	ub::vector<T> newton_step;

	bool flag, resized;

	interval<T> deltat;
	ub::vector< interval<T> > result;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;

	int ret_val;
	interval<T> end2;
	int restart;

	bool save_mode, save_uh, save_rh;

	m = std::numeric_limits<T>::epsilon();
	for (i=0; i<n; i++) {
		m_tmp = norm(init(i)) * std::numeric_limits<T>::epsilon();
		if (m_tmp > m) m = m_tmp;
		m_tmp = rad(init(i)) * TOL1;
		if (m_tmp > m) m = m_tmp;
	}
	tolerance = m;

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
	for (j=0; j<order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< interval<T> >::use_history() = true;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) {
			y(i) = integrate(y(i));
			// 定数関数などの特殊な状況に備えてorderを設定しておく
			y(i) = setorder(y(i), j+1);
		}
		x = init + y;
	}

	if (autostep) {
		// 高次の項から順に見て、最初に見付かった2つの非0項を使う。
		radius = 0.;
		n_rad = 0;
		for (j = order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m_tmp = norm(x(i).v(j));
				m_tmp = mid(x(i).v(j));
				if (m_tmp < 0.) m_tmp = -m_tmp;
				if (m_tmp > m) m = m_tmp;
			}
			if (m == 0.) continue;
			radius_tmp = std::pow((double)m, 1./j);
			// std::cout << j << " " << m << " " << radius_tmp << "\n";
			if (radius_tmp > radius) radius = radius_tmp;
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./order) / radius;
	}

	psa< interval<T> >::mode() = 2;

	restart = 0;
	resized = false;

	while (true) {
		if (autostep) {
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
		t = setorder(torg, order);

		w = f(z, t);

		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, order);
		}
		w = init + w;

		newton_step.resize(n);
		for (i=0; i<n; i++) {
			newton_step(i) = norm(w(i).v(order) - z(i).v(order));
		}
		make_candidate(newton_step);
		for (i=0; i<n; i++) {
			z(i).v(order) += newton_step(i) * interval<T>(-1., 1.);
		}

		if (autostep && resized == false) {
			resized = true;
			m = 0.;
			for (i=0; i<n; i++) {
				m_tmp = rad(eval(z(i), deltat)) - rad(init(i));
				if (m_tmp > m) m = m_tmp;
			}
			m = m / tolerance;
			#if 0
			if (m > TOL2) m = TOL2;
			if (m < TOL3) m = TOL3;
			#endif
			radius /= std::pow((double)m, 1. / order);
			continue;
		}

		w = f(z, t);
		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, order);
		}
		w = init + w;

		flag = true;
		for (i=0; i<n; i++) {
			flag = flag && subset(w(i).v(order), z(i).v(order));
		}
		if (flag) break;

		if (!autostep || restart >= RESTART_MAX) {
			ret_val = 0;
			break;
		}
		radius *= 0.5;
		restart++;
	}

	if (ret_val != 0) {
		for (j=0; j<iter_max; j++) {
			z = w;
			w = f(z, t);
			for (i=0; i<n; i++) {
				temp = integrate(w(i));
				w(i) = setorder(temp, order);
			}
			w = init + w;
			for (i=0; i<n; i++) {
				w(i).v(order) = intersect(w(i).v(order), z(i).v(order));
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
odelong(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0) {

	ub::vector< interval<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	while (1) {
		t1 = end;

		r = ode(f, x, t, t1, order, true, iter_max);
		if (r == 0) {
			if (ret_val == 1) {
				init = x;
				end = t;
			}
			return ret_val;
		}
		ret_val = 1;
		if (verbose == 1) {
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
