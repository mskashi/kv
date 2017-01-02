/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_NV_HPP
#define ODE_NV_HPP

// ODE (not verified)

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/psa.hpp>
#include <kv/ode-param.hpp>


#ifndef ODE_FAST
#define ODE_FAST 1
#endif


namespace ub = boost::numeric::ublas;


namespace kv{


template <class T, class F>
void
ode_nv(F f, ub::vector<T>& init, const T& start, T& end, ode_param<T> p = ode_param<T>()) {
	int n = init.size();
	int i, j;

	ub::vector< psa<T> > x, y;
	psa<T> torg;
	psa<T> t;

	T deltat;
	ub::vector<T> result;

	T m, m_tmp;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;

	bool save_mode, save_uh, save_rh;

	m = p.epsilon;
	for (i=0; i<n; i++) {
		using std::abs;
		m_tmp = abs(init(i)) * p.epsilon;
		if (m_tmp > m) m = m_tmp;
	}
	tolerance = m;

	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa<T>::mode();
	save_uh = psa<T>::use_history();
	save_rh = psa<T>::record_history();
	psa<T>::mode() = 1;
	psa<T>::use_history() = false;
	psa<T>::record_history() = false;
	#if ODE_FAST == 1
	psa<T>::record_history() = true;
	psa<T>::history().clear();
	#endif
	for (j=0; j<p.order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa<T>::use_history() = true;
		if (j == p.order - 1) psa<T>::record_history() = false;
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
		for (j = p.order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				m_tmp = (x(i).v(j) >= 0.) ? x(i).v(j) : -x(i).v(j);
				if (m_tmp > m) m = m_tmp;
			}
			if (m == 0.) continue;
			radius_tmp = std::pow((double)m, 1./j);
			if (radius_tmp > radius) radius = radius_tmp;
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
	}

	deltat = end - start;

	if (p.autostep && radius < deltat) {
		end = start + radius;
		deltat = end - start;
	}

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(x(i), deltat);
	}

	init = result;

	psa<T>::mode() = save_mode;
	psa<T>::use_history() = save_uh;
	psa<T>::record_history() = save_rh;
}

template <class T, class F>
void
odelong_nv(F f, ub::vector<T>& init, const T& start, const T& end, ode_param<T> p = ode_param<T>()) {

	ub::vector<T> x;
	T t, t1;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;
		if (t == t1) break;

		ode_nv(f, x, t, t1, p);
		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << x << "\n";
		}
		t = t1;
	}

	init = x;
}

} // namespace kv

#endif // ODE_NV_HPP
