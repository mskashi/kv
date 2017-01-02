/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_AUTODIF_NV_HPP
#define ODE_AUTODIF_NV_HPP

// ODE (input and output : autodif type) (not verified)

#include <iostream>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/psa.hpp>
#include <kv/autodif.hpp>
#include <kv/ode-param.hpp>


#ifndef ODE_FAST
#define ODE_FAST 1
#endif


namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
void
ode_nv(F f, ub::vector< autodif<T> >& init, const T& start, T& end, ode_param<T> p = ode_param<T>()) {
	int n = init.size();
	int i, j, k, km;

	ub::vector< psa< autodif<T> > > x, y;
	psa< autodif<T> > torg;
	psa< autodif<T> > t;

	T deltat;
	ub::vector< autodif<T> > result, new_init;

	T m, m_tmp;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;

	ub::matrix<T> save;

	bool save_mode, save_uh, save_rh;

	new_init = autodif<T>::compress(init, save);

	m = p.epsilon;
	for (i=0; i<n; i++) {
		using std::abs;
		m_tmp = abs(new_init(i).v) * p.epsilon;
		if (m_tmp > m) m = m_tmp;
		new_init(i).d.resize(n);
		for (j=0; j<n; j++) {
			m_tmp = abs(new_init(i).d(j)) * p.epsilon;
			if (m_tmp > m) m = m_tmp;
		}
	}
	tolerance = m;

	x = new_init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< autodif<T> >::mode();
	save_uh = psa< autodif<T> >::use_history();
	save_rh = psa< autodif<T> >::record_history();
	psa< autodif<T> >::mode() = 1;
	psa< autodif<T> >::use_history() = false;
	psa< autodif<T> >::record_history() = false;
	#if ODE_FAST == 1
	psa< autodif<T> >::record_history() = true;
	psa< autodif<T> >::history().clear();
	#endif
	for (j=0; j<p.order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< autodif<T> >::use_history() = true;
		if (j == p.order - 1) psa< autodif<T> >::record_history() = false;
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
		for (j = p.order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				m_tmp = (x(i).v(j).v >= 0.) ? x(i).v(j).v : -x(i).v(j).v;
				if (m_tmp > m) m = m_tmp;
				// 微分項は考慮しない手もある。
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					m_tmp = (x(i).v(j).d(k) >= 0) ? x(i).v(j).d(k) : -x(i).v(j).d(k);
					if (m_tmp > m) m = m_tmp;
				}
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
		result(i) = eval(x(i), (autodif<T>)deltat);
		result(i).d.resize(n);
		result(i) = autodif<T>::expand(result(i), save);
	}

	init = result;

	psa< autodif<T> >::mode() = save_mode;
	psa< autodif<T> >::use_history() = save_uh;
	psa< autodif<T> >::record_history() = save_rh;
}


template <class T, class F>
void
odelong_nv(F f, ub::vector< autodif<T> >& init, const T& start, const T& end, ode_param<T> p = ode_param<T>()) {

	ub::vector< autodif<T> > x;
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

#endif // ODE_AUTODIF_NV_HPP
