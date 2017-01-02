#ifndef ODE_NV_HPP
#define ODE_NV_HPP

// ODE (not verified)

#include <iostream>
#include <cmath>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "psa.hpp"

namespace ub = boost::numeric::ublas;


namespace kv{


template <class T, class F>
void
ode_nv(F f, ub::vector<T>& init, const T& start, T& end, int order, bool autostep = true) {
	int n = init.size();
	int i, j;

	ub::vector< psa<T> > x, y;
	psa<T> torg;
	psa<T> t;

	T deltat;
	ub::vector<T> result;

	T m, m_tmp;

	T radius, radius_tmp;
	int n_rad;

	bool save_mode, save_uh, save_rh;


	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa<T>::mode();
	save_uh = psa<T>::use_history();
	save_rh = psa<T>::record_history();
	psa<T>::mode() = 1;
	psa<T>::use_history() = false;
	psa<T>::record_history() = false;
	#ifdef ODE_FAST
	psa<T>::record_history() = true;
	psa<T>::history().clear();
	#endif
	for (j=0; j<order; j++) {
		#ifdef ODE_FAST
		if (j == 1) psa<T>::use_history() = true;
		if (j == order - 1) psa<T>::record_history() = false;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) {
			y(i) = integrate(y(i));
			y(i) = setorder(y(i), j+1);
		}
		x = init + y;
	}

	if (autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = order; j>=1; j--) {
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
#ifdef FORCE_STEP
		radius = FORCE_STEP;
#else
#ifdef EXP_STEP
		radius = std::exp(-EXP_STEP) / radius;
#else
		radius = std::pow((double)std::numeric_limits<T>::epsilon()*order/(order-1), 1./order) / radius;
#endif
#endif
	}

	deltat = end - start;

	if (autostep && radius < deltat) {
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
odelong_nv(F f, ub::vector<T>& init, const T& start, const T& end, int order, int verbose = 0) {

	ub::vector<T> x;
	T t, t1;

	x = init;
	t = start;
	while (1) {
		t1 = end;
		if (t == t1) break;

		ode_nv(f, x, t, t1, order, true);
		if (verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << x << "\n";
		}
		t = t1;
	}

	init = x;
}

} // namespace kv

#endif // ODE_NV_HPP
