#ifndef ODE_HPP
#define ODE_HPP

// ODE

#include <iostream>
#include <cmath>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "interval.hpp"
#include "rdouble.hpp"
#include "interval-vector.hpp"
#include "make-candidate.hpp"
#include "psa.hpp"

// autostep==trueで失敗した場合のrestartの回数の上限
#ifndef RESTART_MAX
#define RESTART_MAX 10
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

	bool flag;

	interval<T> deltat;
	ub::vector< interval<T> > result;

	T radius, radius_tmp;
	int n_rad;

	int ret_val = 2;
	interval<T> end2 = end;
	int restart;

	bool save_mode, save_uh, save_rh;


	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< interval<T> >::mode();
	save_uh = psa< interval<T> >::use_history();
	save_rh = psa< interval<T> >::record_history();
	psa< interval<T> >::mode() = 1;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;
	#ifdef ODE_FAST
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	for (j=0; j<order; j++) {
		#ifdef ODE_FAST
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

#ifdef FORCE_STEP
		radius = FORCE_STEP;
#else
#ifdef EXP_STEP
		radius = std::exp(-EXP_STEP) / radius;
#else
		radius = std::pow((double)std::numeric_limits<T>::epsilon()*order/(order-1), 1./order) / radius;
#endif
#endif
		// std::cout << "radius: " << radius << "\n";
	}

	deltat = end2 - start;

	if (autostep && radius < deltat.lower()) {
		// deltatを幅無し区間にする方針
		// end2 = start + radius;
		// deltat = radius;

		// end2を幅無し区間にする方針
		end2 = mid(start + radius);
		deltat = end2 - start;
		if (deltat.lower() <= 0.) {
			ret_val = 0;
			goto last;
		}
		ret_val = 1;
	}

	psa< interval<T> >::mode() = 2;

	restart = 0;

	label:

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

	if (!flag) {
		if (!autostep || restart >= RESTART_MAX) {
			ret_val = 0;
			goto last;
		}
		end2 = mid((start + end2) * 0.5);
		deltat = end2 - start;
		ret_val = 1;
		if (deltat.lower() <= 0.) {
			ret_val = 0;
			goto last;
		}
		restart++;
		goto label;
	}

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

	last:

	if (ret_val != 0) {
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
