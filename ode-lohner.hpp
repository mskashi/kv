#ifndef ODE_LOHNER_HPP
#define ODE_LOHNER_HPP

// ODE by Lohner's method

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
#include "autodif.hpp"

// autostep==trueで失敗した場合のrestartの回数の上限
#ifndef RESTART_MAX
#define RESTART_MAX 10
#endif

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
int
ode_lohner(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, int order, bool autostep = true, int iter_max = 2, ub::vector< psa< interval<T> > >* result_psa = NULL) {
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

	ub::vector< interval<T> > V, V2;
	interval<T> tste;


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
	for (j=0; j<order-1; j++) {
		#ifdef ODE_FAST
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

	if (autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = order-1; j>=1; j--) {
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

	restart = 0;

	label:

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
		// std::cout << "restart\n";
		goto label;
	}

	for (j=0; j<iter_max; j++) {
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
	#ifdef ODE_FAST
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif

	for (j=0; j<order; j++) {
		#ifdef ODE_FAST
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
		for (j=0; j<order; j++) {
			z(i).v(j) = x(i).v(j);
		}
	}

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(z(i), deltat);
	}

	last:

	if (ret_val != 0) {
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
ode_lohner(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, int order, bool autostep = true, int iter_max = 2, ub::vector< psa< autodif< interval<T> > > >* result_psa = NULL) {
	int n = init.size();
	int i, j, k, km;

	ub::vector< psa< autodif< interval<T> > > > x, y;
	psa< autodif< interval<T> > > torg;
	psa< autodif< interval<T> > > t;

	ub::vector< psa< autodif< interval<T> > > > z, w;

	psa< autodif< interval<T> > > temp;
	T m, m_tmp;
	ub::vector<T> newton_step;

	bool flag;

	interval<T> deltat;
	ub::vector< autodif< interval<T> > > result, new_init;

	T radius, radius_tmp;
	int n_rad;

	int ret_val = 2;
	interval<T> end2 = end;
	int restart;

	ub::matrix< interval<T> > save;

	bool save_mode, save_uh, save_rh;

	ub::vector< autodif< interval<T> > > V, V2;
	autodif< interval<T> > tste;


	new_init = autodif< interval<T> >::compress(init, save);

	x = new_init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< autodif< interval<T> > >::mode();
	save_uh = psa< autodif< interval<T> > >::use_history();
	save_rh = psa< autodif< interval<T> > >::record_history();
	psa< autodif< interval<T> > >::mode() = 1;
	psa< autodif< interval<T> > >::use_history() = false;
	psa< autodif< interval<T> > >::record_history() = false;
	#ifdef ODE_FAST
	psa< autodif< interval<T> > >::record_history() = true;
	psa< autodif< interval<T> > >::history().clear();
	#endif
	for (j=0; j<order-1; j++) {
		#ifdef ODE_FAST
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

	if (autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = order-1; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m_tmp = norm(x(i).v(j));
				m_tmp = mid(x(i).v(j).v);
				if (m_tmp < 0.) m_tmp = -m_tmp;
				if (m_tmp > m) m = m_tmp;
				#ifdef IGNORE_DIF_PART
				#else
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					m_tmp = mid(x(i).v(j).d(k));
					if (m_tmp > m) m = m_tmp;
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

	restart = 0;

	label:

	tste = interval<T>(start.lower(), end2.upper());

	V = f(new_init, tste) * interval<T>(0., deltat.upper());
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
	V2 = new_init + f(V, tste) * interval<T>(0., deltat.upper());

	flag = true;
	for (i=0; i<n; i++) {
		flag = flag && subset(V2(i).v, V(i).v);
		V2(i).d.resize(n);
		for (j=0; j<n; j++) {
			flag = flag && subset(V2(i).d(j), V(i).d(j));
		}
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
		// std::cout << "restart\n";
		goto label;
	}

	for (j=0; j<iter_max; j++) {
		V = V2;
		V2 = new_init + f(V, tste) * interval<T>(0., deltat.upper());
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
	#ifdef ODE_FAST
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif

	for (j=0; j<order; j++) {
		#ifdef ODE_FAST
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
		for (j=0; j<order; j++) {
			z(i).v(j) = x(i).v(j);
		}
	}

	for (i=0; i<n; i++) {
		for (j=0; j<=order; j++) {
			z(i).v(j).d.resize(n);
			z(i).v(j) = autodif< interval<T> >::expand(z(i).v(j), save);
		}
	}

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(z(i), (autodif< interval<T> >)deltat);
	}

	last:

	if (ret_val != 0) {
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
odelong_lohner(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0) {

	ub::vector< interval<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	while (1) {
		t1 = end;

		r = ode_lohner(f, x, t, t1, order, true, iter_max);
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


template <class T, class F>
int
odelong_lohner(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0) {

	ub::vector< autodif< interval<T> > > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	while (1) {
		t1 = end;

		r = ode_lohner(f, x, t, t1, order, true, iter_max);
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

#endif // ODE_LOHNER_HPP
