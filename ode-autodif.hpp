#ifndef ODE_AUTODIF_HPP
#define ODE_AUTODIF_HPP

// ODE (input and output : autodif type)

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

#ifndef RESTART_MAX
#define RESTART_MAX 10
#endif

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
int
ode(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, int order, bool autostep = true, int iter_max = 2, ub::vector< psa< autodif< interval<T> > > >* result_psa = NULL) {
	int n = init.size();
	int i, j, k, km;

	ub::vector< psa< autodif< interval<T> > > > x, y;
	psa< autodif< interval<T> > > torg;
	psa< autodif< interval<T> > > t;

	ub::vector< psa< autodif< interval<T> > > > z, w;
	autodif< interval<T> > wmz;

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
	for (j=0; j<order; j++) {
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
		for (j = order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m_tmp = norm(x(i).v(j).v);
				m_tmp = mid(x(i).v(j).v);
				if (m_tmp < 0.) m_tmp = -m_tmp;
				if (m_tmp > m) m = m_tmp;
				// 微分項は考慮しない手もある。
				#ifdef IGNORE_DIF_PART
				#else
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					// m_tmp = norm(x(i).v(j).d(k));
					m_tmp = mid(x(i).v(j).d(k));
					if (m_tmp > m) m = m_tmp;
				}
				#endif
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

	deltat = end2 - start;

	if (autostep && radius < deltat.lower()) {
		end2 = mid(start + radius);
		deltat = end2 - start;
		if (deltat.lower() <= 0.) {
			ret_val = 0;
			goto last;
		}
		ret_val = 1;
	}

	psa< autodif< interval<T> > >::mode() = 2;

	restart = 0;

	label:

	psa< autodif< interval<T> > >::domain() = interval<T>(0., deltat.upper());

	z = x;
	t = setorder(torg, order);

	w = f(z, t);

	for (i=0; i<n; i++) {
		temp = integrate(w(i));
		w(i) = setorder(temp, order);
	}
	w = new_init + w;

	newton_step.resize(n + n * n);
	k = 0;
	for (i=0; i<n; i++) {
		wmz = w(i).v(order) - z(i).v(order);
		newton_step(k++) = norm(wmz.v);
		km = wmz.d.size();
		for (j=0; j<km; j++) {
			newton_step(k++) = norm(wmz.d(j));
		}
		for (j=km; j<n; j++) newton_step(k++) = 0.;
	}
	make_candidate(newton_step);
	k = 0;
	for (i=0; i<n; i++) {
		z(i).v(order).v += newton_step(k++) * interval<T>(-1., 1.);
		z(i).v(order).d.resize(n);
		for (j=0; j<n; j++) {
			z(i).v(order).d(j) += newton_step(k++) * interval<T>(-1., 1.);
		}
	}

	w = f(z, t);
	for (i=0; i<n; i++) {
		temp = integrate(w(i));
		w(i) = setorder(temp, order);
	}
	w = new_init + w;

	flag = true;
	for (i=0; i<n; i++) {
		flag = flag && subset(w(i).v(order).v, z(i).v(order).v);
		w(i).v(order).d.resize(n);
		for (j=0; j<n; j++) {
			flag = flag && subset(w(i).v(order).d(j), z(i).v(order).d(j));
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
		goto label;
	}

	for (k=0; k<iter_max; k++) {
		z = w;
		w = f(z, t);
		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, order);
		}
		w = new_init + w;
		for (i=0; i<n; i++) {
			w(i).v(order).v = intersect(w(i).v(order).v, z(i).v(order).v);
			w(i).v(order).d.resize(n);
			for (j=0; j<n; j++) {
				w(i).v(order).d(j) = intersect(w(i).v(order).d(j), z(i).v(order).d(j));
			}
		}
	}

	for (i=0; i<n; i++) {
		for (j=0; j<=order; j++) {
			w(i).v(j).d.resize(n);
			w(i).v(j) = autodif< interval<T> >::expand(w(i).v(j), save);
		}
	}

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(w(i), (autodif< interval<T> >)deltat);
	}

	last:

	if (ret_val != 0) {
		init = result;
		if (ret_val == 1) end = end2;
		if (result_psa != NULL) *result_psa = w;
	}
	psa< autodif< interval<T> > >::mode() = save_mode;
	psa< autodif< interval<T> > >::use_history() = save_uh;
	psa< autodif< interval<T> > >::record_history() = save_rh;

	return ret_val;
}


template <class T, class F>
int
odelong(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0) {

	ub::vector< autodif < interval<T> > > x;
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

#endif // ODE_AUTODIF_HPP
