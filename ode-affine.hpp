#ifndef ODE_AFFINE_HPP
#define ODE_AFFINE_HPP

// ODE using Affine

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
#include "affine.hpp"

#ifndef RESTART_MAX
#define RESTART_MAX 10
#endif

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
int
ode_affine(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, int order, bool autostep = true, int iter_max = 2, int ep_reduce = 0, int ep_reduce_limit = 0)
{
	int n = init.size();
	int i, j;

	ub::vector< psa< affine<T> > > x, y;
	psa< affine<T> > torg;
	psa< affine<T> > t;

	ub::vector< psa< affine<T> > > z, w;

	psa< affine<T> > temp;
	T m, m_tmp;
	ub::vector<T> newton_step;

	bool flag;

	affine<T> s1, s2;
	interval<T> s2i;
	
	ub::vector< affine<T> > s1_save;
	ub::vector< interval<T> > s2i_save;

	interval<T> deltat;
	ub::vector< affine<T> > result;

	T radius, radius_tmp;
	int n_rad;

	int ret_val = 2;
	interval<T> end2 = end;
	int restart;

	int maxnum_save = affine<T>::maxnum();

	bool save_mode, save_uh, save_rh;


	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< affine<T> >::mode();
	save_uh = psa< affine<T> >::use_history();
	save_rh = psa< affine<T> >::record_history();
	psa< affine<T> >::mode() = 1;
	psa< affine<T> >::use_history() = false;
	psa< affine<T> >::record_history() = false;
	#ifdef ODE_FAST
	psa< affine<T> >::record_history() = true;
	psa< affine<T> >::history().clear();
	#endif
	for (j=0; j<order; j++) {
		#ifdef ODE_FAST
		if (j == 1) psa< affine<T> >::use_history() = true;
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
				// m_tmp = norm(to_interval(x(i).v(j)));
				m_tmp = x(i).v(j).get_mid();
				if (m_tmp < 0.) m_tmp = -m_tmp;
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

	psa< affine<T> >::mode() = 2;

	restart = 0;

	label:

	psa< affine<T> >::domain() = interval<T>(0., deltat.upper());

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
		newton_step(i) = norm(to_interval(w(i).v(order) - z(i).v(order)));
	}
	make_candidate(newton_step);

	s1_save.resize(n);
	s2i_save.resize(n);

	for (i=0; i<n; i++) {
		z(i).v(order) += (affine<T>)(newton_step(i) * interval<T>(-1., 1.));
#ifdef ODE_AFFINE_SIMPLE
		s2i = to_interval(z(i).v(order));
		s2i_save(i) = s2i;
		z(i).v(order) = (affine<T>)s2i;
#else 
		split(z(i).v(order), maxnum_save, s1, s2);
		s2i = to_interval(s2);
		s1_save(i) = s1;
		s2i_save(i) = s2i;
		z(i).v(order) = append(s1, (affine<T>)s2i);
#endif
	}

	w = f(z, t);
	for (i=0; i<n; i++) {
		temp = integrate(w(i));
		w(i) = setorder(temp, order);
	}
	w = init + w;

	flag = true;
	for (i=0; i<n; i++) {
#ifdef ODE_AFFINE_SIMPLE
		s2i = to_interval(w(i).v(order));
		flag = flag && subset(s2i, s2i_save(i));
		s2i_save(i) = s2i;
		w(i).v(order) = (affine<T>)s2i;
#else
		split(w(i).v(order), maxnum_save, s1, s2);
		s2i = to_interval(s2);
		flag = flag && subset(to_interval(s1 - s1_save(i)) + s2i, s2i_save(i));
		s1_save(i) = s1;
		s2i_save(i) = s2i;
		w(i).v(order) = append(s1, (affine<T>)s2i);
#endif
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
		goto label;
	}

	for (j=0; j<iter_max; j++) {
		w = f(w, t);
		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, order);
		}
		w = init + w;

		for (i=0; i<n; i++) {
#ifdef ODE_AFFINE_SIMPLE
			s2i = to_interval(w(i).v(order));
			s2i = intersect(s2i, s2i_save(i));
			s2i_save(i) = s2i;
			w(i).v(order) = (affine<T>)s2i;
#else
			split(w(i).v(order), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			// s2i = intersect(to_interval(s1 - s1_save(i)) + s2i, s2i_save(i));
			s2i = intersect(to_interval(s1 - s1_save(i)) + s2i_save(i), s2i);
			s1_save(i) = s1;
			s2i_save(i) = s2i;
			w(i).v(order) = append(s1, (affine<T>)s2i);
#endif
		}
	}

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(w(i), (affine<T>)deltat);
		if (ep_reduce == 0) {
			split(result(i), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			s1_save(i) = s1;
			s2i_save(i) = s2i;
		}
	}

	if (ep_reduce == 0) {
		affine<T>::maxnum() = maxnum_save;
		for (i=0; i<n; i++) {
			s1_save(i).resize();
			result(i) = append(s1_save(i), (affine<T>)s2i_save(i));
		}
	} else {
		epsilon_reduce(result, ep_reduce, ep_reduce_limit);
	}

	last:

	if (ret_val != 0) {
		init = result;
		if (ret_val == 1) end = end2;
	}

	psa< affine<T> >::mode() = save_mode;
	psa< affine<T> >::use_history() = save_uh;
	psa< affine<T> >::record_history() = save_rh;

	return ret_val;
}

template <class T, class F>
int
odelong_affine(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0, int ep_reduce = 0, int ep_reduce_limit = 0)
{
	int s = init.size();
	ub::vector< affine<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	while (1) {
		t1 = end;

		r = ode_affine(f, x, t, t1, order, true, iter_max, ep_reduce, ep_reduce_limit);
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
			std::cout << to_interval(x) << "\n";
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
odelong_affine(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, int order, int iter_max = 2, int verbose = 0, int ep_reduce = 0, int ep_reduce_limit = 0)
{
	int s = init.size();
	int i;
	ub::vector< affine<T> > x;
	int r;
	interval<T> end2 = end;

	affine<T>::maxnum() = 0;
	x = init;

	r = odelong_affine(f, x, start, end2, order, iter_max, verbose, ep_reduce, ep_reduce_limit);
	if (r == 0) return 0;

	for (i=0; i<s; i++) init(i) = to_interval(x(i));
	if (r == 1) end = end2;

	return r;
}

} // namespace kv

#endif // ODE_AFFINE_HPP
