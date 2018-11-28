/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_AFFINE_HPP
#define ODE_AFFINE_HPP

// ODE using Affine

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
#include <kv/affine.hpp>
#include <kv/ode-param.hpp>
#include <kv/ode-callback.hpp>


#ifndef ODE_FAST
#define ODE_FAST 1
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
ode_affine(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL)
{
	int n = init.size();
	int i, j;

	ub::vector< psa< affine<T> > > x, y;
	psa< affine<T> > torg;
	psa< affine<T> > t;

	ub::vector< psa< affine<T> > > z, w;

	psa< affine<T> > temp;
	T m;
	ub::vector<T> newton_step;

	bool flag, resized;

	affine<T> s1, s2;
	interval<T> s2i;
	
	ub::vector< affine<T> > s1_save;
	ub::vector< interval<T> > s2i_save;

	interval<T> deltat;
	ub::vector< affine<T> > result;

	T radius, radius_tmp;
	T tolerance;
	int n_rad;
	#if ODE_RESTART_RATIO == 1
	T max_ratio;
	#endif

	int ret_val;
	interval<T> end2;
	int restart;

	int maxnum_save = affine<T>::maxnum();

	bool save_mode, save_uh, save_rh;


	m = 1.;
	for (i=0; i<n; i++) {
		m = std::max(m, norm(to_interval(init(i))));
	}
	tolerance = m * p.epsilon;

	x = init;
	torg.v.resize(2);
	torg.v(0) = start; torg.v(1) = 1.;

	save_mode = psa< affine<T> >::mode();
	save_uh = psa< affine<T> >::use_history();
	save_rh = psa< affine<T> >::record_history();
	psa< affine<T> >::mode() = 1;
	psa< affine<T> >::use_history() = false;
	psa< affine<T> >::record_history() = false;
	#if ODE_FAST == 1
	psa< affine<T> >::record_history() = true;
	psa< affine<T> >::history().clear();
	#endif
	for (j=0; j<p.order; j++) {
		#if ODE_FAST == 1
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

	if (p.autostep) {
		radius = 0.;
		n_rad = 0;
		for (j = p.order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				#if ODE_COEF_MID == 1
				using std::abs;
				m = std::max(m, abs(x(i).v(j).get_mid()));
				#else
				m = std::max(m, norm(to_interval(x(i).v(j))));
				#endif
			}
			if (m == 0.) continue;
			radius = std::max(radius, (T)std::pow((double)m, 1./j));
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
	}

	psa< affine<T> >::mode() = 2;

	restart = 0;
	// disable resize because resize algorithm is not suitable
	//  for ode-affine
	// resized = false;
	resized = true;

	while (true) {
		if (p.autostep) {
			end2 = mid(start + radius);
			if (end2 >= end.lower()) {
				end2 = end;
				radius = mid(end2 - start);
				ret_val = 2;
			} else {
				ret_val = 1;
			}
		} else {
			end2 = end;
			ret_val = 2;
		}
		deltat = end2 - start;

		psa< affine<T> >::domain() = interval<T>(0., deltat.upper());

		z = x;
		t = setorder(torg, p.order);

		try { 
			w = f(z, t);
		}
		catch (std::domain_error& e) {
			 if (p.autostep && restart < p.restart_max) {
				psa< affine<T> >::use_history() = false;
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
				throw std::domain_error("ode_affine: evaluation error");
			}
		}

		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, p.order);
		}
		w = init + w;

		newton_step.resize(n);
		for (i=0; i<n; i++) {
			newton_step(i) = norm(to_interval(w(i).v(p.order) - z(i).v(p.order)));
		}
		make_candidate(newton_step);

		s1_save.resize(n);
		s2i_save.resize(n);

		for (i=0; i<n; i++) {
			z(i).v(p.order) += (affine<T>)(newton_step(i) * interval<T>(-1., 1.));
			#ifdef ODE_AFFINE_SIMPLE
			s2i = to_interval(z(i).v(p.order));
			s2i_save(i) = s2i;
			z(i).v(p.order) = (affine<T>)s2i;
			#else 
			split(z(i).v(p.order), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			s1_save(i) = s1;
			s2i_save(i) = s2i;
			z(i).v(p.order) = append(s1, (affine<T>)s2i);
			#endif
		}

		if (p.autostep && ret_val != 2 && resized == false) {
			resized = true;
			m = (std::numeric_limits<T>::min)();
			for (i=0; i<n; i++) {
				m = std::max(m, rad(eval(z(i), (affine<T>)deltat)) - rad(init(i)));
			}
			m = m / tolerance;
			radius_tmp = radius / std::pow((double)m, 1. / p.order);
			if (radius_tmp >= radius && restart > 0) {
				// do nothing, not continue
			} else {
				radius = radius_tmp;
				continue;
			}
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
			#ifdef ODE_AFFINE_SIMPLE
			s2i = to_interval(w(i).v(p.order));
			#if ODE_RESTART_RATIO == 1
			max_ratio = std::max(max_ratio, width(s2i) / width(s2i_save(i)));
			#endif
			flag = flag && subset(s2i, s2i_save(i));
			s2i_save(i) = s2i;
			w(i).v(p.order) = (affine<T>)s2i;
			#else
			split(w(i).v(p.order), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			#if ODE_RESTART_RATIO == 1
			max_ratio = std::max(max_ratio, width(to_interval(s1 - s1_save(i)) + s2i) / width(s2i_save(i)));
			#endif
			flag = flag && subset(to_interval(s1 - s1_save(i)) + s2i, s2i_save(i));
			s1_save(i) = s1;
			s2i_save(i) = s2i;
			w(i).v(p.order) = append(s1, (affine<T>)s2i);
			#endif
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
			w = f(w, t);
			for (i=0; i<n; i++) {
				temp = integrate(w(i));
				w(i) = setorder(temp, p.order);
			}
			w = init + w;

			for (i=0; i<n; i++) {
				#ifdef ODE_AFFINE_SIMPLE
				s2i = to_interval(w(i).v(p.order));
				s2i = intersect(s2i, s2i_save(i));
				s2i_save(i) = s2i;
				w(i).v(p.order) = (affine<T>)s2i;
				#else
				split(w(i).v(p.order), maxnum_save, s1, s2);
				s2i = to_interval(s2);
				// s2i = intersect(to_interval(s1 - s1_save(i)) + s2i, s2i_save(i));
				s2i = intersect(to_interval(s1 - s1_save(i)) + s2i_save(i), s2i);
				s1_save(i) = s1;
				s2i_save(i) = s2i;
				w(i).v(p.order) = append(s1, (affine<T>)s2i);
				#endif
			}
		}

		if (result_psa != NULL) {
			// convert w to interval and store it to *result_psa
			(*result_psa).resize(n);
			for (i=0; i<n; i++) {
				(*result_psa)(i).v.resize(w(i).v.size());
				for (j=0; j<w(i).v.size(); j++) {
					(*result_psa)(i).v(j) = to_interval(w(i).v(j));
				}
			}
		}

		result.resize(n);
		for (i=0; i<n; i++) {
			result(i) = eval(w(i), (affine<T>)deltat);
			if (p.ep_reduce == 0) {
				split(result(i), maxnum_save, s1, s2);
				s2i = to_interval(s2);
				s1_save(i) = s1;
				s2i_save(i) = s2i;
			}
		}

		if (p.ep_reduce == 0) {
			affine<T>::maxnum() = maxnum_save;
			for (i=0; i<n; i++) {
				s1_save(i).resize();
				result(i) = append(s1_save(i), (affine<T>)s2i_save(i));
			}
		} else {
			epsilon_reduce(result, p.ep_reduce, p.ep_reduce_limit);
		}

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
odelong_affine(
	F f,
	ub::vector< affine<T> >& init,
	const interval<T>& start,
	interval<T>& end,
	ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>()
) {
	int s = init.size();
	ub::vector< affine<T> > x, x1;
	interval<T> t, t1;
	int ret_ode;
	int ret_val = 0;
	bool ret_callback;

	ub::vector< psa< interval<T> > > result_tmp;


	x = init;
	t = start;
	p.set_autostep(true);

	while (1) {
		x1 = x;
		t1 = end;

		ret_ode = ode_affine(f, x1, t, t1, p, &result_tmp);
		if (ret_ode == 0) {
			if (ret_val == 1) {
				init = x1;
				end = t;
			}
			return ret_val;
		}
		ret_val = 1;
		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << to_interval(x1) << "\n";
		}

		callback(t, t1, to_interval(x), to_interval(x1), result_tmp);
		
		if (ret_callback == false) {
			init = x1;
			end = t1;
			return 3;
		}

		if (ret_ode == 2) {
			init = x1;
			return 2;
		}

		t = t1;
		x = x1;
	}
}

template <class T, class F>
int
odelong_affine(
	F f,
	ub::vector< interval<T> >& init,
	const interval<T>& start,
	interval<T>& end,
	ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>()
) {
	int s = init.size();
	int i;
	ub::vector< affine<T> > x;
	int maxnum_save;
	int r;

	maxnum_save = affine<T>::maxnum();
	affine<T>::maxnum() = 0;

	x = init;

	r = odelong_affine(f, x, start, end, p, callback);

	affine<T>::maxnum() = maxnum_save;

	if (r == 0) return 0;

	for (i=0; i<s; i++) init(i) = to_interval(x(i));

	return r;
}

} // namespace kv

#endif // ODE_AFFINE_HPP
