/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_AUTODIF_HPP
#define ODE_AUTODIF_HPP

// ODE (input and output : autodif type)

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
ode(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< autodif< interval<T> > > >* result_psa = NULL) {
	int n = init.size();
	int i, j, k, km;

	ub::vector< psa< autodif< interval<T> > > > x, y;
	psa< autodif< interval<T> > > torg;
	psa< autodif< interval<T> > > t;

	ub::vector< psa< autodif< interval<T> > > > z, w;
	autodif< interval<T> > wmz;
	autodif< interval<T> > evalz;

	psa< autodif< interval<T> > > temp;
	T m;
	ub::vector<T> newton_step;

	bool flag, resized;

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
	for (j=0; j<p.order; j++) {
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
		for (j = p.order; j>=1; j--) {
			m = 0.;
			for (i=0; i<n; i++) {
				// m = std::max(m, norm(x(i).v(j).v));
				using std::abs;
				m = std::max(m, abs(mid(x(i).v(j).v)));
				#ifdef IGNORE_DIF_PART
				#else
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					// m = std::max(m, norm(x(i).v(j).d(k)));
					using std::abs;
					m = std::max(m, abs(mid(x(i).v(j).d(k))));
				}
				#endif
			}
			if (m == 0.) continue;
			radius_tmp = std::pow((double)m, 1./j);
			if (radius_tmp > radius) radius = radius_tmp;
			n_rad++;
			if (n_rad == 2) break;
		}
		radius = std::pow((double)tolerance, 1./p.order) / radius;
	}

	psa< autodif< interval<T> > >::mode() = 2;

	restart = 0;
	resized = false;

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

		psa< autodif< interval<T> > >::domain() = interval<T>(0., deltat.upper());

		z = x;
		t = setorder(torg, p.order);

		try {
			w = f(z, t);
		}
		catch (std::domain_error& e) {
			if (restart < p.restart_max) {
				psa< autodif< interval<T> > >::use_history() = false;
				radius *= 0.5;
				restart++;
				continue;
                        } else {
				throw std::domain_error("ode: evaluation error");
                        }

		}

		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, p.order);
		}
		w = new_init + w;

		newton_step.resize(n + n * n);
		k = 0;
		for (i=0; i<n; i++) {
			wmz = w(i).v(p.order) - z(i).v(p.order);
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
			z(i).v(p.order).v += newton_step(k++) * interval<T>(-1., 1.);
			z(i).v(p.order).d.resize(n);
			for (j=0; j<n; j++) {
				z(i).v(p.order).d(j) += newton_step(k++) * interval<T>(-1., 1.);
			}
		}

		if (p.autostep && resized == false) {
			resized = true;
			m = 0.;
			for (i=0; i<n; i++) {
				evalz = eval(z(i), autodif< interval<T> >(deltat));
				m = std::max(m, rad(evalz.v) - rad(new_init(i).v));
				evalz.d.resize(n);
				for (j=0; j<n; j++) {
					m = std::max(m, rad(evalz.d(j)) - rad(new_init(i).d(j)));
				}
			}
			m = m / tolerance;
			if (restart > 0) {
				radius /= std::max(1., std::pow((double)m, 1. / p.order));
			} else {
				radius /= std::pow((double)m, 1. / p.order);
			}
			continue;
		}

		w = f(z, t);
		for (i=0; i<n; i++) {
			temp = integrate(w(i));
			w(i) = setorder(temp, p.order);
		}
		w = new_init + w;

		flag = true;
		for (i=0; i<n; i++) {
			flag = flag && subset(w(i).v(p.order).v, z(i).v(p.order).v);
			w(i).v(p.order).d.resize(n);
			for (j=0; j<n; j++) {
				flag = flag && subset(w(i).v(p.order).d(j), z(i).v(p.order).d(j));
			}
		}
		if (flag) break;

		if (!p.autostep || restart >= p.restart_max) {
			ret_val = 0;
			break;
		}
		radius *= 0.5;
		restart++;
	}

	if (ret_val != 0) {
		for (k=0; k<p.iteration; k++) {
			z = w;
			w = f(z, t);
			for (i=0; i<n; i++) {
				temp = integrate(w(i));
				w(i) = setorder(temp, p.order);
			}
			w = new_init + w;
			for (i=0; i<n; i++) {
				w(i).v(p.order).v = intersect(w(i).v(p.order).v, z(i).v(p.order).v);
				w(i).v(p.order).d.resize(n);
				for (j=0; j<n; j++) {
					w(i).v(p.order).d(j) = intersect(w(i).v(p.order).d(j), z(i).v(p.order).d(j));
				}
			}
		}

		for (i=0; i<n; i++) {
			for (j=0; j<=p.order; j++) {
				w(i).v(j).d.resize(n);
				w(i).v(j) = autodif< interval<T> >::expand(w(i).v(j), save);
			}
		}

		result.resize(n);
		for (i=0; i<n; i++) {
			result(i) = eval(w(i), (autodif< interval<T> >)deltat);
		}

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
odelong(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>()) {

	ub::vector< autodif < interval<T> > > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;

		r = ode(f, x, t, t1, p);
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

#endif // ODE_AUTODIF_HPP
