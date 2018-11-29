/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_AUTODIF_HPP
#define ODE_AUTODIF_HPP

// ODE (input and output : autodif type)

#include <kv/ode.hpp>
#include <kv/autodif.hpp>


#ifndef ODE_AUTODIF_NEW
#define ODE_AUTODIF_NEW 1
#endif


namespace kv {

namespace ub = boost::numeric::ublas;


#if ODE_AUTODIF_NEW == 1

template <class F, class T> struct MakeVariationalEq {
	F f;
	ub::vector< psa<T> > solution;
	int s, s2;

	MakeVariationalEq(F f, ub::vector< psa<T> > solution) : f(f), solution(solution) {
		s = solution.size();
		s2 = s * s;
	}

	ub::vector< psa<T> > operator() (const ub::vector< psa<T> >& x, const psa<T>& t){
		ub::matrix< psa<T> > x2(s, s);
		ub::vector< psa<T> > y(s2);

		ub::vector< psa<T> > solution2(s);
		psa<T> t2;

		ub::vector< psa<T> > rv;
		ub::matrix< psa<T> > rm;

		int i, j, k;
		int order, tmp;

		order = 0;
		k = 0;
		for (i=0; i<s; i++) {
			for (j=0; j<s; j++) {
				tmp = x(k).v.size() - 1;
				if (tmp > order) order = tmp;
				x2(i, j) = x(k);
				k++;
			}
		}

		for (i=0; i<s; i++) {
			solution2(i) = setorder(solution(i), order);
		}
		t2 = setorder(t, order);

		autodif< psa<T> >::split(f(autodif< psa<T> >::init(solution2), autodif< psa<T> >(t2)), rv, rm);

		rm = prod(rm, x2);

		k = 0;
		for (i=0; i<s; i++) {
			for (j=0; j<s; j++) {
				y(k++) = rm(i, j);
			}
		}

		return y;
	}
};

template <class T, class F>
int
ode(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL) {
	int n = init.size();
	int i, j, k;
	int r, ret_val;

	ub::vector< autodif< interval<T> > > result;

	ub::vector< interval<T> > Iv, Fv;
	ub::matrix< interval<T> > Id, Fd;
	ub::matrix< interval<T> > fdI;
	ub::vector< interval<T> > fdI_tmp;
	ub::vector< psa< interval<T> > > result_tmp;

	interval<T> end2 = end;

	kv::autodif< interval<T> >::split(init, Iv, Id);

	Fv = Iv;
	r = ode(f, Fv, start, end2, p, &result_tmp);
	if (r == 0) return 0;
	ret_val = r;

	if (result_psa != NULL) {
		*result_psa = result_tmp;
	}

	MakeVariationalEq< F, interval<T> > g(f, result_tmp);

	fdI_tmp.resize(n * n);
	k = 0;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i == j) fdI_tmp(k) = 1.;
			else fdI_tmp(k) = 0.;
			k++;
		}
	}

	ode_param<T> p2 = p;

	p2.set_autostep(true);
	p2.set_epsilon(std::numeric_limits<T>::infinity());

	r = ode(g, fdI_tmp, start, end2, p2);
	if (r == 0) return 0;
	if (r == 1) {
		if (p.autostep == false) return 0;
		ret_val = 1;
	}

	fdI.resize(n, n);
	k = 0;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			fdI(i, j) = fdI_tmp(k++);
		}
	}

	Fd = prod(fdI, Id);

	result.resize(n);
	int s2 = Fd.size2();
	for (i=0; i<n; i++) {
		result(i).v = Fv(i);
		result(i).d.resize(s2);
		for (j=0; j<s2; j++) {
			result(i).d(j) = Fd(i, j);
		}
	}

	init = result;
	if (ret_val == 1) end = end2;

	return ret_val;
}


#else // ODE_AUTODIF_NEW != 1


#ifndef ODE_FAST
#define ODE_FAST 1
#endif

#ifndef ODE_RESTART_RATIO
#define ODE_RESTART_RATIO 1
#endif

#ifndef ODE_COEF_MID
#define ODE_CORF_MID 0
#endif


template <class T, class F>
int
ode(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL) {
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
	#if ODE_RESTART_RATIO == 1
	T max_ratio;
	#endif

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
				#if ODE_COEF_MID == 1
				using std::abs;
				m = std::max(m, abs(mid(x(i).v(j).v)));
				#else
				m = std::max(m, norm(x(i).v(j).v));
				#endif
				#ifdef IGNORE_DIF_PART
				#else
				km = x(i).v(j).d.size();
				for (k=0; k<km; k++) {
					#if ODE_COEF_MID == 1
					using std::abs;
					m = std::max(m, abs(mid(x(i).v(j).d(k))));
					#else
					m = std::max(m, norm(x(i).v(j).d(k)));
					#endif
				}
				#endif
			}
			if (m == 0.) continue;
			radius = std::max(radius, (T)std::pow((double)m, 1./j));
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

		psa< autodif< interval<T> > >::domain() = interval<T>(0., deltat.upper());

		z = x;
		t = setorder(torg, p.order);

		try {
			w = f(z, t);
		}
		catch (std::domain_error& e) {
			if (p.autostep && restart < p.restart_max) {
				psa< autodif< interval<T> > >::use_history() = false;
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
			z(i).v(p.order).d.resize(n, true);
			for (j=0; j<n; j++) {
				z(i).v(p.order).d(j) += newton_step(k++) * interval<T>(-1., 1.);
			}
		}

		if (p.autostep && ret_val != 2 && resized == false) {
			resized = true;
			m = (std::numeric_limits<T>::min)();
			for (i=0; i<n; i++) {
				evalz = eval(z(i), autodif< interval<T> >(deltat));
				m = std::max(m, rad(evalz.v) - rad(new_init(i).v));
				evalz.d.resize(n);
				for (j=0; j<n; j++) {
					m = std::max(m, rad(evalz.d(j)) - rad(new_init(i).d(j)));
				}
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
		w = new_init + w;

		flag = true;
		#if ODE_RESTART_RATIO == 1
		max_ratio = 0.;
		#endif
		for (i=0; i<n; i++) {
			#if ODE_RESTART_RATIO == 1
			max_ratio = std::max(max_ratio, width(w(i).v(p.order).v) / width(z(i).v(p.order).v));
			#endif
			flag = flag && subset(w(i).v(p.order).v, z(i).v(p.order).v);
			w(i).v(p.order).d.resize(n);
			for (j=0; j<n; j++) {
				#if ODE_RESTART_RATIO == 1
				max_ratio = std::max(max_ratio, width(w(i).v(p.order).d(j)) / width(z(i).v(p.order).d(j)));
				#endif
				flag = flag && subset(w(i).v(p.order).d(j), z(i).v(p.order).d(j));
			}
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
		if (result_psa != NULL) {
			// store w to *result_psa without autodif information
			(*result_psa).resize(n);
			for (i=0; i<n; i++) {
				(*result_psa)(i).v.resize(w(i).v.size());
				for (j=0; j<w(i).v.size(); j++) {
					(*result_psa)(i).v(j) = w(i).v(j).v;
				}
			}
		}
	}

	psa< autodif< interval<T> > >::mode() = save_mode;
	psa< autodif< interval<T> > >::use_history() = save_uh;
	psa< autodif< interval<T> > >::record_history() = save_rh;

	return ret_val;
}

#endif


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
