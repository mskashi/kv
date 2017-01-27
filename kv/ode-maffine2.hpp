/*
 * Copyright (c) 2013-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_MAFFINE2_HPP
#define ODE_MAFFINE2_HPP

// ODE using Affine and Mean Value Form (fast)

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/psa.hpp>
#include <kv/affine.hpp>
#include <kv/ode.hpp>
#include <kv/ode-autodif.hpp>
#include <kv/ode-param.hpp>
#include <kv/ode-callback.hpp>


#ifndef ODE_FAST
#define ODE_FAST 1
#endif


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T, class F>
void
ode_onlytype1(F f, ub::vector< interval<T> >& init, const interval<T>& start, const interval<T>& end, int order) {
	int n = init.size();
	int i, j;

	ub::vector< psa< interval<T> > > x, y;
	psa< interval<T> > torg;
	psa< interval<T> > t;

	ub::vector< interval<T> > result;

	interval<T> deltat;

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
	#if ODE_FAST == 1
	psa< interval<T> >::record_history() = true;
	psa< interval<T> >::history().clear();
	#endif
	for (j=0; j<order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< interval<T> >::use_history() = true;
		if (j == order - 1) psa< interval<T> >::record_history() = false;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) y(i) = integrate(y(i));
		x = init + y;
	}

	deltat = end - start;

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(x(i), deltat);
	}

	init = result;

	psa< interval<T> >::mode() = save_mode;
	psa< interval<T> >::use_history() = save_uh;
	psa< interval<T> >::record_history() = save_rh;
}


template <class T, class F>
void
ode_onlytype1(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, const interval<T>& end, int order) {
	int n = init.size();
	int i, j, k;

	ub::vector< psa< autodif< interval<T> > > > x, y;
	psa< autodif< interval<T> > > torg;
	psa< autodif< interval<T> > > t;

	ub::vector< autodif< interval<T> > > result;

	interval<T> deltat;

	bool save_mode, save_uh, save_rh;


	x = init;

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
	for (j=0; j<order; j++) {
		#if ODE_FAST == 1
		if (j == 1) psa< autodif< interval<T> > >::use_history() = true;
		if (j == order - 1) psa< autodif< interval<T> > >::record_history() = false;
		#endif
		t = setorder(torg, j);
		y = f(x, t);
		for (i=0; i<n; i++) y(i) = integrate(y(i));
		x = init + y;
	}

	deltat = end - start;

	result.resize(n);
	for (i=0; i<n; i++) {
		result(i) = eval(x(i), (autodif< interval<T> >)deltat);
	}

	init = result;

	psa< autodif< interval<T> > >::mode() = save_mode;
	psa< autodif< interval<T> > >::use_history() = save_uh;
	psa< autodif< interval<T> > >::record_history() = save_rh;
}


template <class T, class F>
int
ode_maffine2(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::vector< psa< interval<T> > >* result_psa = NULL)
{
	int n = init.size();
	int i, j;

	ub::vector< interval<T> > c;
	ub::vector< interval<T> > fc;
	ub::vector< interval<T> > I, Idummy, I2;
	ub::vector< autodif< interval<T> > > Iad;

	ub::vector< interval<T> > result_i;
	ub::matrix< interval<T> > result_d;

	ub::vector< affine<T> > result;

	int maxnum_save;

	affine<T> s1, s2;
	interval<T> s2i;
	ub::vector< affine<T> > s1_save;
	ub::vector< interval<T> > s2i_save;

	interval<T> deltat, deltat_n;
	ub::vector< psa< interval<T> > > psa_result;

	int r;

	interval<T> end2 = end;


	I.resize(n);
	c.resize(n);
	for (i=0; i<n; i++) {
		I(i) = to_interval(init(i));
		c(i) = mid(I(i));
	}

	Idummy = I;
	r = ode(f, Idummy, start, end2, p, &psa_result);
	if (r == 0) return 0;

	if (result_psa != NULL) {
		*result_psa = psa_result;
	}

	deltat = end2 - start;
	deltat_n = 1.;
	for (i=0; i<p.order; i++) deltat_n *= deltat;

	I2.resize(n);
	for (i=0; i<n; i++) {
		I2(i) = psa_result(i).v(p.order) * deltat_n;
	}

	Iad = autodif< interval<T> >::init(I);
	// NOTICE: below must be autodif version
	ode_onlytype1(f, Iad, start, end2, p.order-1);

	fc = c;
	ode_onlytype1(f, fc, start, end2, p.order-1);

	autodif< interval<T> >::split(Iad, result_i, result_d);

	if (p.ep_reduce == 0) {
		maxnum_save = affine<T>::maxnum();
	}

	result = I2 + fc + prod(result_d, init - c);

	if (p.ep_reduce == 0) {
		s1_save.resize(n);
		s2i_save.resize(n);
		for (i=0; i<n; i++) {
			split(result(i), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			s1_save(i) = s1;
			s2i_save(i) = s2i;
		}

		affine<T>::maxnum() = maxnum_save;
		for (i=0; i<n; i++) {
			s1_save(i).resize();
			result(i) = append(s1_save(i), (affine<T>)s2i_save(i));
		}
	} else {
		epsilon_reduce(result, p.ep_reduce, p.ep_reduce_limit);
	}

	init = result;
	if (r == 1) end = end2;

	return r;
}

template <class T, class F>
int
odelong_maffine2(
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

		ret_ode = ode_maffine2(f, x1, t, t1, p, &result_tmp);
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

		ret_callback = callback(t, t1, to_interval(x), to_interval(x1), result_tmp);

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
odelong_maffine2(
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

	r = odelong_maffine2(f, x, start, end, p, callback);

	affine<T>::maxnum() = maxnum_save;

	if (r == 0) return 0;

	for (i=0; i<s; i++) init(i) = to_interval(x(i));

	return r;
}

} // namespace kv

#endif // ODE_MAFFINE2_HPP
