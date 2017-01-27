/*
 * Copyright (c) 2016-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_MAFFINE_HPP
#define ODE_MAFFINE_HPP

// ODE using Affine and Mean Value Form (new)

#include <iostream>
#include <list>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
// #include <boost/tuple/tuple.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/affine.hpp>
#include <kv/ode.hpp>
#include <kv/ode-autodif.hpp>
#include <kv/ode-param.hpp>
#include <kv/ode-callback.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;


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
ode_maffine(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>() , ub::matrix< interval<T> >* mat = NULL, ub::vector< psa< interval<T> > >* result_psa = NULL)
{
	int n = init.size();
	int i, j, k;

	ub::vector< interval<T> > c;
	ub::vector< interval<T> > fc;
	ub::vector< interval<T> > I;
	ub::vector< interval<T> > fI;
	ub::matrix< interval<T> > fdI;
	ub::vector< interval<T> > fdI_tmp;

	ub::vector< affine<T> > result;

	ub::vector< psa< interval<T> > > result_tmp;


	int maxnum_save;

	affine<T> s1, s2;
	interval<T> s2i;
	ub::vector< affine<T> > s1_save;
	ub::vector< interval<T> > s2i_save;

	int r;

	int ret_val;
	interval<T> end2 = end;

	ode_param<T> p2 = p;
	p2.set_autostep(false);

	I.resize(n);
	c.resize(n);
	for (i=0; i<n; i++) {
		I(i) = to_interval(init(i));
		c(i) = mid(I(i));
	}

	fI = I;
	r = ode(f, fI, start, end2, p, &result_tmp);
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

	#if 0
	// fixed stepsize
	r = ode(g, fdI_tmp, start, end2, p2);
	if (r == 0) {
		if (p.autostep == false) return 0;
		// autostep
		r = ode(g, fdI_tmp, start, end2, p);
		if (r == 0) return 0;
		ret_val = 1;
	}
	#endif
	r = ode(g, fdI_tmp, start, end2, p2);
	if (r == 0) {
		if (p.autostep == false) return 0;
		for (i=0; i<p.restart_max; i++) {
			end2 = start + (end2 - start) * 0.5;
			r = ode(g, fdI_tmp, start, end2, p2);
			if (r != 0) {
				ret_val = 1;
				break;
			}
		}
		if (r == 0) return 0;
	}

	fdI.resize(n, n);
	k = 0;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			fdI(i, j) = fdI_tmp(k++);
		}
	}

	fc = c;
	while (true) {
		// fixed stepsize
		r = ode(f, fc, start, end2, p2);
		if (r != 0) break;
		p2.order++;
		if (p.verbose == 1) {
			std::cout << "ode_maffine: increase order: " << p2.order << "\n";
		}
	}

	if (p.ep_reduce == 0) {
		maxnum_save = affine<T>::maxnum();
	}

	result = fc + prod(fdI, init - c);

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
	if (ret_val == 1) end = end2;
	if (mat != NULL) *mat = fdI;

	return ret_val;
}

template <class T, class F>
int
odelong_maffine(
	F f,
	ub::vector< affine<T> >& init,
	const interval<T>& start,
	interval<T>& end,
	ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>(),
	ub::matrix< interval<T> >* mat = NULL
) {

	int s = init.size();
	ub::vector< affine<T> > x, x1;
	interval<T> t, t1;
	int ret_ode;
	ub::matrix< interval<T> > M, M_tmp;
	ub::matrix< interval<T> >* M_p;
	int ret_val = 0;
	bool ret_callback;

	ub::vector< psa< interval<T> > > result_tmp;


	if (mat == NULL) {
		M_p = NULL;
	} else {
		M_p = &M_tmp;
		M = ub::identity_matrix< interval<T> >(s);
	}

	x = init;
	t = start;
	p.set_autostep(true);

	while (true) {
		x1 = x;
		t1 = end;

		ret_ode = ode_maffine(f, x1, t, t1, p, M_p, &result_tmp);
		if (ret_ode == 0) {
			if (ret_val == 1) {
				init = x1;
				if (mat != NULL) *mat = M;
				end = t;
			}
			return ret_val;
		}
		ret_val = 1;
		if (mat != NULL) M = prod(M_tmp, M);
		#if 0
		if (result_psa != NULL) {
			(*result_psa).push_back(boost::make_tuple(t, t1, result_tmp));
		}
		#endif
		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << to_interval(x1) << "\n";
		}

		ret_callback = callback(t, t1, to_interval(x), to_interval(x1), result_tmp);

		if (ret_callback == false) {
			init = x1;
			if (mat != NULL) *mat = M;
			end = t1;
			return 3;
		}

		if (ret_ode == 2) {
			init = x1;
			if (mat != NULL) *mat = M;
			return 2;
		}

		t = t1;
		x = x1;
	}
}

template <class T, class F>
int
odelong_maffine(
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

	r = odelong_maffine(f, x, start, end, p, callback);

	affine<T>::maxnum() = maxnum_save;

	if (r == 0) return 0;

	for (i=0; i<s; i++) init(i) = to_interval(x(i));

	return r;
}

template <class T, class F>
int
odelong_maffine(
	F f,
	ub::vector< autodif< interval<T> > >& init,
	const interval<T>& start,
	interval<T>& end,
	ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>()
) {
	int s = init.size();
	int i, j;
	ub::vector< interval<T> > xi;
	ub::vector< affine<T> > x;
	ub::matrix< interval<T> > M, M_tmp;
	int maxnum_save;
	int r;

	autodif< interval<T> >::split(init, xi, M);
	int s2 = M.size2();

	maxnum_save = affine<T>::maxnum();
	affine<T>::maxnum() = 0;

	x = xi;

	r = odelong_maffine(f, x, start, end, p, callback, &M_tmp);

	affine<T>::maxnum() = maxnum_save;

	if (r == 0) return 0;

	M = prod(M_tmp, M);

	for (i=0; i<s; i++) {
		init(i).v = to_interval(x(i));
		init(i).d.resize(s2);
		for (j=0; j<s2; j++) {
			init(i).d(j) = M(i, j);
		}
	}
	
	return r;
}

} // namespace kv

#endif // ODE_MAFFINE_HPP
