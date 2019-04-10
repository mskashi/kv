/*
 * Copyright (c) 2013-2019 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_MAFFINE_HPP
#define ODE_MAFFINE_HPP

//
// ODE using Affine and Mean Value Form
//
//  (2018/11/28) ode-maffine0 and ode-maffine are integrated by
//   porting maffine's algorithm to ode-autodif.hpp .
//
//  use -DODE_AUTODIF_NEW=0 to come back the behaviour of maffine0.
//

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


template <class T, class F>
int
ode_maffine(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>() , ub::matrix< interval<T> >* mat = NULL, ub::vector< psa< interval<T> > >* result_psa = NULL)
{
	int n = init.size();
	int i, j;

	ub::vector< interval<T> > c;
	ub::vector< interval<T> > fc;
	ub::vector< interval<T> > I;
	ub::vector< autodif< interval<T> > > Iad;

	ub::vector< interval<T> > result_i;
	ub::matrix< interval<T> > result_d;

	ub::vector< affine<T> > result;


	int maxnum_save;

	int r;

	int ret_val;
	interval<T> end2 = end;

	I.resize(n);
	c.resize(n);
	for (i=0; i<n; i++) {
		I(i) = to_interval(init(i));
		c(i) = mid(I(i));
	}

	Iad = autodif< interval<T> >::init(I);
	// NOTICE: below must be autodif version of ode
	r = ode(f, Iad, start, end2, p, result_psa);
	if (r == 0) return 0;
	ret_val = r;
	autodif< interval<T> >::split(Iad, result_i, result_d);

	fc = c;
	// Step size should be same as above ode call.
	// Because above ode call is with autodif and interval input and
	// below ode call is without autodif and point input,
	// below ode call is supposed to be easier to succeed than above.
	// If below ode call fails, force success by increasing order.
	ode_param<T> p2 = p;
	p2.set_autostep(false);
	while (true) {
		r = ode(f, fc, start, end2, p2);
		if (r != 0) break;
		p2.order++;
		if (p.verbose == 1) {
			std::cout << "ode_maffine: increase order: " << p.order << "\n";
		}
	}

	if (p.ep_reduce == 0) {
		maxnum_save = affine<T>::maxnum();
	}

	result = fc + prod(result_d, init - c);

	if (p.ep_reduce == 0) {
		epsilon_reduce2(result, maxnum_save);
	} else {
		epsilon_reduce(result, p.ep_reduce, p.ep_reduce_limit);
	}

	init = result;
	if (ret_val == 1) end = end2;
	if (mat != NULL) *mat = result_d;

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
