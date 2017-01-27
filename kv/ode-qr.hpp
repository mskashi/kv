/*
 * Copyright (c) 2013-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_QR_HPP
#define ODE_QR_HPP

// ODE using QR Decomposition

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/qr.hpp>
#include <kv/vleq.hpp>
#include <kv/ode.hpp>
#include <kv/ode-autodif.hpp>
#include <kv/ode-param.hpp>
#include <kv/ode-callback.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T, class F>
int
odelong_qr(
	F f,
	ub::vector< interval<T> >& init,
	const interval<T>& start,
	interval<T>& end,
	ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>(),
	ub::matrix< interval<T> >* mat = NULL
) {
	int s = init.size();
	int i, j;

	ub::vector< interval<T> > c;
	ub::vector< interval<T> > fc;
	ub::vector< autodif< interval<T> > > Iad;

	ub::vector< interval<T> > result_i;
	ub::matrix< interval<T> > result_d;

	ub::vector< interval<T> > x, x1;
	interval<T> t, t1;
	ub::matrix< interval<T> > M;
	int ret_ode, ret_ode2;
	int ret_val = 0;
	bool bo;
	bool ret_callback;

	ub::matrix<T> Q, Q2, R, Q2t;
	ub::matrix< interval<T> > AQ, QAQ, Q2i;
	ub::vector< interval<T> > y, y1, y2, tmp;

	ub::vector< psa< interval<T> > > result_psa;
	ub::vector< psa< autodif< interval<T> > > > result_tmp;


	if (mat != NULL) {
		M = ub::identity_matrix< interval<T> >(s);
	}

	t = start;
	x = init;

	c = mid(x);
	y = x - c;
	Q = ub::identity_matrix<T>(s);

	ode_param<T> p2;

	while (1) {
		x1 = x;
		t1 = end;

		Iad = autodif< interval<T> >::init(x1);
		p2 = p;
		p2.set_autostep(true);
		// NOTICE: below must be autodif version of ode
		ret_ode = ode(f, Iad, t, t1, p2, &result_tmp);
		if (ret_ode == 0) break;

		fc = c;
		// Step size should be same as above ode call.
		// Because above ode call is with autodif and interval input and
		// below ode call is without autodif and point input,
		// below ode call is supposed to be easier to succeed than above.
		// If below ode call fails, force success by increasing order.
		p2 = p;
		p2.set_autostep(false);
		while (1) {
			ret_ode2 = ode(f, fc, t, t1, p2);
			if (ret_ode2 != 0) break;
			p2.order++;
			std::cout << "increase order: " << p2.order << "\n";
		}

		autodif< interval<T> >::split(Iad, result_i, result_d);

		#if 0
		// centering result_d
		fc += prod(result_d - mid(result_d), x - c);
		result_d =  mid(result_d);
		#endif

		AQ = prod(result_d, Q);
		bo = qr(mid(AQ), Q2, R);
		if (bo == false) break;
		Q2i = Q2;
		Q2t = trans(Q2);
		// bo = vleq(Q2i, AQ, QAQ);
		bo = vleq(Q2i, AQ, QAQ, &Q2t);
		if (bo == false) break;
		Q2i = Q2;
		y1 = prod(QAQ, y);
		c = mid(fc);
		tmp = fc - c;
		// bo = vleq(Q2i, tmp, y2);
		bo = vleq(Q2i, tmp, y2, &Q2t);
		if (bo == false) break;
		y = y1 + y2;
		x1 = prod(Q2, y) + c;

		// below seems to have some efficiency.
		// we comment out below because we have not study it
		// theoretically yet.

		// x1 = intersect(x1, result_i);

		Q = Q2;

		ret_val = 1;

		if (mat != NULL) M = prod(result_d, M);

		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << x1 << "\n";
		}

		// store result_tmp to result_psa without autodif information
		result_psa.resize(s);
		for (i=0; i<s; i++) {
			result_psa(i).v.resize(result_tmp(i).v.size());
			for (j=0; j<result_tmp(i).v.size(); j++) {
				result_psa(i).v(j) = result_tmp(i).v(j).v;
			}
		}
		ret_callback = callback(t, t1, x, x1, result_psa);

		if (ret_callback == false) {
			ret_val = 3;
			break;
		}

		if (ret_ode == 2) {
			ret_val = 2;
			break;
		}

		t = t1;
		x = x1;
	}

	if (ret_val >= 1) {
		init = x1;
		if (mat != NULL) *mat = M;
	}
	if (ret_val == 1) {
		end = t;
	}
	if (ret_val == 3) {
		end = t1;
	}

	return ret_val;
}


template <class T, class F>
int
odelong_qr(
	F f,
	ub::vector< autodif< interval<T> > >& init,
	const interval<T>& start,
	interval<T>& end, ode_param<T> p = ode_param<T>(),
	const ode_callback<T>& callback = ode_callback<T>()
) {
	int s = init.size();
	int i, j;
	ub::vector< interval<T> > x;
	ub::matrix< interval<T> > M, M_tmp;
	int r;

	autodif< interval<T> >::split(init, x, M);
	int s2 = M.size2();

	r = odelong_qr(f, x, start, end, p, callback, &M_tmp);

	if (r == 0) return 0;

	M = prod(M_tmp, M);

	for (i=0; i<s; i++) {
		init(i).v = x(i);
		init(i).d.resize(s2);
		for (j=0; j<s2; j++) {
			init(i).d(j) = M(i, j);
		}
	}
	
	return r;
}

} // namespace kv

#endif // ODE_QR_HPP
