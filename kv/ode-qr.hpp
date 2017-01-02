/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
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

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T, class F>
int
odelong_qr(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>(), ub::matrix< interval<T> >* mat = NULL)
{
	int s = init.size();
	int i, j;

	ub::vector< interval<T> > c;
	ub::vector< interval<T> > fc;
	ub::vector< autodif< interval<T> > > Iad;

	ub::vector< interval<T> > result_i;
	ub::matrix< interval<T> > result_d;

	ub::vector< interval<T> > x;
	interval<T> t, t1;
	ub::matrix< interval<T> > M;
	int r, r2;
	int ret_val = 0;
	bool bo;

	ub::matrix<T> Q, Q2, R, Q2t;
	ub::matrix< interval<T> > AQ, QAQ, Q2i;
	ub::vector< interval<T> > y, y1, y2, tmp;

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
		t1 = end;
		Iad = autodif< interval<T> >::init(x);
		p2 = p;
		p2.set_autostep(true);
		// 自動的にautodif対応のものが呼ばれるはず
		r = ode(f, Iad, t, t1, p2);
		if (r == 0) break;

		fc = c;
		// step幅は上のodeと同じでないとまずい。
		// 上が通れば条件の緩いこちらは通るかと思ったが、
		// 通らない例があった。本質的には通るはずなので、次数を
		// 増やしてでも無理矢理通す。
		p2 = p;
		p2.set_autostep(false);
		while (1) {
			r2 = ode(f, fc, t, t1, p2);
			if (r2 != 0) break;
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
		bo = vleq(Q2i, AQ, QAQ, Q2t);
		if (bo == false) break;
		Q2i = Q2;
		y1 = prod(QAQ, y);
		c = mid(fc);
		tmp = fc - c;
		// bo = vleq(Q2i, tmp, y2);
		bo = vleq(Q2i, tmp, y2, Q2t);
		if (bo == false) break;
		y = y1 + y2;
		x = prod(Q2, y) + c;
		Q = Q2;

		ret_val = 1;

		if (mat != NULL) M = prod(result_d, M);

		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << x << "\n";
		}

		if (r == 2) {
			ret_val = 2;
			break;
		}

		t = t1;
	}

	if (ret_val >= 1) {
		init = x;
		if (mat != NULL) *mat = M;
	}
	if (ret_val == 1) {
		end = t;
	}

	return ret_val;
}


template <class T, class F>
int
odelong_qr(F f, ub::vector< autodif< interval<T> > >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>())
{
	int s = init.size();
	int i, j;
	ub::vector< interval<T> > x;
	ub::matrix< interval<T> > M, M_tmp;
	int r;
	interval<T> end2 = end;

	autodif< interval<T> >::split(init, x, M);
	int s2 = M.size2();

	r = odelong_qr(f, x, start, end2, p, &M_tmp);

	if (r == 0) return 0;

	M = prod(M_tmp, M);

	for (i=0; i<s; i++) {
		init(i).v = x(i);
		init(i).d.resize(s2);
		for (j=0; j<s2; j++) {
			init(i).d(j) = M(i, j);
		}
	}
	
	if (r == 1) end = end2;

	return r;
}

} // namespace kv

#endif // ODE_QR_HPP
