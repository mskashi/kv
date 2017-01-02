#ifndef POINCAREMAP_HPP
#define POINCAREMAP_HPP

#include "strobomap.hpp"

namespace ub = boost::numeric::ublas;

namespace kv {

template <class F1, class F2, class TT> class PoincareMap {
	public:
	F1 f1;
	F2 f2;
	TT start;
	int order;
	int iter_max;
	int verbose;
	int ep_reduce;
	int ep_reduce_limit;

	PoincareMap(F1 f1_v, F2 f2_v, TT start_v, int order_v, int iter_max_v = 2, int verbose_v = 0, int ep_reduce_v = 0, int ep_reduce_limit_v = 0)
	: f1(f1_v), f2(f2_v), start(start_v), order(order_v), iter_max(iter_max_v), verbose(verbose_v), ep_reduce(ep_reduce_v), ep_reduce_limit(ep_reduce_limit_v) {
        }

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
		int s = x.size();
		int i;
		ub::vector<T> x2(s-1);
		ub::vector<T> y;
		ub::vector<T> r(s);
		T t;

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1);

		StroboMap<F1, TT> st(f1, start, TT(t), order, iter_max, verbose, ep_reduce, ep_reduce_limit);

		y = st(x2);

		for (i=0; i<s-1; i++) r(i) = y(i) - x(i);
		r(s-1) = f2(x);

		return r;
	}

	template <class T> ub::vector< autodif<T> > operator() (const ub::vector< autodif<T> >& x) {
		int s = x.size();
		int i;
		ub::vector< autodif<T> > x2(s-1);
		ub::vector< autodif<T> > y;
		ub::vector< autodif<T> > r(s);
		T t;
		ub::vector<T> y2(s-1);
		ub::vector<T> dxdt(s-1);

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1).v;

		StroboMap<F1, TT> st(f1, start, TT(t), order, iter_max, verbose, ep_reduce, ep_reduce_limit);

		y = st(x2);

		/*
		  解の終了時刻tに関する微分が出来ないので、後付けしている。
		  dx/dt = f(x,t)なので、解をtで微分したものは方程式の
		  右辺に解と終了時刻を代入すると得られる。
		  pをそれに関して微分している変数とすると、
		  dx/dp = (dx/dt) * (dt/dp)
		*/
		for (i=0; i<s-1; i++) y2(i) = y(i).v;
		dxdt = f1(y2, t);

		for (i=0; i<s-1; i++) {
			y(i).d += dxdt(i) * x(s-1).d;
		}

		for (i=0; i<s-1; i++) r(i) = y(i) - x(i);
		r(s-1) = f2(x);

		return r;
	}
};

} // namespace kv

#endif // POINCAREMAP_HPP
