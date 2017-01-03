/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef POINCAREMAP_HPP
#define POINCAREMAP_HPP

#include <kv/strobomap.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template <class F1, class F2, class T> class PoincareMap {
	public:
	F1 f1;
	F2 f2;
	interval<T> start;
	ode_param<T> p;

	PoincareMap(F1 f1, F2 f2, interval<T> start, ode_param<T> p = ode_param<T>())
	: f1(f1), f2(f2), start(start), p(p) {}

	ub::vector<T> operator() (const ub::vector<T>& x) {
		int s = x.size();
		int i;
		ub::vector<T> x2(s-1);
		ub::vector<T> y;
		ub::vector<T> r(s);
		interval<T> t;

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1);

		StroboMap<F1, T> st(f1, start, t, p);

		y = st(x2);

		for (i=0; i<s-1; i++) r(i) = y(i) - x(i);
		r(s-1) = f2(x);

		return r;
	}

	ub::vector< autodif<T> > operator() (const ub::vector< autodif<T> >& x) {
		int s = x.size();
		int i;
		ub::vector< autodif<T> > x2(s-1);
		ub::vector< autodif<T> > y;
		ub::vector< autodif<T> > r(s);
		interval<T> t;
		ub::vector<T> y2(s-1);
		ub::vector<T> dxdt(s-1);

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1).v;

		StroboMap<F1, T> st(f1, start, t, p);

		y = st(x2);

		/*
		  adding "derivative of solution w.r.t. end time" afterwards
		  using dx/dt = f(x,t) and dx/dp = (dx/dt) * (dt/dp)
		*/
		for (i=0; i<s-1; i++) y2(i) = y(i).v;
		dxdt = f1(y2, mid(t));

		for (i=0; i<s-1; i++) {
			y(i).d += dxdt(i) * x(s-1).d;
		}

		for (i=0; i<s-1; i++) r(i) = y(i) - x(i);
		r(s-1) = f2(x);

		return r;
	}

	ub::vector< interval<T> > operator() (const ub::vector< interval<T> >& x) {
		int s = x.size();
		int i;
		ub::vector< interval<T> > x2(s-1);
		ub::vector< interval<T> > y;
		ub::vector< interval<T> > r(s);
		interval<T> t;

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1);

		StroboMap<F1, T> st(f1, start, t, p);

		y = st(x2);

		for (i=0; i<s-1; i++) r(i) = y(i) - x(i);
		r(s-1) = f2(x);

		return r;
	}

	ub::vector< autodif< interval<T> > > operator() (const ub::vector< autodif< interval<T> > >& x) {
		int s = x.size();
		int i;
		ub::vector< autodif< interval<T> > > x2(s-1);
		ub::vector< autodif< interval<T> > > y;
		ub::vector< autodif< interval<T> > > r(s);
		interval<T> t;
		ub::vector< interval<T> > y2(s-1);
		ub::vector< interval<T> > dxdt(s-1);

		for (i=0; i<s-1; i++) x2(i) = x(i);
		t = x(s-1).v;

		StroboMap<F1, T> st(f1, start, t, p);

		y = st(x2);

		/*
		  adding "derivative of solution w.r.t. end time" afterwards
		  using dx/dt = f(x,t) and dx/dp = (dx/dt) * (dt/dp)
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
