/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DOUBLEINTEGRAL_HPP
#define DOUBLEINTEGRAL_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>

namespace kv {

template <class T, class F>
interval<T>
doubleintegral(F f, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int order, int n, bool triangle = false) {
	interval<T> r1, r2, step1, step2, c1, c2, result, w;
	psa< psa< interval<T> > > x1, x2, y;
	psa< interval<T> > z;
	psa< interval<T> > upper_bound;
	int i, j;
	bool save_mode1, save_uh1, save_rh1;
	bool save_mode2, save_uh2, save_rh2;


	step1 = (end1 - start1) / (T)n;
	step2 = (end2 - start2) / (T)n;
	r1 = step1 / 2.;
	r2 = step2 / 2.;

	x1.v.resize(order+1);
	for (i=0; i<=order; i++) {
		x1.v(i).v.resize(order+1);
		for (j=0; j<=order; j++) {
			x1.v(i).v(j) = 0.;
		}
	}
	x1.v(1).v(0) = 1.;

	x2.v.resize(order+1);
	for (i=0; i<=order; i++) {
		x2.v(i).v.resize(order+1);
		for (j=0; j<=order; j++) {
			x2.v(i).v(j) = 0.;
		}
	}
	x2.v(0).v(1) = 1.;

	save_mode1 = psa< psa< interval<T> > >::mode();
	save_uh1 = psa< psa< interval<T> > >::use_history();
	save_rh1 = psa< psa< interval<T> > >::record_history();
	save_mode2 = psa< interval<T> >::mode();
	save_uh2 = psa< interval<T> >::use_history();
	save_rh2 = psa< interval<T> >::record_history();
	psa< psa< interval<T> > >::mode() = 2;
	psa< psa< interval<T> > >::use_history() = false;
	psa< psa< interval<T> > >::record_history() = false;
	psa< interval<T> >::mode() = 2;
	psa< interval<T> >::use_history() = false;
	psa< interval<T> >::record_history() = false;

	psa< psa< interval<T> > >::domain() = interval<T>(-r1.upper(), r1.upper());
	psa< interval<T> >::domain() = interval<T>(-r2.upper(), r2.upper());

	if (triangle) {
		upper_bound.v.resize(order+1);
		upper_bound.v(1) = -(end1 - start1) / (end2 - start2);
	}

	result = 0.;

	for (i=0; i<n; i++) {
		c2 = start2 + (T)i * step2 + r2;
		x2.v(0).v(0) = c2;
		for (j=0; j<n; j++) {
			if (triangle && i+j > n-1) continue;
			c1 = start1 + (T)j * step1 + r1;
			x1.v(0).v(0) = c1;
			y = integrate(f(x1, x2));
			if (triangle && i+j == n-1) {
				z = integrate(eval(y, upper_bound) - eval(y, (psa< interval<T> >)(-r1)));
			} else {
				z = integrate(eval(y, (psa< interval<T> >)r1) - eval(y, (psa< interval<T> >)(-r1)));
			}
			w = eval(z, r2) - eval(z, -r2);
			result += w;
		}
	}

	psa< psa< interval<T> > >::mode() = save_mode1;
	psa< psa< interval<T> > >::use_history() = save_uh1;
	psa< psa< interval<T> > >::record_history() = save_rh1;
	psa< interval<T> >::mode() = save_mode2;
	psa< interval<T> >::use_history() = save_uh2;
	psa< interval<T> >::record_history() = save_rh2;

	return result;
}

template <class F, class TT> class DoubleIntegralTriangleConv {
	public:
	F f;
	TT x1, y1, x2x1, x3x1, y2y1, y3y1, J;

	DoubleIntegralTriangleConv(F f_v, TT x1_v, TT y1_v, TT x2_v, TT y2_v, TT x3_v, TT y3_v) :
	f(f_v), x1(x1_v), y1(y1_v), x2x1(x2_v-x1_v), x3x1(x3_v-x1_v), y2y1(y2_v-y1_v), y3y1(y3_v-y1_v) {
		J = x2x1 * y3y1 - x3x1 * y2y1;
		// J = (J >= 0.)? J : -J;
		J = abs(J);
	}

	template <class T> T operator() (T u, T v) {
		return f(x1 + x2x1*u + x3x1*v, y1 + y2y1*u + y3y1*v) * J;
	}
};

template <class T, class F>
interval<T>
doubleintegral_triangle(F f, interval<T> x1, interval<T> y1, interval<T> x2, interval<T> y2, interval<T> x3, interval<T> y3, int order, int n) {
	DoubleIntegralTriangleConv< F, interval<T> > g(f, x1, y1, x2, y2, x3, y3);
	return doubleintegral(g, (interval<T>)0., (interval<T>)1., (interval<T>)0., (interval<T>)1., order, n, true);
}

} // namespace kv

#endif // DOUBLEINTEGRAL_HPP
