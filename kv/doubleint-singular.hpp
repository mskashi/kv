/*
 * Copyright (c) 2019 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DOUBLEINT_SINGULAR_HPP
#define DOUBLEINT_SINGULAR_HPP

#include <kv/defint-singular.hpp>

namespace kv {

/*
 * helper functions for double integral with singular edge
 */

// Reverse Time with argument 1 : convert f(x, y) to f(-x, y)
template <class F> class Doubleint_Reverse1 {
	F f;
	public:
	Doubleint_Reverse1(F f) : f(f) {}
	template <class T> T operator()(const T& x, const T& y) {
                return f(-x, y);
        }
};

// Reverse Time with argument 2 : convert f(x, y) to f(x, -y)
template <class F> class Doubleint_Reverse2 {
	F f;
	public:
	Doubleint_Reverse2(F f) : f(f) {}
	template <class T> T operator()(const T& x, const T& y) {
                return f(x, -y);
        }
};

// Reverse Time with both arguments : convert f(x, y) to f(-x, -y)
template <class F> class Doubleint_Reverse12 {
	F f;
	public:
	Doubleint_Reverse12(F f) : f(f) {}
	template <class T> T operator()(const T& x, const T& y) {
                return f(-x, -y);
        }
};

/*
 * double integral with singular edge
 *
 *   multiplicity1 multiplicity2
 *      not zero     not zero    : singular with edge start1=0, start2=0
 *      not zero        0        : singular with edge start1=0
 *         0         not zero    : singular with edge start2=0
 */          

template <class T, class F1, class F2>
interval<T>
doubleintegral_power3
(F1 f, F2 g, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int order, interval<T> power, int multiplicity1 = 1, int multiplicity2 = 1) {
	interval<T> step1, step2, result;
	psa< psa< interval<T> > > x1, x2, y;
	psa< interval<T> > z;
	int i, j;
	interval<T> tx, tp;
	bool save_mode1, save_uh1, save_rh1;
	bool save_mode2, save_uh2, save_rh2;


	step1 = end1 - start1;
	step2 = end2 - start2;

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

	psa< psa< interval<T> > >::domain() = interval<T>(0., step1.upper());
	psa< interval<T> >::domain() = interval<T>(0., step2.upper());

	x2.v(0).v(0) = start2;
	x1.v(0).v(0) = start1;

	y = f(x1, x2);
	y = setorder(y, order);
	y = div_tn(y, multiplicity1);
	for (i=0; i<y.v.size(); i++) {
		y.v(i) = setorder(y.v(i), order);
		y.v(i) = div_tn(y.v(i), multiplicity2);
	}

	y = pow(y, power) * g(x1, x2);
	
	z = 0.;
	tp = power * multiplicity1 + 1.;
	tx = pow(step1, tp);
	for (i=0; i<y.v.size(); i++) {
		z += y.v(i) * tx / tp;
		tp += 1.;
		tx *= step1;
	}

	result = 0.;
	tp = power * multiplicity2 + 1.;
	tx = pow(step2, tp);
	for (i=0; i<y.v.size(); i++) {
		result += z.v(i) * tx / tp;
		tp += 1.;
		tx *= step2;
	}

	psa< psa< interval<T> > >::mode() = save_mode1;
	psa< psa< interval<T> > >::use_history() = save_uh1;
	psa< psa< interval<T> > >::record_history() = save_rh1;
	psa< interval<T> >::mode() = save_mode2;
	psa< interval<T> >::use_history() = save_uh2;
	psa< interval<T> >::record_history() = save_rh2;

	return result;
}

/*
 * doubleintegral_power3 with reversed 1st argument
 */

template <class T, class F1, class F2>
interval<T>
doubleintegral_power3_r1
(F1 f, F2 g, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int order, interval<T> power, int multiplicity1 = 1, int multiplicity2 = 1) {
	return doubleintegral_power3(Doubleint_Reverse1<F1>(f), Doubleint_Reverse1<F2>(g), -end1, -start1, start2, end2, order, power, multiplicity1, multiplicity2);
}

/*
 * doubleintegral_power3 with reversed 2nd argument
 */

template <class T, class F1, class F2>
interval<T>
doubleintegral_power3_r2
(F1 f, F2 g, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int order, interval<T> power, int multiplicity1 = 1, int multiplicity2 = 1) {
	return doubleintegral_power3(Doubleint_Reverse2<F1>(f), Doubleint_Reverse2<F2>(g), start1, end1, -end2, -start2, order, power, multiplicity1, multiplicity2);
}

/*
 * doubleintegral_power3 with reversed 1st and 2nd arguments
 */

template <class T, class F1, class F2>
interval<T>
doubleintegral_power3_r12
(F1 f, F2 g, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int order, interval<T> power, int multiplicity1 = 1, int multiplicity2 = 1) {
	return doubleintegral_power3(Doubleint_Reverse12<F1>(f), Doubleint_Reverse12<F2>(g), -end1, -start1, -end2, -start2, order, power, multiplicity1, multiplicity2);
}


/*
 *  internal classes for doubleintegral_singular_point
 */

struct Func_f_for_singular_point {
	template <class T> T operator() (const T& r, const T& s) {
		return r;
	}
};

template <class F, class G, class TT>
struct Func_g_for_singular_point {
	F f;
	G g;
	TT x0, y0, x1x0, y1y0, x2x1, y2y1, p;
	int order;

	Func_g_for_singular_point(F f, G g, TT x0, TT y0, TT x1, TT y1, TT x2, TT y2, TT p, int order) : f(f), g(g), x0(x0), y0(y0), x1x0(x1-x0), y1y0(y1-y0), x2x1(x2-x1), y2y1(y2-y1), p(p), order(order) {}
	template <class T> T operator() (const T& r, const T& s) {
		T x, y;
		x = x0 + x1x0 * r + x2x1 * r * s;
		y = y0 + y1y0 * r + y2y1 * r * s;
		return pow(div_tn(f(x, y), order), p) * g(x, y);
	}
};

/*
 *  doubleintegral of (f^power g) over triangle (x0,y0)-(x1,y1)-(x2,y2)
 *  with singular point (x0,y0) 
 *
 *  multiplicity: multiplicity of zero of f at (x0,y0)
 */

template <class T, class F1, class F2>
interval<T>
doubleintegral_singular_point
(F1 f, F2 g, interval<T> x0, interval<T> y0, interval<T> x1, interval<T> y1, interval<T> x2, interval<T>y2, int order, interval<T> power, int multiplicity = 1, int div = 4)
{

	interval<T> det, r;
	int i;

	using std::abs;
	det = abs((x1-x0) * (y2-y1) - (x2-x1) * (y1-y0));

	r = 0.;
	for (i=0; i<div; i++) {

		r += doubleintegral_power3(
			Func_f_for_singular_point(),
			Func_g_for_singular_point<F1,F2,interval<T>>(f, g, x0, y0, x1, y1, x2, y2, power, multiplicity),
			interval<T>(0.),
			interval<T>(1.),
			interval<T>(1. / div * i),
			interval<T>(1. / div * (i+1)),
			order, 
			multiplicity * power + 1,
			1,
			0);
	}

	return r * det;
}

} // namespace kv

#endif // DOUBLEINT_SINGULAR_HPP
