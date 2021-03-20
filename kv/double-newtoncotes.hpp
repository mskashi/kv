/*
 * Copyright (c) 2021 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DOUBLE_NEWTONCOTES_HPP
#define DOUBLE_NEWTONCOTES_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/highderiv.hpp>
#include <kv/optimize.hpp>

namespace kv {

namespace ub = boost::numeric::ublas;

template <class F>
struct Two_Arg_To_Vec {
	F f;
	Two_Arg_To_Vec(F f) : f(f) {}
	template <class T> T operator()(const ub::vector<T>& x) {
		return f(x(0), x(1));
	}
};

template<class F, class T>
interval<T> double_newtoncotes(F f, interval<T> start1, interval<T> end1, interval<T> start2, interval<T> end2, int n, int m1, int m2, int errorterm_div = 1)
{
	int i, j;
	interval<T> x, y, w1, w2, h1, h2, k1, k2, r, tmp;

	interval<T> coeff[8][8] = {
		{interval<T>(0)}, //dummy
		{interval<T>(1)/2, interval<T>(1)},
		{interval<T>(1)/3, interval<T>(4)/3, interval<T>(2)/3},
		{interval<T>(3)/8, interval<T>(9)/8, interval<T>(9)/8, interval<T>(3)/4},
		{interval<T>(14)/45, interval<T>(64)/45, interval<T>(8)/15, interval<T>(64)/45, interval<T>(28)/45},
		{interval<T>(95)/288, interval<T>(125)/96, interval<T>(125)/144, interval<T>(125)/144, interval<T>(125)/96, interval<T>(95)/144},
		{interval<T>(41)/140, interval<T>(54)/35, interval<T>(27)/140, interval<T>(68)/35, interval<T>(27)/140, interval<T>(54)/35, interval<T>(41)/70},
		{interval<T>(5257)/17280, interval<T>(25039)/17280, interval<T>(343)/640, interval<T>(20923)/17280, interval<T>(20923)/17280, interval<T>(343)/640, interval<T>(25039)/17280, interval<T>(5257)/8640}
	};
	interval<T> errterm[8] = {
		interval<T>(0), //dummy
		interval<T>(-1)/12,
		interval<T>(-1)/180,
		interval<T>(-3)/80,
		interval<T>(-2)/945,
		interval<T>(-275)/12096,
		interval<T>(-3)/2800,
		interval<T>(-8183)/518400
	};
	int errord[8] = {
		0, // dummy
		2,
		4,
		4,
		6,
		6,
		8,
		8
	};
	
	if (n < 1 || n > 7) throw std::domain_error("double_newtoncotes: 1<=n<=7 required.");

	w1 = end1 - start1;
	w2 = end2 - start2;

	// calculate interval higher derivative of f

	Two_Arg_To_Vec<F> fv(f);

	PartialDeriv< Two_Arg_To_Vec<F> > f_xh(fv, 0, errord[n]);
	PartialDeriv< Two_Arg_To_Vec<F> > f_yh(fv, 1, errord[n]);

	ub::vector< interval<T> > v(2);
	v(0) = interval<T>::hull(start1, end1);
	v(1) = interval<T>::hull(start2, end2);

	interval<T> hd_x, hd_y;
	hd_x = interval<T>::hull(kv::minimize_value(v, f_xh, errorterm_div), kv::maximize_value(v, f_xh, errorterm_div));
	hd_y = interval<T>::hull(kv::minimize_value(v, f_yh, errorterm_div), kv::maximize_value(v, f_yh, errorterm_div));

	//  Estimate best m1,m2 which will minimize the error

	if (m1 == 0 || m2 == 0) {
		int tm1, tm2; // temporary m1,m2 for parameter estimation
		tm1 = tm2 = 10;
		while (tm1 % n != 0) tm1++;
		while (tm2 % n != 0) tm2++;
		T tx, ty, th1, th2, a, b1, b2;
		th1 = mid(w1) / tm1;
		th2 = mid(w2) / tm2;

		r = 0;
		for (i=0; i<=tm1; i++) {
			tx = mid(start1) + i * th1;
			if (i == 0 || i == m1) k1 = coeff[n][0];
			else if (i % n == 0) k1 = coeff[n][n];
			else k1 = coeff[n][i % n];
			for (j=0; j<=tm2; j++) {
				ty = mid(start2) + j * th2;
				if (j == 0 || j == m2) k2 = coeff[n][0];
				else if (j % n == 0) k2 = coeff[n][n];
				else k2 = coeff[n][j % n];
				r += mid(k1 * k2) * f(x, y);
			}
		}
		r *= th1 * th2;
		// error is estimated as a*m1*m2 + b1/m1^p + b2/m2^p
		a = width(r) / (tm1 * tm2);
		tmp = hd_x * w1 * w2 * pow(th1, errord[n]) * errterm[n];
		b1 = width(tmp) * pow(tm1, errord[n]);
		tmp = hd_y * w1 * w2 * pow(th2, errord[n]) * errterm[n];
		b2 = width(tmp) * pow(tm2, errord[n]);
		// add small value to avoid error
		a += std::numeric_limits<T>::epsilon();
		b1 += std::numeric_limits<T>::epsilon();
		b2 += std::numeric_limits<T>::epsilon();
		// estimation of best m1, m2
		T m1m2;
		m1m2 = pow(errord[n] * errord[n] * b1 * b2 / (a * a), T(1)/(errord[n] + 2));
		m1 = (int)ceil(sqrt(m1m2 * pow(b1/b2, T(1)/errord[n])));
		m2 = (int)ceil(sqrt(m1m2 * pow(b2/b1, T(1)/errord[n])));
		#ifdef NEWTONCOTES_SHOW_ESTIMATED_M
		std::cout << "estimated m1,m2: " << m1 << ", " << m2 << "\n";
		#endif
	}

	while (m1 % n != 0) m1++;
	while (m2 % n != 0) m2++;

	h1 = w1 / m1;
	h2 = w2 / m2;

	r = 0;
	for (i=0; i<=m1; i++) {
		x = start1 + i * h1;
		if (i == 0 || i == m1) k1 = coeff[n][0];
		else if (i % n == 0) k1 = coeff[n][n];
		else k1 = coeff[n][i % n];
		for (j=0; j<=m2; j++) {
			y = start2 + j * h2;
			if (j == 0 || j == m2) k2 = coeff[n][0];
			else if (j % n == 0) k2 = coeff[n][n];
			else k2 = coeff[n][j % n];
			r += (k1 * k2) * f(x, y);
		}
	}
	r *= h1 * h2;

	r += hd_x * w1 * w2 * pow(h1, errord[n]) * errterm[n];
	r += hd_y * w1 * w2 * pow(h2, errord[n]) * errterm[n];

	return r;
}

}; // namespace kv

#endif // DOUBLE_NEWTONCOTES_HPP
