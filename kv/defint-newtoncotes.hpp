/*
 * Copyright (c) 2021 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DEFINT_NEWTONCOTES_HPP
#define DEFINT_NEWTONCOTES_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/highderiv.hpp>
#include <kv/optimize.hpp>

namespace kv {

template<class F, class T>
interval<T> defint_newtoncotes(F f, interval<T> start, interval<T> end, int n, int m, int errorterm_div = 1)
{
	int i, j;
	interval<T> x, y, w, h, k, r, tmp;

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
	
	if (n < 1 || n > 7) throw std::domain_error("defint_newtoncotes: 1<=n<=7 required.");

	w = end - start;

	// calculate interval higher derivative of f

	interval<T> v, hd;
	HighDeriv<F> fh(f, errord[n]);
	v = interval<T>::hull(start, end);
	hd = interval<T>::hull(kv::minimize_value(v, fh, errorterm_div), kv::maximize_value(v, fh, errorterm_div));

	//  Estimate best m which will minimize the error

	if (m == 0) {
		int tm; // temporary m for parameter estimation
		tm = 10;
		while (tm % n != 0) tm++;
		T tx, th, a, b;
		th = mid(w) / tm;
		r = 0;
		for (i=0; i<=tm; i++) {
			tx = mid(start) + i * th;
			if (i == 0 || i == m) k = coeff[n][0];
			else if (i % n == 0) k = coeff[n][n];
			else k = coeff[n][i % n];
			r += mid(k) * f(x);
		}
		r *= th;
		// error is estimated as a*m + b/m^p
		a = width(r) / tm;
		tmp = hd * w * pow(th, errord[n]) * errterm[n];
		b = width(tmp) * pow(tm, errord[n]);
                // add small value to avoid error
		a += std::numeric_limits<T>::epsilon();
		b += std::numeric_limits<T>::epsilon();
		// estimation of best m
		m = (int)ceil(pow(errord[n] * b / a, T(1)/(errord[n] + 1)));
		#ifdef NEWTONCOTES_SHOW_ESTIMATED_M
		std::cout << "estimated m: " << m << "\n";
		#endif
	}

	while (m % n != 0) m++;

	h = w / m;

	r = 0;
	for (i=0; i<=m; i++) {
		x = start + i * h;
		if (i == 0 || i == m) k = coeff[n][0];
		else if (i % n == 0) k = coeff[n][n];
		else k = coeff[n][i % n];
		r += k * f(x);
	}
	r *= h;

	r += hd * w * pow(h, errord[n]) * errterm[n];

	return r;
}

}; // namespace kv

#endif // DEFINT_NEWTONCOTES_HPP
