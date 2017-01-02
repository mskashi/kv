/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RKF45_HPP
#define RKF45_HPP

// Runge-Kutta-Fehlberg Method (RKF45)

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ub = boost::numeric::ublas;

namespace kv {

template <class T, class F>
void
rkf45(F f, ub::vector<T>& init, T start, T end, T init_h, T tol, int verbose = 0) {

	int i;
	int size = init.size();

	ub::vector<T> x, x4, x5, k1, k2, k3, k4, k5, k6, tmp;
	bool flag = false;

	T t = start;

	T h = init_h;
	T h_min = init_h / 64;
	T h_max = init_h * 64;

	T err, scale;

	x = init;

	while (1) {
		if (t + h >= end) {
			h = end - t;
			flag = true;
		}

		k1 = h * f(x, t);

		// k2 = h * f(x + k1/4., t + h/4.);
		tmp = x + k1/4.;
		k2 = h * f(tmp, t + h/4.);

		// k3 = h * f(x + (3./32.)*k1 + (9./32.)*k2, t + (3./8.)*h);
		tmp = x + (3./32.)*k1 + (9./32.)*k2;
		k3 = h * f(tmp, t + (3./8.)*h);

		// k4 = h * f(x + (1932./2197.)*k1 - (7200./2197.)*k2 + (7296./2197.)*k3, t + (12./13.)*h);
		tmp = x + (1932./2197.)*k1 - (7200./2197.)*k2 + (7296./2197.)*k3;
		k4 = h * f(tmp, t + (12./13.)*h);

		// k5 = h * f(x + (439./216.)*k1 - 8.*k2 + (3680./513.)*k3 - (845./4104.)*k4, t + h);
		tmp = x + (439./216.)*k1 - 8.*k2 + (3680./513.)*k3 - (845./4104.)*k4;
		k5 = h * f(tmp, t + h);

		// k6 = h * f(x - (8./27.)*k1 + 2.*k2 - (3544./2565.)*k3 + (1859./4104.)*k4 - (11./40.)*k5, t + h/2.);
		tmp = x - (8./27.)*k1 + 2.*k2 - (3544./2565.)*k3 + (1859./4104.)*k4 - (11./40.)*k5;
		k6 = h * f(tmp, t + h/2.);

		x4 = x + (25./216.)*k1 + (1408./2565.)*k3 + (2197./4101.)*k4 - k5/5.;
		x5 = x + (16./135.)*k1 + (6656./12825.)*k3 + (28561./56430.)*k4 - (9./50.)*k5 + (2./55.)*k6;

		x = x4; // x = x5; RKF54
		t += h;

		if (verbose == 1) {
			std::cout << "t: " << t << "\n";
			std::cout << x << "\n";
		}

		if (flag) break;

		err = norm_2(x5 - x4);
		scale = sqrt(sqrt(tol * h / (2. * err)));

		if (scale < 0.25) scale = 0.25;
		if (scale > 4.) scale = 4.;
		h *= scale;
		if (h < h_min) h = h_min;
		if (h > h_max) h = h_max;
	}

	init = x;
}

} // namespace kv

#endif // RKF45_HPP
