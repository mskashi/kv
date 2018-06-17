/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RKF45_HPP
#define RKF45_HPP

// Runge-Kutta-Fehlberg Method (RKF45)

#include <iostream>
#include <algorithm>
#include <kv/ode-param.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template <class T, class F>
void
rkf45(F f, ub::vector<T>& init, T start, T end, ode_param<T> p = ode_param<T>()) {
	ub::vector<T> x, x4, x5, dx1, dx2, dx3, dx4, dx5, dx6, tmp;
	bool flag = false;

	T t = start;

	T h = sqrt(sqrt(p.epsilon));

	T err, scale;

	x = init;

	while (true) {
		if (t + h >= end) {
			h = end - t;
			flag = true;
		}

		dx1 = h * f(x, t);

		// dx2 = h * f(x + dx1/4, t + h/4);
		tmp = x + dx1/4;
		dx2 = h * f(tmp, t + h/4);

		// dx3 = h * f(x + 3/32*dx1 + 9/32*dx2, t + 3/8*h);
		tmp = x + T(3)/32*dx1 + T(9)/32*dx2;
		dx3 = h * f(tmp, t + T(3)/8*h);

		// dx4 = h * f(x + 1932/2197*dx1 - 7200/2197*dx2 + 7296/2197*dx3, t + 12/13*h);
		tmp = x + T(1932)/2197*dx1 - T(7200)/2197*dx2 + T(7296)/2197*dx3;
		dx4 = h * f(tmp, t + T(12)/13*h);

		// dx5 = h * f(x + 439/216*dx1 - 8*dx2 + 3680/513*dx3 - 845/4104*dx4, t + h);
		tmp = x + T(439)/216*dx1 - 8*dx2 + T(3680)/513*dx3 - T(845)/4104*dx4;
		dx5 = h * f(tmp, t + h);

		// dx6 = h * f(x - 8/27*dx1 + 2*dx2 - 3544/2565*dx3 + 1859/4104*dx4 - 11/40*dx5, t + h/2);
		tmp = x - T(8)/27*dx1 + 2*dx2 - T(3544)/2565*dx3 + T(1859)/4104*dx4 - T(11)/40*dx5;
		dx6 = h * f(tmp, t + h/2);

		x4 = x + T(25)/216*dx1 + T(1408)/2565*dx3 + T(2197)/4104*dx4 - dx5/5;
		x5 = x + T(16)/135*dx1 + T(6656)/12825*dx3 + T(28561)/56430*dx4 - T(9)/50*dx5 + T(2)/55*dx6;

		x = x4; // x = x5; RKF54
		t += h;

		if (p.verbose == 1) {
			std::cout << "t: " << t << "\n";
			std::cout << x << "\n";
		}

		if (flag) break;

		err = norm_2(x5 - x4);
		scale = sqrt(sqrt(std::max(norm_2(x), 1.) * p.epsilon / err));

		if (scale < 0.25) scale = 0.25;
		if (scale > 4.) scale = 4.;
		h *= scale;
	}

	init = x;
}

} // namespace kv

#endif // RKF45_HPP
