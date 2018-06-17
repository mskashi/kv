/*
 * Copyright (c) 2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RKF78_HPP
#define RKF78_HPP

// Runge-Kutta-Fehlberg Method (RKF78)

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
rkf78(F f, ub::vector<T>& init, T start, T end, ode_param<T> p = ode_param<T>()) {
	ub::vector<T> x, x7, x8, dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10, dx11, dx12, dx13, tmp;
	bool flag = false;

	T t = start;

	T h = pow(p.epsilon, 1/7.);

	T err, scale;

	x = init;

	while (true) {
		if (t + h >= end) {
			h = end - t;
			flag = true;
		}

		dx1 = h * f(x, t);

		// dx2 = h * f(x + 2/27*dx1, t + 2/27*h);
		tmp = x + T(2)/27*dx1;
		dx2 = h * f(tmp, t + T(2)/27*h);

		// dx3 = h * f(x + dx1/36 + dx2/12, t + h/9);
		tmp = x + dx1/36 + dx2/12;
		dx3 = h * f(tmp, t + h/9);

		// dx4 = h * f(x + dx1/24 + dx3/8, t + h/6);
		tmp = x + dx1/24 + dx3/8;
		dx4 = h * f(tmp, t + h/6);

		// dx5 = h * f(x + 5/12*dx1 - 25/16*dx3 + 25/16*dx4, t + 5/12*h);
		tmp = x + T(5)/12*dx1 - T(25)/16*dx3 + T(25)/16*dx4;
		dx5 = h * f(tmp, t + T(5)/12*h);

		// dx6 = h * f(x + dx1/20 + dx4/4 + dx5/5, t + h/2);
		tmp = x + dx1/20 + dx4/4 + dx5/5;
		dx6 = h * f(tmp, t + h/2);

		// dx7 = h * f(x - 25/108*dx1 + 125/108*dx4 - 65/27*dx5 + 125/54*dx6, t + 5/6*h);
		tmp = x - T(25)/108*dx1 + T(125)/108*dx4 - T(65)/27*dx5 + T(125)/54*dx6;
		dx7 = h * f(tmp, t + T(5)/6*h);

		// dx8 = h * f(x + 31/300*dx1 + 61/225*dx5 - 2/9*dx6 + 13/900*dx7, t + h/6);
		tmp = x + T(31)/300*dx1 + T(61)/225*dx5 - T(2)/9*dx6 + T(13)/900*dx7;
		dx8 = h * f(tmp, t + h/6);

		// dx9 = h * f(x + 2*dx1 - 53/6*dx4 + 704/45*dx5 - 107/9*dx6 + 67/90*dx7 + 3*dx8, t + 2/3*h);
		tmp = x + 2*dx1 - T(53)/6*dx4 + T(704)/45*dx5 - T(107)/9*dx6 + T(67)/90*dx7 + 3*dx8;
		dx9 = h * f(tmp, t + T(2)/3*h);

		// dx10 = h * f(x - 91/108*dx1 + 23/108*dx4 - 976/135*dx5 + 311/54*dx6 - 19/60*dx7 + 17/6*dx8 - dx9/12, t + h/3);
		tmp = x - T(91)/108*dx1 + T(23)/108*dx4 - T(976)/135*dx5 + T(311)/54*dx6 - T(19)/60*dx7 + T(17)/6*dx8 - dx9/12;
		dx10 = h * f(tmp, t + h/3);

		// dx11 = h * f(x + 2383/4100*dx1 - 341/164*dx4 + 4496/1025*dx5 - 301/82*dx6 + 2133/4100*dx7 + 45/82*dx8 + 45/164*dx9 + 18/41*dx10, t + h);
		tmp = x + T(2383)/4100*dx1 - T(341)/164*dx4 + T(4496)/1025*dx5 - T(301)/82*dx6 + T(2133)/4100*dx7 + T(45)/82*dx8 + T(45)/164*dx9 + T(18)/41*dx10;
		dx11 = h * f(tmp, t + h);

		// dx12 = h * f(x + 3/205*dx1 - 6/41*dx6 - 3/205*dx7 - 3/41*dx8 + 3/41*dx9 + 6/41*dx10, t);
		tmp = x + T(3)/205*dx1 - T(6)/41*dx6 - T(3)/205*dx7 - T(3)/41*dx8 + T(3)/41*dx9 + T(6)/41*dx10;
		dx12 = h * f(tmp, t);

		/*
		  Coefficient of dx10 is 19/41 on Hairer's text "Solving
		  Ordinary Differential Equations I". 12/41 is correct.
 		*/
		// dx13 = h * f(x - 1777/4100*dx1 - 341/164*dx4 + 4496/1025*dx5 - 289/82*dx6 + 2193/4100*dx7 + 51/82*dx8 + 33/164*dx9 + 12/41*dx10 + dx12, t + h);
		tmp = x - T(1777)/4100*dx1 - T(341)/164*dx4 + T(4496)/1025*dx5 - T(289)/82*dx6 + T(2193)/4100*dx7 + T(51)/82*dx8 + T(33)/164*dx9 + T(12)/41*dx10 + dx12;
		dx13 = h * f(tmp, t + h);


		x7 = x + T(41)/840*dx1 + T(34)/105*dx6 + T(9)/35*dx7 + T(9)/35*dx8 + T(9)/280*dx9 + T(9)/280*dx10 + T(41)/840*dx11;
		x8 = x + T(34)/105*dx6 + T(9)/35*dx7 + T(9)/35*dx8 + T(9)/280*dx9 + T(9)/280*dx10 + T(41)/840*dx12 + T(41)/840*dx13;

		x = x7; // x = x8; RKF87
		t += h;

		if (p.verbose == 1) {
			std::cout << "t: " << t << "\n";
			std::cout << x << "\n";
		}

		if (flag) break;

		err = norm_2(x8 - x7);
		scale = pow(std::max(norm_2(x), 1.) * p.epsilon / err, 1/7.);

		if (scale < 0.25) scale = 0.25;
		if (scale > 4.) scale = 4.;
		h *= scale;
	}

	init = x;
}

} // namespace kv

#endif // RKF78_HPP
