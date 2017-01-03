/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RK_HPP
#define RK_HPP

// Runge-Kutta

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace kv {

namespace ub = boost::numeric::ublas;

template <class T, class F>
void
rk(F f, ub::vector<T>& init, T start, T end) {

	ub::vector<T> x, dx1, dx2, dx3, dx4, x1, x2, x3;

	T t = start;
	T dt = end - start;

	x = init;
	dx1 = f(x, t) * dt;
	x1 = x + dx1 / 2.;
	dx2 = f(x1, t + dt/2.) * dt;
	x2 = x + dx2 / 2.;
	dx3 = f(x2, t + dt/2.) * dt;
	x3 = x + dx3;
	dx4 = f(x3, t + dt) * dt;

	init = x + dx1 / 6. + dx2 / 3. + dx3 / 3. + dx4 / 6.;
}

} // namespace kv

#endif // RK_HPP
