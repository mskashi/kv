/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef JOINTRANGE_HPP
#define JOINTRANGE_HPP

#include <boost/numeric/ublas/vector.hpp>

#include <kv/affine.hpp>
#include <kv/gnuplot.hpp>

namespace ub = boost::numeric::ublas;


namespace kv {

template<class T>
void jointrange(const affine<T>& x, const affine<T>& y, const gnuplot& g, int t=1)
{
	int n, i, j, s;
	T tmp;
	T dx, dy;
	ub::vector<T> vx, vy;
	ub::matrix<int> signcache;

	n = std::max(x.a.size(), y.a.size()) - 1;

#if AFFINE_SIMPLE >= 1
	vx.resize(n+3);
	vy.resize(n+3);
#else
	vx.resize(n+1);
	vy.resize(n+1);
#endif

	for (i=0; i<=n; i++) {
		vx(i) = x.a(i);
		vy(i) = y.a(i);
	}

#if AFFINE_SIMPLE >= 1
	vx(n+1) = x.er;
	vy(n+1) = 0.;
	vx(n+2) = 0.;
	vy(n+2) = y.er;
	n += 2;
#endif

	signcache.resize(n+1, n+1);
	for (i=1; i<=n; i++) {
		for (j=i+1; j<=n; j++) {
			tmp = vy(i) * vx(j) - vx(i) * vy(j);
			if (tmp >= 0.) s = 1;
			else s = -1;
			signcache(i, j) = s;
			signcache(j, i) = -s;
		}
	}

	for (i=1; i<=n; i++) {
		dx = 0.;
		dy = 0.;
		for (j=1; j<=n; j++) {
			if (j == i) continue;
			// tmp = vy(i) * vx(j) - vx(i) * vy(j);
			// if (tmp >= 0.) s = 1;
			// else s = -1;
			s = signcache(i, j);
			dx += s * vx(j);
			dy += s * vy(j);
		}
		g.line(vx(0)+dx+vx(i), vy(0)+dy+vy(i), 
		       vx(0)+dx-vx(i), vy(0)+dy-vy(i), t);
		g.line(vx(0)-dx+vx(i), vy(0)-dy+vy(i), 
		       vx(0)-dx-vx(i), vy(0)-dy-vy(i), t);
	}
}

} // namespace kv

#endif // JOINTRANGE_HPP
