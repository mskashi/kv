/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODESCALE_HPP
#define ODESCAPE_HPP

#include <boost/numeric/ublas/vector.hpp>

/*
 *  伊理, 藤野: 数値計算の常識, pp.95-96 の方法
 */

namespace ub = boost::numeric::ublas;

namespace kv {

template <class F> class Odescale {
	public:
	F f;
	Odescale(F f_v) : f(f_v) {}
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y, tx, ty;
		T g0;
		int s = x.size();
		int i;

		y.resize(s);
		tx.resize(s-1);

		for (i=0; i<s-1; i++) {
			tx(i) = x(i);
		}

		ty = f(tx, x(s-1));

		g0 = 1.;
		for (i=0; i<s-1; i++) {
			g0 += ty(i) * ty(i);
		}
		g0 = 1. / sqrt(g0);

		for (i=0; i<s-1; i++) {
			y(i) = ty(i) * g0;
		}
		y(s-1) = g0;

		return y;
	}
};

} // namespace kv

#endif // ODESCALE_HPP
