/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODESCALE_HPP
#define ODESCALE_HPP

#include <boost/numeric/ublas/vector.hpp>

/*
 *  from "Iri, Fujino: Suuchi Keisan no Joushiki, pp.95-96"
 */

namespace kv {

namespace ub = boost::numeric::ublas;

template <class F> class Odescale {
	public:
	F f;
	Odescale(F f) : f(f) {}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
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
