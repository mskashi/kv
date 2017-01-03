/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef BESSEL_HPP
#define BESSEL_HPP

#include <kv/defint.hpp>

namespace kv {

template <class TT> struct Bessel_f {
	int n;
	TT x;

	Bessel_f(int n, TT x) : n(n), x(x) {
	}

	template <class T> T operator() (const T& s) {
		return cos(x * sin(s) - n * s);
	}
};

template <class T> interval<T> Bessel(int n, interval<T> x) {
	static const interval<T> p = kv::constants< interval<T> >::pi();

	return defint_autostep(kv::Bessel_f< interval<T> >(n, x), interval<T>(0.), p, 12) / p;
}

} // namespace kv

#endif // BESSEL_HPP
