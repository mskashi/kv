/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef BESSEL_HPP
#define BESSEL_HPP

#include <kv/defint.hpp>
#include <kv/defint-singular.hpp>
#include <kv/gamma.hpp>

namespace kv {


// Bessel function J_n(x) for integer n


template <class TT> struct Besselj_i {
	int n;
	TT x;

	Besselj_i(int n, TT x) : n(n), x(x) {
	}

	template <class T> T operator() (const T& s) {
		return cos(x * sin(s) - n * s);
	}
};

template <class T> interval<T> besselj(int n, interval<T> x) {
	static const interval<T> p = kv::constants< interval<T> >::pi();

	return defint_autostep(Besselj_i< interval<T> >(n, x), interval<T>(0.), p, 12) / p;
}


// Bessel function J_nu(x) for non-integer nu
//  nu > -0.5 is required.


template <class TT> struct Besselj_ni {
	TT nu;
	TT x;

	Besselj_ni(TT nu, TT x) : nu(nu), x(x) {
	}

	template <class T> T operator() (const T& s) {
		return cos(x * cos(s)) * pow(sin(s), 2. * nu);
	}
};

struct Besselj_ni_f {
	template <class T> T operator() (const T& s) {
		return sin(s);
	}
};


template <class TT> struct Besselj_ni_g {
	TT x;

	Besselj_ni_g(TT x) : x(x) {
	}

	template <class T> T operator() (const T& s) {
		return cos(x * cos(s));
	}
};

template <class T> interval<T> besselj(interval<T> nu, interval<T> x) {
	static const interval<T> p = constants< interval<T> >::pi();
	interval<T> tmp;

	tmp =  defint_power3_autostep(Besselj_ni< interval<T> >(nu, x), Besselj_ni_f(), Besselj_ni_g< interval<T> >(x), interval<T>(0.), p/2, 12, 2 * nu);
	tmp += defint_power3_autostep_r(Besselj_ni< interval<T> >(nu, x), Besselj_ni_f(), Besselj_ni_g< interval<T> >(x), p/2, p, 12, 2 * nu);

	return tmp * pow(0.5 * x, nu) / sqrt(p) / gamma(nu + 0.5);
}

} // namespace kv

#endif // BESSEL_HPP
