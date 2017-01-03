/*
 * Copyright (c) 2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

/*
 * calculate Gaussian hypergeometric function 2F1(a,b,c;z)
 */

#ifndef HYPERGEOM_HPP
#define HYPERGEOM_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/geoseries.hpp>
#include <algorithm>
#include <limits>

#if 0
#include <cmath>
#include <kv/defint-singular.hpp>
#include <kv/beta.hpp>
#endif

namespace kv {

#if 0
/*
 * calculate Gaussian hypergeometric function 2F1(a,b,c;z)
 * using Euler type integral formula
 *  c > b > 0, |z| < 1 is required.
 */

template <class TT> struct HyperGeom {
	TT a, b, c, z;
	HyperGeom(TT a, TT b, TT c, TT z) : a(a), b(b), c(c), z(z) {}

	template <class T> T operator() (const T& t) {
		return pow(t, T(b)-1) * pow(1-t, T(c)-T(b)-1) * pow(1-T(z)*t, -T(a));
	}
};

struct HyperGeom_s1 {
	HyperGeom_s1() {}

	template <class T> T operator() (const T& t) {
		return t;
	}
};

template <class TT> struct HyperGeom_s2 {
	TT a, b, c, z;
	HyperGeom_s2(TT a, TT b, TT c, TT z) : a(a), b(b), c(c), z(z) {}

	template <class T> T operator() (const T& t) {
		return pow(1-t, T(c)-T(b)-1) * pow(1-T(z)*t, -T(a));
	}
};

struct HyperGeom_s3 {
	HyperGeom_s3() {}

	template <class T> T operator() (const T& t) {
		return 1-t;
	}
};

template <class TT> struct HyperGeom_s4 {
	TT a, b, z;
	HyperGeom_s4(TT a, TT b, TT z) : a(a), b(b), z(z) {}

	template <class T> T operator() (const T& t) {
		return pow(t, T(b)-1) * pow(1-T(z)*t, -T(a));
	}
};

#if !defined(HYPERGEOM_ORDER)
#define HYPERGEOM_ORDER 18
#endif

template <class T>
interval<T> hypergeom(const interval<T>& a, const interval<T>& b, const interval<T>& c, const interval<T>& z) {
	interval<T> tmp;

	tmp = defint_power3_autostep(HyperGeom< interval<T> >(a, b, c, z), HyperGeom_s1(), HyperGeom_s2< interval<T> >(a, b, c, z), interval<T>(0.), interval<T>(0.5), HYPERGEOM_ORDER, interval<T>(b - 1));
	tmp += defint_power3_autostep_r(HyperGeom< interval<T> >(a, b, c, z), HyperGeom_s3(), HyperGeom_s4< interval<T> >(a, b, z), interval<T>(0.5), interval<T>(1.), HYPERGEOM_ORDER, interval<T>(c - b - 1));

	return tmp / beta(b, c-b);
}

#endif

template <class T>
interval<T> hypergeom(const interval<T>& a, const interval<T>& b, const interval<T>& c, const interval<T>& z) {
	interval<T> result, ta, tb, tc, td, tr, ratio, residual, tmp;
	T offset;

	result = 1.;
	tr = 1.;
	ta = a;
	tb = b;
	tc = c;
	td = 1;

	while (true) {
		tr *= ta;
		tr *= tb;
		tr /= tc;
		tr /= td;
		tr *= z;
		offset = std::min(std::min(std::min(ta.lower(), tb.lower()), tc.lower()), td.lower());
		if (offset > 0) {
			tmp = 1 / interval<T>(offset, std::numeric_limits<T>::infinity());
			ratio = z * (1 + (ta - offset) * tmp) * (1 + (tb - offset) * tmp)
				/ (1 + (tc - offset) * tmp) / (1 + (td - offset) * tmp);
			if (mag(ratio) < 1) {
				residual = tr * geoseries(ratio);
				if (width(residual) < mag(result) * std::numeric_limits<T>::epsilon()) {
					result += residual;
					break;
				}
			}
		}
		result += tr;
		ta += 1;
		tb += 1;
		tc += 1;
		td += 1;
	}

	return result;
}

} // namespace kv

#endif // HYPERGEOM_HPP
