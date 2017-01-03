/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef AIRY_HPP
#define AIRY_HPP

// Airy function

#include <iostream>
#include <kv/ode-maffine.hpp>
#include <kv/gamma.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;


struct Airy_p {
	template <class T> ub::vector<T> operator() (const ub::vector<T> & x, const T& t) {
		ub::vector<T> y(2);
		y(0) = x(1);
		y(1) = t * x(0);
		return y;
	}
};

struct Airy_m {
	template <class T> ub::vector<T> operator() (const ub::vector<T> & x, const T& t) {
		ub::vector<T> y(2);
		y(0) = -x(1);
		y(1) = t * x(0);
		return y;
	}
};

/*
 * calculate Airy function using verified ODE solver
 *  mode == 0: calculate Ai
 *  mode == 1: calculate Bi
 *  d == true: calculate derivative
 */

template <class T> interval<T> airy(const interval<T>& x, int mode = 0, bool d = false)
{
	static const interval<T> g23 = gamma(interval<T>(2.) / 3.);
	static const interval<T> g13 = gamma(interval<T>(1.) / 3.);

	static const interval<T> ai0 = 1. / (pow(interval<T>(3.), interval<T>(2.) / 3.) * g23);
	static const interval<T> aid0 = -1. / (pow(interval<T>(3.), interval<T>(1.) / 3.) * g13);
	static const interval<T> bi0 = 1. / (pow(interval<T>(3.), interval<T>(1.) / 6.) * g23);
	static const interval<T> bid0 = pow(interval<T>(3.), interval<T>(1.) / 6.) / g13;

	ub::vector< interval<T> > ix(2), tx;
	interval<T> end, result;
	int r;
	bool flag;
	int suffix;

	if (mode == 0)  {
		ix(0) = ai0;
		ix(1) = aid0;
	} else {
		ix(0) = bi0;
		ix(1) = bid0;
	}

	if (d == false) {
		suffix = 0;
	} else {
		suffix = 1;
	}

	if (x.upper() == 0. && x.lower() == 0.) {
		return ix(suffix);
	}

	flag = false;

	if (x.upper() > 0) {
		tx = ix;
		end = x;
		if (end.lower() < 0.) end.lower() = 0.;
		r = odelong_maffine(Airy_p(), tx, kv::interval<T>(0.), end);
		if (r != 2) {
			throw std::domain_error("airy(): cannot calculate verified solution.");
		}
		result = tx(suffix);
		flag = true;
	}

	if (x.lower() < 0.) {
		tx = ix;
		end = x;
		if (end.upper() > 0.) end.upper() = 0.;
		end = -end;
		r = odelong_maffine(Airy_m(), tx, kv::interval<T>(0.), end);
		if (r != 2) {
			throw std::domain_error("airy(): cannot calculate verified solution.");
		}
		if (flag) {
			result = interval<T>::hull(result, tx(suffix));
		} else {
			result = tx(suffix);
		}
	}

	return result;
}

template <class T> interval<T> airy_Ai(const interval<T>& x)
{
	return airy(x, 0, false);
}

template <class T> interval<T> airy_Bi(const interval<T>& x)
{
	return airy(x, 1, false);
}

template <class T> interval<T> airy_Ai_d(const interval<T>& x)
{
	return airy(x, 0, true);
}

template <class T> interval<T> airy_Bi_d(const interval<T>& x)
{
	return airy(x, 1, true);
}

} // namespace kv

#endif // AIRY_HPP
