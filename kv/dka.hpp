/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DKA_HPP
#define DKA_HPP

// Durand Kerner Aberth algorithm

#include "kv/complex.hpp"
#include "kv/psa.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/constants/constants.hpp>

namespace ub = boost::numeric::ublas;

namespace kv {

template <class T>
ub::vector< complex<T> > dka(const psa< complex<T> >& p, T tor = 1e2)
{
	int i, j;
	int s = p.v.size();
	int n = s - 1;
	ub::vector< complex<T> > a(s);
	complex<T> c;
	psa<T> b;
	T r, tmp, m1, m2;
	psa<T> db;
	ub::vector< complex<T> > x(n);
	T pi = boost::math::constants::pi<T>();
	complex<T> f, df, tmp2;
	bool flag;


	// make equation b(r) = 0

	for (i=0; i<s; i++) a(i) = p.v(i);

	c = -a(n-1) / (a(n) * n);
	for (i=1; i<=n; i++) {
		for (j = n; j >= i; j--) {
			a(j-1) += a(j) * c;
		}
	}

	b.v.resize(s);
	for (i=0; i<n-1; i++) b.v(i) = -abs(a(i));
	b.v(n-1) = 0.;
	b.v(n) = abs(a(n));

	// calculate initial radius for Newton's method
	// r = 2^n such that g(2^n) > 0 and g(2^(n-1)) < 0

	r = 1.;
	while (true) {
		if (eval(b, r) > 0.) break;
		r *= 2.;
	}

	// db(r) = b'(r)

	db.v.resize(n);
	for (i=1; i<s; i++) db.v(i-1) = b.v(i) * i;

	// calculate radius by Newton's method for equation b(r) = 0

	while (true) {
		tmp = eval(b, r) / eval(db, r);
		r -= tmp;
		using std::abs;
		if (abs(tmp) < abs(r) * std::numeric_limits<T>::epsilon() * tor) break;
	}

	// set Aberth's initial values
	for (i=0; i<n; i++) {
		x(i) = c + r * exp((2. * i * pi / n + pi / (2. * n)) * complex<T>::i());
	}

	// Durand Kerner algorithm

	while (true) {
		m1 = 0.;
		m2 = 0.;
		flag = false;
		for (i=0; i<n; i++) {
			f = eval(p, x(i));
			df = p.v(n);
			for (j=0; j<n; j++) {
				if (j == i) continue;
				df *= x(i) - x(j);
			}
			tmp2 = f / df;
			if (tmp2.real() != tmp2.real() || tmp2.imag() != tmp2.imag()) { // isnan
				flag = true; break;
			}
			if (abs(tmp2) > m1) m1 = abs(tmp2);
			if (abs(x(i)) > m2) m2 = abs(x(i));
			x(i) -= tmp2;
		}
		if (flag) break;
		if (m1 < m2 * std::numeric_limits<T>::epsilon() * tor) break;
	}

	return x;
}

} // namespace kv

#endif // DKA_HPP
