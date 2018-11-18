/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DKA_HPP
#define DKA_HPP

#include <limits>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <kv/constants.hpp>

namespace kv {

namespace ub = boost::numeric::ublas;

template <class T>
static T inline eval_polynomial (const ub::vector<T>& p, const T& x)
{
	int i;
	T r;
	int s = p.size();

	r = p(s-1);
	for (i=s-2; i>=0; i--) {
		r = r * x + p(i);
	}

	return r;
}

// Durand Kerner Aberth algorithm

template <class T>
bool dka(const ub::vector< kv::complex<T> >& p, ub::vector< kv::complex<T> >& x, T epsilon = std::numeric_limits<T>::epsilon())
{
	int i, j;
	int s = p.size();
	int n = s - 1;
	ub::vector< kv::complex<T> > a(s);
	kv::complex<T> c;
	ub::vector<T> b;
	T r, tmp, norm1, norm2, norm3;
	ub::vector<T> db;
	ub::vector< kv::complex<T> > dx(n);
	T pi = kv::constants<T>::pi();
	kv::complex<T> f, df;

	using std::abs;
	using std::pow;

	if (p(n).real() * p(n).real() + p(n).imag() * p(n).imag() == 0.) {
		return false;
	}

	// make equation b(r) = 0

	for (i=0; i<s; i++) a(i) = p(i);

	c = -a(n-1) / (a(n) * n);
	for (i=1; i<=n; i++) {
		for (j = n; j >= i; j--) {
			a(j-1) += a(j) * c;
		}
	}

	b.resize(s);
	for (i=0; i<n-1; i++) b(i) = -abs(a(i));
	b(n-1) = 0.;
	b(n) = abs(a(n));
	// prepare for the case where b(0)...b(n-2) are 0.
	b(0) = std::min(b(0), -std::numeric_limits<T>::epsilon());

	// initial guess for b(r) = 0
	r = 0.;
	for (i=0; i<n-1; i++) {
		tmp = pow(n * abs(b[i]/b[n]), 1.0/(n-i));
		if (tmp > r) r = tmp;
	}

	// db(r) = b'(r)
	db.resize(n);
	for (i=1; i<s; i++) db(i-1) = b(i) * i;

	// calculate radius by Newton's method for equation b(r) = 0
	while (true) {
		tmp = eval_polynomial(b, r) / eval_polynomial(db, r);
		r -= tmp;
		if (abs(tmp) < n * abs(r) * epsilon) break;
	}

	// set Aberth's initial values

	x.resize(n);
	for (i=0; i<n; i++) {
		x(i) = c + r * exp((2. * i * pi / n + pi / (2. * n)) * kv::complex<T>::i());
	}

	// max(|pi|)
	norm3 = 0.;
	for (i=0; i<s; i++) {
		norm3 = std::max(norm3, abs(p(i)));
	}

	// Durand Kerner algorithm

	while (true) {
		for (i=0; i<n; i++) {
			f = eval_polynomial(p, x(i));
			df = p(n);
			for (j=0; j<n; j++) {
				if (j == i) continue;
				df *= x(i) - x(j);
			}
			if (df.real() * df.real() + df.imag() * df.imag() == 0.) {
				dx(i) = 0.;
			} else {
				dx(i) = f / df;
			}
		}

		x -= dx;

		norm1 = 1.;
		norm2 = 0.;
		for (i=0; i<n; i++) {
			norm1 = std::max(norm1, abs(x(i)));
			norm2 = std::max(norm2, abs(dx(i)));
		}
		if (norm2 <= n * norm1 * epsilon) break;
		if (abs(f) <= n * norm3 * epsilon) break;
	}

	return true;
}

// error estimation using Smith's theorem

template <class T>
ub::vector< kv::complex< interval<T> > > smith_error(const ub::vector< kv::complex< interval<T> > >& p, const ub::vector< kv::complex<T> >& x)
{
	int i, j;
	int s = p.size();
	int n = s - 1;
	ub::vector< kv::complex< interval<T> > > x2, x3(n);
	kv::complex< interval<T> > f, df;
	T err;
	ub::vector<int> unified(n);
	bool flag;

	x2 = x;

	for (i=0; i<n; i++) {
		f = eval_polynomial(p, x2(i));
		df = p(n);
		for (j=0; j<n; j++) {
			if (j == i) continue;
			df *= x2(i) - x2(j);
		}
		if (zero_in(df.real() * df.real() + df.imag() * df.imag())) {
			err = std::numeric_limits<T>::infinity();
		} else {
			err = abs(n * f / df).upper();
		}
		x3(i).real() = x2(i).real() + err * interval<T>(-1., 1.);
		x3(i).imag() = x2(i).imag() + err * interval<T>(-1., 1.);
	}

	/* unify overlapping solutions */

	for (i=0; i<n; i++) unified(i) = -1; /* not unified */

	while (true) {
		flag = true;
		for (i=0; i<n-1; i++) {
			if (unified(i) != -1) continue;
			for (j=i; j<n; j++) {
				if (unified(j) != -1) continue;
				if (overlap(x3(i).real(), x3(j).real()) && overlap(x3(i).imag(), x3(j).imag())) {
					x3(i).real() = interval<T>::hull(x3(i).real(), x3(j).real());
					x3(i).imag() = interval<T>::hull(x3(i).imag(), x3(j).imag());
					unified(j) = i; /* x3(j) is unified to x3(i) */
					flag = false;
				}
			}
		}
		if (flag) break;
	}

	for (i=0; i<n; i++) {
		if (unified(i) != -1) x3(i) = x3(unified(i));
	}

	return x3;
}

// verified Durand Kerner Aberth

template <class T>
bool vdka(const ub::vector< kv::complex< interval<T> > >& p, ub::vector< kv::complex< interval<T> > >& result, T epsilon = std::numeric_limits<T>::epsilon())
{
	int i;
	int s = p.size();
	int n = s - 1;
	ub::vector< kv::complex< T > > p2(s), x2;

	if (zero_in(p(n).real() * p(n).real() + p(n).imag() * p(n).imag())) {
		return false;
	}

	for (i=0; i<s; i++) {
		p2(i).real() = mid(p(i).real());
		p2(i).imag() = mid(p(i).imag());
	}

	dka(p2, x2, epsilon);

	result = smith_error(p, x2);
	return true;
}

} // namespace kv

#endif // DKA_HPP
