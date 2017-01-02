/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef CARDANO_HPP
#define CARDANO_HPP

// Solve polynomial equations of degree 3/4
// by Cardano/Ferrari's methods

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/complex.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace ub = boost::numeric::ublas;

namespace kv {

template <class T> bool cardano(const ub::vector< complex< interval<T> > > &in, ub::vector< complex< interval<T> > >& out)
{
	complex< interval<T> > a, b, c, p, q, w, m, n, nt, nn, tmp;
	int i;
	bool flag;

	if (in.size() != 4) return false;
	if (zero_in(in(3).real()) && zero_in(in(3).imag())) return false;

	a = in(2) / in(3);
	b = in(1) / in(3);
	c = in(0) / in(3);

	p = b - pow(a, 2) / 3.;
	q = 2 * pow(a, 3) /27. - a * b / 3. + c;
	w = complex< interval<T> >(-1., sqrt(interval<T>(3.))) / 2.;

	tmp = sqrt(pow(q, 2)/4. + pow(p, 3) / 27.) ;
	m = pow(-q / 2. + tmp, 1. / interval<T>(3.));
	n = pow(-q / 2. - tmp, 1. / interval<T>(3.));

	nt = n;
	flag = false;

	for (i=0; i<3; i++) {
		tmp = m * nt + p/3.;
		if (zero_in(tmp.real()) && zero_in(tmp.imag())) {
			if (flag == false) {
				nn = nt;
				flag = true;
			} else {
				nn.real() = interval<T>::hull(nn.real(), nt.real());
				nn.imag() = interval<T>::hull(nn.imag(), nt.imag());
			}
		}
		nt = w * nt;
	}

	out.resize(3);
	out(0) = m + nn - a / 3.;
	out(1) = w * m + w * w * nn - a / 3.;
	out(2) = w * w * m + w * nn - a / 3.;

	return true;
}

template <class T> bool ferrari(const ub::vector< complex< interval<T> > > &in, ub::vector< complex< interval<T> > >& out) 
{
	complex< interval<T> > a, b, c, d, p, q, r, tmp;
	ub::vector< complex< interval<T> > > ci, co;
	complex< interval<T> > l, m, n, nt, nn;
	int i;
	bool flag;

	if (in.size() != 5) return false;
	if (zero_in(in(4).real()) && zero_in(in(3).imag())) return false;

	a = in(3) / in(4);
	b = in(2) / in(4);
	c = in(1) / in(4);
	d = in(0) / in(4);

	p = b - pow(a, 2) * 3. / 8.;
	q = c - b * a / 2. + pow(a, 3) / 8.;
	r = d - c * a / 4. + b * pow(a, 2) / 16. - pow(a, 4) * 3. / 256.;

	ci.resize(4);
	ci(3) = 1.;
	ci(2) = -p / 2.;
	ci(1) = -r;
	ci(0) = r * p / 2. - pow(q, 2) / 8.;

	cardano(ci, co);

	l = co(0);
	m = sqrt(2. * l - p);
	n = sqrt(pow(l, 2) - r);

	nt = n;
	flag = false;

	for (i=0; i<2; i++) {
		tmp = 2. * m * nt + q;
		if (zero_in(tmp.real()) && zero_in(tmp.imag())) {
			if (flag == false) {
				nn = nt;
				flag = true;
			} else {
				nn.real() = interval<T>::hull(nn.real(), nt.real());
				nn.imag() = interval<T>::hull(nn.imag(), nt.imag());
			}
		}
		nt = -nt;
	}

	out.resize(4);
	out(0) = (m + sqrt(pow(m, 2) - 4. * (l - nn))) / 2. - a / 4.;
	out(1) = (m - sqrt(pow(m, 2) - 4. * (l - nn))) / 2. - a / 4.;
	out(2) = (-m + sqrt(pow(m, 2) - 4. * (l + nn))) / 2. - a / 4.;
	out(3) = (-m - sqrt(pow(m, 2) - 4. * (l + nn))) / 2. - a / 4.;

	return true;
}

} // namespace kv

#endif // CARDANO_HPP
