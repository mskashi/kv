/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef QR_HPP
#define QR_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

// QR decomposition using Gram-Schmit orthogonalization

namespace ub = boost::numeric::ublas;


namespace kv {


template <class T> bool qr(const ub::matrix<T>& in, ub::matrix<T>& q, ub::matrix<T>& r)
{
	int i, j, k;
	int n;
	ub::vector<T> v;
	T tmp;

	n = in.size1();
	if (n != in.size2()) return false;

	q.resize(n, n);
	r.resize(n, n);
	v.resize(n);

	for (i=0; i<n; i++) {
		for (k=0; k<n; k++) v(k) = in(k, i);
		for (j=0; j<i; j++) {
			tmp = 0.;
			for (k=0; k<n; k++) {
				tmp += in(k, i) * q(k, j);
			}
			for (k=0; k<n; k++) {
				v(k) -= tmp * q(k, j);
			}
			r(j, i) = tmp;
		}
		tmp = norm_2(v);
		if (tmp == 0.) return false;
		for (k=0; k<n; k++) q(k, i) = v(k) / tmp;
		r(i, i) = tmp;
		for (j=i+1; j<n; j++) r(j, i) = 0.;
	}

	return true;
}

} // namespace kv

#endif // QR_HPP
