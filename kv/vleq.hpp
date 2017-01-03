/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef VLEQ_HPP
#define VLEQ_HPP

#include <cstddef>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <kv/matrix-inversion.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

// verified linear equation solver
// for matrix-rhs
// set r if you already have approximate inverse of a

template <class T>
bool vleq(
	const ub::matrix< interval<T> >& a,
	const ub::matrix< interval<T> >& b,
	ub::matrix<interval<T> >& x,
	const ub::matrix<T>* r = NULL )
{
	int i, j;
	int s1 = b.size1();
	int s2 = b.size2();
	ub::matrix<T> R, E;
	ub::matrix< interval<T> > EmRA;
	ub::vector< interval<T> > xtmp(s1), btmp(s1), rtmp(s1);
	bool bo;
	T norm1, norm2, norm3, err;

	if (r == NULL) {
		bo = invert(mid(a), R);
		if (bo == false) return false;
	} else {
		R = *r;
	}

	E = ub::identity_matrix<T>(s1);

	EmRA = E - prod(R, a);
	norm1 = max_norm(EmRA);
	rop<T>::begin();
	norm1 = rop<T>::sub_down(T(1.), norm1);
	rop<T>::end();

	if (norm1 <= 0.) return false;

	x = prod(R, mid(b));

	norm2 = max_norm(R);

	for (i=0; i<s2; i++) {
		for (j=0; j<s1; j++) {
			xtmp(j) = x(j, i);
			btmp(j) = b(j, i);
		}
		rtmp = prod(a, xtmp) - btmp;
		norm3 = max_norm(rtmp);
		rop<T>::begin();
		err = rop<T>::div_up(rop<T>::mul_up(norm2, norm3), norm1);
		rop<T>::end();
		for (j=0; j<s1; j++) {
			x(j, i) += interval<T>(-err, err);
		}
	}

	return true;
}

// verified linear equation solver
// for vector-rhs
// set r if you already have approximate inverse of a

template <class T>
bool vleq(
	const ub::matrix< interval<T> >& a,
	const ub::vector< interval<T> >& b,
	ub::vector<interval<T> >& x,
	const ub::matrix<T>* r = NULL )
{
	int s = b.size();
	int i;
	ub::matrix< interval<T> > b1, x1;
	bool bo;

	b1.resize(s, 1);
	x1.resize(s, 1);
	x.resize(s);

	for (i=0; i<s; i++) b1(i, 0) = b(i);

	bo = vleq(a, b1, x1, r);
	if (bo == false) return false;

	for (i=0; i<s; i++) x(i) = x1(i, 0);

	return true;
}

} // namespace kv

#endif // VLEQ_HPP
