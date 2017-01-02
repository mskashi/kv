/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef INTERVAL_VECTOR_HPP
#define INTERVAL_VECTOR_HPP

// utilities for interval vector/matrix

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/interval.hpp>

namespace ub = boost::numeric::ublas;

namespace kv {

template <class T> inline ub::vector<T> mid (const ub::vector< interval<T> >& I) {
	int i;
	int s = I.size();
	ub::vector<T> r(s);

	for (i=0; i<s; i++) {
		r(i) = mid(I(i));
	}

	return r;
}

template <class T> inline ub::matrix<T> mid (const ub::matrix< interval<T> >& I) {
	int i, j;
	int s1 = I.size1();
	int s2 = I.size2();
	ub::matrix<T> r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i,j) = mid(I(i,j));
		}
	}

	return r;
}


template <class T> inline ub::vector<T> rad (const ub::vector< interval<T> >& I) {
	int i;
	int s = I.size();
	ub::vector<T> r(s);

	for (i=0; i<s; i++) {
		r(i) = rad(I(i));
	}

	return r;
}

template <class T> inline ub::matrix<T> rad (const ub::matrix< interval<T> >& I) {
	int i, j;
	int s1 = I.size1();
	int s2 = I.size2();
	ub::matrix<T> r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i,j) = rad(I(i,j));
		}
	}

	return r;
}


template <class T> inline ub::vector<T> mag (const ub::vector< interval<T> >& I) {
	int i;
	int s = I.size();
	ub::vector<T> r(s);

	for (i=0; i<s; i++) {
		r(i) = mag(I(i));
	}

	return r;
}

template <class T> inline ub::matrix<T> mag (const ub::matrix< interval<T> >& I) {
	int i, j;
	int s1 = I.size1();
	int s2 = I.size2();
	ub::matrix<T> r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i,j) = mag(I(i,j));
		}
	}

	return r;
}

template <class T> inline ub::vector<T> mag (const ub::vector< T >& v) {
	int i;
	int s = v.size();
	ub::vector<T> r(s);

	for (i=0; i<s; i++) {
		r(i) = mag(v(i));
	}

	return r;
}

template <class T> inline ub::matrix<T> mag (const ub::matrix< T >& m) {
	int i, j;
	int s1 = m.size1();
	int s2 = m.size2();
	ub::matrix<T> r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i,j) = mag(m(i,j));
		}
	}

	return r;
}



template <class T> inline ub::vector<T> mig (const ub::vector< interval<T> >& I) {
	int i;
	int s = I.size();
	ub::vector<T> r(s);

	for (i=0; i<s; i++) {
		r(i) = mig(I(i));
	}

	return r;
}

template <class T> inline ub::matrix<T> mig (const ub::matrix< interval<T> >& I) {
	int i, j;
	int s1 = I.size1();
	int s2 = I.size2();
	ub::matrix<T> r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i,j) = mig(I(i,j));
		}
	}

	return r;
}


template <class T> inline bool zero_in (const ub::vector< interval<T> >& I) {
	int i;
	int s = I.size();

	for (i=0; i<s; i++) {
		if (!zero_in(I(i))) return false;
	}

	return true;
}

template <class T> inline bool subset (const ub::vector< interval<T> >& I, const ub::vector< interval<T> > & J) {
	int i;
	int s = I.size();

	for (i=0; i<s; i++) {
		if (!subset(I(i), J(i))) return false;
	}

	return true;
}

template <class T> inline bool proper_subset (const ub::vector< interval<T> >& I, const ub::vector< interval<T> > & J) {
	int i;
	int s = I.size();

	for (i=0; i<s; i++) {
		if (!proper_subset(I(i), J(i))) return false;
	}

	return true;
}

template <class T> inline bool overlap (const ub::vector< interval<T> >& I, const ub::vector< interval<T> > & J) {
	int i;
	int s = I.size();

	for (i=0; i<s; i++) {
		if (!overlap(I(i), J(i))) return false;
	}

	return true;
}

template <class T> inline ub::vector< interval<T> > intersect (const ub::vector< interval<T> >& I, const ub::vector< interval<T> > & J) {
	int i;
	int s = I.size();
	ub::vector< interval<T> > r(s);

	for (i=0; i<s; i++) {
		r(i) = intersect(I(i), J(i));
	}

	return r;
}

template <class T> inline ub::matrix< interval<T> > intersect (const ub::matrix< interval<T> >& I, const ub::matrix< interval<T> > & J) {
	int i, j;
	int s1 = I.size1();
	int s2 = I.size2();
	ub::matrix< interval<T> > r(s1, s2);

	for (i=0; i<s1; i++) {
		for (j=0; j<s2; j++) {
			r(i, j) = intersect(I(i, j), J(i, j));
		}
	}

	return r;
}

template <class T> inline T max_norm (const ub::vector<T>& x) {
	int i;
	int s = x.size();
	T tmp;
	T r = 0.;

	for (i=0; i<s; i++) {
		// tmp = abs(x(i));
		tmp = (x(i) >= 0.) ? x(i) : -x(i);
		if (tmp > r) r = tmp;
	}

	return r;
}

template <class T> inline T max_norm (const ub::matrix<T>& x) {
	int i, j;
	int s1 = x.size1();
	int s2 = x.size2();
	T tmp, tmp2;
	T r = 0;

	for (i=0; i<s1; i++) {
		tmp = 0.;
		rop<T>::begin();
		for (j=0; j<s2; j++) {
			tmp2 = (x(i, j) >= 0.) ? x(i, j) : -x(i, j);
			tmp = rop<T>::add_up(tmp, tmp2);
		}
		rop<T>::end();
		if (tmp > r) r = tmp;
	}

	return r;
}

template <class T> inline T max_norm (const ub::vector< interval<T> >& x) {
	int i;
	int s = x.size();
	T tmp;
	T r = 0.;

	for (i=0; i<s; i++) {
		tmp = norm(x(i));
		if (tmp > r) r = tmp;
	}

	return r;
}

template <class T> inline T max_norm (const ub::matrix< interval<T> >& x) {
	int i, j;
	int s1 = x.size1();
	int s2 = x.size2();
	T tmp;
	T r = 0;

	for (i=0; i<s1; i++) {
		tmp = 0.;
		rop<T>::begin();
		for (j=0; j<s2; j++) {
			tmp = rop<T>::add_up(tmp, norm(x(i, j)));
		}
		rop<T>::end();
		if (tmp > r) r = tmp;
	}

	return r;
}

} // namespace kv

#endif // INTERVAL_VECTOR_HPP
