/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef LP_HPP
#define LP_HPP

#include <stdexcept>
#include <list>
#include <limits>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template <class T> inline T lp_minimize(ub::vector<T>& objfunc, std::list< ub::vector<T> >& constraints)
{
	int i, j;
	int m = constraints.size();
	int lp_n = objfunc.size() - 1;
	int n = lp_n + m;
	int pivot_i, pivot_j;
	T tmp, minmax;
	typename std::list< ub::vector<T> >::iterator p;

	ub::matrix<T> a(m+1, n+1);

	a(0, 0) = objfunc(0);

	for (i=1; i<n+1; i++) {
		if (i <= lp_n) a(0, i) = -objfunc(i);
		else a(0, i) = 0.;
	}

	p = constraints.begin();
	for (i=0; i<m; i++) {
		if ((*p)[0] > 0) {
			throw std::domain_error("lp_minimize: constraints sign error");
		}
		a(i+1, 0) = -(*p)[0];
		for (j=1; j<=lp_n; j++) {
			if (j < (*p).size()) a(i+1, j) = (*p)[j];
			else a(i+1, j) = 0.;
		}
		for (j=0; j<m; j++) {
			a(i+1, lp_n+j+1) = (i == j) ? 1 : 0;
		}
		p++;
	}

	ub::vector<bool> isbasic(n+1);
	for (i=0; i<lp_n+1; i++) isbasic(i) = false;
	for (i=0; i<m; i++) isbasic(lp_n+i+1) = true;

	ub::vector<int> basic(m+1);
	for (i=0; i<m; i++) basic(i+1) = lp_n+i+1;
	
	while (true) {
		pivot_j = -1;
		minmax = 0.;
		for (i=1; i<n+1; i++) {
			if (isbasic(i)) continue;
			if (a(0, i) > minmax) {
				pivot_j = i;
				minmax = a(0, i);
			}
		}
		if (pivot_j == -1) break;

		pivot_i = -1;
		minmax = std::numeric_limits<T>::max();
		for (i=1; i<m+1; i++) {
			if (a(i, pivot_j) <= 0.) continue;
			tmp = a(i, 0) / a(i, pivot_j);
			if (tmp < minmax) {
				minmax = tmp;
				pivot_i = i;
			}
		}
		if (pivot_i == -1) {
			throw std::domain_error("lp_minimize: no optimal solution");
		}

		isbasic(basic(pivot_i)) = false;
		isbasic(pivot_j) = true;
		basic(pivot_i) = pivot_j;

		tmp = a(pivot_i, pivot_j);
		for (i=0; i<n+1; i++) {
			if (isbasic(i)) {
				if (i == pivot_j) a(pivot_i, i) = 1.;
				continue;
			}
			a(pivot_i, i) /= tmp;
		}

		for (i=0; i<m+1; i++) {
			if (i == pivot_i) continue;
			tmp = a(i, pivot_j);
			for (j=0; j<n+1; j++) {
				if (isbasic(j)) {
					if (j == pivot_j) a(i, j) = 0.;
					continue;
				}
				a(i, j) -= a(pivot_i, j) * tmp;
			}
		}
	}

	return a(0, 0);
}

template <class T> inline T lp_minimize_verified(ub::vector<T>& objfunc, std::list< ub::vector<T> >& constraints, int round)
{
	int i, j;
	int m = constraints.size();
	int lp_n = objfunc.size() - 1;
	int n = lp_n + m;
	int pivot_i, pivot_j;
	T tmp, minmax;
	typename std::list< ub::vector<T> >::iterator p;
	interval<T> Itmp, Imin;

	ub::matrix<T> a(m+1, n+1);
	ub::matrix< interval<T> > Ia(m+1, n+1);

	a(0, 0) = objfunc(0);

	for (i=1; i<n+1; i++) {
		if (i <= lp_n) a(0, i) = -objfunc(i);
		else a(0, i) = 0.;
	}

	p = constraints.begin();
	for (i=0; i<m; i++) {
		if ((*p)[0] > 0) {
			throw std::domain_error("lp_minimize: constraints sign error");
		}
		a(i+1, 0) = -(*p)[0];
		for (j=1; j<=lp_n; j++) {
			if (j < (*p).size()) a(i+1, j) = (*p)[j];
			else a(i+1, j) = 0.;
		}
		for (j=0; j<m; j++) {
			a(i+1, lp_n+j+1) = (i == j) ? 1 : 0;
		}
		p++;
	}

	ub::vector<bool> isbasic(n+1);
	for (i=0; i<lp_n+1; i++) isbasic(i) = false;
	for (i=0; i<m; i++) isbasic(lp_n+i+1) = true;

	ub::vector<int> basic(m+1);
	for (i=0; i<m; i++) basic(i+1) = lp_n+i+1;
	
	while (true) {
		pivot_j = -1;
		minmax = 0.;
		for (i=1; i<n+1; i++) {
			if (isbasic(i)) continue;
			if (a(0, i) > minmax) {
				pivot_j = i;
				minmax = a(0, i);
			}
		}
		if (pivot_j == -1) break;

		pivot_i = -1;
		minmax = std::numeric_limits<T>::max();
		for (i=1; i<m+1; i++) {
			if (a(i, pivot_j) <= 0.) continue;
			tmp = a(i, 0) / a(i, pivot_j);
			if (tmp < minmax) {
				minmax = tmp;
				pivot_i = i;
			}
		}
		if (pivot_i == -1) {
			throw std::domain_error("lp_minimize: no optimal solution");
		}

		Imin = (interval<T>)a(pivot_i, 0) / a(pivot_i, pivot_j);

		for (i=1; i<m+1; i++) {
			if (a(i, pivot_j) <= 0.) continue;
			if (i == pivot_i) {continue;}
			if (round == -1) {
				while (true) {
					Itmp = (interval<T>)a(i, 0) / a(i, pivot_j);
					if (Imin.upper() < Itmp.lower()) {
						break;
					}
					// succ
					a(i, 0) = ((interval<T>)a(i, 0) + std::numeric_limits<T>::denorm_min()).upper();
				}
			} else {
				Itmp = (interval<T>)a(i, 0) / a(i, pivot_j);
				while (true) {
					if (Imin.upper() < Itmp.lower()) {
						break;
					}
					// pred
					a(pivot_i, 0) = ((interval<T>)a(pivot_i, 0) - std::numeric_limits<T>::denorm_min()).lower();

					Imin = (interval<T>)a(pivot_i, 0) / a(pivot_i, pivot_j);
				}
			}
		}

		isbasic(basic(pivot_i)) = false;
		isbasic(pivot_j) = true;
		basic(pivot_i) = pivot_j;

		Itmp = a(pivot_i, pivot_j);
		for (i=0; i<n+1; i++) {
			if (isbasic(i)) {
				if (i == pivot_j) a(pivot_i, i) = 1.;
				continue;
			}
			Ia(pivot_i, i) = a(pivot_i, i) / Itmp;
		}

		for (i=0; i<m+1; i++) {
			if (i == pivot_i) continue;
			Itmp = a(i, pivot_j);
			for (j=0; j<n+1; j++) {
				if (isbasic(j)) {
					if (j == pivot_j) a(i, j) = 0.;
					continue;
				}
				Ia(i, j) = a(i, j) - Ia(pivot_i, j) * Itmp;
			}
		}

		for (i=0; i<n+1; i++) {
			if (isbasic(i)) continue;
			if ((i == 0) ^ (round == 1))  {
				a(pivot_i, i) = Ia(pivot_i, i).upper();
			} else {
				a(pivot_i, i) = Ia(pivot_i, i).lower();
			}
		}
		for (i=0; i<m+1; i++) {
			if (i == pivot_i) continue;
			for (j=0; j<n+1; j++) {
				if (isbasic(j)) continue;
				if ((j == 0) ^ (round == 1) ^ (i == 0))  {
					a(i, j) = Ia(i, j).upper();
				} else {
					a(i, j) = Ia(i, j).lower();
				}
			}
		}
	}

	return a(0, 0);
}

} // namespace kv

#endif // LP_HPP
