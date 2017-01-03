/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef OPTIMIZE_HPP
#define OPTIMIZE_HPP

#include <iostream>
#include <list>
#include <limits>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/autodif.hpp>


// 0: Do not use TRIM
// 1: use TRIM and use slow trim algorithm
// 2: use TRIM and use fast bat slightly inefficient trim algorithm
// 3: use TRIM and use new trim algorithm

#ifndef USE_TRIM
#define USE_TRIM 0
#endif


namespace kv {

namespace ub = boost::numeric::ublas;


template <class T, class F>
std::list< ub::vector< interval<T> > >
optimize(const ub::vector< interval<T> >& init, F f, T limit, bool unify = true, int verbose = 0)
{
	std::list< ub::vector< interval<T> > > targets;
	targets.push_back(init);
	return optimize_list(targets, f, limit, unify, verbose);
}

template <class T, class F>
std::list< ub::vector< interval<T> > >
optimize_list(std::list< ub::vector< interval<T> > > targets, F f, T limit, bool unify = true, int verbose = 0)
{
	int s = (targets.front()).size();
	ub::vector< interval<T> > I, C, I1, I2, IR, fdi, C2;
	interval<T> fc, fi, mvf, fc2; 
	T tmp, tmp2;
	std::list< ub::vector< interval<T> > > solutions;
	typename std::list< ub::vector< interval<T> > >::iterator p;
	int i, j, k, mi;
	bool flag;
	interval<T> A, B, J, J2, Itmp; 
#if USE_TRIM == 3
	ub::vector< interval<T> > A0, A1, A2; // for new trim algorithm
#endif // USE_TRIM == 3

	C2.resize(s);

	T delta = std::numeric_limits<T>::max();

	while (!targets.empty()) {
		I = targets.front();
		targets.pop_front();

		try {
			fi = f(I);
		}
		catch (std::domain_error& e) {
			goto label;
		}

		if (fi.lower() > delta) {
			continue;
		}

		C = mid(I);
		try {
			fc = f(C);
			autodif< interval<T> >::split(f(autodif< interval<T> >::init(I)), fi, fdi);
		}
		catch (std::domain_error& e) {
			goto label;
		}

		mvf = fc + inner_prod(fdi, I - C);
		if (mvf.lower() > delta) {
			continue;
		}

		// update delta at C
		tmp = fc.upper();
		if (tmp < delta) delta = tmp;

		// C2 is likely to give small value
		for (i=0; i<s; i++) {
			tmp = mid(fdi(i));
			if (tmp > 0.) tmp2 = I(i).lower();
			else if (tmp < 0.) tmp2 = I(i).upper();
			else tmp2 = mid(I(i));
			C2(i).assign(tmp2, tmp2);
		}

		// update delta at C2
		try {
			fc2 = f(C2);
		}
		catch (std::domain_error& e) {
			goto label;
		}
		tmp = fc2.upper();
		if (tmp < delta) delta = tmp;

#if USE_TRIM >= 1
		// interval shrinking

#if USE_TRIM == 3
		// prepare for new trim algorithm
		A0.resize(s);
		A1.resize(s);
		A2.resize(s);
		for (j=0; j<s; j++) {
			A0(j) = fdi(j) * (I(j)-C(j));
		}
		Itmp = 0.;
		for (j=0; j<s; j++) {
			A1(j) = Itmp;
			Itmp += A0(j);
		}
		Itmp = 0.;
		for (j=s-1; j>=0; j--) {
			A2(j) = Itmp;
			Itmp += A0(j);
		}
#endif // USE_TRIM == 3

		IR = I;
		flag = false; // non-existence in I turns out or not
		for (j=0; j<s; j++) {
			B = fdi(j);

#if USE_TRIM == 1
			// calculate A simply
			// simple but slow
			A = 0.;
			for (k=0; k<s; k++) {
				if (k == j) continue;
				A += fdi(k) * (I(k)-C(k));
			}
			A += fc;
#endif // USE_TRIM == 1
#if USE_TRIM == 2
			// old trim algorithm
			// calculate back A from mvf
			A = mvf;
			Itmp = B * (I(j)-C(j));
			rop<T>::begin();
			tmp = rop<T>::sub_down(A.lower(), Itmp.lower());
			tmp2 = rop<T>::sub_up(A.upper(), Itmp.upper());
			rop<T>::end();
			A.assign(tmp, tmp2);
#endif // USE_TRIM == 2
#if USE_TRIM == 3
			// new trim algorithm
			A = fc + A1(j) + A2(j);
#endif // USE_TRIM == 3

			A -= delta;

			if (zero_in(B)) {
#if USE_ZERODIVIDE >= 1
				bool bdummy;
				if (rad(B) <= 0) continue;
				J = C(j) - division_part1(A, B, bdummy);
				J2 = C(j) - division_part2(A, B);
				if (overlap(IR(j), J)) {
					if (overlap(IR(j), J2)) {
#if USE_ZERODIVIDE == 2
						if (overlap(J, J2)) continue;
						// interval division
						I1 = IR;
						I2 = IR;
						I1(j) = intersect(IR(j), J);
						I2(j) = intersect(IR(j), J2);
						targets.push_back(I1);
						targets.push_back(I2);
						flag = true;
						break;
#else
						continue;
#endif
					} else {
						if (!subset(IR(j), J)) {
							IR(j) = intersect(IR(j), J);
						}
					}
				} else {
					if (overlap(IR(j), J2)) {
						if (!subset(IR(j), J2)) {
							IR(j) = intersect(IR(j), J2);
						}
					} else {
						flag = true;
						break;
					}
				}
#else
				continue;
#endif
			} else {
				J = C(j) - A/B;
				if (overlap(IR(j), J)) {
					if (!subset(IR(j), J)) {
						IR(j) = intersect(IR(j), J);
					}
				} else {
					flag = true;
					break;
				}
			}
		}

		// non-existence in I turns out
		if (flag == true) {
			continue;
		}

		I = IR;

#endif // USE_TRIM >= 1

		label:;

		tmp2 = 0.;
		for (i=0; i<s; i++) {
			tmp = width(I(i));
			if (tmp > tmp2) {
				tmp2 = tmp; mi = i;
			}
		}

		if (tmp2 < limit) {
			if (verbose >= 1) {
				std::cout << I << "\n";
			}
			if (unify) {
				while (true) {
					flag = false;
                                        p = solutions.begin();
					while (p != solutions.end()) {
						if (overlap(*p, I)) {
							I = hull(I, *p);
							p = solutions.erase(p);
							flag = true;
							continue;
						}
						p++;
					}
					if (flag == false) break;
				}
			}
			solutions.push_back(I);
			continue;
		}

		tmp = mid(I(mi));
		I1 = I; I2 = I;
		I1(mi).assign(I1(mi).lower(), tmp);
		I2(mi).assign(tmp, I2(mi).upper());
		targets.push_back(I1);
		targets.push_back(I2);
	}

	if (verbose >= 1) {
		std::cout << delta << "\n";
	}

	return solutions;
}

} // namespace kv

#endif // OPTIMIZE_HPP
