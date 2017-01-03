/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ALLSOL_HPP
#define ALLSOL_HPP

#include <iostream>
#include <list>
#include <limits>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
// #include <boost/random.hpp>
#include <kv/matrix-inversion.hpp>
#include <kv/autodif.hpp>


#ifndef EDGE_RATIO
#define EDGE_RATIO 0.9
#endif

#ifndef USE_TRIM
#define USE_TRIM 1
#endif

#ifndef USE_ZERODIVIDE
#define USE_ZERODIVIDE 2
#endif

#ifndef USE_MULTITRIM
#define USE_MULTITRIM 0
#endif

#ifndef USE_AFFINEABS
#define USE_AFFINEABS 0
#endif

#ifndef USE_MIGSKIP
#define USE_MIGSKIP 1
#endif

#ifndef RECOVER_RATIO
#define RECOVER_RATIO 0.1
#endif

#ifndef USE_SLOWDIVIDE
#define USE_SLOWDIVIDE 1
#endif

#ifndef USE_FI
#define USE_FI 1
#endif

#ifndef USE_SUPERLINEAR
#define USE_SUPERLINEAR 1
#endif

#ifndef UNIFY_REST
#define UNIFY_REST 1
#endif

#ifndef EXISTENCE_RATIO
#define EXISTENCE_RATIO 0.8
#endif

#ifndef ITER_STOP_RATIO
#define ITER_STOP_RATIO 0.9
#endif


namespace kv {

namespace bn = boost::numeric;
namespace ub = boost::numeric::ublas;

namespace allsol_sub {

// return index of I_i which has maximum width

template <class T> int search_maxwidth (const ub::vector< interval<T> >& I) {
	int s = I.size();
	int i, mi;
	T m, tmp;

	m = 0.;
	for (i=0; i<s; i++) {
		tmp = width(I(i));
		if (tmp > m) {
			m = tmp; mi = i;
		}
	}

	return mi;
}

// return max width(I_i) / width(J_i)

template <class T> T widthratio_max (const ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J) {
	int s = I.size();
	int i;
	T tmp, r;

	r = 0.;

	for (i=0; i<s; i++) {
		tmp = width(I(i)) / width(J(i));
		if (tmp > r) r = tmp;
	}

	return r;
}

// return min width(I_i) / width(J_i)

template <class T> T widthratio_min (const ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J) {
	int s = I.size();
	int i;
	T tmp, r;

	r = std::numeric_limits<T>::max();

	for (i=0; i<s; i++) {
		tmp = width(I(i)) / width(J(i));
		if (tmp < r) r = tmp;
	}

	return r;
}

// J: original interval, I: shrinked interval
// If width(I(i)) < ratio * width(J(i)) then inflate I(i)
//  until width(I(i)) = ratio * width(J(i))
// Inflation outside J(i) is "not allowed".

template <class T> void recovery_inflation (ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J, T ratio) {
	int s = I.size();
	int i;
	T tmp, tmp2, l, u;

	for (i=0; i<s; i++) {
		if ( width(I(i)) < ratio * width(J(i)) ) {
			tmp = mid(I(i));
			tmp2 = rad(J(i)) * ratio;
			l = tmp - tmp2;
			u = tmp + tmp2;
			if (l < J(i).lower()) {
				u += (J(i).lower() - l);
				l = J(i).lower();
			} else if (u > J(i).upper()) {
				l -= (u - J(i).upper());
				u = J(i).upper();
			}
			// to be sure that new I(i) includes original I(i)
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}

// J: original interval, I: shrinked interval
// If width(I(i)) < ratio * width(J(i)) then inflate I(i)
//  until width(I(i)) = ratio * width(J(i))
// Inflation outside J(i) is "allowed".

template <class T> void recovery_inflation2 (ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J, T ratio) {
	int s = I.size();
	int i;
	T tmp, tmp2, l, u;

	for (i=0; i<s; i++) {
		if ( width(I(i)) < ratio * width(J(i)) ) {
			tmp = mid(I(i));
			tmp2 = rad(J(i)) * ratio;
			l = tmp - tmp2;
			u = tmp + tmp2;
			// to be sure that new I(i) includes original I(i)
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}

#if UNIFY_REST == 1
// convex hull of interval vector
template <class T> ub::vector< interval<T> > iv_hull (const ub::vector< interval<T> >& x, const ub::vector< interval<T> >& y) {
	int i;
	int s = x.size();
	ub::vector< interval<T> > r;

	r.resize(s);
	for (i=0; i<s; i++) {
		r(i) = interval<T>::hull(x(i), y(i));
	}

	return r;
}
#endif // UNIFY_REST == 1

// **not used**
// return index of division such that norm of M seems to be smallest
// under the scaled norm u = rad(I:divided)

template <class T> int search_optimal_divide (const ub::matrix< interval<T> >& M, const ub::vector< interval<T> >& I) {
	int n = I.size();
	ub::vector<T> u(n), s(n);
	ub::matrix<T> m(n,n);
	int i, j, p, mp;
	T tmp, tmp2, tmp3;

	for (i=0; i<n; i++) u(i) = width(I(i));
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			m(i,j) = norm(M(i,j));
		}
	}
	for (i=0; i<n; i++) {
		tmp = 0.;
		for (j=0; j<n; j++) {
			tmp += m(i,j) * u(j);
		}
		s(i) = tmp / u(i);
	}

	tmp3 = std::numeric_limits<T>::max();
	for (p=0; p<n; p++) {
		tmp2 = 0.;
		for (i=0; i<n; i++) {
			tmp = s(i) - m(i,p)*u(p)/u(i)/2.;
			if (i == p) tmp *= 2.;
			if (tmp > tmp2) tmp2 = tmp;
		}
		if (tmp2 < tmp3) {
			tmp3 = tmp2;
			mp = p;
		}
	}

	return mp;
}

// **not used**
// return index such that width(K(i)) is relatively smallest
// compared with width(I(i))

template <class T> int search_optimal_divide2 (const ub::vector< interval<T> >& K, const ub::vector< interval<T> >& I) {
	int n = I.size();
	T tmp, tmp2;
	int i, r;

	tmp2 = std::numeric_limits<T>::max();
	for (i=0; i<n; i++) {
		tmp =  width(K(i))/ width(I(i));
		if (tmp < tmp2) {
			tmp2 = tmp;
			r = i;
		}
	}

	return r;
}

} // namespace allsol_sub


// find all solution of f in I

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol (
F f,
const ub::vector< interval<T> >& I,
int verbose=1,
T giveup = T(0.),
std::list< ub::vector < interval<T> > >* rest=NULL
)
{
	std::list< ub::vector < interval<T> > > targets;
	targets.push_back(I);
	return allsol_list(f, targets, verbose, giveup, rest);
}


// find all solution of f in targets (list of intervals)

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol_list (
F f,
std::list< ub::vector< interval<T> > > targets,
int verbose=1,
T giveup = T(0.),
std::list< ub::vector < interval<T> > >* rest=NULL
)
{
	int s = (targets.front()).size();
	std::list< ub::vector< interval<T> > > solutions, solutions_big;
	int count_ne_test = 0;
	int count_ex_test = 0;
	int count_unknown = targets.size();
	int count_ne = 0;
	int count_ex = 0;
	int count_giveup = 0;

	#pragma omp parallel
	{

	ub::vector< interval<T> > I, fc, fi, C, CK, K, mvf, I1, I2, IR, Iorg, g;
	ub::vector<T> v;
	ub::matrix< interval<T> > fdi, M, L2;
	ub::matrix<T> L, R, E;
	typename std::list< ub::vector< interval<T> > >::iterator p, p2;
	int i, j, k, mi;
	T tmp, tmp2;
	bool r, M_calculated, flag, flag2, flag3;
	interval<T> A, B, J, J2, Itmp;

	// boost::variate_generator<boost::mt19937, boost::uniform_int<> > rand (boost::mt19937(time(0)), boost::uniform_int<>(0, s-1));

	E = ub::identity_matrix<T>(s);

	while (true) {
		if (verbose >= 2) {
			#pragma omp critical (cout)
			{
			std::cout << "ne_test: " << count_ne_test << ", ex_test: " << count_ex_test << ", unknown: " << count_unknown << ", ne: " << count_ne << ", ex: " << count_ex << ", giveup: " << count_giveup << "    \r" << std::flush;
			}
		}

		#ifdef _OPENMP

		int iflag = 0;
		#pragma omp critical (targets)
		{
		if (count_unknown == 0) iflag = 2;
		else {
			if (targets.empty()) {
				iflag = 1;
			} else {
				I = targets.front();
				targets.pop_front();
			}
		}
		}
		if (iflag == 2)  break;
		if (iflag == 1) continue;

		#else // _OPENMP

		if (targets.empty()) break;
		I = targets.front();
		targets.pop_front();

		#endif // _OPENMP

		Iorg = I;

		// non-existence test

		#pragma omp atomic
		count_ne_test++;

#if USE_FI == 1
		try {
			fi = f(I);
		}
		catch (std::domain_error& e) {
			goto label;
		}

		if (!zero_in(fi)) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}
#endif

		try {
			autodif< interval<T> >::split(f(autodif< interval<T> >::init(I)), fi, fdi);
		}
		catch (std::domain_error& e) {
			goto label;
		}

#if USE_FI != 1
		if (!zero_in(fi)) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}
#endif

		C = mid(I);
		try {
			fc = f(C);
		}
		catch (std::domain_error& e) {
			goto label;
		}

		mvf = fc + prod(fdi, I - C);
		if (!zero_in(mvf)) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}

#if USE_TRIM == 1
		// interval shrinking
		IR = I;
		flag = false; // non-existence in I turns out or not
		flag2 = false; // shrinking of I occurs or not
#if USE_ZERODIVIDE == 2
		flag3 = false; // division of I occurs or not
#endif
		for (i=0; i<s; i++) {
			for (j=0; j<s; j++) {
				B = fdi(i, j);

				// calculate back A from mvf
				A = mvf(i);
				Itmp = B * (I(j)-C(j));
				rop<T>::begin();
				tmp = rop<T>::sub_down(A.lower(), Itmp.lower());
				tmp2 = rop<T>::sub_up(A.upper(), Itmp.upper());
				rop<T>::end();
				A.assign(tmp, tmp2);
#if 0
				// calculate A simply
				// below is simple but slow
				A = 0.;
				for (k=0; k<s; k++) {
					if (k == j) continue;
					A += fdi(i, k) * (I(k)-C(k));
				}
				A += fc(i);
#endif

				if (zero_in(B)) {
#if USE_ZERODIVIDE >= 1
					bool bdummy;
					if (rad(B) <= 0.) continue;
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
							allsol_sub::recovery_inflation(I1, Iorg, (T)RECOVER_RATIO);
							allsol_sub::recovery_inflation(I2, Iorg, (T)RECOVER_RATIO);
							#pragma omp critical (targets)
							{
							targets.push_back(I1);
							targets.push_back(I2);
							count_unknown += 1;
							}
							flag3 = true;
							break;
#else
							continue;
#endif
						} else {
							if (!subset(IR(j), J)) {
								IR(j) = intersect(IR(j), J);
								flag2 = true;
							}
						}
					} else {
						if (overlap(IR(j), J2)) {
							if (!subset(IR(j), J2)) {
								IR(j) = intersect(IR(j), J2);
								flag2 = true;
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
							flag2 = true;
						}
					} else {
						flag = true;
						break;
					}
				}
			} // loop j
			if (flag == true) break;
#if USE_ZERODIVIDE == 2
			if (flag3 == true) break;
#endif
		} // loop i

#if USE_ZERODIVIDE == 2
		// interval division occurs
		if (flag3 == true) continue;
#endif

		// non-existence in I turns out
		if (flag == true) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}

		if (flag2 == true) {
			allsol_sub::recovery_inflation(IR, Iorg, (T)RECOVER_RATIO);
#if USE_MULTITRIM == 1
			// if radius of I is smaller than half on the 
			// "division scheduled" index by interval shrinking,
			// skip the existence test and do interval
			// shrinking again.
			mi = allsol_sub::search_maxwidth(I);
			if (rad(IR(mi)) <= 0.5 * rad(I(mi))) {
				#pragma omp critical (targets)
				{
				targets.push_back(IR);
				}
				continue;
			}
#endif

			// renew I, C and f(C)
			// do not renew f'(I) because of calculation cost
			I = IR;
			C = mid(I);
			try {
				fc = f(C);
			}
			catch (std::domain_error& e) {
				goto label;
			}

			// re-check mvf
			mvf = fc + prod(fdi, I - C);
			if (!zero_in(mvf)) {
				#pragma omp atomic
				count_ne++;
				#pragma omp critical (targets)
				{
				count_unknown--;
				}
				continue;
			}
		}

#endif // USE_TRIM == 1

#if USE_AFFINEABS == 1
		g.resize(s);
		for (i=0; i<s; i++) {
			g(i) = fc(i) / rad(mvf(i));
		}
		L2 = mid(fdi);
		I1 = fc + prod(fdi - L2, I - C);
		I2 = prod(g, L2);
		if (inner_prod(mag(g), mig(I1)) - inner_prod(mag(I2), rad(I)) > 0.) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}
#endif

		// existence test

		M_calculated = false;

		// skip existence test using mig
#if USE_MIGSKIP >= 1
#if USE_MIGSKIP == 1
		// use "loose" condition
		// interval<T> is for avoiding compile error
		mvf = fc + prod(mig(fdi), interval<T>(1. + EDGE_RATIO) * (I - C));
#else
		// use "strict" condition
		mvf = fc + prod(mig(fdi), I - C);
#endif
		if (!zero_in(mvf)) goto label;
#endif

		L = mid(fdi);

#if 0
		// **not used**
		// below is same as mig?
		I1 = - fc + prod(L-fdi, I-C);
		I2 = prod(L, I-C);
		if (!subset(I1, I2)) goto label;
#endif

#if 0
		// **not used**
		// prod(mag(L) - mag(L - fdi), rad(I));
		L2 = L - fdi;
		v = prod(mag(L) - mag(L2), rad(I));

		flag = true;
		for (i=0; i<s; i++) {
			if (abs(fc(i)) > v(i)) {flag=false; break;}
		}
		if (flag == false) goto label;
#endif

		#pragma omp atomic
		count_ex_test++;

		r = invert(L, R);
		if (!r) goto label;

		M = E - prod(R, fdi);
		M_calculated = true;
		CK = C - prod(R, fc);
		K = CK +  prod(M, I - C);
		if (!overlap(K, I)) {
			#pragma omp atomic
			count_ne++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}

#if 0
		// **not used**
		// same condition as above !overlap(K, I)
		I1 = Rfc + prod(Rfdi, I - C);
		if (!zero_in(I1)) {
			count_ne++;
			continue;
		}
#endif

		if (proper_subset(K, I) && allsol_sub::widthratio_max(K, I) < EXISTENCE_RATIO ) {
			#pragma omp critical (solutions)
			{
			// check whether the solution is already found or not
			flag = true;
			p = solutions.begin();
			p2 = solutions_big.begin();
			while (p != solutions.end()) {
				if (overlap(K, *p)) {
					if (subset(K, *p2)||subset(*p, I)) {
						flag = false;
						break;
					}
					while (true) {
						C = mid(K);
						I1 = C - prod(R, f(C)) + prod(M, K - C);
						K = intersect(K, I1);
						if (subset(K, *p2)) {
							flag2 = true;
							break;
						}
						if (!overlap(K, *p)) {
							flag2 = false;
							break;
						}
					}
					if (flag2 == true) {
						flag = false;
						break;
					} else {
						/* never reach? */
						std::cout << "two overlap intervals includes different solutions";
					}
				}
				p++;
				p2++;
			}
			if (flag) { // new solution found
				if (verbose >= 1) {
					#pragma omp critical (cout)
					{
					std::cout << I << "(ex)\n";
					}
				}
				solutions_big.push_back(I);
				// iterative refinement
				while (1) {
					C = mid(K);
					#if USE_SUPERLINEAR == 1
					autodif< interval<T> >::split(f(autodif< interval<T> >::init(K)), fi, fdi);
					L = mid(fdi);
					r = invert(L, R);
					M = E - prod(R, fdi);
					#endif
					I1 = C - prod(R, f(C)) + prod(M, K - C);
					I1 = intersect(K, I1);
					tmp = allsol_sub::widthratio_min(I1, K);
					K = I1;
					if (tmp > ITER_STOP_RATIO) break;
				}
				solutions.push_back(K);
				count_ex++;
				if (verbose >= 1) {
					#pragma omp critical (cout)
					{
					std::cout << K << "(ex:improved)\n";
					}
				}
			}
			} // pragma omp critical (solutions)
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}

		// check the case that solution may exist near boundary.
		// If so, use K as next interval
		if (allsol_sub::widthratio_max(K, I) < EDGE_RATIO) {
			allsol_sub::recovery_inflation2(K, Iorg, (T)RECOVER_RATIO);
			#pragma omp critical (targets)
			{
			targets.push_back(K);
			}
			continue;
		}


		I = intersect(I, K);
		allsol_sub::recovery_inflation(I, Iorg, (T)RECOVER_RATIO);

		label:

		// divide interval

		// if (M_calculated) mi = search_optimal_divide(M, I);
		// if (M_calculated) mi = search_optimal_divide2(K, I);
		// else mi = allsol_sub::search_maxwidth(I);
		// mi = rand();

#if USE_SLOWDIVIDE == 1
		if (allsol_sub::widthratio_max(I, Iorg) <= 0.5) {
			#pragma omp critical (targets)
			{
			targets.push_back(I);
			}
			continue;
		}
#endif

#ifdef USE_SLOWDIVIDE_OLD
		
		// if radius of I is smaller than half on the 
		// "division scheduled" index by interval shrinking,
		// skip division.

		mi = allsol_sub::search_maxwidth(Iorg);
		if (rad(I(mi)) <= 0.5 * rad(Iorg(mi))) {
			#pragma omp critical (targets)
			{
			targets.push_back(I);
			}
			continue;
		}
#else
		mi = allsol_sub::search_maxwidth(I);
#endif

		tmp = mid(I(mi));
		if (width(I(mi)) < giveup || tmp == I(mi).lower() || tmp == I(mi).upper()) {
			if (verbose >= 2) {
				std::cout << "too small interval (may be multiple root?):\n" << I << "\n";
			}
			if (rest != NULL) {
				#pragma omp critical (rest)
				{
				#if UNIFY_REST == 1
				I1 = I;
				while (true) {
					flag = false;
					p = (*rest).begin();
					while (p != (*rest).end()) {
						if (overlap(*p, I1)) {
							I1 = allsol_sub::iv_hull(I1, *p);
							p = (*rest).erase(p);
							// #pragma omp atomic
							// count_giveup--;
							flag = true;
							continue;
						}
						p++;
					}
					if (flag == false) break;
				}
				(*rest).push_back(I1);
				#else // UNIFY_REST == 1
				(*rest).push_back(I);
				#endif // UNIFY_REST == 1
				}
			}
			#pragma omp atomic
			count_giveup++;
			#pragma omp critical (targets)
			{
			count_unknown--;
			}
			continue;
		}

		I1 = I; I2 = I;
		I1(mi).assign(I1(mi).lower(), tmp);
		I2(mi).assign(tmp, I2(mi).upper());
		#pragma omp critical (targets)
		{
		targets.push_back(I1);
		targets.push_back(I2);
		count_unknown += 1;
		}
	}

	} // pragma omp parallel

	if (verbose >= 1) {
			std::cout << "ne_test: " << count_ne_test << ", ex_test: " << count_ex_test << ", ne: " << count_ne << ", ex: " << count_ex << ", giveup: " << count_giveup << "    \n";
	}

	return solutions;
}

} // namespace kv

#endif // ALLSOL_HPP
