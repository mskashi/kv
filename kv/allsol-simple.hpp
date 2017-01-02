/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ALLSOL_HPP
#define ALLSOL_HPP

#include <iostream>
#include <list>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/matrix-inversion.hpp>
#include <kv/autodif.hpp>

namespace ub = boost::numeric::ublas;

namespace kv {


// 最大幅の成分を与える

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

// 区間ベクトルの幅の比を成分毎に計算し、その最大値

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

// 区間ベクトルの幅の比を成分毎に計算し、その最小値

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


// 全解探索 中身はallsol_list

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol (const ub::vector< interval<T> >& init, F f, int verbose=1)
{
	std::list< ub::vector < interval<T> > > targets;
	targets.push_back(init);
	return allsol_list(targets, f, verbose);
}


// 初期区間をリストで複数与えるallsol

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol_list (std::list< ub::vector< interval<T> > > targets, F f, int verbose=1)
{
	int s = (targets.front()).size();
	ub::vector< interval<T> > I, fc, fi, C, CK, K, mvf, I1, I2;
	ub::matrix< interval<T> > fdi, M;
	ub::matrix<T> L, R, E;
	std::list< ub::vector< interval<T> > > solutions, solutions_big;
	typename std::list< ub::vector< interval<T> > >::iterator p, p2;
	int i, j, k, mi;
	T tmp;
	bool r, flag, flag2;
	int count_ne_test = 0;
	int count_ex_test = 0;
	int count_unknown = targets.size();
	int count_ne = 0;
	int count_ex = 0;

	E = ub::identity_matrix<T>(s);

	while (!targets.empty()) {
		if (verbose >= 2) {
			std::cout << "ne_test: " << count_ne_test << ", ex_test: " << count_ex_test << ", unknown: " << count_unknown << ", ne: " << count_ne << ", ex: " << count_ex << "    \r" << std::flush;
		}

		I = targets.front();
		targets.pop_front();
		count_unknown--;

		// 以下、非存在テスト

		count_ne_test++;

		try {
			fi = f(I);
		}
		catch (std::range_error& e) {
			goto label;
		}

		if (!zero_in(fi)) {
			count_ne++;
			continue;
		}

		C = mid(I);
		try {
			fc = f(C);
			autodif< interval<T> >::split(f(autodif< interval<T> >::init(I)), fi, fdi);
		}
		catch (std::range_error& e) {
			goto label;
		}

		mvf = fc + prod(fdi, I - C);
		if (!zero_in(mvf)) {
			count_ne++;
			continue;
		}

		// 以下、存在テスト

		L = mid(fdi);

		count_ex_test++;

		r = invert(L, R);
		if (!r) goto label;

		M = E - prod(R, fdi);
		CK = C - prod(R, fc);
		K = CK +  prod(M, I - C);
		if (!overlap(K, I)) {
			count_ne++;
			continue;
		}
		if (proper_subset(K, I)) {
			// 以下、既に見付かっている解との重複チェック
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
			if (flag) { // 新しい解が見付かった
				if (verbose >= 1) std::cout << I << "(ex)\n";
				solutions_big.push_back(I);
				// 反復改良
				while (1) {
					C = mid(K);
					I1 = C - prod(R, f(C)) + prod(M, K - C);
					I1 = intersect(K, I1);
					tmp = widthratio_min(I1, K);
					K = I1;
					if (tmp > 0.9) break;
				}
				solutions.push_back(K);
				count_ex++;
				if (verbose >= 1) std::cout << K << "(ex:improved)\n";
			}
			continue;
		}

		// 解が境界またはその近くに有るかも知れないケース
		// Kをそのまま次回のテスト区間に使う。
		if (widthratio_max(K, I) < 0.9) {
			targets.push_back(K);
			count_unknown++;
			continue;
		}

		I = intersect(I, K);

		label:

		// どの成分で分割するか選択

		mi = search_maxwidth(I);

		tmp = mid(I(mi));
		if (tmp == I(mi).lower() || tmp == I(mi).upper()) {
			std::cout << "too small interval (may be multiple root?):\n" << I << "\n";
			continue;
		}

		I1 = I; I2 = I;
		I1(mi).assign(I1(mi).lower(), tmp);
		I2(mi).assign(tmp, I2(mi).upper());
		targets.push_back(I1);
		targets.push_back(I2);
		count_unknown += 2;
	}

	if (verbose >= 1) {
			std::cout << "ne_test: " << count_ne_test << ", ex_test: " << count_ex_test << ", ne: " << count_ne << ", ex: " << count_ex << "    \n";
	}

	return solutions;
}

} // namespace kv

#endif // ALLSOL_HPP
