/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ALLSOL_AFFINE_HPP
#define ALLSOL_AFFINE_HPP

#include <iostream>
#include <list>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/interval-vector.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/matrix-inversion.hpp>
#include <kv/affine.hpp>
#include <kv/lp.hpp>

namespace ub = boost::numeric::ublas;


#ifndef USE_TRIM
#define USE_TRIM 1
#endif

#ifndef USE_SUMABS
#define USE_SUMABS 0
#endif

#ifndef USE_LP
#define USE_LP 0
#endif

#ifndef USE_FI
#define USE_FI 1
#endif


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
allsol_affine(const ub::vector< interval<T> >& init, F f, int verbose=1)
{
	std::list< ub::vector < interval<T> > > targets;
	targets.push_back(init);
	return allsol_list_affine(targets, f, verbose);
}


// 初期区間をリストで複数与えるallsol

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol_list_affine(std::list< ub::vector< interval<T> > > targets, F f, int verbose=1)
{
	int s = (targets.front()).size();
	ub::vector< interval<T> > I, fc, fi, C, CK, K, mvf, I1, I2;
	ub::matrix< interval<T> > fdi, M;
	ub::matrix<T> L, R, E;
	std::list< ub::vector< interval<T> > > solutions, solutions_big;
	typename std::list< ub::vector< interval<T> > >::iterator p, p2;
	int i, j, k, mi;
	T tmp;
	bool r, M_calculated, flag, flag2;
	int count_ne_test = 0;
	int count_ex_test = 0;
	int count_unknown = targets.size();
	int count_ne = 0;
	int count_ex = 0;
	ub::vector< affine<T> > ax, ay, ak;
	ub::vector< interval<T> > IR;
	T trim_tmp;
	interval<T> trim_I;
	affine<T> sum_abs;
	ub::vector< ub::vector<T> > ay2;
	interval<T> I_tmp;
	T lp_tmp;
	ub::vector< interval<T> > objfunc_I;
	ub::vector<T> objfunc;
	std::list< ub::vector<T> > constraints;
	ub::vector<T> c_tmp;
	int ep_size, lp_size;
	bool lp_flag;
	T lp_return;

	E = ub::identity_matrix<T>(s);
	L.resize(s,s);

	while (!targets.empty()) {
		if (verbose >= 2) {
			std::cout << "ne_test: " << count_ne_test << ", ex_test: " << count_ex_test << ", unknown: " << count_unknown << ", ne: " << count_ne << ", ex: " << count_ex << "    \r" << std::flush;
		}

		I = targets.front();
		targets.pop_front();
		count_unknown--;

		// 以下、非存在テスト

		count_ne_test++;

#if USE_FI == 1
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
#endif

		affine<T>::maxnum() = 0;
		ax = I;
		try {
			ay = f(ax);
		}
		catch (std::range_error& e) {
			goto label;
		}

		fi = to_interval(ay);
		if (!zero_in(fi)) {
			count_ne++;
			continue;
		}

#if USE_SUMABS == 1
		sum_abs = 0.;
		for (i=0; i<s; i++) {
			sum_abs += abs(ay(i));
		}
		if (to_interval(sum_abs).lower() > 0.) {
			count_ne++;
			continue;
		}
#endif

#if USE_LP == 1
#if USE_LP_FAST == 1
		ep_size = 2 * s;
		lp_size = 1 + ep_size * 2;
		ay2.resize(s + ep_size);
		for (i=0; i<s; i++) {
			ay2(i).resize(lp_size);
			for (j=s + 1; j < lp_size; j++) {
				ay2(i)(j) = 0.;
			}
			I_tmp = ay(i).get_mid();
			for (j=1; j <= s; j++) {
				tmp = ay(i).get_coef(j);
				ay2(i)(j) = tmp;
				I_tmp -= tmp;
			}
			lp_tmp = ay(i).get_err();
			rop<T>::begin();
			for (j=s+1; j <= affine<T>::maxnum(); j++) {
				tmp = ay(i).get_coef(j);
				if (tmp < 0.) tmp = -tmp;
				lp_tmp = rop<T>::add_up(lp_tmp, tmp);
			}
			rop<T>::end();
			I_tmp -= lp_tmp;
			rop<T>::begin();
			lp_tmp = rop<T>::add_up(lp_tmp, rad(I_tmp));
			rop<T>::end();
			ay2(i)(s + 1 + i) = lp_tmp;
			ay2(i)(0) = I_tmp.lower();
			if (ay2(i)(0) > 0.) {
				for (j=0; j<lp_size; j++) {
					ay2(i)(j) = -ay2(i)(j);
				}
			}
		}
#else
		ep_size = affine<T>::maxnum() + s;
		lp_size = 1 + ep_size * 2;
		ay2.resize(s + ep_size);
		for (i=0; i<s; i++) {
			ay2(i).resize(lp_size);
			for (j=affine<T>::maxnum() + 1; j < lp_size; j++) {
				ay2(i)(j) = 0.;
			}
			I_tmp = ay(i).get_mid();
			for (j=1; j <= affine<T>::maxnum(); j++) {
				tmp = ay(i).get_coef(j);
				ay2(i)(j) = tmp;
				I_tmp -= tmp;
			}
			lp_tmp = ay(i).get_err();
			I_tmp -= lp_tmp;
			rop<T>::begin();
			lp_tmp = rop<T>::add_up(lp_tmp, rad(I_tmp));
			rop<T>::end();
			ay2(i)(affine<T>::maxnum() + 1 + i) = lp_tmp;
			ay2(i)(0) = I_tmp.lower();
			if (ay2(i)(0) > 0.) {
				for (j=0; j<lp_size; j++) {
					ay2(i)(j) = -ay2(i)(j);
				}
			}
		}
#endif
		for (i=0; i<ep_size; i++) {
			ay2(s + i).resize(lp_size);
			ay2(s + i)(0) = -2.;
			for (j=1; j<lp_size; j++) {
				ay2(s + i)(j) = 0.;
			}
			ay2(s + i)(1 + i) = 1.;
			ay2(s + i)(ep_size + 1 + i) = 1.;
		}
		objfunc_I.resize(lp_size);
		for (j=0; j<lp_size; j++) objfunc_I(j) = 0.;
		for (i=0; i<s + ep_size; i++) {
			for (j=0; j<lp_size; j++) {
				objfunc_I(j) -= ay2(i)(j);
			}
		}
		objfunc.resize(lp_size);
		for (j=0; j<lp_size; j++) {
			objfunc(j) = objfunc_I(j).lower();
		}

		lp_flag = true;
		constraints.clear();
		for (i=0; i<s + ep_size; i++) {
			c_tmp.resize(lp_size);
			for (j=0; j<lp_size; j++) {
				c_tmp(j) = ay2(i)(j);
			}
			/*
			if (c_tmp(0) >= 0.) {
				lp_flag = false;
			}
			*/
			constraints.push_back(c_tmp);
		}

		if (lp_flag) {
			lp_return = lp_minimize_verified(objfunc, constraints, -1);
			if (lp_return > 0.) {
				count_ne++;
				continue;
			}
		}
		
#endif

#if USE_TRIM == 1
		IR = I;
		flag = false;
		flag2 = false;
		for (i=0; i<s; i++) {
			trim_tmp = rad(ay(i));
			for (j=0; j<s; j++) {
				tmp = ay(i).get_coef(j+1);
				if (tmp == 0.) continue;
				if (tmp < 0.) tmp = -tmp;
				rop<T>::begin();
				tmp = rop<T>::sub_up(trim_tmp, tmp);
				rop<T>::end();
				trim_I = ax(j).get_mid() - (ay(i).get_mid() + tmp * interval<T>(-1., 1.)) * ax(j).get_coef(j+1) / ay(i).get_coef(j+1);
				if (overlap(IR(j), trim_I)) {
					if (!subset(IR(j), trim_I)) {
						IR(j) = intersect(IR(j), trim_I);
						flag2 = true;
					}
				} else {
					flag = true;
					break;
				}
			}
			if (flag == true) break;
		}

		if (flag == true) {
			count_ne++;
			continue;
		}
#endif

		// 以下、存在テスト

		for (i=0; i<s; i++) {
			for (j=0; j<s; j++) {
				L(i, j) = ay(i).get_coef(j+1) / ax(j).get_coef(j+1);
			}
		}

		count_ex_test++;

		r = invert(L, R);
		if (!r) goto label;

		ak = ax - prod(R, ay);
		K = to_interval(ak);

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
						affine<T>::maxnum() = 0;
						ax = K;
						ay = f(ax);
						#if 0
						for (i=0; i<s; i++) {
							for (j=0; j<s; j++) {
								L(i, j) = ay(i).get_coef(j+1) / ax(j).get_coef(j+1);
							}
						}
						invert(L, R);
						#endif
						ak = ax - prod(R, ay);
						I1 = to_interval(ak);
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
					affine<T>::maxnum() = 0;
					ax = K;
					ay = f(ax);
					#if 0
					for (i=0; i<s; i++) {
						for (j=0; j<s; j++) {
							L(i, j) = ay(i).get_coef(j+1) / ax(j).get_coef(j+1);
						}
					}
					r = invert(L, R);
					#endif
					ak = ax - prod(R, ay);
					I1 = to_interval(ak);
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

#if USE_TRIM == 1
		if (!overlap(IR, K)) {
			count_ne++;
			continue;
		} else {
			I = intersect(IR, K);
		}
#else
		I = intersect(I, K);
#endif

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

#endif // ALLSOL_AFFINE_HPP
