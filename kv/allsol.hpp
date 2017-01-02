/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
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

namespace bn = boost::numeric;
namespace ub = boost::numeric::ublas;


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

// Jが元、Iが縮小後として、成分毎に、縮小比が小さすぎる(ratio以下)な
// 成分があったら、ratioの比になるようにIを膨らませる。
// 膨らませるときは、Jの外にははみ出ないようにする。

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
			// 数学的厳密性のため念のため
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}

// はみ出ても良いversion

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
			// 数学的厳密性のため念のため
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}

#if 0
// 渡されたJより比率が悪くなるのを禁止するversion
// rateに意味は無い。
template <class T> void recovery_inflation (ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J, T rate) {
	int s = I.size();
	int i;
	T tmp, tmp2, l, u, m, limit;

	u = 0.;
	l = std::numeric_limits<T>::max();
	m = 0.;
	for (i=0; i<s; i++) {
		tmp = width(J(i));
		if (tmp > u) u = tmp;
		if (tmp < l) l = tmp;
		tmp = width(I(i));
		if (tmp > m) m = tmp;
	}
	limit = m / u * l;

	for (i=0; i<s; i++) {
		if ( width(I(i)) < limit ) {
			tmp = mid(I(i));
			tmp2 = limit / 2.;
			l = tmp - tmp2;
			u = tmp + tmp2;
			if (l < J(i).lower()) {
				u += (J(i).lower() - l);
				l = J(i).lower();
			} else if (u > J(i).upper()) {
				l -= (u - J(i).upper());
				u = J(i).upper();
			}
			// 数学的厳密性のため念のため
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}

template <class T> void recovery_inflation2 (ub::vector< interval<T> >& I, const ub::vector< interval<T> >& J, T rate) {
	int s = I.size();
	int i;
	T tmp, tmp2, l, u, m, limit;

	u = 0.;
	l = std::numeric_limits<T>::max();
	m = 0.;
	for (i=0; i<s; i++) {
		tmp = width(J(i));
		if (tmp > u) u = tmp;
		if (tmp < l) l = tmp;
		tmp = width(I(i));
		if (tmp > m) m = tmp;
	}
	limit = m / u * l;

	for (i=0; i<s; i++) {
		if ( width(I(i)) < limit ) {
			tmp = mid(I(i));
			tmp2 = limit / 2.;
			l = tmp - tmp2;
			u = tmp + tmp2;
			// 数学的厳密性のため念のため
			I(i) = interval<T>::hull(I(i), interval<T>(l, u));
		}
	}
}
#endif

// 分割後のIの半径をscaled normとしたとき、
// Mのノルムが最も小さく見えるような分割方向を返す。

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

// K(I)とIの幅を比べて相対的に最もKが小さいような成分を返す。

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


// 全解探索 中身はallsol_list

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol (const ub::vector< interval<T> >& init,
F f,
int verbose=1,
T giveup = T(0.),
std::list< ub::vector < interval<T> > >* rest=NULL
)
{
	std::list< ub::vector < interval<T> > > targets;
	targets.push_back(init);
	return allsol_list(targets, f, verbose, giveup, rest);
}


// 初期区間をリストで複数与えるallsol

template <class T, class F> std::list< ub::vector< interval<T> > >
allsol_list (
std::list< ub::vector< interval<T> > > targets,
F f,
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

		// 以下、非存在テスト

		#pragma omp atomic
		count_ne_test++;

#if USE_FI == 1
		try {
			fi = f(I);
		}
		catch (std::range_error& e) {
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
		catch (std::range_error& e) {
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
		catch (std::range_error& e) {
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
		// 区間縮小
		IR = I;
		flag = false; // 解の非存在が言えたか
		flag2 = false; // Iの削減が発生したか
#if USE_ZERODIVIDE == 2
		flag3 = false; // Iの分割が発生したか
#endif
		for (i=0; i<s; i++) {
			for (j=0; j<s; j++) {
				B = fdi(i, j);

				// Aをmvfから逆算して計算
				A = mvf(i);
				Itmp = B * (I(j)-C(j));
				rop<T>::begin();
				tmp = rop<T>::sub_down(A.lower(), Itmp.lower());
				tmp2 = rop<T>::sub_up(A.upper(), Itmp.upper());
				rop<T>::end();
				A.assign(tmp, tmp2);
#if 0
				// Aを真面目に計算すると遅い。
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
							// 区間分割
							I1 = IR;
							I2 = IR;
							I1(j) = intersect(IR(j), J);
							I2(j) = intersect(IR(j), J2);
							recovery_inflation(I1, Iorg, (T)RECOVER_RATIO);
							recovery_inflation(I2, Iorg, (T)RECOVER_RATIO);
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
		// 区間縮小で区間分割が発生
		if (flag3 == true) continue;
#endif

		// 区間縮小で非存在が判明
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
			recovery_inflation(IR, Iorg, (T)RECOVER_RATIO);
#if USE_MULTITRIM == 1
			// 区間縮小が「分割予定成分」に関して半分以下なら
			// 存在チェックをせずに再度区間縮小する。
			mi = search_maxwidth(I);
			if (rad(IR(mi)) <= 0.5 * rad(I(mi))) {
				#pragma omp critical (targets)
				{
				targets.push_back(IR);
				}
				continue;
			}
#endif

			// Iを更新し、Cとf(C)も更新してしまう。
			// f'(I)はコストが高いと思われるので更新せず。
			I = IR;
			C = mid(I);
			try {
				fc = f(C);
			}
			catch (std::range_error& e) {
				goto label;
			}

			// 一応mvfをもう一度
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

		// 以下、存在テスト

		M_calculated = false;

		// migを用いて存在テストをスキップ
#if USE_MIGSKIP >= 1
#if USE_MIGSKIP == 1
		// Inflationが行われなくなるのを防ぐためにやや条件を緩く
		// interval<T> が無いと場合によってはcompile error
		mvf = fc + prod(mig(fdi), interval<T>(1. + EDGE_RATIO) * (I - C));
#else
		// 手加減無し
		mvf = fc + prod(mig(fdi), I - C);
#endif
		if (!zero_in(mvf)) goto label;
#endif

		L = mid(fdi);

#if 0
		// migと同等?
		I1 = - fc + prod(L-fdi, I-C);
		I2 = prod(L, I-C);
		if (!subset(I1, I2)) goto label;
#endif

#if 0
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
		// 実験したところ、上の!overlap(K, I)と 
		// 全く同じ条件なので意味無し。
		I1 = Rfc + prod(Rfdi, I - C);
		if (!zero_in(I1)) {
			count_ne++;
			continue;
		}
#endif

		if (proper_subset(K, I)) {
			#pragma omp critical (solutions)
			{
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
				if (verbose >= 1) {
					#pragma omp critical (cout)
					{
					std::cout << I << "(ex)\n";
					}
				}
				solutions_big.push_back(I);
				// 反復改良
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
					tmp = widthratio_min(I1, K);
					K = I1;
					if (tmp > 0.9) break;
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

		// 解が境界またはその近くに有るかも知れないケース
		// Kをそのまま次回のテスト区間に使う。
		if (widthratio_max(K, I) < EDGE_RATIO) {
			recovery_inflation2(K, Iorg, (T)RECOVER_RATIO);
			#pragma omp critical (targets)
			{
			targets.push_back(K);
			}
			continue;
		}


		I = intersect(I, K);
		recovery_inflation(I, Iorg, (T)RECOVER_RATIO);

		label:

		// どの成分で分割するか選択

		// if (M_calculated) mi = search_optimal_divide(M, I);
		// if (M_calculated) mi = search_optimal_divide2(K, I);
		// else mi = search_maxwidth(I);
		// mi = rand();

#if USE_SLOWDIVIDE == 1
		if (widthratio_max(I, Iorg) <= 0.5) {
			#pragma omp critical (targets)
			{
			targets.push_back(I);
			}
			continue;
		}
#endif

#ifdef USE_SLOWDIVIDE_OLD
		// 区間縮小が「分割予定成分」に関して半分以下なら
		// 分割しない。そうでなければ普通に分割。

		mi = search_maxwidth(Iorg);
		if (rad(I(mi)) <= 0.5 * rad(Iorg(mi))) {
			#pragma omp critical (targets)
			{
			targets.push_back(I);
			}
			continue;
		}
#else
		mi = search_maxwidth(I);
#endif

		tmp = mid(I(mi));
		if (width(I(mi)) < giveup || tmp == I(mi).lower() || tmp == I(mi).upper()) {
			if (verbose >= 2) {
				std::cout << "too small interval (may be multiple root?):\n" << I << "\n";
			}
			if (rest != NULL) {
				#pragma omp critical (rest)
				{
				(*rest).push_back(I);
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
