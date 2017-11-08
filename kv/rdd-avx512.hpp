/*
 * Copyright (c) 2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_AVX512_HPP
#define RDD_AVX512_HPP

#include <cmath>
#include <limits>
#include <immintrin.h>


#ifndef DD_NEW_SQRT
#define DD_NEW_SQRT 1
#endif


namespace kv {

template <> struct rop <dd> {

	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		volatile __m128d a1, b1, x1, y1;
		x = a * b;
		a1 = _mm_set_sd(a);
		b1 = _mm_set_sd(b);
		x1 = _mm_set_sd(x);
		y1 = _mm_fmsub_round_sd(a1, b1, x1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&y, y1);
	}

	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		volatile __m128d a1, b1, x1, y1;
		x = a * b;
		a1 = _mm_set_sd(a);
		b1 = _mm_set_sd(b);
		x1 = _mm_set_sd(x);
		y1 = _mm_fmsub_round_sd(a1, b1, x1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&y, y1);
	}

	static dd add_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mx, my, mz;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		mx = _mm_set_sd(x.a2);
		my = _mm_set_sd(y.a2);
		mz = _mm_set_sd(z2);
		mx = _mm_add_round_sd(mx, my, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		mz = _mm_add_round_sd(mz, mx, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mx, my, mz;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		mx = _mm_set_sd(x.a2);
		my = _mm_set_sd(y.a2);
		mz = _mm_set_sd(z2);
		mx = _mm_add_round_sd(mx, my, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		mz = _mm_add_round_sd(mz, mx, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mx, my, mz;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		mx = _mm_set_sd(x.a2);
		my = _mm_set_sd(y.a2);
		mz = _mm_set_sd(z2);
		mx = _mm_sub_round_sd(mx, my, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		mz = _mm_add_round_sd(mz, mx, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mx, my, mz;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		mx = _mm_set_sd(x.a2);
		my = _mm_set_sd(y.a2);
		mz = _mm_set_sd(z2);
		mx = _mm_sub_round_sd(mx, my, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		mz = _mm_add_round_sd(mz, mx, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mz2, mxa1, mxa2, mya1, mya2;

		twoproduct_up(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		mz2 = _mm_set_sd(z2);
		mxa1 = _mm_set_sd(x.a1);
		mxa2 = _mm_set_sd(x.a2);
		mya1 = _mm_set_sd(y.a1);
		mya2 = _mm_set_sd(y.a2);

		mz2 = _mm_fmadd_round_sd(mxa1, mya2, mz2, _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
		mz2 = _mm_fmadd_round_sd(mxa2, mya1, mz2, _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);
		mz2 = _mm_fmadd_round_sd(mxa2, mya2, mz2, _MM_FROUND_TO_POS_INF |_MM_FROUND_NO_EXC);

		_mm_store_sd(&z2, mz2);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mz2, mxa1, mxa2, mya1, mya2;

		twoproduct_down(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		mz2 = _mm_set_sd(z2);
		mxa1 = _mm_set_sd(x.a1);
		mxa2 = _mm_set_sd(x.a2);
		mya1 = _mm_set_sd(y.a1);
		mya2 = _mm_set_sd(y.a2);

		mz2 = _mm_fmadd_round_sd(mxa1, mya2, mz2, _MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC);
		mz2 = _mm_fmadd_round_sd(mxa2, mya1, mz2, _MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC);
		mz2 = _mm_fmadd_round_sd(mxa2, mya2, mz2, _MM_FROUND_TO_NEG_INF |_MM_FROUND_NO_EXC);

		_mm_store_sd(&z2, mz2);

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mz1, mz2, mz3, mz4, mxa1, mxa2, mya1, mya2, mtmp;

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		// if (z1 != z1) return std::numeric_limits<dd>::infinity();
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y.a1 >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);

				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1 * 0.5);
				mxa2 = _mm_set_sd(x.a2 * 0.5);
				mya1 = _mm_set_sd(y.a1 * 0.5);
				mya2 = _mm_set_sd(y.a2 * 0.5);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			} else {
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1);
				mxa2 = _mm_set_sd(x.a2);
				mya1 = _mm_set_sd(y.a1);
				mya2 = _mm_set_sd(y.a2);

				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			}
		} else {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1 * 0.5);
				mxa2 = _mm_set_sd(x.a2 * 0.5);
				mya1 = _mm_set_sd(y.a1 * 0.5);
				mya2 = _mm_set_sd(y.a2 * 0.5);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			} else {
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1);
				mxa2 = _mm_set_sd(x.a2);
				mya1 = _mm_set_sd(y.a1);
				mya2 = _mm_set_sd(y.a2);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			}
		}

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd div_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile __m128d mz1, mz2, mz3, mz4, mxa1, mxa2, mya1, mya2, mtmp;

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		// if (z1 != z1) return std::numeric_limits<dd>::infinity();
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y.a1 >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);

				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1 * 0.5);
				mxa2 = _mm_set_sd(x.a2 * 0.5);
				mya1 = _mm_set_sd(y.a1 * 0.5);
				mya2 = _mm_set_sd(y.a2 * 0.5);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			} else {
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1);
				mxa2 = _mm_set_sd(x.a2);
				mya1 = _mm_set_sd(y.a1);
				mya2 = _mm_set_sd(y.a2);

				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			}
		} else {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1 * 0.5);
				mxa2 = _mm_set_sd(x.a2 * 0.5);
				mya1 = _mm_set_sd(y.a1 * 0.5);
				mya2 = _mm_set_sd(y.a2 * 0.5);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			} else {
				mz1 = _mm_set_sd(z1);
				mz3 = _mm_set_sd(z3);
				mz4 = _mm_set_sd(z4);
				mxa1 = _mm_set_sd(x.a1);
				mxa2 = _mm_set_sd(x.a2);
				mya1 = _mm_set_sd(y.a1);
				mya2 = _mm_set_sd(y.a2);
				
				// z2 = (((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4;
				mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_fnmadd_round_sd(mz1, mya2, mz2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

				_mm_store_sd(&z2, mz2);

				if (z2 > 0.) {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				} else {
					mtmp = _mm_add_round_sd(mya1, mya2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				}
				mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
				_mm_store_sd(&z2, mz2);
			}
		}

		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

#if DD_NEW_SQRT == 1
	static dd sqrt_up(const dd& x) {
		double z1, z2, z3, z4;
		volatile __m128d mz1, mz2, mz3, mz4, mxa1, mxa2, mtmp;

		#if 0
		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
		}
		#endif

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		z1 = std::sqrt(x.a1);
		twoproduct_up(-z1, z1, z3, z4);

		mz1 = _mm_set_sd(z1);
		mz3 = _mm_set_sd(z3);
		mz4 = _mm_set_sd(z4);
		mxa1 = _mm_set_sd(x.a1);
		mxa2 = _mm_set_sd(x.a2);

		// z2 = (z3 + x.a1) + x.a2 + z4;
		mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz2);
		if (z2 > 0.) {
			// tmp = std::sqrt(x.a1 + x.a2) + z1;
			mtmp = _mm_add_round_sd(mxa1, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_sqrt_round_sd(mtmp, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_add_round_sd(mtmp, mz1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		} else {
			mtmp = _mm_add_round_sd(mxa1, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_sqrt_round_sd(mtmp, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_add_round_sd(mtmp, mz1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		}
		mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz2);
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sqrt_down(const dd& x) {
		double z1, z2, z3, z4;
		volatile __m128d mz1, mz2, mz3, mz4, mxa1, mxa2, mtmp;

		#if 0
		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
		}
		#endif

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		z1 = std::sqrt(x.a1);
		twoproduct_down(-z1, z1, z3, z4);

		mz1 = _mm_set_sd(z1);
		mz3 = _mm_set_sd(z3);
		mz4 = _mm_set_sd(z4);
		mxa1 = _mm_set_sd(x.a1);
		mxa2 = _mm_set_sd(x.a2);

		// z2 = (z3 + x.a1) + x.a2 + z4;
		mz2 = _mm_add_round_sd(mz3, mxa1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		mz2 = _mm_add_round_sd(mz2, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		mz2 = _mm_add_round_sd(mz2, mz4, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz2);
		if (z2 > 0.) {
			// tmp = std::sqrt(x.a1 + x.a2) + z1;
			mtmp = _mm_add_round_sd(mxa1, mxa2, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_sqrt_round_sd(mtmp, mtmp, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_add_round_sd(mtmp, mz1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		} else {
			mtmp = _mm_add_round_sd(mxa1, mxa2, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_sqrt_round_sd(mtmp, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			mtmp = _mm_add_round_sd(mtmp, mz1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		}
		mz2 = _mm_div_round_sd(mz2, mtmp, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&z2, mz2);
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}
#else
	static dd sqrt_up(const dd& x) {
		dd r, r2;

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_up(x, r);

		if (r > r2) return r;
		else return r2;
	}

	static dd sqrt_down(const dd& x) {
		dd r, r2;

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_down(x, r);

		if (r > r2) return r2;
		else return r;
	}
#endif

	static void begin() {
	}

	static void end() {
	}

	static void print_up(const dd& x, std::ostream& s) {
		char format;
		if (s.flags() & s.scientific) {
			if (s.flags() & s.fixed) {
				format = 'g';
			} else {
				format = 'e';
			}
		} else {
			if (s.flags() & s.fixed) {
				format = 'f';
			} else {
				format = 'g';
			}
		}
		s << conv_dd::ddtostring(x.a1, x.a2, s.precision(), format, 1);
	}

	static void print_down(const dd& x, std::ostream& s) {
		char format;
		if (s.flags() & s.scientific) {
			if (s.flags() & s.fixed) {
				format = 'g';
			} else {
				format = 'e';
			}
		} else {
			if (s.flags() & s.fixed) {
				format = 'f';
			} else {
				format = 'g';
			}
		}
		s << conv_dd::ddtostring(x.a1, x.a2, s.precision(), format, -1);
	}

	static dd fromstring_up(const std::string& s) {
		double x1, x2;
		conv_dd::stringtodd(s, x1, x2, 1);
		return dd(x1, x2);
	}

	static dd fromstring_down(const std::string& s) {
		double x1, x2;
		conv_dd::stringtodd(s, x1, x2, -1);
		return dd(x1, x2);
	}
};

} // namespace kv

#endif // RDD_AVX512_HPP
