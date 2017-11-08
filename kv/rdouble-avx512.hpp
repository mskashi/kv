/*
 * Copyright (c) 2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDOUBLE_AVX512_HPP
#define RDOUBLE_AVX512_HPP

#include <iostream>
#include <string>
#include <kv/conv-double.hpp>

#include <immintrin.h>


namespace kv {

template <> struct rop <double> {

	static double add_up(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_add_round_sd(x1, y1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double add_down(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_add_round_sd(x1, y1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double sub_up(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_sub_round_sd(x1, y1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double sub_down(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_sub_round_sd(x1, y1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double mul_up(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_mul_round_sd(x1, y1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double mul_down(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_mul_round_sd(x1, y1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double div_up(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_div_round_sd(x1, y1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double div_down(const double& x, const double& y) {
		volatile __m128d x1, y1, z1;
		double r;
		x1 = _mm_set_sd(x);
		y1 = _mm_set_sd(y);
		z1 = _mm_div_round_sd(x1, y1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double sqrt_up(const double& x) {
		volatile __m128d x1, z1;
		double r;
		x1 = _mm_set_sd(x);
		z1 = _mm_sqrt_round_sd(x1, x1, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static double sqrt_down(const double& x) {
		volatile __m128d x1, z1;
		double r;
		x1 = _mm_set_sd(x);
		z1 = _mm_sqrt_round_sd(x1, x1, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		_mm_store_sd(&r, z1);
		return r;
	}

	static void begin() {
	}

	static void end() {
	}

	static void print_up(const double& x, std::ostream& s) {
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
		s << conv_double::dtostring(x, s.precision(), format, 1);
	}

	static void print_down(const double& x, std::ostream& s) {
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
		s << conv_double::dtostring(x, s.precision(), format, -1);
	}

	static double fromstring_up(const std::string& s) {
		return conv_double::stringtod(s, 1);
	}

	static double fromstring_down(const std::string& s) {
		return conv_double::stringtod(s, -1);
	}
};

} // namespace kv

#endif // RDOUBLE_AVX512_HPP
