/*
 * Copyright (c) 2021-2024 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RFLOAT64X_HPP
#define RFLOAT64X_HPP

#include <cmath>

#if !defined(__HAVE_FLOAT64X) || __GNUC__ >= 13
#error "_Float64x is not available on this compiler"
#endif

#ifdef KV_FASTROUND
#error "KV_FASTROUND is not available for rounding mode change of _Float64x"
#endif 

#ifdef KV_NOHWROUND
#error "rounding emulation for _Float64x is not yet supported"
#endif 

#include <iostream>
#include <string>
#include <kv/conv-float64x.hpp>
#include <kv/hwround.hpp>


namespace kv {

template <> struct rop <_Float64x> {

	static _Float64x add_up(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = x, y1 = y;
		r = x1 + y1;
		return r;
	}

	static _Float64x add_down(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = -x, y1 = -y;
		r = x1 + y1;
		r = -r;
		return r;
	}

	static _Float64x sub_up(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = x, y1 = y;
		r = x1 - y1;
		return r;
	}

	static _Float64x sub_down(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = -x, y1 = -y;
		r = x1 - y1;
		r = -r;
		return r;
	}

	static _Float64x mul_up(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = x, y1 = y;
		r = x1 * y1;
		return r;
	}

	static _Float64x mul_down(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = -x, y1 = y;
		r = x1 * y1;
		r = -r;
		return r;
	}

	static _Float64x div_up(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = x, y1 = y;
		r = x1 / y1;
		return r;
	}

	static _Float64x div_down(const _Float64x& x, const _Float64x& y) {
		volatile _Float64x r, x1 = -x, y1 = y;
		r = x1 / y1;
		r = -r;
		return r;
	}

	static _Float64x sqrt_up(const _Float64x& x) {
		volatile _Float64x r, x1 = x;
		r = std::sqrt(x1);
		return r;
	}

	static _Float64x sqrt_down(const _Float64x& x) {
		volatile _Float64x r, x1 = x;
		hwround::rounddown();
		r = std::sqrt(x1);
		hwround::roundup();

		return r;
	}

	static void begin() {
		hwround::roundup();
	}

	static void end() {
		hwround::roundnear();
	}

	static void print_up(const _Float64x& x, std::ostream& s) {
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
		s << conv_float64x::float64xtostring(x, s.precision(), format, 1);
	}

	static void print_down(const _Float64x& x, std::ostream& s) {
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
		s << conv_float64x::float64xtostring(x, s.precision(), format, -1);
	}

	static _Float64x fromstring_up(const std::string& s) {
		return conv_float64x::stringtofloat64x(s, 1);
	}

	static _Float64x fromstring_down(const std::string& s) {
		return conv_float64x::stringtofloat64x(s, -1);
	}
};

template <> struct constants< interval<_Float64x> > {
	static interval<_Float64x> pi() {
		static const interval<_Float64x> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<_Float64x> e() {
		static const interval<_Float64x> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<_Float64x> ln2() {
		static const interval<_Float64x> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}
	static interval<_Float64x> str(const std::string& s) {
		return interval<_Float64x>(s, s);
	}
	static interval<_Float64x> str(const std::string& s1, const std::string& s2) {
		return interval<_Float64x>(s1, s2);
	}
};

} // namespace kv

#endif // RFLOAT64X_HPP
