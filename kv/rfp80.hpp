/*
 * Copyright (c) 2021-2026 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RFP80_HPP
#define RFP80_HPP

#include <kv/fp80.hpp>

#ifdef KV_HAVE_FP80

#ifdef KV_FASTROUND
#error "KV_FASTROUND is not available for rounding mode change of fp80"
#endif 

#ifdef KV_NOHWROUND
#error "rounding emulation for fp80 is not yet supported"
#endif 

#include <iostream>
#include <cmath>
#include <string>
#include <kv/conv-fp80.hpp>
#include <kv/hwround.hpp>


namespace kv {

template <> struct rop <fp80> {

	static fp80 add_up(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = x, y1 = y;
		r = x1 + y1;
		return r;
	}

	static fp80 add_down(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = -x, y1 = -y;
		r = x1 + y1;
		r = -r;
		return r;
	}

	static fp80 sub_up(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = x, y1 = y;
		r = x1 - y1;
		return r;
	}

	static fp80 sub_down(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = -x, y1 = -y;
		r = x1 - y1;
		r = -r;
		return r;
	}

	static fp80 mul_up(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = x, y1 = y;
		r = x1 * y1;
		return r;
	}

	static fp80 mul_down(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = -x, y1 = y;
		r = x1 * y1;
		r = -r;
		return r;
	}

	static fp80 div_up(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = x, y1 = y;
		r = x1 / y1;
		return r;
	}

	static fp80 div_down(const fp80& x, const fp80& y) {
		volatile fp80 r, x1 = -x, y1 = y;
		r = x1 / y1;
		r = -r;
		return r;
	}

	static fp80 sqrt_up(const fp80& x) {
		volatile fp80 r, x1 = x;
		r = std::sqrt(x1);
		return r;
	}

	static fp80 sqrt_down(const fp80& x) {
		volatile fp80 r, x1 = x;
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

	static void print_up(const fp80& x, std::ostream& s) {
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
		s << conv_fp80::fp80tostring(x, s.precision(), format, 1);
	}

	static void print_down(const fp80& x, std::ostream& s) {
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
		s << conv_fp80::fp80tostring(x, s.precision(), format, -1);
	}

	static fp80 fromstring_up(const std::string& s) {
		return conv_fp80::stringtofp80(s, 1);
	}

	static fp80 fromstring_down(const std::string& s) {
		return conv_fp80::stringtofp80(s, -1);
	}
};

template <> struct constants< interval<fp80> > {
	static interval<fp80> pi() {
		static const interval<fp80> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<fp80> e() {
		static const interval<fp80> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<fp80> ln2() {
		static const interval<fp80> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}
	static interval<fp80> str(const std::string& s) {
		return interval<fp80>(s, s);
	}
	static interval<fp80> str(const std::string& s1, const std::string& s2) {
		return interval<fp80>(s1, s2);
	}
};

} // namespace kv

#endif // KV_HAVE_FP80

#endif // RFP80_HPP
