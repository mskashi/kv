/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDOUBLE_HWROUND_HPP
#define RDOUBLE_HWROUND_HPP

#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <kv/hwround.hpp>
#include <kv/conv-double.hpp>

namespace kv {

template <> struct rop <double> {

	static double add_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		r = x1 + y1;
		// if (r != r) return std::numeric_limits<double>::infinity();
		return r;
	}

	static double add_down(const double& x, const double& y) {
		volatile double r, x1 = -x, y1 = -y;
		r = x1 + y1;
		r = -r;
		// if (r != r) return -std::numeric_limits<double>::infinity();
		return r;
	}

	static double sub_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		r = x1 - y1;
		// if (r != r) return std::numeric_limits<double>::infinity();
		return r;
	}

	static double sub_down(const double& x, const double& y) {
		volatile double r, x1 = -x, y1 = -y;
		r = x1 - y1;
		r = -r;
		// if (r != r) return -std::numeric_limits<double>::infinity();
		return r;
	}

	static double mul_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		r = x1 * y1;
		// if (r != r) return std::numeric_limits<double>::infinity();
		return r;
	}

	static double mul_down(const double& x, const double& y) {
		volatile double r, x1 = -x, y1 = y;
		r = x1 * y1;
		r = -r;
		// if (r != r) return -std::numeric_limits<double>::infinity();
		return r;
	}

	static double div_up(const double& x, const double& y) {
		volatile double r, x1 = x, y1 = y;
		r = x1 / y1;
		// if (r != r) return std::numeric_limits<double>::infinity();
		return r;
	}

	static double div_down(const double& x, const double& y) {
		volatile double r, x1 = -x, y1 = y;
		r = x1 / y1;
		r = -r;
		// if (r != r) return -std::numeric_limits<double>::infinity();
		return r;
	}

	static double sqrt_up(const double& x) {
		volatile double r, x1 = x;
		r = sqrt(x1);
		return r;
	}

	static double sqrt_down(const double& x) {
		volatile double r, x1 = x;
		hwround::rounddown();
		r = sqrt(x1);
		hwround::roundup();

		#if 0
		volatile double r, x1 = x, tmp;
		tmp = -sqrt(x1);
		r = -(x1 / tmp);
		#endif

		return r;
	}

	static void begin() {
		hwround::roundup();
	}

	static void end() {
		hwround::roundnear();
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

#endif // RDOUBLE_HWROUND_HPP
