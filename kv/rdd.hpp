/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_HPP
#define RDD_HPP

#include <kv/hwround.hpp>


namespace kv {

template <> struct rop <dd> {

	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);

		double na, nb, a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;

		x = a * b;
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
		if (std::fabs(a) > th) {
			na = a * c1;
			nb = b * c2;
		} else if (std::fabs(b) > th) {
			na = a * c2;
			nb = b * c1;
		} else {
			na = a;
			nb = b;
		}
		dd::split(na, a1, a2);
		dd::split(nb, b1, b2);
		hwround::roundup();
		v1 = x; v3 = a1; v4 = a2; v5 =b1; v6 = b2; 
		v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) + v4 * v6;
		y = v2;
		hwround::roundnear();
	}

	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);

		double na, nb, a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;

		x = a * b;
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
		if (std::fabs(a) > th) {
			na = a * c1;
			nb = b * c2;
		} else if (std::fabs(b) > th) {
			na = a * c2;
			nb = b * c1;
		} else {
			na = a;
			nb = b;
		}
		dd::split(na, a1, a2);
		dd::split(nb, b1, b2);
		hwround::rounddown();
		v1 = x; v3 = a1; v4 = a2; v5 =b1; v6 = b2; 
		v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) + v4 * v6;
		y = v2;
		hwround::roundnear();
	}

	static dd add_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		if (z1 != z1) return std::numeric_limits<dd>::infinity();

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);

	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (z1 != z1) return -std::numeric_limits<dd>::infinity();

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		if (z1 != z1) return std::numeric_limits<dd>::infinity();

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (z1 != z1) return -std::numeric_limits<dd>::infinity();

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_up(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		if (z1 != z1) return std::numeric_limits<dd>::infinity();

		hwround::roundup();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);


		return dd(z3, z4);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_down(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (z1 != z1) return -std::numeric_limits<dd>::infinity();

		hwround::rounddown();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		if (z1 != z1) return std::numeric_limits<dd>::infinity();
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			hwround::rounddown();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			hwround::roundup();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		} else {
			twoproduct_down(-z1, y.a1, z3, z4);
			hwround::roundup();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			hwround::rounddown();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		}
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd div_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (z1 != z1) return -std::numeric_limits<dd>::infinity();
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			hwround::roundup();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			hwround::rounddown();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		} else {
			twoproduct_up(-z1, y.a1, z3, z4);
			hwround::rounddown();
			v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
			tmp =  v7 + v8;
			hwround::roundup();
			v2 = ((((v3 + v5) + (-v1) * v8) + v6) + v4) / tmp;
			z2 = v2;
		}
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sqrt_up(const dd& x) {
		dd r, r2;

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_up(x, r);

		if (r > r2) return r;
		else return r2;
	}

	static dd sqrt_down(const dd& x) {
		dd r, r2;

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_down(x, r);

		if (r > r2) return r2;
		else return r;
	}

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

#endif // RDD_HPP
