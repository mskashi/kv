/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_NOHWROUND_HPP
#define RDD_NOHWROUND_HPP

#include <kv/rdouble-nohwround.hpp>

namespace kv {

template <> struct rop <dd> {

	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		static const double th = std::ldexp(1., 996);
		static const double c1 = std::ldexp(1., -28);
		static const double c2 = std::ldexp(1., 28);
		static const double th2 = std::ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;

		x = a * b;
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
		if (std::fabs(x) > th2) {
			// y = ((((a1 * 0.5) * b1 - (x * 0.5)) * 2. + a2 * b1) + a1 * b2) + a2 * b2;
			y = rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::sub_up(rop<double>::mul_up(a1 * 0.5, b1), x * 0.5) * 2., rop<double>::mul_up(a2, b1)), rop<double>::mul_up(a1, b2)), rop<double>::mul_up(a2, b2));
		} else {
			// y = (((a1 * b1 - x) + a2 * b1) + a1 * b2) + a2 * b2;
			y = rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::sub_up(rop<double>::mul_up(a1, b1), x), rop<double>::mul_up(a2, b1)), rop<double>::mul_up(a1, b2)), rop<double>::mul_up(a2, b2));
		}
	}

	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		static const double th = std::ldexp(1., 996);
		static const double c1 = std::ldexp(1., -28);
		static const double c2 = std::ldexp(1., 28);
		static const double th2 = std::ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;

		x = a * b;
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
		if (std::fabs(x) > th2) {
			// y = ((((a1 * 0.5) * b1 - (x * 0.5)) * 2. + a2 * b1) + a1 * b2) + a2 * b2;
			y = rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::sub_down(rop<double>::mul_down(a1 * 0.5, b1), x * 0.5) * 2., rop<double>::mul_down(a2, b1)), rop<double>::mul_down(a1, b2)), rop<double>::mul_down(a2, b2));
		} else {
			// y = (((a1 * b1 - x) + a2 * b1) + a1 * b2) + a2 * b2;
			y = rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::sub_down(rop<double>::mul_down(a1, b1), x), rop<double>::mul_down(a2, b1)), rop<double>::mul_down(a1, b2)), rop<double>::mul_down(a2, b2));
		}
	}

	static dd add_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		// z2 += x.a2 + y.a2;
		z2 = rop<double>::add_up(z2, rop<double>::add_up(x.a2, y.a2));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);

	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		// z2 += x.a2 + y.a2;
		z2 = rop<double>::add_down(z2, rop<double>::add_down(x.a2, y.a2));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		// z2 += x.a2 - y.a2;
		z2 = rop<double>::add_up(z2, rop<double>::sub_up(x.a2, y.a2));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		// z2 += x.a2 - y.a2;
		z2 = rop<double>::add_down(z2, rop<double>::sub_down(x.a2, y.a2));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		twoproduct_up(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}

		// z2 += x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		z2 = rop<double>::add_up(z2, rop<double>::add_up(rop<double>::add_up(rop<double>::mul_up(x.a1, y.a2), rop<double>::mul_up(x.a2, y.a1)), rop<double>::mul_up(x.a2, y.a2)));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		twoproduct_down(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		// z2 += x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		z2 = rop<double>::add_down(z2, rop<double>::add_down(rop<double>::add_down(rop<double>::mul_down(x.a1, y.a2), rop<double>::mul_down(x.a2, y.a1)), rop<double>::mul_down(x.a2, y.a2)));
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);
				// z2 = ((((z3 + x.a1 * 0.5) + (-z1) * (y.a2 * 0.5)) + x.a2 * 0.5) + z4) / (y.a1 * 0.5 + y.a2 * 0.5);
				z2 = rop<double>::div_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(z3, x.a1 * 0.5), rop<double>::mul_up(-z1, y.a2 * 0.5)), x.a2 * 0.5), z4), rop<double>::add_down(y.a1 * 0.5, y.a2 * 0.5));
			} else {
				// z2 = ((((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
				z2 = rop<double>::div_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(z3, x.a1), rop<double>::mul_up(-z1, y.a2)), x.a2), z4), rop<double>::add_down(y.a1, y.a2));
			}
		} else {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);
				// z2 = ((((z3 + x.a1 * 0.5) + (-z1) * (y.a2 * 0.5)) + x.a2 * 0.5) + z4) / (y.a1 * 0.5 + y.a2 * 0.5);
				z2 = rop<double>::div_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(z3, x.a1 * 0.5), rop<double>::mul_down(-z1, y.a2 * 0.5)), x.a2 * 0.5), z4), rop<double>::add_up(y.a1 * 0.5, y.a2 * 0.5));
			} else {
				// z2 = ((((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
				z2 = rop<double>::div_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(z3, x.a1), rop<double>::mul_down(-z1, y.a2)), x.a2), z4), rop<double>::add_up(y.a1, y.a2));
			}
		}
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
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		if (y >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);
				// z2 = ((((z3 + x.a1 * 0.5) + (-z1) * (y.a2 * 0.5)) + x.a2 * 0.5) + z4) / (y.a1 * 0.5 + y.a2 * 0.5);
				z2 = rop<double>::div_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(z3, x.a1 * 0.5), rop<double>::mul_down(-z1, y.a2 * 0.5)), x.a2 * 0.5), z4), rop<double>::add_up(y.a1 * 0.5, y.a2 * 0.5));
			} else {
				// z2 = ((((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
				z2 = rop<double>::div_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(rop<double>::add_down(z3, x.a1), rop<double>::mul_down(-z1, y.a2)), x.a2), z4), rop<double>::add_up(y.a1, y.a2));
			}
		} else {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);
				// z2 = ((((z3 + x.a1 * 0.5) + (-z1) * (y.a2 * 0.5)) + x.a2 * 0.5) + z4) / (y.a1 * 0.5 + y.a2 * 0.5);
				z2 = rop<double>::div_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(z3, x.a1 * 0.5), rop<double>::mul_up(-z1, y.a2 * 0.5)), x.a2 * 0.5), z4), rop<double>::add_down(y.a1 * 0.5, y.a2 * 0.5));
			} else {
				// z2 = ((((z3 + x.a1) + (-z1) * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
				z2 = rop<double>::div_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(rop<double>::add_up(z3, x.a1), rop<double>::mul_up(-z1, y.a2)), x.a2), z4), rop<double>::add_down(y.a1, y.a2));
			}
		}
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sqrt_up(const dd& x) {
		dd r, r2;

		if (x == 0.) return dd(0.);
		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		r2 = div_up(x, r);

		if (r > r2) return r;
		else return r2;
	}

	static dd sqrt_down(const dd& x) {
		dd r, r2;

		if (x == 0.) return dd(0.);
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

#endif // RDD_NOHWROUND_HPP
