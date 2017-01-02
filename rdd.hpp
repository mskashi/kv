#ifndef RDD_HPP
#define RDD_HPP

#include "hwround.hpp"


namespace kv {

template <> struct rop <dd> {

	static void safe_split(const double& a, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);

		if (-th <= a && a <= th) {
			dd::split(a, x, y);
		} else {
			double a1 = a * c1;
			dd::split(a1, x, y);
			x *= c2;
			y *= c2;
		}
	}

	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		double a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;
		x = a * b;
		safe_split(a, a1, a2);
		safe_split(b, b1, b2);
		hwround::roundup();
		v1 = x; v3 = a1; v4 = a2; v5 =b1; v6 = b2; 
		v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) + v4 * v6;
		y = v2;
		hwround::roundnear();
	}

	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		double a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;
		x = a * b;
		safe_split(a, a1, a2);
		safe_split(b, b1, b2);
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
		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		#endif

		return dd(z3, z4);

	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, y.a1, z1, z2);
		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}
		#endif

		return dd(z3, z4);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);
		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		#endif

		return dd(z3, z4);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		dd::twosum(x.a1, -y.a1, z1, z2);
		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}
		#endif

		return dd(z3, z4);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_up(x.a1, y.a1, z1, z2);
		hwround::roundup();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		#endif

		return dd(z3, z4);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		twoproduct_down(x.a1, y.a1, z1, z2);
		hwround::rounddown();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}
		#endif

		return dd(z3, z4);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

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

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return -(std::numeric_limits<dd>::max)();
		}
		#endif

		return dd(z3, z4);
	}

	static dd div_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		z1 = x.a1 / y.a1;

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

		#if defined(DD_INFINITY)
		if (z3 == std::numeric_limits<double>::infinity()) {
			return (std::numeric_limits<dd>::max)();
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}
		#endif

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

	static void finish() {
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
		s << dd::ddtostring(x, s.precision(), format, 1);
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
		s << dd::ddtostring(x, s.precision(), format, -1);
	}

	static dd fromstring_up(const std::string& s) {
		return dd::stringtodd(s, 1);
	}

	static dd fromstring_down(const std::string& s) {
		return dd::stringtodd(s, -1);
	}
};

} // namespace kv

#endif // RDD_HPP
