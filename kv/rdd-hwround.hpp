/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_HWROUND_HPP
#define RDD_HWROUND_HPP

#include <cmath>
#include <limits>
#include <kv/hwround.hpp>


#ifndef DD_NEW_SQRT
#define DD_NEW_SQRT 1
#endif

#ifndef KV_USE_TPFMA
#define KV_USE_TPFMA 0
#endif


namespace kv {

template <> struct rop <dd> {

#if KV_USE_TPFMA == 1
	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		volatile double v1, v2, v3, v4;
		x = a * b;
		hwround::roundup();
		v1 = a; v2 = b; v3 = x;
		// y = fma(a, b, -x);
		v4 = fma(v1, v2, -v3);
		y = v4;
		hwround::roundnear();
	}
#else // KV_USE_TPFMA
	static void twoproduct_up(const double& a, const double& b, double& x, double& y) {
		static const double th = std::ldexp(1., 996);
		static const double c1 = std::ldexp(1., -28);
		static const double c2 = std::ldexp(1., 28);
		static const double th2 = std::ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;

		x = a * b;
		#if 0
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
		#endif
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
		v1 = x; v3 = a1; v4 = a2; v5 = b1; v6 = b2; 
		if (std::fabs(x) > th2) {
			v2 = ((((v3 * 0.5) * v5 - (v1 * 0.5)) * 2. + v4 * v5) + v3 * v6) + v4 * v6;
		} else {
			v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) + v4 * v6;
		}
		y = v2;
		hwround::roundnear();
	}
#endif // KV_USE_TPFMA

#if KV_USE_TPFMA == 1
	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		volatile double v1, v2, v3, v4;
		x = a * b;
		hwround::rounddown();
		v1 = a; v2 = b; v3 = x;
		// y = fma(a, b, -x);
		v4 = fma(v1, v2, -v3);
		y = v4;
		hwround::roundnear();
	}
#else // KV_USE_TPFMA
	static void twoproduct_down(const double& a, const double& b, double& x, double& y) {
		static const double th = std::ldexp(1., 996);
		static const double c1 = std::ldexp(1., -28);
		static const double c2 = std::ldexp(1., 28);
		static const double th2 = std::ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;
		volatile double v1, v2, v3, v4, v5, v6;

		x = a * b;
		#if 0
		if (std::fabs(x) == std::numeric_limits<double>::infinity()) {
			y = 0.;
			return;
		}
		#endif
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
		v1 = x; v3 = a1; v4 = a2; v5 = b1; v6 = b2; 
		if (std::fabs(x) > th2) {
			v2 = ((((v3 * 0.5) * v5 - (v1 * 0.5)) * 2. + v4 * v5) + v3 * v6) + v4 * v6;
		} else {
			v2 = (((v3 * v5 - v1) + v4 * v5) + v3 * v6) + v4 * v6;
		}
		y = v2;
		hwround::roundnear();
	}
#endif // KV_USE_TPFMA

	static dd add_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 + y.a1, 0.);
		}

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			z1 = -(std::numeric_limits<dd>::max)().a1;
			z2 = -(std::numeric_limits<dd>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			z3 = -(std::numeric_limits<dd>::max)().a1;
			z4 = -(std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd add_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 + y.a1, 0.);
		}

		dd::twosum(x.a1, y.a1, z1, z2);

		if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == std::numeric_limits<double>::infinity()) {
			z1 = (std::numeric_limits<dd>::max)().a1;
			z2 = (std::numeric_limits<dd>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == std::numeric_limits<double>::infinity()) {
			z3 = (std::numeric_limits<dd>::max)().a1;
			z4 = (std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd sub_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 - y.a1, 0.);
		}

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			z1 = -(std::numeric_limits<dd>::max)().a1;
			z2 = -(std::numeric_limits<dd>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			z3 = -(std::numeric_limits<dd>::max)().a1;
			z4 = -(std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd sub_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 - y.a1, 0.);
		}

		dd::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == std::numeric_limits<double>::infinity()) {
			z1 = (std::numeric_limits<dd>::max)().a1;
			z2 = (std::numeric_limits<dd>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == std::numeric_limits<double>::infinity()) {
			z3 = (std::numeric_limits<dd>::max)().a1;
			z4 = (std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd mul_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 * y.a1, 0.);
		}

		twoproduct_up(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			z1 = -(std::numeric_limits<dd>::max)().a1;
			z2 = -(std::numeric_limits<dd>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			z3 = -(std::numeric_limits<dd>::max)().a1;
			z4 = -(std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd mul_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 * y.a1, 0.);
		}

		twoproduct_down(x.a1, y.a1, z1, z2);

		if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == std::numeric_limits<double>::infinity()) {
			z1 = (std::numeric_limits<dd>::max)().a1;
			z2 = (std::numeric_limits<dd>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == std::numeric_limits<double>::infinity()) {
			z3 = (std::numeric_limits<dd>::max)().a1;
			z4 = (std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd div_up(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 / y.a1, 0.);
		}

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == -std::numeric_limits<double>::infinity()) {
			z1 = -(std::numeric_limits<dd>::max)().a1;
			z2 = -(std::numeric_limits<dd>::max)().a2;
		}

		if (y.a1 >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);
				hwround::roundup();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v5 *= 0.5; v6 *= 0.5; v7 *= 0.5; v8 *= 0.5;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					hwround::rounddown();
					tmp =  v7 + v8;
					hwround::roundup();
				} else {
					tmp =  v7 + v8;
				}
				hwround::roundup();
				v2 /= tmp;
				z2 = v2;
			} else {
				hwround::roundup();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					hwround::rounddown();
					tmp =  v7 + v8;
					hwround::roundup();
				} else {
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			}
		} else {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);
				hwround::rounddown();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v5 *= 0.5; v6 *= 0.5; v7 *= 0.5; v8 *= 0.5;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					tmp =  v7 + v8;
					hwround::roundup();
				} else {
					hwround::roundup();
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			} else {
				hwround::rounddown();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					tmp =  v7 + v8;
					hwround::roundup();
				} else {
					hwround::roundup();
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			}
		}
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == -std::numeric_limits<double>::infinity()) {
			z3 = -(std::numeric_limits<dd>::max)().a1;
			z4 = -(std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

	static dd div_down(const dd& x, const dd& y) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6, v7, v8;
		volatile double tmp;

		if (std::fabs(x.a1) == std::numeric_limits<double>::infinity() || std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(x.a1 / y.a1, 0.);
		}

		z1 = x.a1 / y.a1;

		if (z1 == -std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		} else if (z1 == std::numeric_limits<double>::infinity()) {
			z1 = (std::numeric_limits<dd>::max)().a1;
			z2 = (std::numeric_limits<dd>::max)().a2;
		}

		if (y.a1 >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_down(-z1, y.a1 * 0.5, z3, z4);
				hwround::rounddown();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v5 *= 0.5; v6 *= 0.5; v7 *= 0.5; v8 *= 0.5;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					hwround::roundup();
					tmp =  v7 + v8;
					hwround::rounddown();
				} else {
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			} else {
				hwround::rounddown();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					hwround::roundup();
					tmp =  v7 + v8;
					hwround::rounddown();
				} else {
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			}
		} else {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
				twoproduct_up(-z1, y.a1 * 0.5, z3, z4);
				hwround::roundup();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v5 *= 0.5; v6 *= 0.5; v7 *= 0.5; v8 *= 0.5;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					tmp =  v7 + v8;
					hwround::rounddown();
				} else {
					hwround::rounddown();
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			} else {
				hwround::roundup();
				v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2; v7 = y.a1; v8 = y.a2;
				v2 = (((v3 + v5) + (-v1) * v8) + v6) + v4;
				if (v2 > 0.) {
					tmp =  v7 + v8;
					hwround::rounddown();
				} else {
					hwround::rounddown();
					tmp =  v7 + v8;
				}
				v2 /= tmp;
				z2 = v2;
			}
		}
		hwround::roundnear();

		dd::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		} else if (z3 == std::numeric_limits<double>::infinity()) {
			z3 = (std::numeric_limits<dd>::max)().a1;
			z4 = (std::numeric_limits<dd>::max)().a2;
		}

		return dd(z3, z4);
	}

#if DD_NEW_SQRT == 1
	static dd sqrt_up(const dd& x) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6;
		volatile double tmp;

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
		// z2 = ((z3 + x.a1) + x.a2 + z4) / (sqrt(x.a1 + x.a2) + z1);
		hwround::roundup();
		v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2;
		v2 = (v3 + v5) + v6 + v4;
		if (v2 > 0.) {
			hwround::rounddown();
			tmp = std::sqrt(v5 + v6) + v1;
			hwround::roundup();
		} else {
			tmp = std::sqrt(v5 + v6) + v1;
		}
		v2 /= tmp;
		z2 = v2;
		hwround::roundnear();
		dd::twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	static dd sqrt_down(const dd& x) {
		double z1, z2, z3, z4;
		volatile double v1, v2, v3, v4, v5, v6;
		volatile double tmp;

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
		// z2 = ((z3 + x.a1) + x.a2 + z4) / (sqrt(x.a1 + x.a2) + z1);
		hwround::rounddown();
		v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2;
		v2 = (v3 + v5) + v6 + v4;
		if (v2 > 0.) {
			hwround::roundup();
			tmp = std::sqrt(v5 + v6) + v1;
			hwround::rounddown();
		} else {
			tmp = std::sqrt(v5 + v6) + v1;
		}
		v2 /= tmp;
		z2 = v2;
		hwround::roundnear();
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

#endif // RDD_HWROUND_HPP
