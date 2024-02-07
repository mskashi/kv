/*
 * Copyright (c) 2022-2024 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDDX_HPP
#define RDDX_HPP

#include <cmath>

#if !defined(__HAVE_FLOAT64X) || __GNUC__ >= 13
#error "_Float64x is not available on this compiler"
#endif

#ifdef KV_FASTROUND
#error "KV_FASTROUND is not available for rounding mode change of _Float64x"
#endif 

#ifdef KV_NOHWROUND
#error "rounding emulation for ddx is not yet supported"
#endif 

#include <limits>
#include <kv/hwround.hpp>


namespace kv {

template <> struct rop <ddx> {

	static void twoproduct_up(const _Float64x& a, const _Float64x& b, _Float64x& x, _Float64x& y) {
		static const _Float64x th = std::ldexp((_Float64x)1., 16351);
		static const _Float64x c1 = std::ldexp((_Float64x)1., -33);
		static const _Float64x c2 = std::ldexp((_Float64x)1., 33);
		static const _Float64x th2 = std::ldexp((_Float64x)1., 16383);

		_Float64x na, nb, a1, a2, b1, b2;
		volatile _Float64x v1, v2, v3, v4, v5, v6;

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
		ddx::split(na, a1, a2);
		ddx::split(nb, b1, b2);
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

	static void twoproduct_down(const _Float64x& a, const _Float64x& b, _Float64x& x, _Float64x& y) {
		static const _Float64x th = std::ldexp((_Float64x)1., 16351);
		static const _Float64x c1 = std::ldexp((_Float64x)1., -33);
		static const _Float64x c2 = std::ldexp((_Float64x)1., 33);
		static const _Float64x th2 = std::ldexp((_Float64x)1., 16383);

		_Float64x na, nb, a1, a2, b1, b2;
		volatile _Float64x v1, v2, v3, v4, v5, v6;

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
		ddx::split(na, a1, a2);
		ddx::split(nb, b1, b2);
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

	static ddx add_up(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 + y.a1, 0.);
		}

		ddx::twosum(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			z1 = -(std::numeric_limits<ddx>::max)().a1;
			z2 = -(std::numeric_limits<ddx>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			z3 = -(std::numeric_limits<ddx>::max)().a1;
			z4 = -(std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx add_down(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 + y.a1, 0.);
		}

		ddx::twosum(x.a1, y.a1, z1, z2);

		if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			z1 = (std::numeric_limits<ddx>::max)().a1;
			z2 = (std::numeric_limits<ddx>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 + v3;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			z3 = (std::numeric_limits<ddx>::max)().a1;
			z4 = (std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx sub_up(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 - y.a1, 0.);
		}

		ddx::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			z1 = -(std::numeric_limits<ddx>::max)().a1;
			z2 = -(std::numeric_limits<ddx>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			z3 = -(std::numeric_limits<ddx>::max)().a1;
			z4 = -(std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx sub_down(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 - y.a1, 0.);
		}

		ddx::twosum(x.a1, -y.a1, z1, z2);

		if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			z1 = (std::numeric_limits<ddx>::max)().a1;
			z2 = (std::numeric_limits<ddx>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a2; v3 = y.a2;
		v1 += v2 - v3;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			z3 = (std::numeric_limits<ddx>::max)().a1;
			z4 = (std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx mul_up(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 * y.a1, 0.);
		}

		twoproduct_up(x.a1, y.a1, z1, z2);

		if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			z1 = -(std::numeric_limits<ddx>::max)().a1;
			z2 = -(std::numeric_limits<ddx>::max)().a2;
		}

		hwround::roundup();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			z3 = -(std::numeric_limits<ddx>::max)().a1;
			z4 = -(std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx mul_down(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 * y.a1, 0.);
		}

		twoproduct_down(x.a1, y.a1, z1, z2);

		if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			z1 = (std::numeric_limits<ddx>::max)().a1;
			z2 = (std::numeric_limits<ddx>::max)().a2;
		}

		hwround::rounddown();
		v1 = z2; v2 = x.a1; v3 = x.a2; v4 = y.a1; v5 = y.a2;
		v1 += v2 * v5 + v3 * v4 + v3 * v5;
		z2 = v1;
		hwround::roundnear();

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			z3 = (std::numeric_limits<ddx>::max)().a1;
			z4 = (std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx div_up(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5, v6, v7, v8;
		volatile _Float64x tmp;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 / y.a1, 0.);
		}

		z1 = x.a1 / y.a1;

		if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			z1 = -(std::numeric_limits<ddx>::max)().a1;
			z2 = -(std::numeric_limits<ddx>::max)().a2;
		}

		if (y.a1 >= 0.) {
			twoproduct_up(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
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
			if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
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

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			z3 = -(std::numeric_limits<ddx>::max)().a1;
			z4 = -(std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx div_down(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5, v6, v7, v8;
		volatile _Float64x tmp;

		if (std::fabs(x.a1) == std::numeric_limits<_Float64x>::infinity() || std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1 / y.a1, 0.);
		}

		z1 = x.a1 / y.a1;

		if (z1 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		} else if (z1 == std::numeric_limits<_Float64x>::infinity()) {
			z1 = (std::numeric_limits<ddx>::max)().a1;
			z2 = (std::numeric_limits<ddx>::max)().a2;
		}

		if (y.a1 >= 0.) {
			twoproduct_down(-z1, y.a1, z3, z4);
			if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
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
			if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
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

		ddx::twosum(z1, z2, z3, z4);

		if (z3 == -std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		} else if (z3 == std::numeric_limits<_Float64x>::infinity()) {
			z3 = (std::numeric_limits<ddx>::max)().a1;
			z4 = (std::numeric_limits<ddx>::max)().a2;
		}

		return ddx(z3, z4);
	}

	static ddx sqrt_up(const ddx& x) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5, v6;
		volatile _Float64x tmp;

		#if 0
		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
		}
		#endif

		if (x == 0.) return ddx(0.);
		if (x.a1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1, 0.);
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
			if (v5 == (std::numeric_limits<_Float64x>::max)()) {
				tmp = std::sqrt(v5*0.25 + v6*0.25)*2 + v1;
			} else {
				tmp = std::sqrt(v5 + v6) + v1;
			}
		}
		v2 /= tmp;
		z2 = v2;
		hwround::roundnear();
		ddx::twosum(z1, z2, z3, z4);

		return ddx(z3, z4);
	}

	static ddx sqrt_down(const ddx& x) {
		_Float64x z1, z2, z3, z4;
		volatile _Float64x v1, v2, v3, v4, v5, v6;
		volatile _Float64x tmp;

		#if 0
		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
		}
		#endif

		if (x == 0.) return ddx(0.);
		if (x.a1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1, 0.);
		}

		z1 = std::sqrt(x.a1);
		twoproduct_down(-z1, z1, z3, z4);
		// z2 = ((z3 + x.a1) + x.a2 + z4) / (sqrt(x.a1 + x.a2) + z1);
		hwround::rounddown();
		v1 = z1; v3 = z3; v4 = z4; v5 = x.a1; v6 = x.a2;
		v2 = (v3 + v5) + v6 + v4;
		if (v2 > 0.) {
			hwround::roundup();
			if (v5 == (std::numeric_limits<_Float64x>::max)()) {
				tmp = std::sqrt(v5*0.25 + v6*0.25)*2 + v1;
			} else {
				tmp = std::sqrt(v5 + v6) + v1;
			}
			hwround::rounddown();
		} else {
			tmp = std::sqrt(v5 + v6) + v1;
		}
		v2 /= tmp;
		z2 = v2;
		hwround::roundnear();
		ddx::twosum(z1, z2, z3, z4);

		return ddx(z3, z4);
	}

	static void begin() {
	}

	static void end() {
	}

	static void print_up(const ddx& x, std::ostream& s) {
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
		s << conv_ddx::ddxtostring(x.a1, x.a2, s.precision(), format, 1);
	}

	static void print_down(const ddx& x, std::ostream& s) {
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
		s << conv_ddx::ddxtostring(x.a1, x.a2, s.precision(), format, -1);
	}

	static ddx fromstring_up(const std::string& s) {
		_Float64x x1, x2;
		conv_ddx::stringtoddx(s, x1, x2, 1);
		return ddx(x1, x2);
	}

	static ddx fromstring_down(const std::string& s) {
		_Float64x x1, x2;
		conv_ddx::stringtoddx(s, x1, x2, -1);
		return ddx(x1, x2);
	}
};

template <> struct constants< interval<ddx> > {
	static interval<ddx> pi() {
		static const interval<ddx> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<ddx> e() {
		static const interval<ddx> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<ddx> ln2() {
		static const interval<ddx> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}

	static interval<ddx> str(const std::string& s) {
		return interval<ddx>(s, s);
	}

	static interval<ddx> str(const std::string& s1, const std::string& s2) {
		return interval<ddx>(s1, s2);
	}
};

} // namespace kv

#endif // RDDX_HPP
