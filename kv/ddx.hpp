/*
 * Copyright (c) 2021-2024 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DDX_HPP
#define DDX_HPP

#include <cmath>

#if !defined(__HAVE_FLOAT64X) || __GNUC__ >= 13
#error "_Float64x is not available on this compiler"
#endif

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

#include <kv/convert.hpp>
#include <kv/constants.hpp>
#include <kv/conv-ddx.hpp>


#ifndef DD_FASTMULT
#define DD_FASTMULT 0
#endif


namespace kv {


class ddx;

template <class C> struct convertible<C, ddx> {
	static const bool value = boost::is_arithmetic<C>::value || boost::is_same<C, ddx>::value || boost::is_convertible<C, std::string>::value;
};
template <class C> struct acceptable_n<C, ddx> {
	static const bool value = boost::is_arithmetic<C>::value;
};
template <class C> struct acceptable_s<C, ddx> {
	static const bool value = boost::is_convertible<C, std::string>::value;
};


class ddx {
	public:

	_Float64x a1;
	_Float64x a2;


	static void fasttwosum(const _Float64x& a, const _Float64x& b, _Float64x& x, _Float64x& y) {
		_Float64x tmp;
		x = a + b;
		tmp = x - a;
		y = b - tmp;
	}

	static void twosum(const _Float64x& a, const _Float64x& b, _Float64x& x, _Float64x& y) {
		_Float64x tmp;

		x = a + b;
		if (std::fabs(a) > std::fabs(b)) {
			tmp = x - a;
			y = b - tmp;
		} else {
			tmp = x - b;
			y = a - tmp;
		}
	}

	static void split(const _Float64x& a, _Float64x& x, _Float64x& y) {
		static const _Float64x sigma = std::ldexp((_Float64x)1., 32) + 1;
		_Float64x tmp;

		tmp = a * sigma;
		x = tmp - (tmp - a);
		y = a - x;
	}

	static void twoproduct(const _Float64x& a, const _Float64x& b, _Float64x& x, _Float64x& y) {
		static const _Float64x th = std::ldexp((_Float64x)1., 16351);
		static const _Float64x c1 = std::ldexp((_Float64x)1., -33);
		static const _Float64x c2 = std::ldexp((_Float64x)1., 33);
		static const _Float64x th2 = std::ldexp((_Float64x)1., 16383);

		_Float64x na, nb, a1, a2, b1, b2;

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
		split(na, a1, a2);
		split(nb, b1, b2);
		if (std::fabs(x) > th2) {
			y = ((((a1 * 0.5) * b1 - (x * 0.5)) * 2 + a2 * b1) + a1 * b2) + a2 * b2;
		} else {
			y = (((a1 * b1 - x) + a2 * b1) + a1 * b2) + a2 * b2;
		}
	}

	ddx() {
		a1 = 0.;
		a2 = 0.;
	}

	template <class C> explicit ddx(const C& x, typename boost::enable_if_c< acceptable_n<C, ddx>::value >::type* =0) {
		a1 = x;
		a2 = 0.;
	}

	template <class C1, class C2> ddx(const C1& x, const C2& y, typename boost::enable_if_c< acceptable_n<C1, ddx>::value && acceptable_n<C2, ddx>::value >::type* =0) {
		a1 = x;
		a2 = y;
	}

	template <class C> explicit ddx(const C& x, typename boost::enable_if_c< acceptable_s<C, ddx>::value >::type* =0) {
		conv_ddx::stringtoddx(x, a1, a2);
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx& >::type operator=(const C& x) {
		a1 = x;
		a2 = 0.;
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx& >::type operator=(const C& x) {
		conv_ddx::stringtoddx(x, a1, a2);
		return *this;
	}

	friend ddx operator+(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x.a1, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x.a2 + y.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator+(const ddx& x, const C& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x.a1, y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator+(const ddx& x, const C& y) {
		return x + ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator+(const C& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += y.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator+(const C& x, const ddx& y) {
		return ddx(x) + y;
	}

	friend ddx& operator+=(ddx& x, const ddx& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx& >::type operator+=(ddx& x, const C& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx& >::type operator+=(ddx& x, const C& y) {
		x = x + ddx(y);
		return x;
	}

	friend ddx operator-(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x.a1, -y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x.a2 - y.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator-(const ddx& x, const C& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x.a1, -y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator-(const ddx& x, const C& y) {
		return x - ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator-(const C& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twosum(x, -y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 -= y.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator-(const C& x, const ddx& y) {
		return ddx(x) - y;
	}

	friend ddx& operator-=(ddx& x, const ddx& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx& >::type operator-=(ddx& x, const C& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx& >::type operator-=(ddx& x, const C& y) {
		x = x - ddx(y);
		return x;
	}

	friend ddx operator-(const ddx& x) {
		ddx r;

		r.a1 = - x.a1;
		r.a2 = - x.a2;

		return r;
	}

	friend ddx operator*(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twoproduct(x.a1, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}

		// x.a2 * y.a2 is very small but sometimes important
		#if DD_FASTMULT == 1
		z2 += x.a1 * y.a2 + x.a2 * y.a1;
		#else
		z2 += x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		#endif

		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator*(const ddx& x, const C& y) {
		_Float64x z1, z2, z3, z4;

		twoproduct(x.a1, y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x.a2 * y;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator*(const ddx& x, const C& y) {
		return x * ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator*(const C& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		twoproduct(x, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		z2 += x * y.a2;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator*(const C& x, const ddx& y) {
		return ddx(x) * y;
	}

	friend ddx& operator*=(ddx& x, const ddx& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx& >::type operator*=(ddx& x, const C& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx& >::type operator*=(ddx& x, const C& y) {
		x = x * ddx(y);
		return x;
	}

	friend ddx operator/(const ddx& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		z1 = x.a1 / y.a1;
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		if (std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}

		twoproduct(-z1, y.a1, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			twoproduct(-z1, y.a1 * 0.5, z3, z4);
			z2 = ((((z3 + (x.a1 * 0.5)) - z1 * (y.a2 * 0.5)) + (x.a2 * 0.5)) + z4) / (y.a1 * 0.5);
		} else {
			// z2 = ((((z3 + x.a1) - z1 * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
			z2 = ((((z3 + x.a1) - z1 * y.a2) + x.a2) + z4) / y.a1;
		}
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator/(const ddx& x, const C& y) {
		_Float64x z1, z2, z3, z4;

		z1 = x.a1 / y;
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		if (std::fabs((_Float64x)y) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}

		twoproduct(-z1, y, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			twoproduct(-z1, y * 0.5, z3, z4);
			z2 = (((z3 + (x.a1 * 0.5)) + (x.a2 * 0.5)) + z4) / (y * 0.5);
		} else {
			z2 = (((z3 + x.a1) + x.a2) + z4) / y;
		}
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator/(const ddx& x, const C& y) {
		return x / ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type operator/(const C& x, const ddx& y) {
		_Float64x z1, z2, z3, z4;

		z1 = x / y.a1;
		if (std::fabs(z1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}
		if (std::fabs(y.a1) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z1, 0.);
		}

		twoproduct(-z1, y.a1, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			z2 = (((z3 + (x * 0.5)) - z1 * (y.a2 * 0.5)) + z4) / (y.a1 * 0.5);
		} else {
			z2 = (((z3 + x) - z1 * y.a2) + z4) / y.a1;
		}
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(z3, 0.);
		}

		return ddx(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type operator/(const C& x, const ddx& y) {
		return ddx(x) / y;
	}

	friend ddx& operator/=(ddx& x, const ddx& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx& >::type operator/=(ddx& x, const C& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx& >::type operator/=(ddx& x, const C& y) {
		x = x / ddx(y);
		return x;
	}

	friend std::ostream& operator<<(std::ostream& s, const ddx& x) {
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
		s << conv_ddx::ddxtostring(x.a1, x.a2, s.precision(), format, 0);
		return s;
	}

	friend ddx sqrt(const ddx& x) {
		_Float64x z1, z2, z3, z4;

		if (x < 0.) {
			throw std::domain_error("ddx: sqrt of negative value");
                }

		if (x == 0.) return ddx(0.);
		if (x.a1 == std::numeric_limits<_Float64x>::infinity()) {
			return ddx(x.a1, 0.);
		}

		z1 = std::sqrt(x.a1);
		twoproduct(-z1, z1, z3, z4);
		z2 = ((z3 + x.a1) + x.a2 + z4) / (2 * z1);
		twosum(z1, z2, z3, z4);

		return ddx(z3, z4);
	}

	friend ddx abs(const ddx& x) {
		if (x.a1 >= 0.) {
			return x;
		} else {
			return ddx(-x.a1, -x.a2);
		}
	}

	friend bool operator<(const ddx& x, const ddx& y) {
		if (x.a1 < y.a1) return true;
		else if (x.a1 > y.a1) return false;
		else return (x.a2 < y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator<(const ddx& x, const C& y) {
		if (x.a1 < y) return true;
		else if (x.a1 > y) return false;
		else return (x.a2 < 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator<(const ddx& x, const C& y) {
		return x < ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator<(const C& x, const ddx& y) {
		if (x < y.a1) return true;
		else if (x > y.a1) return false;
		else return (0. < y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator<(const C& x, const ddx& y) {
		return ddx(x) < y;
	}

	friend bool operator<=(const ddx& x, const ddx& y) {
		if (x.a1 < y.a1) return true;
		else if (x.a1 > y.a1) return false;
		else return (x.a2 <= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator<=(const ddx& x, const C& y) {
		if (x.a1 < y) return true;
		else if (x.a1 > y) return false;
		else return (x.a2 <= 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator<=(const ddx& x, const C& y) {
		return x <= ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator<=(const C& x, const ddx& y) {
		if (x < y.a1) return true;
		else if (x > y.a1) return false;
		else return (0. <= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator<=(const C& x, const ddx& y) {
		return ddx(x) <= y;
	}

	friend bool operator>(const ddx& x, const ddx& y) {
		if (x.a1 > y.a1) return true;
		else if (x.a1 < y.a1) return false;
		else return (x.a2 > y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator>(const ddx& x, const C& y) {
		if (x.a1 > y) return true;
		else if (x.a1 < y) return false;
		else return (x.a2 > 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator>(const ddx& x, const C& y) {
		return x > ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator>(const C& x, const ddx& y) {
		if (x > y.a1) return true;
		else if (x < y.a1) return false;
		else return (0. > y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator>(const C& x, const ddx& y) {
		return ddx(x) > y;
	}

	friend bool operator>=(const ddx& x, const ddx& y) {
		if (x.a1 > y.a1) return true;
		else if (x.a1 < y.a1) return false;
		else return (x.a2 >= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator>=(const ddx& x, const C& y) {
		if (x.a1 > y) return true;
		else if (x.a1 < y) return false;
		else return (x.a2 >= 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator>=(const ddx& x, const C& y) {
		return x >= ddx(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator>=(const C& x, const ddx& y) {
		if (x > y.a1) return true;
		else if (x < y.a1) return false;
		else return (0. >= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, bool >::type operator>=(const C& x, const ddx& y) {
		return ddx(x) >= y;
	}

	friend bool operator==(const ddx& x, const ddx& y) {
		return (x.a1 == y.a1) && (x.a2 == y.a2);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator==(const ddx& x, const C& y) {
		return (x.a1 == y) && (x.a2 == 0.);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator==(const C& x, const ddx& y) {
		return (x == y.a1) && (0. == y.a2);
	}

	friend bool operator!=(const ddx& x, const ddx& y) {
		return (x.a1 != y.a1) || (x.a2 != y.a2);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator!=(const ddx& x, const C& y) {
		return (x.a1 != y) || (x.a2 != 0.);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, bool >::type operator!=(const C& x, const ddx& y) {
		return (x != y.a1) || (0. != y.a2);
	}

	friend ddx floor(const ddx& x) {
		_Float64x z1, z2, z3, z4;
		z1 = std::floor(x.a1);
		if (z1 != x.a1) {
			return ddx(z1, 0.);
		} else {
			z2 = std::floor(x.a2);
			twosum(z1, z2, z3, z4);
			return ddx(z3, z4);
		}
	}

	friend ddx ceil(const ddx& x) {
		_Float64x z1, z2, z3, z4;
		z1 = std::ceil(x.a1);
		if (z1 != x.a1) {
			return ddx(z1, 0.);
		} else {
			z2 = std::ceil(x.a2);
			twosum(z1, z2, z3, z4);
			return ddx(z3, z4);
		}
	}

	friend ddx frexp(const ddx& x, int* m) {
		_Float64x z1, z2;
		z1 = std::frexp(x.a1, m);
		z2 = std::ldexp(x.a2, -(*m));
		if ((z1 == 0.5 && z2 < 0) || (z1 == -0.5 && z2 > 0)) {
			z1 *= 2;
			z2 *= 2;
			(*m)--;
		}
		return ddx(z1, z2);
	}

	friend ddx ldexp(const ddx& x, int m) {
		_Float64x z1, z2;
		z1 = std::ldexp(x.a1, m);
		z2 = std::ldexp(x.a2, m);
		return ddx(z1, z2);
	}

	operator int() const {
		return (int)a1;
	}

	operator double() const {
		return (double)a1;
	}

	operator _Float64x() const {
		return (_Float64x)a1;
	}


	static ddx pi() {
		static const ddx tmp(
			"3.1415926535897932384626433832795028841971693993751"
		);
		return tmp;
	}

	static ddx e() {
		static const ddx tmp(
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static ddx ln2() {
		static const ddx tmp(
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}

	friend ddx pow(const ddx& x, int y) {
		ddx r, xp;
		int a, tmp;

		if (y == 0) return ddx(1.);

		a = (y >= 0) ? y : -y;

		tmp = a;
		r = 1.;
		xp = x;
		while (tmp != 0) {
			if (tmp % 2 != 0) {
				r *= xp;
			}
			tmp /= 2;
			xp = xp * xp;
		}

		if (y < 0) {
			r = 1. / r;
		}

		return r;
	}

	friend ddx pow(const ddx& x, const ddx& y) {
		return exp(y * log(x));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value && ! boost::is_integral<C>::value, ddx >::type pow(const ddx& x, const C& y) {
                return pow(x, ddx(y));
        }

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type pow(const ddx& x, const C& y) {
                return pow(x, ddx(y));
        }

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, ddx>::value, ddx >::type pow(const C& x, const ddx& y) {
                return pow(ddx(x), y);
        }

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, ddx>::value, ddx >::type pow(const C& x, const ddx& y) {
                return pow(ddx(x), y);
        }

	// power by integer (i is passed by double)
	static ddx ipower(const ddx& x, double i) {
		double tmp;
		ddx xp = x;
		ddx r(1.);

		if (i != i) return (ddx)i; // NaN check

		while (i != 0.) {
			i *= 0.5;
			using std::floor;
			tmp = floor(i);
			if (tmp != i) {
				i = tmp;
				r *= xp;
			}
			xp = xp * xp;
		}

		return r;
	}

	friend ddx exp(const ddx& x) {
		ddx x_i, x_f, tmp;
		ddx r, y;
		int i;

		if (x == ddx(std::numeric_limits<_Float64x>::infinity(), 0.)) {
			return ddx(std::numeric_limits<_Float64x>::infinity(), 0.);
		}
		if (x == -ddx(std::numeric_limits<_Float64x>::infinity(), 0.)) {
			return (ddx)0.;
		}

		if (x >= 0.) {
			x_i = floor(x);
			x_f = x - x_i;
			if (x_f >= 0.5) {
				x_f -= 1.;
				x_i += 1.;
			}
		} else {
			x_i = -floor(-x);
			x_f = x - x_i;
			if (x_f <= -0.5) {
				x_f += 1.;
				x_i -= 1.;
			}
		}

		r = 1.;
		y = 1.;
		for (i=1;  i<=29 ; i++) {
			y *= x_f;
			y /= i;
			r += y;
		}

		if (x_i >= 0.) {
			// r *= pow(constants<dd>::e(), (int)x_i);
			r *= ipower(e(), (double)x_i);
		} else {
			// r /= pow(constants<dd>::e(), -(int)x_i);
			r /= ipower(e(), -(double)x_i);
		}

		return r;
	}

	static ddx expm1_origin(const ddx& x) {
		ddx r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=29 ; i++) {
			y *= x;
			y /= i;
			r += y;
		}

		return r;
	}

	friend ddx expm1(const ddx& x) {
		if (x >= -0.5 && x <= 0.5) {
			return expm1_origin(x);
		} else {
			return exp(x) - 1.;
		}
	}

	friend ddx log(const ddx& x) {
		ddx x2, x2m1;
		ddx r;
		ddx xn;
		ddx sqrt2 = sqrt(ddx(2.));
		int p_i;
		ddx p;
		int i;

		if (x < 0.) {
			throw std::domain_error("ddx: log of negative value");
                }

		if (x == ddx(std::numeric_limits<_Float64x>::infinity(), 0.)) {
			return ddx(std::numeric_limits<_Float64x>::infinity(), 0.);
		}
		if (x == 0.) {
			return -ddx(std::numeric_limits<_Float64x>::infinity(), 0.);
		}

		x2 = frexp(x, &p_i);
		p = p_i;

		while (x2 > 4. * std::sqrt(2.) - 4.) {
			x2 *= 0.5;
			p += 1.;
		}
		while (x2 > 4. - 2. * std::sqrt(2.)) {
			x2 /= sqrt2;
			p += 0.5;
		}
		while (x2 < 2. - std::sqrt(2.)) {
			x2 *= 2.;
			p -= 1.;
		}
		while (x2 < 2. * std::sqrt(2.) - 2.) {
			x2 *= sqrt2;
			p -= 0.5;
		}

		x2m1 = x2 - 1.;
		r = 0.;
		xn = -1.;
		for (i=1; i<=54; i++) {
			xn = -xn * x2m1; 
			r += xn / i;
		}

		r += ln2() * p;

		return r;
	}

	static ddx log1p_origin(const ddx& x) {
		ddx r;
		ddx xn;
		int i;

		r = 0.;
		xn = -1.;
		for (i=1; i<=54 ; i++) {
			xn = -xn * x; 
			r += xn / i;
		}

		return r;
	}

	friend ddx log1p(const ddx& x) {

		if (x >= -(3. - 2. * std::sqrt(2.)) && x <= 3. - 2. * std::sqrt(2.)) {
			return log1p_origin(x);
		} else {
			return log(x + 1.);
		}
	}

	static ddx sin_origin(const ddx& I) {
		ddx r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=32 ; i++) {
			y *= I;
			y /= i;
			if (i % 2 != 0) {
				if (i % 4 == 1) {
					r += y;
				} else {
					r -= y;
				}
			}
		}

		return r;
	}

	static ddx cos_origin(const ddx& I) {
		ddx r, y;
		int i;

		r = 1.;
		y = 1.;
		for (i=1; i<=32 ; i++) {
			y *= I;
			y /= i;
			if (i % 2 == 0) {
				if (i % 4 == 0) {
					r += y;
				} else {
					r -= y;
				}
			}
		}

		return r;
	}

	friend ddx sin(ddx I) {
		const ddx p = pi();
		int n;

		if (I <= -p || I >= p) {
			using std::floor;
			n = floor(I / (p * 2.) + 0.5);
			I -= n * p * 2.;
		}

		if (I <= -p * 3. / 4.) {
			return -sin_origin(I + p);
		}
		if (I <= -p * 0.5) {
			return -cos_origin(-p * 0.5 - I);
		}
		if (I <= -p * 0.25) {
			return -cos_origin(I + p * 0.5);
		}
		if (I <= 0.) {
			return -sin_origin(-I);
		}
		if (I <= p * 0.25) {
			return sin_origin(I);
		}
		if (I <= p * 0.5) {
			return cos_origin(p * 0.5 - I);
		}
		if (I <= p * 3. / 4.) {
			return cos_origin(I - p * 0.5);
		}
		return sin_origin(p - I);
	}

	friend ddx cos(ddx I) {
		const ddx p = pi();
		int n;

		if (I <= -p || I >= p) {
			using std::floor;
			n = floor(I / (p * 2.) + 0.5);
			I -= n * p * 2.;
		}

		if (I <= -p * 3. / 4.) {
			return -cos_origin(I + p);
		}
		if (I <= -p * 0.5) {
			return -sin_origin(-p * 0.5 - I);
		}
		if (I <= -p * 0.25) {
			return sin_origin(I + p * 0.5);
		}
		if (I <= 0.) {
			return cos_origin(-I);
		}
		if (I <= p * 0.25) {
			return cos_origin(I);
		}
		if (I <= p * 0.5) {
			return sin_origin(p * 0.5 - I);
		}
		if (I <= p * 3. / 4.) {
			return -sin_origin(I - p * 0.5);
		}
		return -cos_origin(p - I);
	}

	friend ddx tan(const ddx& x) {
		return sin(x) / cos(x);
	}

	static ddx atan_origin(const ddx& I) {
		ddx r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=96 ; i++) {
			y *= I;
			if (i % 2 != 0) {
				if (i % 4 == 1) {
					r += y / i;
				} else {
					r -= y / i;
				}
			}
		}

		return r;
	}

	friend ddx atan(const ddx& x) {
		const ddx p = pi();

		if (x < -(std::sqrt(2.) + 1.)) {
			return -p * 0.5 - atan_origin(1. / x);
		}
		if (x < -(std::sqrt(2.) - 1.)) {
			return -p * 0.25 + atan_origin((1. + x)/(1 - x));
		}
		if (x < (std::sqrt(2.) - 1.)) {
			return atan_origin(x);
		}
		if (x < (std::sqrt(2.) + 1.)) {
			return p * 0.25 + atan_origin((x - 1.)/(x + 1.));
		}
		return p * 0.5 - atan_origin(1. / x);
	}

	friend ddx asin(const ddx& x) {
		const ddx pih = pi() * 0.5;

		if (x == 1) return pih;
		if (x == -1) return -pih;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return atan(x / sqrt(1. - ddx(x) * x));
		} else {
			if (x > 0.) {
				return atan(x / sqrt((1. + ddx(x)) * (1. - x)));
			} else {
				return atan(x / sqrt((1. + x) * (1. - ddx(x))));
			}
		}
	}

	static ddx pih_m_atan(const ddx& I) {
		const ddx p = pi();

		if (I < -(std::sqrt(2.) + 1.)) {
			return p + atan_origin(1. / I);
		}
		if (I < -(std::sqrt(2.) - 1.)) {
			return p * 0.75 - atan_origin((1. + I)/(1 - I));
		}
		if (I < (std::sqrt(2.) - 1.)) {
			return p * 0.5 - atan_origin(I);
		}
		if (I < (std::sqrt(2.) + 1.)) {
			return p * 0.25 - atan_origin((I - 1.)/(I + 1.));
		}
		return atan_origin(1. / I);
	}

	friend ddx acos(const ddx& x) {
		const ddx p = pi();

		if (x == 1.) return ddx(0.);
		if (x == -1.) return p;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return pih_m_atan(x / sqrt(1. - ddx(x) * x));
		} else {
			if (x > 0.) {
				return pih_m_atan(x / sqrt((1. + ddx(x)) * (1. - x)));
			} else {
				return pih_m_atan(x / sqrt((1. + x) * (1. - ddx(x))));
			}
		}
	}

	friend ddx atan2(const ddx& y, const ddx& x) {
		const ddx p = pi();

		if (y <= x && y > -x) {
			return atan(y / x);
		}
		if (y > x && y > -x) {
			return p * 0.5 - atan(x / y);
		}
		if (y > x && y <= -x) {
			if (y >= 0.) {
				return p + atan(y / x);
			} else {
				return -p + atan(y / x);
			}
		}
		return -p * 0.5 - atan(x / y);
	}

	static ddx sinh_origin(const ddx& x) {
		ddx r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=29 ; i++) {
			y *= x;
			y /= i;
			if (i % 2 != 0) {
				r += y;
			}
		}

		return r;
	}

	friend ddx sinh(const ddx& x) {
		ddx tmp;
		if (x >= -0.5 && x <= 0.5) {
			return sinh_origin(x);
		} else if (x > 0.) {
			tmp = exp(x);
			return (tmp - 1. / tmp) * 0.5;
		} else {
			tmp = exp(-x);
			return (1. / tmp - tmp) * 0.5;
		}
	}

	friend ddx cosh(const ddx& x) {
		ddx tmp;
		if (x >= 0.) {
			tmp = exp(x);
			return (tmp + 1. / tmp) * 0.5;
		} else {
			tmp = exp(-x);
			return (1. / tmp + tmp) * 0.5;
		}
	}

	friend ddx tanh(const ddx& x) {
		if (x > 0.5) {
			return 1. - 2. / (1. + exp(2. * x));
		} else if (x < -0.5) {
			return 2. / (1. + exp(-2. * x)) - 1.;
		} else {
			return sinh_origin(x) / cosh(x);
		}
	}

	friend ddx asinh(const ddx& x) {
		if (x < -0.5) {
			return -log(-x + sqrt(1. + ddx(x) * x));
		} else if (x <= 0.5) {
			return log1p((1. + x / (1. + sqrt(1. + ddx(x) * x))) * x);
		} else {
			return log(x + sqrt(1. + ddx(x) * x));
		}
	}

	friend ddx acosh(const ddx& x) {
		if (x == 1.) {
			return ddx(0.);
		} else if (x <= 1.5) {
			ddx y(x - 1.);
			return log1p(y + sqrt(y * (ddx(x) + 1.)));
		} else {
			return log(x + sqrt(ddx(x) * x - 1.));
		}
	}

	friend ddx atanh(const ddx& x) {
		if (x == -1.) return -ddx(std::numeric_limits<_Float64x>::infinity(), 0.);
		if (x == 1.) return ddx(std::numeric_limits<_Float64x>::infinity(), 0.);
		if (x < -0.5) {
			return 0.5 * log((1. + x) / (1. - ddx(x)));
		} else if (x <= 0.5) {
			return 0.5 * log1p(2. * x / (1. - ddx(x)));
		} else {
			return 0.5 * log((1. + ddx(x)) / (1. - x));
		}
	}
};

} // namespace kv

namespace std {
template <> class numeric_limits<kv::ddx> {
	public:

	static kv::ddx epsilon() {
		static const kv::ddx tmp(std::ldexp((_Float64x)1., -127));
		return tmp;
	}
	static kv::ddx infinity() {
		static const _Float64x tmp = numeric_limits<_Float64x>::infinity();
		static const kv::ddx tmp2(tmp, 0.);
		return tmp2;
	}

	static kv::ddx max() {
		static const _Float64x tmp = numeric_limits<_Float64x>::max();
		static const kv::ddx tmp2(tmp, std::ldexp(tmp, -65));
		return tmp2;
	}

	static kv::ddx min() {
		static const _Float64x tmp = numeric_limits<_Float64x>::min();
		static const kv::ddx tmp2(tmp, 0.);
		return tmp2;
	}

	static kv::ddx denorm_min() {
		static const _Float64x tmp = numeric_limits<_Float64x>::denorm_min();
		static const kv::ddx tmp2(tmp, 0.);
		return tmp2;
	}

	static const int digits = 128;
	static const int digits10 = (digits-1)*0.30103; // ln(2)/ln(10)
	static const int radix = 2;
	static const int min_exponent = -16381;
	static const int min_exponent10 = -4931;
	static const int max_exponent = 16384;
	static const int max_exponent10 = 4932;
};
} // namespace std

namespace kv {
template <> struct constants<ddx> {
	static ddx pi() {
		return ddx::pi();
	}

	static ddx e() {
		return ddx::e();
	}

	static ddx ln2() {
		return ddx::ln2();
	}

	static ddx str(const std::string& s) {
		return ddx(s);
	}
};

} // namespace kv

#endif // DDX_HPP
