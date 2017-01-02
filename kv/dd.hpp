/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DD_HPP
#define DD_HPP

#include <iostream>
#include <limits>
#include <cmath>
#include <string>

#include <kv/convert.hpp>
#include <kv/constants.hpp>
#include <kv/conv-dd.hpp>
#include <kv/fpu53.hpp>


namespace kv {


class dd;

template <class C> struct convertible<C, dd> {
	static const bool value = boost::is_arithmetic<C>::value || boost::is_same<C, dd>::value || boost::is_convertible<C, std::string>::value;
};
template <class C> struct acceptable_n<C, dd> {
	static const bool value = boost::is_arithmetic<C>::value;
};
template <class C> struct acceptable_s<C, dd> {
	static const bool value = boost::is_convertible<C, std::string>::value;
};


class dd {
	public:

	double a1;
	double a2;


	static void fasttwosum(const double& a, const double& b, double& x, double& y) {
		double tmp;
		x = a + b;
		tmp = x - a;
		y = b - tmp;
	}

	static void twosum(const double& a, const double& b, double& x, double& y) {
		double tmp;

		x = a + b;
		if (std::fabs(a) > std::fabs(b)) {
			tmp = x - a;
			y = b - tmp;
		} else {
			tmp = x - b;
			y = a - tmp;
		}
	}

	static void split(const double& a, double& x, double& y) {
		static const double sigma = (double)((1L << 27) + 1);
		double tmp;

		tmp = a * sigma;
		x = tmp - (tmp - a);
		y = a - x;
	}

	static void twoproduct(const double& a, const double& b, double& x, double& y) {
		static const double th = ldexp(1., 996);
		static const double c1 = ldexp(1., -28);
		static const double c2 = ldexp(1., 28);
		static const double th2 = ldexp(1., 1023);

		double na, nb, a1, a2, b1, b2;

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
			y = a2 * b2 - ((((x * 0.5) - (a1 * 0.5)  * b1) * 2. - a2 * b1) - a1 * b2);
		} else {
			y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
		}
	}

	dd() {
		a1 = 0.;
		a2 = 0.;
	}

	template <class C> explicit dd(const C& x, typename boost::enable_if_c< acceptable_n<C, dd>::value >::type* =0) {
		a1 = x;
		a2 = 0.;
	}

	template <class C1, class C2> dd(const C1& x, const C2& y, typename boost::enable_if_c< acceptable_n<C1, dd>::value && acceptable_n<C2, dd>::value >::type* =0) {
		a1 = x;
		a2 = y;
	}

	template <class C> explicit dd(const C& x, typename boost::enable_if_c< acceptable_s<C, dd>::value >::type* =0) {
		conv_dd::stringtodd(x, a1, a2);
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, dd>::value, dd& >::type operator=(const C& x) {
		a1 = x;
		a2 = 0.;
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_s<C, dd>::value, dd& >::type operator=(const C& x) {
		conv_dd::stringtodd(x, a1, a2);
		return *this;
	}

	friend dd operator+(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		twosum(x.a1, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2 + y.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator+(const dd& x, const C& y) {
		double z1, z2, z3, z4;

		twosum(x.a1, y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator+(const dd& x, const C& y) {
		return x + dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator+(const C& x, const dd& y) {
		double z1, z2, z3, z4;

		twosum(x, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += y.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator+(const C& x, const dd& y) {
		return dd(x) + y;
	}

	friend dd& operator+=(dd& x, const dd& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd& >::type operator+=(dd& x, const C& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd& >::type operator+=(dd& x, const C& y) {
		x = x + dd(y);
		return x;
	}

	friend dd operator-(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		twosum(x.a1, -y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2 - y.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator-(const dd& x, const C& y) {
		double z1, z2, z3, z4;

		twosum(x.a1, -y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator-(const dd& x, const C& y) {
		return x - dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator-(const C& x, const dd& y) {
		double z1, z2, z3, z4;

		twosum(x, -y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 -= y.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator-(const C& x, const dd& y) {
		return dd(x) - y;
	}

	friend dd& operator-=(dd& x, const dd& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd& >::type operator-=(dd& x, const C& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd& >::type operator-=(dd& x, const C& y) {
		x = x - dd(y);
		return x;
	}

	friend dd operator-(const dd& x) {
		dd r;

		r.a1 = - x.a1;
		r.a2 = - x.a2;

		return r;
	}

	friend dd operator*(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		twoproduct(x.a1, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		#if 0
		// boost LU successfully run on Linux -m32
		volatile double v, v2;
		v = z2;
		v2 = x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		v += v2;
		z2 = v; 
		#endif
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator*(const dd& x, const C& y) {
		double z1, z2, z3, z4;
		volatile double v;

		twoproduct(x.a1, y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2 * y;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator*(const dd& x, const C& y) {
		return x * dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator*(const C& x, const dd& y) {
		double z1, z2, z3, z4;

		twoproduct(x, y.a1, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x * y.a2;
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator*(const C& x, const dd& y) {
		return dd(x) * y;
	}

	friend dd& operator*=(dd& x, const dd& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd& >::type operator*=(dd& x, const C& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd& >::type operator*=(dd& x, const C& y) {
		x = x * dd(y);
		return x;
	}

	friend dd operator/(const dd& x, const dd& y) {
		double z1, z2, z3, z4;

		z1 = x.a1 / y.a1;
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		twoproduct(-z1, y.a1, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			twoproduct(-z1, y.a1 * 0.5, z3, z4);
			z2 = ((((z3 + (x.a1 * 0.5)) - z1 * (y.a2 * 0.5)) + (x.a2 * 0.5)) + z4) / (y.a1 * 0.5);
		} else {
			// z2 = ((((z3 + x.a1) - z1 * y.a2) + x.a2) + z4) / (y.a1 + y.a2);
			z2 = ((((z3 + x.a1) - z1 * y.a2) + x.a2) + z4) / y.a1;
		}
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator/(const dd& x, const C& y) {
		double z1, z2, z3, z4;

		z1 = x.a1 / y;
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (std::fabs(y) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		twoproduct(-z1, y, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			twoproduct(-z1, y * 0.5, z3, z4);
			z2 = (((z3 + (x.a1 * 0.5)) + (x.a2 * 0.5)) + z4) / (y * 0.5);
		} else {
			z2 = (((z3 + x.a1) + x.a2) + z4) / y;
		}
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator/(const dd& x, const C& y) {
		return x / dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator/(const C& x, const dd& y) {
		double z1, z2, z3, z4;

		z1 = x / y.a1;
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (std::fabs(y.a1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}

		twoproduct(-z1, y.a1, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			z2 = (((z3 + (x * 0.5)) - z1 * (y.a2 * 0.5)) + z4) / (y.a1 * 0.5);
		} else {
			z2 = (((z3 + x) - z1 * y.a2) + z4) / y.a1;
		}
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type operator/(const C& x, const dd& y) {
		return dd(x) / y;
	}

	friend dd& operator/=(dd& x, const dd& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd& >::type operator/=(dd& x, const C& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd& >::type operator/=(dd& x, const C& y) {
		x = x / dd(y);
		return x;
	}

	friend std::ostream& operator<<(std::ostream& s, const dd& x) {
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
		s << conv_dd::ddtostring(x.a1, x.a2, s.precision(), format, 0);
		return s;
	}

	friend dd sqrt(const dd& x) {
		dd r;

		if (x == 0.) return dd(0.);
		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		return r;
	}

	friend dd abs(const dd& x) {
		if (x.a1 >= 0.) {
			return x;
		} else {
			return dd(-x.a1, -x.a2);
		}
	}

	friend bool operator<(const dd& x, const dd& y) {
		if (x.a1 < y.a1) return true;
		else if (x.a1 > y.a1) return false;
		else return (x.a2 < y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator<(const dd& x, const C& y) {
		if (x.a1 < y) return true;
		else if (x.a1 > y) return false;
		else return (x.a2 < 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator<(const dd& x, const C& y) {
		return x < dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator<(const C& x, const dd& y) {
		if (x < y.a1) return true;
		else if (x > y.a1) return false;
		else return (0. < y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator<(const C& x, const dd& y) {
		return dd(x) < y;
	}

	friend bool operator<=(const dd& x, const dd& y) {
		if (x.a1 < y.a1) return true;
		else if (x.a1 > y.a1) return false;
		else return (x.a2 <= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator<=(const dd& x, const C& y) {
		if (x.a1 < y) return true;
		else if (x.a1 > y) return false;
		else return (x.a2 <= 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator<=(const dd& x, const C& y) {
		return x <= dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator<=(const C& x, const dd& y) {
		if (x < y.a1) return true;
		else if (x > y.a1) return false;
		else return (0. <= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator<=(const C& x, const dd& y) {
		return dd(x) <= y;
	}

	friend bool operator>(const dd& x, const dd& y) {
		if (x.a1 > y.a1) return true;
		else if (x.a1 < y.a1) return false;
		else return (x.a2 > y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator>(const dd& x, const C& y) {
		if (x.a1 > y) return true;
		else if (x.a1 < y) return false;
		else return (x.a2 > 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator>(const dd& x, const C& y) {
		return x > dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator>(const C& x, const dd& y) {
		if (x > y.a1) return true;
		else if (x < y.a1) return false;
		else return (0. > y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator>(const C& x, const dd& y) {
		return dd(x) > y;
	}

	friend bool operator>=(const dd& x, const dd& y) {
		if (x.a1 > y.a1) return true;
		else if (x.a1 < y.a1) return false;
		else return (x.a2 >= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator>=(const dd& x, const C& y) {
		if (x.a1 > y) return true;
		else if (x.a1 < y) return false;
		else return (x.a2 >= 0.);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator>=(const dd& x, const C& y) {
		return x >= dd(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator>=(const C& x, const dd& y) {
		if (x > y.a1) return true;
		else if (x < y.a1) return false;
		else return (0. >= y.a2);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, bool >::type operator>=(const C& x, const dd& y) {
		return dd(x) >= y;
	}

	friend bool operator==(const dd& x, const dd& y) {
		return (x.a1 == y.a1) && (x.a2 == y.a2);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator==(const dd& x, const C& y) {
		return (x.a1 == y) && (x.a2 == 0.);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator==(const C& x, const dd& y) {
		return (x == y.a1) && (0. == y.a2);
	}

	friend bool operator!=(const dd& x, const dd& y) {
		return (x.a1 != y.a1) || (x.a2 != y.a2);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator!=(const dd& x, const C& y) {
		return (x.a1 != y) || (x.a2 != 0.);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, bool >::type operator!=(const C& x, const dd& y) {
		return (x != y.a1) || (0. != y.a2);
	}

	friend dd floor(const dd& x) {
		double z1, z2, z3, z4;
		z1 = std::floor(x.a1);
		if (z1 != x.a1) {
			return dd(z1, 0.);
		} else {
			z2 = std::floor(x.a2);
			twosum(z1, z2, z3, z4);
			return dd(z3, z4);
		}
	}

	friend dd frexp(const dd& x, int* m) {
		double z1, z2;
		z1 = std::frexp(x.a1, m);
		if (*m == 0) {
			return x;
		}
		z2 = ldexp(x.a2, -(*m));
		return dd(z1, z2);
	}

	operator int() const {
		return (int)a1;
	}

	operator double() const {
		return (double)a1;
	}
};

} // namespace kv

namespace std {
template <> class numeric_limits<kv::dd> {
	public:

	static kv::dd epsilon() {
		static const kv::dd tmp(ldexp(1., -105));
		return tmp;
	}
	static kv::dd infinity() {
		static const double tmp = numeric_limits<double>::infinity();
		static const kv::dd tmp2(tmp, 0.);
		return tmp2;
	}

	static kv::dd max() {
		static const double tmp = numeric_limits<double>::max();
		static const kv::dd tmp2(tmp, ldexp(tmp, -54));
		return tmp2;
	}

	static kv::dd min() {
		static const double tmp = numeric_limits<double>::min();
		static const kv::dd tmp2(tmp, 0.);
		return tmp2;
	}

	static const int digits = 106;
	static const int digits10 = 31;
	static const int radix = 2;
	static const int min_exponent = -1021;
	static const int min_exponent10 = -307;
	static const int max_exponent = 1024;
	static const int max_exponent10 = 308;
};
} // namespace std

#endif // DD_HPP
