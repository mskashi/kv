/*
 * Copyright (c) 2013-2018 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef DD_HPP
#define DD_HPP

#include <iostream>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <string>

#include <kv/convert.hpp>
#include <kv/constants.hpp>
#include <kv/conv-dd.hpp>
#include <kv/fpu53.hpp>


#ifndef DD_FASTMULT
#define DD_FASTMULT 0
#endif

#ifndef DD_NEW_SQRT
#define DD_NEW_SQRT 1
#endif

#ifndef KV_USE_TPFMA
#define KV_USE_TPFMA 0
#endif


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

#if KV_USE_TPFMA == 1
	static void twoproduct(const double& a, const double& b, double& x, double& y) {
		x = a * b;
		y = fma(a, b, -x);
	}
#else // KV_USE_TPFMA
	static void twoproduct(const double& a, const double& b, double& x, double& y) {
		static const double th = std::ldexp(1., 996);
		static const double c1 = std::ldexp(1., -28);
		static const double c2 = std::ldexp(1., 28);
		static const double th2 = std::ldexp(1., 1023);

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
			y = ((((a1 * 0.5) * b1 - (x * 0.5)) * 2 + a2 * b1) + a1 * b2) + a2 * b2;
		} else {
			y = (((a1 * b1 - x) + a2 * b1) + a1 * b2) + a2 * b2;
		}
	}
#endif // KV_USE_TPFMA

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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

		// x.a2 * y.a2 is very small but sometimes important
		#if DD_FASTMULT == 1
		z2 += x.a1 * y.a2 + x.a2 * y.a1;
		#else
		z2 += x.a1 * y.a2 + x.a2 * y.a1 + x.a2 * y.a2;
		#endif

		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator*(const dd& x, const C& y) {
		double z1, z2, z3, z4;

		twoproduct(x.a1, y, z1, z2);
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		z2 += x.a2 * y;
		twosum(z1, z2, z3, z4);
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

		return dd(z3, z4);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type operator/(const dd& x, const C& y) {
		double z1, z2, z3, z4;

		z1 = x.a1 / y;
		if (std::fabs(z1) == std::numeric_limits<double>::infinity()) {
			return dd(z1, 0.);
		}
		if (std::fabs((double)y) == std::numeric_limits<double>::infinity()) {
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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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
		if (std::fabs(z3) == std::numeric_limits<double>::infinity()) {
			return dd(z3, 0.);
		}

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

#if DD_NEW_SQRT == 1
	friend dd sqrt(const dd& x) {
		double z1, z2, z3, z4;

		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
                }

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		z1 = std::sqrt(x.a1);
		twoproduct(-z1, z1, z3, z4);
		z2 = ((z3 + x.a1) + x.a2 + z4) / (2 * z1);
		twosum(z1, z2, z3, z4);

		return dd(z3, z4);
	}
#else
	friend dd sqrt(const dd& x) {
		dd r;

		if (x < 0.) {
			throw std::domain_error("dd: sqrt of negative value");
                }

		if (x == 0.) return dd(0.);
		if (x.a1 == std::numeric_limits<double>::infinity()) {
			return dd(x.a1, 0.);
		}

		r = std::sqrt(x.a1);
		r = (r + x / r) * 0.5;
		return r;
	}
#endif

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

	friend dd ceil(const dd& x) {
		double z1, z2, z3, z4;
		z1 = std::ceil(x.a1);
		if (z1 != x.a1) {
			return dd(z1, 0.);
		} else {
			z2 = std::ceil(x.a2);
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
		z2 = std::ldexp(x.a2, -(*m));
		return dd(z1, z2);
	}

	friend dd ldexp(const dd& x, int m) {
		double z1, z2;
		z1 = std::ldexp(x.a1, m);
		z2 = std::ldexp(x.a2, m);
		return dd(z1, z2);
	}

	operator int() const {
		return (int)a1;
	}

	operator double() const {
		return (double)a1;
	}

	static dd pi() {
		static const dd tmp(
			"3.1415926535897932384626433832795028841971693993751"
		);
		return tmp;
	}

	static dd e() {
		static const dd tmp(
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static dd ln2() {
		static const dd tmp(
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}

	friend dd pow(const dd& x, int y) {
		dd r, xp;
		int a, tmp;

		if (y == 0) return dd(1.);

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

	friend dd pow(const dd& x, const dd& y) {
		return exp(y * log(x));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value && ! boost::is_integral<C>::value, dd >::type pow(const dd& x, const C& y) {
                return pow(x, dd(y));
        }

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type pow(const dd& x, const C& y) {
                return pow(x, dd(y));
        }

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, dd>::value, dd >::type pow(const C& x, const dd& y) {
                return pow(dd(x), y);
        }

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, dd>::value, dd >::type pow(const C& x, const dd& y) {
                return pow(dd(x), y);
        }

	// power by integer (i is passed by double)
	static dd ipower(const dd& x, double i) {
		double tmp;
		dd xp = x;
		dd r(1.);

		if (i != i) return (dd)i; // NaN check

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

	friend dd exp(const dd& x) {
		dd x_i, x_f, tmp;
		dd r, y;
		int i;

		if (x == dd(std::numeric_limits<double>::infinity(), 0.)) {
			return dd(std::numeric_limits<double>::infinity(), 0.);
		}
		if (x == -dd(std::numeric_limits<double>::infinity(), 0.)) {
			return (dd)0.;
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
		for (i=1;  i<=25 ; i++) {
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

	static dd expm1_origin(const dd& x) {
		dd r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=25 ; i++) {
			y *= x;
			y /= i;
			r += y;
		}

		return r;
	}

	friend dd expm1(const dd& x) {
		if (x >= -0.5 && x <= 0.5) {
			return expm1_origin(x);
		} else {
			return exp(x) - 1.;
		}
	}

	friend dd log(const dd& x) {
		dd x2, x2m1;
		dd r;
		dd xn;
		dd sqrt2 = sqrt(dd(2.));
		int p_i;
		dd p;
		int i;

		if (x < 0.) {
			throw std::domain_error("dd: log of negative value");
                }

		if (x == dd(std::numeric_limits<double>::infinity(), 0.)) {
			return dd(std::numeric_limits<double>::infinity(), 0.);
		}
		if (x == 0.) {
			return -dd(std::numeric_limits<double>::infinity(), 0.);
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
		for (i=1; i<=45; i++) {
			xn = -xn * x2m1; 
			r += xn / i;
		}

		r += ln2() * p;

		return r;
	}

	static dd log1p_origin(const dd& x) {
		dd r;
		dd xn;
		int i;

		r = 0.;
		xn = -1.;
		for (i=1; i<=45 ; i++) {
			xn = -xn * x; 
			r += xn / i;
		}

		return r;
	}

	friend dd log1p(const dd& x) {

		if (x >= -(3. - 2. * std::sqrt(2.)) && x <= 3. - 2. * std::sqrt(2.)) {
			return log1p_origin(x);
		} else {
			return log(x + 1.);
		}
	}

	static dd sin_origin(const dd& I) {
		dd r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=28 ; i++) {
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

	static dd cos_origin(const dd& I) {
		dd r, y;
		int i;

		r = 1.;
		y = 1.;
		for (i=1; i<=28 ; i++) {
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

	friend dd sin(const dd& I) {
		const dd p = pi();

		if (I >= p) {
			return sin(I - p * 2.);
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

	friend dd cos(const dd& I) {
		const dd p = pi();

		if (I >= p) {
			return cos(I - p * 2.);
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

	friend dd tan(const dd& x) {
		return sin(x) / cos(x);
	}

	static dd atan_origin(const dd& I) {
		dd r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=79 ; i++) {
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

	friend dd atan(const dd& x) {
		const dd p = pi();

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

	friend dd asin(const dd& x) {
		const dd pih = pi() * 0.5;

		if (x == 1) return pih;
		if (x == -1) return -pih;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return atan(x / sqrt(1. - dd(x) * x));
		} else {
			if (x > 0.) {
				return atan(x / sqrt((1. + dd(x)) * (1. - x)));
			} else {
				return atan(x / sqrt((1. + x) * (1. - dd(x))));
			}
		}
	}

	static dd pih_m_atan(const dd& I) {
		const dd p = pi();

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

	friend dd acos(const dd& x) {
		const dd p = pi();

		if (x == 1.) return dd(0.);
		if (x == -1.) return p;
		using std::abs;
		if (abs(x) < std::sqrt(6.) / 3.) {
			return pih_m_atan(x / sqrt(1. - dd(x) * x));
		} else {
			if (x > 0.) {
				return pih_m_atan(x / sqrt((1. + dd(x)) * (1. - x)));
			} else {
				return pih_m_atan(x / sqrt((1. + x) * (1. - dd(x))));
			}
		}
	}

	friend dd atan2(const dd& y, const dd& x) {
		const dd p = pi();

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

	static dd sinh_origin(const dd& x) {
		dd r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1; i<=25 ; i++) {
			y *= x;
			y /= i;
			if (i % 2 != 0) {
				r += y;
			}
		}

		return r;
	}

	friend dd sinh(const dd& x) {
		dd tmp;
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

	friend dd cosh(const dd& x) {
		dd tmp;
		if (x >= 0.) {
			tmp = exp(x);
			return (tmp + 1. / tmp) * 0.5;
		} else {
			tmp = exp(-x);
			return (1. / tmp + tmp) * 0.5;
		}
	}

	friend dd tanh(const dd& x) {
		if (x > 0.5) {
			return 1. - 2. / (1. + exp(2. * x));
		} else if (x < -0.5) {
			return 2. / (1. + exp(-2. * x)) - 1.;
		} else {
			return sinh_origin(x) / cosh(x);
		}
	}

	friend dd asinh(const dd& x) {
		if (x < -0.5) {
			return -log(-x + sqrt(1. + dd(x) * x));
		} else if (x <= 0.5) {
			return log1p((1. + x / (1. + sqrt(1. + dd(x) * x))) * x);
		} else {
			return log(x + sqrt(1. + dd(x) * x));
		}
	}

	friend dd acosh(const dd& x) {
		if (x == 1.) {
			return dd(0.);
		} else if (x <= 1.5) {
			dd y(x - 1.);
			return log1p(y + sqrt(y * (dd(x) + 1.)));
		} else {
			return log(x + sqrt(dd(x) * x - 1.));
		}
	}

	friend dd atanh(const dd& x) {
		if (x == -1.) return -dd(std::numeric_limits<double>::infinity(), 0.);
		if (x == 1.) return dd(std::numeric_limits<double>::infinity(), 0.);
		if (x < -0.5) {
			return 0.5 * log((1. + x) / (1. - dd(x)));
		} else if (x <= 0.5) {
			return 0.5 * log1p(2. * x / (1. - dd(x)));
		} else {
			return 0.5 * log((1. + dd(x)) / (1. - x));
		}
	}
};

} // namespace kv

namespace std {
template <> class numeric_limits<kv::dd> {
	public:

	static kv::dd epsilon() {
		static const kv::dd tmp(std::ldexp(1., -105));
		return tmp;
	}
	static kv::dd infinity() {
		static const double tmp = numeric_limits<double>::infinity();
		static const kv::dd tmp2(tmp, 0.);
		return tmp2;
	}

	static kv::dd max() {
		static const double tmp = numeric_limits<double>::max();
		static const kv::dd tmp2(tmp, std::ldexp(tmp, -54));
		return tmp2;
	}

	static kv::dd min() {
		static const double tmp = numeric_limits<double>::min();
		static const kv::dd tmp2(tmp, 0.);
		return tmp2;
	}

	static kv::dd denorm_min() {
		static const double tmp = numeric_limits<double>::denorm_min();
		static const kv::dd tmp2(tmp, 0.);
		return tmp2;
	}

	static const int digits = 106;
	static const int digits10 = (digits-1)*0.30103; // ln(2)/ln(10)
	static const int radix = 2;
	static const int min_exponent = -1021;
	static const int min_exponent10 = -307;
	static const int max_exponent = 1024;
	static const int max_exponent10 = 308;
};
} // namespace std

namespace kv {
template <> struct constants<dd> {
	static dd pi() {
		return dd::pi();
	}

	static dd e() {
		return dd::e();
	}

	static dd ln2() {
		return dd::ln2();
	}

	static dd str(const std::string& s) {
		return dd(s);
	}
};
} // namespace kv

#endif // DD_HPP
