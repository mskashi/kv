/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef MPFR_HPP
#define MPFR_HPP

#include <iostream>
#include <limits>
#include <string>
#include <kv/convert.hpp>
#include <kv/constants.hpp>
#include <mpfr.h>

namespace kv {

template <int N> class mpfr;

template <class C, int N> struct convertible<C, mpfr<N> > {
	static const bool value = boost::is_arithmetic<C>::value || boost::is_same<C, mpfr<N> >::value || boost::is_convertible<C, std::string>::value;
};
template <class C, int N> struct acceptable_n<C, mpfr<N> > {
	static const bool value = boost::is_arithmetic<C>::value;
};
template <class C, int N> struct acceptable_s<C, mpfr<N> > {
	static const bool value = boost::is_convertible<C, std::string>::value;
};

template <int N = 53> class mpfr {
	public:
	mpfr_t a;

	mpfr () {
		mpfr_init2(a, N);
		mpfr_set_si(a, 0, MPFR_RNDN);
	}

	// copy constructor is needed
	mpfr (const mpfr& x) {
		mpfr_init2(a, N);
		mpfr_set(a, x.a, MPFR_RNDN);
	}

	template <class C> explicit mpfr(const C& x, typename boost::enable_if_c< acceptable_n<C, mpfr>::value >::type* =0) {
		mpfr_init2(a, N);
		mpfr_set_d(a, (double)x, MPFR_RNDN);
	}

	template <class C> explicit mpfr(const C& x, typename boost::enable_if_c< acceptable_s<C, mpfr>::value >::type* =0) {
		mpfr_init2(a, N);
		mpfr_set_str(a, std::string(x).c_str(), 10, MPFR_RNDN);
	}

	~mpfr () {
		mpfr_clear(a);
	}

	// assignment operator must be overloaded
	mpfr& operator=(const mpfr& x) {
		mpfr_set(a, x.a, MPFR_RNDN);
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr& >::type operator=(const C& x) {
		mpfr_set_d(a, (double)x, MPFR_RNDN);
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr& >::type operator=(const C& x) {
		mpfr_set_str(a, x, 10, MPFR_RNDN);
		return *this;
	}

	friend mpfr operator+(const mpfr& x, const mpfr& y) {
		mpfr z;
		mpfr_add(z.a, x.a, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator+(const mpfr& x, const C& y) {
		mpfr z;
		mpfr_add_d(z.a, x.a, (double)y, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator+(const mpfr& x, const C& y) {
		return x + mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator+(const C& x, const mpfr& y) {
		mpfr z;
		mpfr_add_d(z.a, y.a, (double)x, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator+(const C& x, const mpfr& y) {
		return mpfr(x) + y;
	}

	friend mpfr& operator+=(mpfr& x, const mpfr& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr& >::type operator+=(mpfr& x, const C& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr& >::type operator+=(mpfr& x, const C& y) {
		x = x + mpfr(y);
		return x;
	}

	friend mpfr operator-(const mpfr& x, const mpfr& y) {
		mpfr z;
		mpfr_sub(z.a, x.a, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator-(const mpfr& x, const C& y) {
		mpfr z;
		mpfr_sub_d(z.a, x.a, (double)y, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator-(const mpfr& x, const C& y) {
		return x - mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator-(const C& x, const mpfr& y) {
		mpfr z;
		mpfr_d_sub(z.a, (double)x, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator-(const C& x, const mpfr& y) {
		return mpfr(x) - y;
	}

	friend mpfr& operator-=(mpfr& x, const mpfr& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr& >::type operator-=(mpfr& x, const C& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr& >::type operator-=(mpfr& x, const C& y) {
		x = x - mpfr(y);
		return x;
	}

	friend mpfr operator-(const mpfr& x) {
		mpfr r;

		mpfr_neg(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr operator*(const mpfr& x, const mpfr& y) {
		mpfr z;
		mpfr_mul(z.a, x.a, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator*(const mpfr& x, const C& y) {
		mpfr z;
		mpfr_mul_d(z.a, x.a, (double)y, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator*(const mpfr& x, const C& y) {
		return x * mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator*(const C& x, const mpfr& y) {
		mpfr z;
		mpfr_mul_d(z.a, y.a, (double)x, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator*(const C& x, const mpfr& y) {
		return mpfr(x) * y;
	}

	friend mpfr& operator*=(mpfr& x, const mpfr& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr& >::type operator*=(mpfr& x, const C& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr& >::type operator*=(mpfr& x, const C& y) {
		x = x * mpfr(y);
		return x;
	}

	friend mpfr operator/(const mpfr& x, const mpfr& y) {
		mpfr z;
		mpfr_div(z.a, x.a, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator/(const mpfr& x, const C& y) {
		mpfr z;
		mpfr_div_d(z.a, x.a, (double)y, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator/(const mpfr& x, const C& y) {
		return x / mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type operator/(const C& x, const mpfr& y) {
		mpfr z;
		mpfr_d_div(z.a, (double)x, y.a, MPFR_RNDN);
		return z;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr >::type operator/(const C& x, const mpfr& y) {
		return mpfr(x) / y;
	}

	friend mpfr& operator/=(mpfr& x, const mpfr& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr& >::type operator/=(mpfr& x, const C& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, mpfr& >::type operator/=(mpfr& x, const C& y) {
		x = x / mpfr(y);
		return x;
	}

	friend mpfr sqrt(const mpfr& x) {
		mpfr r;

		mpfr_sqrt(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr abs(const mpfr& x) {
		mpfr r;

		mpfr_abs(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend bool operator<(const mpfr& x, const mpfr& y) {
		if (mpfr_cmp(x.a, y.a) < 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator<(const mpfr& x, const C& y) {
		if (mpfr_cmp_d(x.a, (double)y) < 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator<(const mpfr& x, const C& y) {
		return x < mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator<(const C& x, const mpfr& y) {
		if (mpfr_cmp_d(y.a, (double)x) > 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator<(const C& x, const mpfr& y) {
		return mpfr(x) < y;
	}

	friend bool operator<=(const mpfr& x, const mpfr& y) {
		if (mpfr_cmp(x.a, y.a) <= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator<=(const mpfr& x, const C& y) {
		if (mpfr_cmp_d(x.a, (double)y) <= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator<=(const mpfr& x, const C& y) {
		return x <= mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator<=(const C& x, const mpfr& y) {
		if (mpfr_cmp_d(y.a, (double)x) >= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator<=(const C& x, const mpfr& y) {
		return mpfr(x) <= y;
	}

	friend bool operator>(const mpfr& x, const mpfr& y) {
		if (mpfr_cmp(x.a, y.a) > 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator>(const mpfr& x, const C& y) {
		if (mpfr_cmp_d(x.a, (double)y) > 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator>(const mpfr& x, const C& y) {
		return x > mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator>(const C& x, const mpfr& y) {
		if (mpfr_cmp_d(y.a, (double)x) < 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator>(const C& x, const mpfr& y) {
		return mpfr(x) > y;
	}

	friend bool operator>=(const mpfr& x, const mpfr& y) {
		if (mpfr_cmp(x.a, y.a) >= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator>=(const mpfr& x, const C& y) {
		if (mpfr_cmp_d(x.a, (double)y) >= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator>=(const mpfr& x, const C& y) {
		return x >= mpfr(y);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator>=(const C& x, const mpfr& y) {
		if (mpfr_cmp_d(y.a, (double)x) <= 0) return true;
		else return false;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_s<C, mpfr>::value, bool >::type operator>=(const C& x, const mpfr& y) {
		return mpfr(x) >= y;
	}

	friend bool operator==(const mpfr& x, const mpfr& y) {
		return (mpfr_cmp(x.a, y.a) == 0);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator==(const mpfr& x, const C& y) {
		return (mpfr_cmp_d(x.a, (double)y) == 0);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator==(const C& x, const mpfr& y) {
		return (mpfr_cmp_d(y.a, (double)x) == 0);
	}

	friend bool operator!=(const mpfr& x, const mpfr& y) {
		return (mpfr_cmp(x.a, y.a) != 0);
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator!=(const mpfr& x, const C& y) {
		return (mpfr_cmp_d(x.a, (double)y) != 0);
	}
	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, bool >::type operator!=(const C& x, const mpfr& y) {
		return (mpfr_cmp_d(y.a, (double)x) != 0);
	}

	friend mpfr floor(const mpfr& x) {
		mpfr r;

		mpfr_floor(r.a, x.a);

		return r;
	}

	friend mpfr ceil(const mpfr& x) {
		mpfr r;

		mpfr_ceil(r.a, x.a);

		return r;
	}

	friend mpfr frexp(const mpfr& x, int* m) {
		mpfr_exp_t e;
		mpfr r;

		mpfr_frexp(&e, r.a, x.a, MPFR_RNDN);

		*m = e;

		return r;
	}

	friend mpfr ldexp(const mpfr& x, int m) {
		mpfr r;

		mpfr_mul_2si(r.a, x.a, m, MPFR_RNDN);

		return r;
	}

	operator int() const {
		return mpfr_get_si(a, MPFR_RNDZ);
	}

	operator double() const {
		return mpfr_get_d(a, MPFR_RNDN);
	}

	friend std::ostream& operator<<(std::ostream& s, const mpfr& x) {
		char format;
		char fmt[100];
		char *str;

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

		sprintf(fmt, "%%.%dRN%c", (int)s.precision(), format);
		mpfr_asprintf(&str, fmt, x.a);
		s << str;
		mpfr_free_str(str);

		return s;
	}

	friend mpfr pow(const mpfr& x, int y) {
		mpfr r;

		mpfr_pow_si(r.a, x.a, y, MPFR_RNDN);

		return r;
	}

	friend mpfr pow(const mpfr& x, const mpfr& y) {
		mpfr r;

		mpfr_pow(r.a, x.a, y.a, MPFR_RNDN);

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value && ! boost::is_integral<C>::value, mpfr >::type pow(const mpfr& x, const C& y) {
		return pow(x, mpfr(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, mpfr>::value, mpfr >::type pow(const C& x, const mpfr& y) {
		return pow(mpfr(x), y);
	}

	friend mpfr exp(const mpfr& x) {
		mpfr r;

		mpfr_exp(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr expm1(const mpfr& x) {
		mpfr r;

		mpfr_expm1(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr log(const mpfr& x) {
		mpfr r;

		mpfr_log(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr log1p(const mpfr& x) {
		mpfr r;

		mpfr_log1p(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr sin(const mpfr& x) {
		mpfr r;

		mpfr_sin(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr cos(const mpfr& x) {
		mpfr r;

		mpfr_cos(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr tan(const mpfr& x) {
		mpfr r;

		mpfr_tan(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr asin(const mpfr& x) {
		mpfr r;

		mpfr_asin(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr acos(const mpfr& x) {
		mpfr r;

		mpfr_acos(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr atan(const mpfr& x) {
		mpfr r;

		mpfr_atan(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr atan2(const mpfr& y, const mpfr& x) {
		mpfr r;

		mpfr_atan2(r.a, y.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr sinh(const mpfr& x) {
		mpfr r;

		mpfr_sinh(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr cosh(const mpfr& x) {
		mpfr r;

		mpfr_cosh(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr tanh(const mpfr& x) {
		mpfr r;

		mpfr_tanh(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr asinh(const mpfr& x) {
		mpfr r;

		mpfr_asinh(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr acosh(const mpfr& x) {
		mpfr r;

		mpfr_acosh(r.a, x.a, MPFR_RNDN);

		return r;
	}

	friend mpfr atanh(const mpfr& x) {
		mpfr r;

		mpfr_atanh(r.a, x.a, MPFR_RNDN);

		return r;
	}
};

} // namespace kv

namespace std {
template <int N> class numeric_limits< kv::mpfr<N> > {
	public: static kv::mpfr<N> epsilon() {
		kv::mpfr<N> tmp;
		mpfr_set_ui_2exp(tmp.a, 1, 1-N, MPFR_RNDN);
		return tmp;
	}
	static kv::mpfr<N> infinity() {
		kv::mpfr<N> tmp;
		mpfr_set_inf(tmp.a, 1);
		return tmp;
	}

	static kv::mpfr<N> max() {
		kv::mpfr<N> tmp, tmp2;
		mpfr_set_ui_2exp(tmp.a, 1, mpfr_get_emax() - 1, MPFR_RNDN);
		tmp2 = 2 - epsilon();
		return tmp2 * tmp;
	}

	static kv::mpfr<N> min() {
		kv::mpfr<N> tmp;
		mpfr_set_ui_2exp(tmp.a, 1, mpfr_get_emin() - 1, MPFR_RNDN);
		return tmp;
	}
	static kv::mpfr<N> denorm_min() {
		return min();
	}

	static const int digits = N;
	static const int digits10 = (digits-1)*0.30103; // ln(2)/ln(10)
	static const int radix = 2;
	// static const int min_exponent = mpfr_get_emin();
	static const int min_exponent = MPFR_EMIN_DEFAULT;
	static const int min_exponent10 = min_exponent*0.30103;
	// static const int max_exponent = mpfr_get_emax();
	static const int max_exponent = MPFR_EMAX_DEFAULT;
	static const int max_exponent10 = max_exponent*0.30103;
	static const float_denorm_style has_denorm = denorm_absent;
};
} // namespace std

namespace kv {
template <int N> struct constants< mpfr<N> > {
	static mpfr<N> pi() {
		static mpfr<N> tmp(0);
		if (tmp != 0) return tmp;
		mpfr_const_pi(tmp.a, MPFR_RNDN);
		return tmp;
	}

	static mpfr<N> e() {
		static mpfr<N> tmp(0);
		if (tmp != 0) return tmp;
		mpfr<N> one(1);
		mpfr_exp(tmp.a, one.a, MPFR_RNDN);
		return tmp;
	}

	static mpfr<N> ln2() {
		static mpfr<N> tmp(0);
		if (tmp != 0) return tmp;
		mpfr_const_log2(tmp.a, MPFR_RNDN);
		return tmp;
	}

	static mpfr<N> str(const std::string& s) {
		return mpfr<N>(s);
	}
};
} // namespace kv

#endif // MPFR_HPP
