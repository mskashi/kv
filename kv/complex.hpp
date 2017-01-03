/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <iostream>
#include <cmath>

#include <kv/convert.hpp>


namespace kv {


template <class T> class complex;

template <class C, class T> struct convertible<C, complex<T> > {
	static const bool value = convertible<C, T>::value || boost::is_same<C, complex<T> >::value;
};
template <class C, class T> struct acceptable_n<C, complex<T> > {
	static const bool value = convertible<C, T>::value;
};


template <class T> class complex {
	T re;
	T im;

	public:

	typedef T base_type;

	complex() {
		re = 0.;
		im = 0.;
	}

	template <class C> explicit complex(const C& x, typename boost::enable_if_c< acceptable_n<C, complex>::value >::type* =0) {
		re = x;
		im = 0.;
	}

	template <class C> explicit complex(const complex<C>& x, typename boost::enable_if_c< acceptable_n<C, complex>::value >::type* =0) {
		re = x.real();
		im = x.imag();
	}

	template <class C1, class C2> complex(const C1& x, const C2& y, typename boost::enable_if_c< acceptable_n<C1, complex>::value && acceptable_n<C2, complex>::value >::type* =0) {
		re = x;
		im = y;
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator=(const C& x) {
		re = x;
		im = 0.;
		return *this;
	}

	template <class C> typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator=(const complex<C>& x) {
		re = x.real();
		im = x.imag();
		return *this;
	}

	friend complex operator+(const complex& x, const complex& y) {
		complex r;

		r.re = x.re + y.re;
		r.im = x.im + y.im;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator+(const complex& x, const C& y) {
		complex r;

		r.re = x.re + y;
		r.im = x.im;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator+(const C& x, const complex& y) {
		complex r;

		r.re = x + y.re;
		r.im = y.im;

		return r;
	}

	friend complex& operator+=(complex& x, const complex& y) {
		x = x + y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator+=(complex& x, const C& y) {
		x.re += y;
		return x;
	}

	friend complex operator-(const complex& x, const complex& y) {
		complex r;

		r.re = x.re - y.re;
		r.im = x.im - y.im;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator-(const complex& x, const C& y) {
		complex r;

		r.re = x.re - y;
		r.im = x.im;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator-(const C& x, const complex& y) {
		complex r;

		r.re = x - y.re;
		r.im = - y.im;

		return r;
	}

	friend complex& operator-=(complex& x, const complex& y) {
		x = x - y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator-=(complex& x, const C& y) {
		x.re -= y;
		return x;
	}

	friend complex operator-(const complex& x) {
		complex r;

		r.re = - x.re;
		r.im = - x.im;

		return r;
	}

	friend complex operator*(const complex& x, const complex& y) {
		complex r;

		r.re = x.re * y.re - x.im * y.im;
		r.im = x.re * y.im + x.im * y.re;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator*(const complex& x, const C& y) {
		complex r;

		r.re = x.re * y;
		r.im = x.im * y;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator*(const C& x, const complex& y) {
		complex r;

		r.re = x * y.re;
		r.im = x * y.im;

		return r;
	}

	friend complex& operator*=(complex& x, const complex& y) {
		x = x * y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator*=(complex& x, const C& y) {
		x.re *= y;
		x.im *= y;
		return x;
	}

	friend complex operator/(const complex& x, const complex& y) {
		complex r;
		T tmp, tmp2;

		tmp = y.re * y.re + y.im * y.im;
		r.re = (y.re * x.re + y.im * x.im) / tmp;
		r.im = (y.re * x.im - x.re * y.im) / tmp;

		#if 0
		using std::abs;
		if (abs(y.re) > abs(y.im)) {
			tmp2 = y.im / y.re;
			tmp = y.re + y.im * tmp2;
			r.re = (x.re + tmp2 * x.im) / tmp;
			r.im = (x.im - x.re * tmp2) / tmp;
		} else {
			tmp2 = y.re / y.im;
			tmp = y.re * tmp2 + y.im;
			r.re = (tmp2 * x.re + x.im) / tmp;
			r.im = (tmp2 * x.im - x.re) / tmp;
		}
		#endif

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator/(const complex& x, const C& y) {
		complex r;

		r.re = x.re / y;
		r.im = x.im / y;

		return r;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type operator/(const C& x, const complex& y) {
		complex r;
		T tmp;

		tmp = y.re * y.re + y.im * y.im;
		r.re = (y.re * x) / tmp;
		r.im = (- x * y.im) / tmp;

		#if 0
		T tmp2;
		using std::abs;
		if (abs(y.re) > abs(y.im)) {
			tmp2 = y.im / y.re;
			tmp = y.re + y.im * tmp2;
			r.re = x / tmp;
			r.im = (- x * tmp2) / tmp;
		} else {
			tmp2 = y.re / y.im;
			tmp = y.re * tmp2 + y.im;
			r.re = (tmp2 * x) / tmp;
			r.im = (- x) / tmp;
		}
		#endif

		return r;
	}

	friend complex& operator/=(complex& x, const complex& y) {
		x = x / y;
		return x;
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex& >::type operator/=(complex& x, const C& y) {
		x.re /= y;
		x.im /= y;
		return x;
	}


	friend std::ostream& operator<<(std::ostream& s, const complex& x) {
		s << '(' << x.re << ")+(" << x.im << ")i";
		return s;
	}


	const T& real() const {
		return re;
	}
	const T& imag() const {
		return im;
	}
	T& real() {
		return re;
	}
	T& imag() {
		return im;
	}

	static complex i() {
		return complex(T(0.), T(1.));
	}


	friend T abs(const complex& x) {
		using std::sqrt;
		using std::pow;
		return sqrt(pow(x.re, 2) + pow(x.im, 2));
	}

	friend T arg(const complex& x) {
		using std::atan2;
		return atan2(x.im, x.re);
	}

	friend complex conj(const complex& x) {
		return complex(x.re, -x.im);
	}

	friend complex sqrt(const complex& x) {
		T a, r;
		r = sqrt(abs(x));
		a = arg(x) * 0.5;
		using std::sqrt;
		using std::cos;
		using std::sin;
		return complex(r * cos(a), r * sin(a));
	}

	friend complex pow(const complex& x, int y) {
		complex r, xp;
		int tmp;

		if (y == 0) return complex(1.);

		tmp = (y >= 0) ? y : -y;

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

	friend complex pow(const complex& x, const complex& y) {
		return exp(y * log(x));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value && ! boost::is_integral<C>::value, complex >::type pow(const complex& x, const C& y) {
		return pow(x, complex(y));
	}

	template <class C> friend typename boost::enable_if_c< acceptable_n<C, complex>::value, complex >::type pow(const C& x, const complex& y) {
		return pow(complex(x), y);
	}

	friend complex exp(const complex& x) {
		T tmp;
		using std::exp;
		using std::cos;
		using std::sin;
		tmp = exp(x.re);
		return complex(tmp * cos(x.im), tmp * sin(x.im));
	}

	friend complex log(const complex& x) {
		using std::log;
		return complex(log(abs(x)), arg(x));
	}

	friend complex sin(const complex& x) {
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		return complex(sin(x.re) * cosh(x.im), cos(x.re) * sinh(x.im));
	}

	friend complex cos(const complex& x) {
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		return complex(cos(x.re) * cosh(x.im), - sin(x.re) * sinh(x.im));
	}

	friend complex tan(const complex& x) {
		T tx, ty, tmp;
		tx = 2. * x.re;
		ty = 2. * x.im;
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		tmp = cos(tx) + cosh(ty);
		return complex(sin(tx) / tmp, sinh(ty) / tmp);
	}

	friend complex asin(const complex& x) {
		return -i() * log(i() * x + sqrt(1. - x * x)); 
	}

	friend complex acos(const complex& x) {
		return -i() * log(x + i() * sqrt(1. - x * x)); 
	}

	friend complex atan(const complex& x) {
		return i() * 0.5 * log((i() + x) / (i() - x)); 
	}

	friend complex sinh(const complex& x) {
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		return complex(sinh(x.re) * cos(x.im), cosh(x.re) * sin(x.im));
	}

	friend complex cosh(const complex& x) {
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		return complex(cosh(x.re) * cos(x.im), sinh(x.re) * sin(x.im));
	}

	friend complex tanh(const complex& x) {
		T tx, ty, tmp;
		tx = 2. * x.re;
		ty = 2. * x.im;
		using std::sin;
		using std::cos;
		using std::cosh;
		using std::sinh;
		tmp = cosh(tx) + cos(ty);
		return complex(sinh(tx) / tmp, sin(ty) / tmp);
	}

	friend complex asinh(const complex& x) {
		return log(x + sqrt(x * x + 1.)); 
	}

	friend complex acosh(const complex& x) {
		return log(x + sqrt(x * x - 1.)); 
	}

	friend complex atanh(const complex& x) {
		return 0.5 * log((1. + x) / (1. - x)); 
	}
};

} // namespace kv

#endif // COMPLEX_HPP
