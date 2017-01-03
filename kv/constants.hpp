/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <sstream>
#include <string>

namespace kv {

#if 0
template <class T> struct constants {
	private:

	static T atan_forpi(const T& x) {
		T t, r, y;
		int i;

		r = 0.;
		y = 1.;
		for (i=1;  ; i++) {
			y *= x;
			t = y / (T)i;
			if (i % 2 != 0) {
				if (i % 4 == 1) {
					r += t;
				} else {
					r -= t;
				}
			}
			if (abs(t) < std::numeric_limits<T>::epsilon()) {
				break;
			}
		}

		return r;
	}

	public:

	static T pi() {
		/*
		static const T tmp(
			"3.1415926535897932384626433832795028841971693993751"
		);
		return tmp;
		*/
		static T tmp(0.);
		if (tmp != 0.) return tmp;

		tmp = 16. * atan_forpi(1. / T(5.)) - 4. * atan_forpi(1. / T(239.));

		return tmp;
	}

	static T e() {
		/*
		static const T tmp(
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
		*/
		static T tmp(0.);
		if (tmp != 0.) return tmp;

		T y;
		int i;

		tmp = 1.;
		y = 1.;
		for (i=1;  ; i++) {
			y /= (T)i;
			tmp += y;
			if (y < std::numeric_limits<T>::epsilon()) {
				break;
			}
		}

		return tmp;
	}

	static T ln2() {
		/*
		static const T tmp(
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
		*/

		static T tmp(0.);
		if (tmp != 0.) return tmp;

		T y, x2, x2m1, xn, t;
		int i;

		x2 = sqrt(sqrt(T(2.)));
		x2m1 = x2 - 1.;
		tmp = 0.;
		xn = -1.;
		for (i=1;  ; i++) {
			xn = -xn * x2m1; 
			t = xn / (T)i;
			tmp += t;
			using std::abs;
			if (abs(t) < std::numeric_limits<T>::epsilon()) {
				break;
			}
		}

		tmp = tmp * 4.;

		return tmp;
	}

	static T str(const std::string& s) {
		return T(s);
	}

	static T str(const std::string& s1, const std::string& s2) {
		return (T)constants<typename T::base_type>::str(s1, s2);
	}
};
#endif

template <class T> struct constants {
	static T pi() {
		return (T)constants<typename T::base_type>::pi();
	}
	static T e() {
		return (T)constants<typename T::base_type>::e();
	}

	static T ln2() {
		return (T)constants<typename T::base_type>::ln2();
	}

	static T str(const std::string& s) {
		return T(s);
	}

	static T str(const std::string& s1, const std::string& s2) {
		return (T)constants<typename T::base_type>::str(s1, s2);
	}
};

template <> struct constants<double> {
	static double pi() {
		return 3.1415926535897932384626433832795028841971693993751;
	}
	static double e() {
		return 2.7182818284590452353602874713526624977572470937000;
	}
	static double ln2() {
		return 0.69314718055994530941723212145817656807550013436026;
	}
	static double str(const std::string& s) {
		std::istringstream is(s);
		double r;
		is >> r;
		return r;
	}
};

} // namespace kv

#endif // CONSTANTS_HPP
