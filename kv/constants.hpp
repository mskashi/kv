/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <sstream>
#include <string>

namespace kv {

template <class T> struct constants {
	static T pi() {
		static const T tmp(
			"3.1415926535897932384626433832795028841971693993751"
		);
		return tmp;
	}
	static T e() {
		static const T tmp(
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}
	static T ln2() {
		static const T tmp(
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
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
