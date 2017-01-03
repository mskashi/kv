/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDOUBLE_HPP
#define RDOUBLE_HPP

#ifdef KV_NOHWROUND
#include <kv/rdouble-nohwround.hpp>
#else
#include <kv/rdouble-hwround.hpp>
#endif

namespace kv {
template <> struct constants< interval<double> > {
	static interval<double> pi() {
		static const interval<double> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<double> e() {
		static const interval<double> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<double> ln2() {
		static const interval<double> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}
	static interval<double> str(const std::string& s) {
		return interval<double>(s, s);
	}
	static interval<double> str(const std::string& s1, const std::string& s2) {
		return interval<double>(s1, s2);
	}
};
} // namespace kv

#endif // RDOUBLE_HPP
