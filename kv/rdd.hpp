/*
 * Copyright (c) 2013-2015 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef RDD_HPP
#define RDD_HPP

#ifdef KV_NOHWROUND
#include <kv/rdd-nohwround.hpp>
#else
#include <kv/rdd-hwround.hpp>
#endif

namespace kv {
template <> struct constants< interval<dd> > {
	static interval<dd> pi() {
		static const interval<dd> tmp(
			"3.1415926535897932384626433832795028841971693993751",
			"3.1415926535897932384626433832795028841971693993752"
		);
		return tmp;
	}

	static interval<dd> e() {
		static const interval<dd> tmp(
			"2.7182818284590452353602874713526624977572470936999",
			"2.7182818284590452353602874713526624977572470937000"
		);
		return tmp;
	}

	static interval<dd> ln2() {
		static const interval<dd> tmp(
			"0.69314718055994530941723212145817656807550013436025",
			"0.69314718055994530941723212145817656807550013436026"
		);
		return tmp;
	}

	static interval<dd> str(const std::string& s) {
		return interval<dd>(s, s);
	}

	static interval<dd> str(const std::string& s1, const std::string& s2) {
		return interval<dd>(s1, s2);
	}
};
} // namespace kv

#endif // RDD_HPP
