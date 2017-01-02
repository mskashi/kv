/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

namespace kv {
	// C is convertible to T or not
	template <class C, class T> struct convertible {
		static const bool value = boost::is_convertible<C, T>::value;
	};
	// C is acceptable as numeric constant of T
	template <class C, class T> struct acceptable_n;
	// C is acceptable as string constant of T
	template <class C, class T> struct acceptable_s;
}

#endif // CONVERT_HPP
