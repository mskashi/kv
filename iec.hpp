#ifndef IEC_HPP
#define IEC_HPP

#include <boost/type_traits.hpp>

namespace kv {

template <typename From, typename To> struct is_explicitly_convertible {
	template <int> struct dummy_size {};
	typedef char yes[1];
	typedef char no[2];

	// template<typename T1, typename T2> static yes &selector(dummy_size<sizeof(static_cast<T2>(boost::declval<T1>()))> *);
	template<typename T1, typename T2> static yes &selector(dummy_size<sizeof(T2(*(T1*)0))> *);
	template<typename T1, typename T2> static no &selector(...);

	static bool const value = (sizeof(selector<From, To>(0)) == sizeof(yes)) || boost::is_convertible<From, To>::value;
};

} // namespace kv

#endif // IEC_HPP
