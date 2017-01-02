#include <iostream>
#include <boost/type_traits.hpp>

namespace kv {
	template <class C, class T> struct convertible {
		static const bool value = boost::is_convertible<C, T>::value;
	};
	template <class C, class T> struct acceptable_n;
	template <class C, class T> struct acceptable_s;


	template <class T> class interval;
	template <class C, class T> struct convertible<C, interval<T> > {
		static const bool value = convertible<C, T>::value || boost::is_same<C, interval<T> >::value || boost::is_convertible<C, std::string>::value;
	};

	template <class C, class T> struct acceptable_n<C, interval<T> > {
		static const bool value = convertible<C, T>::value && (! boost::is_convertible<C, std::string>::value);
	};

	template <class C, class T> struct acceptable_s<C, interval<T> > {
		static const bool value = boost::is_convertible<C, std::string>::value;
	};

	template <class T> class interval {
		public:
		T inf, sup;
	};


	class dd;
	template <class C> struct convertible<C, dd> {
		static const bool value = boost::is_arithmetic<C>::value || boost::is_same<C, dd>::value || boost::is_convertible<C, std::string>::value;
	};

	template <class C> struct acceptable_n<C, dd> {
		static const bool value = boost::is_arithmetic<C>::value && (! boost::is_convertible<C, std::string>::value);
	};

	template <class C> struct acceptable_s<C, dd> {
		static const bool value = boost::is_convertible<C, std::string>::value;
	};

	class dd {
		public:
		double a1, a2;
	};

}

int main()
{
	std::cout << kv::convertible<int, kv::interval< kv::dd > >::value << "\n";
}
