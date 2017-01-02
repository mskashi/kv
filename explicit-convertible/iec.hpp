template <typename From, typename To> struct is_explicitly_convertible {
	template <int> struct D { };
	typedef char yes[1];
	typedef char no[2];

	template<typename T, typename U> static yes &f(int, D<sizeof T(*(U*)0)>* = 0);
	template<typename T, typename U> static no &f(...);

	static bool const value = (sizeof f<To, From>(0) == sizeof(yes));
};
