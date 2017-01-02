#include <iostream>
#include <boost/type_traits.hpp>
#include "iec.hpp"

class hoge {
	public:
	hoge(const int& x) {
	}
};

class hoge2 {
	public:
	explicit hoge2(const int& x) {
	}
};

int main()
{
	std::cout << is_explicitly_convertible<int, hoge>::value << "\n";
	std::cout << is_explicitly_convertible<int, hoge2>::value << "\n";
	std::cout << is_explicitly_convertible<void *, hoge2>::value << "\n";
	std::cout << boost::is_convertible<int, hoge>::value << "\n";
	std::cout << boost::is_convertible<int, hoge2>::value << "\n";
	std::cout << boost::is_convertible<void*, hoge2>::value << "\n";
}
