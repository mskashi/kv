#include <iostream>
#include <kv/defint-singular.hpp>

typedef kv::interval<double> itv;

/*
 * \int_{-1}^1 cos(x) (x-(-1))^{-0.99} (1-x)^{-0.99} dx
 * from https://twitter.com/HideOgata/status/614639138688512000
 */

struct Func1 {
	template <class T> T operator()(const T& x) {
		static const T c("-0.99");
		return cos(x) * pow(x-(-1), c) * pow(1-x, c);
	}
};

struct Func1_s1 {
	template <class T> T operator()(const T& x) {
		static const T c("-0.99");
		return cos(x) * pow(1-x, c);
	}
};

struct Func1_s2 {
	template <class T> T operator()(const T& x) {
		static const T c("-0.99");
		return cos(x) * pow(x-(-1), c);
	}
};

int main() {
	std::cout.precision(17);

	std::cout <<
	kv::defint_power_autostep(Func1(), Func1_s1(), itv(-1.), itv(0.), 12, itv("-0.99"))
	+ kv::defint_power_autostep_r(Func1(), Func1_s2(), itv(0.), itv(1.), 12, itv("-0.99"))
	<< "\n";
}
