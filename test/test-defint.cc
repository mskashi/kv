#include <iostream>
#include <kv/defint.hpp>

typedef kv::interval<double> itv;


struct Func {
	template <class T> T operator() (const T& x) {
		return 1./x;
	}
};

struct Kahaner10 {
	template <class T> T operator() (const T& x) {
		return 1./(1. + x);
	}
};

int main() {
	std::cout.precision(17);

	std::cout << kv::defint(Func(), (itv)1., (itv)3., 10, 10) << "\n";
	std::cout << kv::defint(Func(), (itv)3., (itv)1., 10, 10) << "\n";
	std::cout << kv::defint_autostep(Func(), (itv)1., (itv)3., 12) << "\n";
	std::cout << kv::defint_autostep(Func(), (itv)3., (itv)1., 12) << "\n";
	std::cout << kv::defint(Kahaner10(), (itv)0., (itv)1., 10, 10) << "\n";
}
