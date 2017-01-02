#include <iostream>
#include <kv/defint.hpp>

typedef kv::interval<double> itvd;


struct Func {
	template <class T> T operator() (T x) {
		return 1./x;
	}
};

struct Kahaner10 {
	template <class T> T operator() (T x) {
		return 1./(1. + x);
	}
};

int main() {
	std::cout.precision(17);

	std::cout << kv::defint(Func(), (itvd)1., (itvd)3., 10, 10) << "\n";
	std::cout << kv::defint(Func(), (itvd)3., (itvd)1., 10, 10) << "\n";
	std::cout << kv::defint_autostep(Func(), (itvd)1., (itvd)3., 12) << "\n";
	std::cout << kv::defint_autostep(Func(), (itvd)3., (itvd)1., 12) << "\n";
	std::cout << kv::defint(Kahaner10(), (itvd)0., (itvd)1., 10, 10) << "\n";
}
