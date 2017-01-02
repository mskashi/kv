#include <iostream>
#include "defint.hpp"

namespace bn = boost::numeric;
namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> T operator() (T x) {
		return 1./x;
	}
};

class Kahaner10 {
	public:
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
