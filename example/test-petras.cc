#include <iostream>
#include <kv/defint.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;

/*
 * An example found in
 *  Knut Petras: "Principles of verified numerical integration".
 * QUADPACK returns wrong solution.
 */

class Petras {
	public:
	template <class T> T operator() (T x) {
		return 5. * sin(x) + (9.*x-4.)*(9*x-8.)*(3*x-4.)*(9.*x-10.)*(kv::constants<double>::pi() - 2.*x)/(1.+(90.*x-110.)*(90.*x-110.)*(90.*x-110.)*(90.*x-110.));
	}
};

int main() {
	std::cout.precision(17);

	std::cout << kv::defint_autostep(Petras(), (itvd)0., kv::constants<double>::pi(), 12) << "\n";
}
