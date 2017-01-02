#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class NeuEX2 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = 9.9*x(0) - 7.6 * x(1) + 7.6;
		y(1) = 12.6 * x(0) - 9.9 * x(1) + 9.9 + x(0)*x(0)*x(0)/7.6;

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	ix.resize(2);
	ix(0) = 1.;
	ix(1) = 1.;
	for (i=0; i<8; i++) ix(0) /= 2;

	std::cout.precision(17);

	// end = std::numeric_limits<double>::infinity();
	// end = 30;
	end = 300;
	r = kv::odelong_maffine(NeuEX2(), ix, itvd(0.), end, kv::ode_param<double>().set_verbose(1).set_order(18));
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		for (i=0; i<ix.size(); i++) {
			std::cout << ix(i) << "\n";
		}
		std::cout << end << "\n";
	}
}
