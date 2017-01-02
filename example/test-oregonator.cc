#include <iostream>
#include <limits>

#include <kv/ode.hpp>
#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Oregonator {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 77.27*(x(1) - x(0)*x(1) + x(0) - 8.375e-06 * x(0)*x(0));
		y(1) = (x(2) - x(0)*x(1) - x(1)) / 77.27;
		y(2) = 0.161 * (x(0) - x(2));

		return y;
	}
};

int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;
	itvd end;

	std::cout.precision(17);

	x.resize(3);
	x(0) = 4.; x(1) = 1.1; x(2) = 4.;

	ix = x;
	end = 4.;
	r = kv::odelong_maffine(Oregonator(), ix, itvd(0.), end, kv::ode_param<double>().set_verbose(1).set_restart_max(10));
	if (r == 0) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
