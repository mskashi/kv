#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>
#include <kv/ode-maffine2.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;

/*
  van der Pol equaion
   x'' - K(1-x^2)x'+x = 0
 */

class VDP {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = 0.01 * (1. - x(0)*x(0)) * x(1) - x(0);

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;
	itvd end;
	kv::ode_param<double> p;

	std::cout.precision(17);

	ix.resize(2);
	ix(0) = itvd(1., 1.01);
	ix(1) = itvd(1., 1.01);
	end = 5000.;

	p.set_verbose(1);
	p.set_order(18);
	p.set_ep_reduce(1000);
	p.set_ep_reduce_limit(1100);
	// r = kv::odelong_maffine2(VDP(), ix, itvd(0.), end, p);
	r = odelong_maffine(VDP(), ix, itvd(0.), end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
