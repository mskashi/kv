#include <iostream>
#include <limits>
#include <kv/ode-lohner.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Lorenz {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};


int main()
{
	ub::vector<itvd> x;
	ub::vector< kv::autodif<itvd> > xd;
	int r;
	itvd end;

	std::cout.precision(17);

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = std::numeric_limits<double>::infinity();
	r = kv::ode_lohner(Lorenz(), x, itvd(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = 1.;
	r = kv::odelong_lohner(Lorenz(), x, itvd(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	xd = kv::autodif<itvd>::init(x);
	end = std::numeric_limits<double>::infinity();
	r = kv::ode_lohner(Lorenz(), xd, itvd(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << xd << "\n";
		std::cout << end << "\n";
	}

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	xd = kv::autodif<itvd>::init(x);
	end = 1.;
	r = kv::odelong_lohner(Lorenz(), xd, itvd(0.), end);

	if (!r) std::cout << "can't calculate verified solution\n";
	else { 
		std::cout << xd << "\n";
		std::cout << end << "\n";
	}
}
