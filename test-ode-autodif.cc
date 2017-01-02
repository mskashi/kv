#include <iostream>
#include <limits>
#include "ode-autodif.hpp"

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};

class Lorenz {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = (T)10. * ( x(1) - x(0) );
		y(1) = (T)28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (T)(-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

int main()
{
	ub::vector<itvd> x;
	ub::vector< kv::autodif<itvd> > xd;
	bool r;

	itvd end;

	end = std::numeric_limits<double>::infinity();

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	xd = kv::autodif<itvd>::init(x);

	std::cout.precision(17);

	r = kv::ode(Lorenz(), xd, itvd(0.), end, 12);

	if (!r) std::cout << "No Solution\n";
	else {
		std::cout << xd << "\n";
		std::cout << end << "\n";
	}

	end = 1.;
	r = kv::odelong(Lorenz(), xd, itvd(0.), end, 12);

	if (!r) std::cout << "No Solution\n";
	else {
		std::cout << xd << "\n";
		std::cout << end << "\n";
	}
}
