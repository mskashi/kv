#include <iostream>
#include <limits>
#include "ode.hpp"

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

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

class QuadTest {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = t * t;
		// y(0) = 0.;

		return y;
	}
};

int main()
{
	ub::vector<itvd> x;
	int r;
	itvd end;

	std::cout.precision(17);

	// x.resize(2);
	// x(0) = 1.; x(1) = 0.;
	// r = ode(Func(), x, itvd(0.), end, 20, true);

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	end = std::numeric_limits<double>::infinity();
	r = kv::ode(Lorenz(), x, itvd(0.), end, 20, true);

	// x.resize(1);
	// x(0) = 0.;
	// r = ode(QuadTest(), x, itvd(0.), end, 20, true);

	if (!r) std::cout << "No Solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	end = 1.;
	r = kv::odelong(Lorenz(), x, itvd(0.), end, 12);

	if (!r) std::cout << "No Solution\n";
	else { 
		std::cout << x << "\n";
		std::cout << end << "\n";
	}
}
