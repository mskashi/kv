#include <iostream>
#include <limits>
#include "ode-nv.hpp"

namespace ub = boost::numeric::ublas;

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

int main()
{
	ub::vector<double> x;

	// x.resize(2);
	// x(0) = 1.; x(1) = 0.;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	double end = std::numeric_limits<double>::infinity();

	// ode_nv(Func(), x, 0., end, 10);
	kv::ode_nv(Lorenz(), x, 0., end, 12);

	std::cout << x << "\n";
	std::cout << end << "\n";

	kv::odelong_nv(Lorenz(), x, 0., 1., 12);
	std::cout << x << "\n";
}
