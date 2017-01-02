#include <iostream>
#include <limits>
#include <kv/ode-nv.hpp>

namespace ub = boost::numeric::ublas;

struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
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

	std::cout.precision(17);

	x.resize(3);

	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	double end = std::numeric_limits<double>::infinity();

	kv::ode_nv(Lorenz(), x, 0., end);

	std::cout << x << "\n";
	std::cout << end << "\n";

	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	kv::odelong_nv(Lorenz(), x, 0., 1.);
	std::cout << x << "\n";
	std::cout << 1. << "\n";
}
