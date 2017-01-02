#include <iostream>
#include "newton.hpp"

namespace ub = boost::numeric::ublas;

class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}
};

int main()
{
	ub::vector<double> x;

	x.resize(2);
	x(0) = 0.6;
	x(1) = 0.6;

	kv::newton(x, Func());
	std::cout  << x << "\n";

	kv::rand_newton(x, Func());
	std::cout  << x << "\n";
}
