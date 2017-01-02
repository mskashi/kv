#include <iostream>
#include <kv/newton.hpp>

namespace ub = boost::numeric::ublas;

struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}
};

int main()
{
	ub::vector<double> x;

	std::cout.precision(17);

	x.resize(2);
	x(0) = 0.6;
	x(1) = 0.6;

	kv::newton(Func(), x);
	std::cout << x << "\n";

	kv::newton_random(Func(), x);
	std::cout << x << "\n";
}
