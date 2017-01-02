#include <iostream>
#include "poincaremap.hpp"
#include "newton.hpp"
#include "kraw-approx.hpp"

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Lorenz {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = -8./3. * x(2) + x(0) * x(1);

		return y;
	}
};

class LorenzPoincareSection {
	public:
	template <class T> T operator() (ub::vector<T> x){
		T y;

		y = x(2) - 27.;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;

	std::cout.precision(17);

	Lorenz lo;
	LorenzPoincareSection lops;
	kv::PoincareMap<Lorenz,LorenzPoincareSection,itvd> lopo(lo, lops, (itvd)0., 12);

	x.resize(4);

	x(0) = -13.7; 
	x(1) = -19.6;
	x(2) = 27.;
	x(3) = 1.56;

	newton(x, lopo);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}

	x(0) = -12.6; 
	x(1) = -16.97;
	x(2) = 27.;
	x(3) = 2.306;

	newton(x, lopo);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}

	x(0) = -13.06; 
	x(1) = -17.99;
	x(2) = 27.;
	x(3) = 3.82;

	// 40次くらいないと通らない。
	lopo = kv::PoincareMap<Lorenz,LorenzPoincareSection,itvd>(lo, lops, (itvd)0., 40);

	newton(x, lopo);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
