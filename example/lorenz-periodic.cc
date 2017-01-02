#include <iostream>
#include <kv/poincaremap.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = -8./3. * x(2) + x(0) * x(1);

		return y;
	}
};

struct LorenzPoincareSection {
	template <class T> T operator() (const ub::vector<T>& x){
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
	kv::PoincareMap<Lorenz,LorenzPoincareSection,double> lopo(lo, lops, (itvd)0.);

	x.resize(4);

	x(0) = -13.7; 
	x(1) = -19.6;
	x(2) = 27.;
	x(3) = 1.56;

	kv::newton(lopo, x);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}

	x(0) = -12.6; 
	x(1) = -16.97;
	x(2) = 27.;
	x(3) = 2.306;

	kv::newton(lopo, x);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}

	x(0) = -13.06; 
	x(1) = -17.99;
	x(2) = 27.;
	x(3) = 3.82;

	// 50次くらいないと通らない。
	lopo = kv::PoincareMap<Lorenz,LorenzPoincareSection,double>(lo, lops, (itvd)0., kv::ode_param<double>().set_order(50));

	kv::newton(lopo, x);

	r = kv::krawczyk_approx(lopo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
