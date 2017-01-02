#include <iostream>
#include <kv/poincaremap.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Craik {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = x(1) * x(2) - x(1);
		y(1) = - x(2) * x(0) + x(2);
		y(2) = x(0) * x(1) - x(0);

		return y;
	}
};

struct CraikPoincareSection {
	template <class T> T operator() (const ub::vector<T>& x){
		T y;

		y = x(0) + x(2) - 1.;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix;
	bool r;

	std::cout.precision(17);

	Craik cr;
	CraikPoincareSection crps;
	kv::PoincareMap<Craik,CraikPoincareSection,double> crpo(cr, crps, (itv)0., kv::ode_param<double>().set_order(30));

	x.resize(4);

	x(0) = 8.043; 
	x(1) = 0.5;
	x(2) = -7.043;
	x(3) = 2.53;

	kv::newton(crpo, x, 1e-8);

	r = kv::krawczyk_approx(crpo, x, ix, 2);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
