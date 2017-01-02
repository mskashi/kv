#include <iostream>
#include <kv/poincaremap.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


template <class TT> struct Rossler {
	TT a;

	Rossler (TT a_v): a(a_v) {}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = -x(1) - x(2);
		y(1) = x(0) + 0.2 * x(1);
		y(2) = 0.2 + x(2) * (x(0) - a);

		return y;
	}
};

struct RosslerPoincareSection {
	template <class T> T operator() (ub::vector<T> x){
		T y;

		y = x(0) - 0.;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;

	std::cout.precision(17);

	Rossler<double> f(2.2);
	RosslerPoincareSection g;
	kv::PoincareMap<Rossler<double>,RosslerPoincareSection,double> h(f, g, (itvd)0.);

	x.resize(4);

	x(0) = 0.; 
	x(1) = -3.92;
	x(2) = 0.064;
	x(3) = 5.73;

	kv::newton(h, x);

	r = kv::krawczyk_approx(h, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}

	f = Rossler<double>(5.7);
	h = kv::PoincareMap<Rossler<double>,RosslerPoincareSection,double>(f, g, (itvd)0.);
	x(0) = 0.; 
	x(1) = -8.38;
	x(2) = 0.0296;
	x(3) = 5.88;

	kv::newton(h, x);

	r = kv::krawczyk_approx(h, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
