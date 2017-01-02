#include <iostream>
#include <kv/poincaremap.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};

struct FuncPoincareSection {
	template <class T> T operator() (const ub::vector<T>& x){
		T y;

		y = x(0) - 0.5;

		return y;
	}
};

/*
  van der Pol equaion
   x'' - K(1-x^2)x'+x = 0
 */

struct VDP {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = 0.01 * (1. - x(0)*x(0)) * x(1) - x(0);

		return y;
	}
};

struct VDPPoincareSection {
	template <class T> T operator() (const ub::vector<T>& x){
		T y;

		y = x(0) - 0.;

		return y;
	}
};


int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix;
	ub::vector< kv::autodif<double> > dx;
	ub::vector< kv::autodif<itv> > dix;
	bool r;

	std::cout.precision(17);


	Func f;

	FuncPoincareSection ps;

	kv::PoincareMap<Func,FuncPoincareSection,double> po(f, ps, (itv)0.);

	x.resize(3);
	x(0) = 0;
	x(1) = 1;
	x(2) = 3;

	std::cout << po(x) << "\n";

	dx = kv::autodif<double>::init(x);
	std::cout << po(dx) << "\n";

	ix = x;
	std::cout << po(ix) << "\n";

	dix = kv::autodif<itv>::init(ix);
	std::cout << po(dix) << "\n";


	VDP vdp;
	VDPPoincareSection vdpps;
	kv::PoincareMap<VDP,VDPPoincareSection,double> vdppo(vdp, vdpps, (itv)0.);

	x(0) = 0.; 
	x(1) = 2.;
	x(2) = 6;
	kv::newton(vdppo, x);

	r = kv::krawczyk_approx(vdppo, x, ix);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
