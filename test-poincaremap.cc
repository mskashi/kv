#include <iostream>
#include "poincaremap.hpp"
// #include "newton.hpp"
#include "kraw-approx.hpp"

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

class FuncPoincareSection {
	public:
	template <class T> T operator() (ub::vector<T> x){
		T y;

		y = x(0) - 0.5;

		return y;
	}
};

/*
  van der Pol equaion
   x'' - K(1-x^2)x'+x = 0
 */

class VDP {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = 0.01 * (1. - x(0)*x(0)) * x(1) - x(0);

		return y;
	}
};

class VDPPoincareSection {
	public:
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
	ub::vector< kv::autodif<double> > dx;
	ub::vector< kv::autodif<itvd> > dix;
	bool r;

	std::cout.precision(17);


	Func f;

	FuncPoincareSection ps;

	kv::PoincareMap<Func,FuncPoincareSection,itvd> po(f, ps, (itvd)0., 12);

	x.resize(3);
	x(0) = 0;
	x(1) = 1;
	x(2) = 3;

	std::cout << po(x) << "\n";

	dx = kv::autodif<double>::init(x);
	std::cout << po(dx) << "\n";

	ix = x;
	std::cout << po(ix) << "\n";

	dix = kv::autodif<itvd>::init(ix);
	std::cout << po(dix) << "\n";


	VDP vdp;
	VDPPoincareSection vdpps;
	kv::PoincareMap<VDP,VDPPoincareSection,itvd> vdppo(vdp, vdpps, (itvd)0., 12);

	x(0) = 0.; 
	x(1) = 2.;
	x(2) = 6;
	// newton(x, vdppo);

	r = krawczyk_approx(vdppo, x, ix, 10);
	if (r) {
		std::cout << "solution found.\n";
		std::cout << ix << "\n";
	}
}
