#include <iostream>
#include <limits>

#include "ode-affine.hpp"
#include "ode-affine-wrapper.hpp"


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

class VdP {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = 10000.* (1. - x(0)*x(0))*x(1) - x(0);

		return y;
	}
};

class Nobi {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = x(0) - x(0)*x(0)*x(0);

		return y;
	}
};

int main()
{
	ub::vector< kv::affine<double> > x;
	bool r;
	int i, j;
	itvd tmp;

	itvd start, end;

	// x.resize(2);
	// x(0) = itvd(1., 1.);
	// x(1) = itvd(0., 0.);

	x.resize(3);
	// x(0) = 15.; x(1) = 15.; x(2) = 36.;
	x(0) = (itvd)15.;
	x(1) = (itvd)15.;
	x(2) = (itvd)36.;
	// x(0) = itvd(15-1e-4, 15+1e-4);
	// x(1) = itvd(15-1e-4, 15+1e-4);
	// x(2) = itvd(36-1e-4, 36+1e-4);

	// x.resize(2);
	// x(0) = (itvd)2.;
	// x(1) = (itvd)0.;

	// x.resize(2);
	// x(0) = (itvd)0.;
	// x(1) = (itvd)4.;
	// x(0) = itvd(0-0.05, 0+0.05);
	// x(1) = itvd(4-0.05, 4+0.05);

	std::cout.precision(17);

	#if 0
	start = 0.;
	while(1) {
		end = std::numeric_limits<double>::infinity();
		// r = ode_affine(Lorenz(), x, start, end, 12);
		r = ode_wrapper(Lorenz(), x, start, end, 12);
		if (!r) {
			std::cout << "No Solution\n";
			break;
		}
		std::cout << "t: " << end << "\n";
		for (j=0; j<x.size(); j++) {
			tmp = to_interval(x(j));
			std::cout << tmp << "\n";
			std::cout << width(tmp) << "\n";
		}
		start = end;
	}
	#endif
	end = std::numeric_limits<double>::infinity();
	r = kv::odelong_affine(Lorenz(), x, start, end, 12, 2, 1);
	r = kv::odelong_wrapper(Lorenz(), x, start, end, 12, 2, 1);
}
