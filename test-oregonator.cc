#include <iostream>
#include <limits>

#include "ode.hpp"
#include "ode-maffine.hpp"
#include "ode-maffine2.hpp"


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

class Oregonator {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = 77.27*(x(1) - x(0)*x(1) + x(0) - 8.375e-06 * x(0)*x(0));
		y(1) = (x(2) - x(0)*x(1) - x(1)) / 77.27;
		y(2) = 0.161 * (x(0) - x(2));

		return y;
	}
};

int main()
{
	int i;
	ub::vector<double> x;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	x.resize(3);
	x(0) = 4.; x(1) = 1.1; x(2) = 4.;

	std::cout.precision(17);

	ix = x;
	end = 4.;
	r = kv::odelong_maffine(Oregonator(), ix, itvd(0.), end, 12, 2, 1);
	// r = odelong_maffine2(Oregonator(), ix, itvd(0.), end, 12, 2, 1);
	// r = odelong(Oregonator(), ix, itvd(0.), end, 12, 2, 1);
	/*
	end = 1.;
	r = ode(Oregonator(), ix, itvd(0.), end, 12);
	std::cout << end << "\n";
	*/
	if (r == 0) {
		std::cout << "No Solution\n";
	} else {
		std::cout << r << "\n";
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
