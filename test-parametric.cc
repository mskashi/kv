#include <iostream>
#include <limits>

#include "ode-maffine.hpp"


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = -itvd(4.9,5.1) * x(0);

		return y;
	}
};

class Func2 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1) * x(0);
		y(1) = 0.;

		return y;
	}
};


class LotkaVolterra {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = itvd(2.99, 3.01) * x(0) * (1. - x(1));
		y(1) = itvd(0.99, 1.01) * x(1) * (x(0) - 1.);

		return y;
	}
};

class LotkaVolterra2 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(4);

		y(0) = x(2) * x(0) * (1. - x(1));
		y(1) = x(3) * x(1) * (x(0) - 1.);
		y(2) = 0.;
		y(3) = 0.;

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;

	itvd end;
	std::cout.precision(17);

	ix.resize(1);
	ix(0) = 1.;

	end = 1.;
	r = kv::odelong_maffine(Func(), ix, itvd(0.), end, 12);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(2);
	ix(0) = 1.;
	ix(1) = -itvd(4.9,5.1);

	end = 1.;
	r = kv::odelong_maffine(Func2(), ix, itvd(0.), end, 12);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(2);
	ix(0) = 1.2;
	ix(1) = 1.1;

	end = 10.;
	r = kv::odelong_maffine(LotkaVolterra(), ix, itvd(0.), end, 12);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(4);
	ix(0) = 1.2;
	ix(1) = 1.1;
	ix(2) = itvd(2.99, 3.01);
	ix(3) = itvd(0.99, 1.01);

	end = 10.;
	r = kv::odelong_maffine(LotkaVolterra2(), ix, itvd(0.), end, 12);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
