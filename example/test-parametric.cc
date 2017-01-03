#include <iostream>
#include <limits>

#include <kv/ode-maffine.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(1);

		y(0) = -itv(4.9,5.1) * x(0);

		return y;
	}
};

struct Func2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1) * x(0);
		y(1) = 0.;

		return y;
	}
};


struct LotkaVolterra {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = itv(2.99, 3.01) * x(0) * (1. - x(1));
		y(1) = itv(0.99, 1.01) * x(1) * (x(0) - 1.);

		return y;
	}
};

struct LotkaVolterra2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
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
	ub::vector<itv> ix;
	bool r;

	itv end;
	std::cout.precision(17);

	ix.resize(1);
	ix(0) = 1.;

	end = 1.;
	r = kv::odelong_maffine(Func(), ix, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(2);
	ix(0) = 1.;
	ix(1) = -itv(4.9,5.1);

	end = 1.;
	r = kv::odelong_maffine(Func2(), ix, itv(0.), end);
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
	r = kv::odelong_maffine(LotkaVolterra(), ix, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}

	ix.resize(4);
	ix(0) = 1.2;
	ix(1) = 1.1;
	ix(2) = itv(2.99, 3.01);
	ix(3) = itv(0.99, 1.01);

	end = 10.;
	r = kv::odelong_maffine(LotkaVolterra2(), ix, itv(0.), end);
	if (!r) {
		std::cout << "No Solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
