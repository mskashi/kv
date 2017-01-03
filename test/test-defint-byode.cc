#include <iostream>
#include <limits>
#include <kv/ode.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;


struct Func {
	template <class T> T operator() (T x) {
		return 1./x;
	}
};

struct Kahaner10 {
	template <class T> T operator() (T x) {
		return 1./(1. + x);
	}
};

template <class F> struct DefintByODE {
	F f;
	DefintByODE(F f) : f(f) {}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t) {
		ub::vector<T> y(1);

		y(0) = f(t);

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itv> ix;
	int r;

	itv end;

	ix.resize(1);
	ix(0) = 0.;

	std::cout.precision(17);

	Func f;
	DefintByODE<Func> g(f);

	end = 3.;
	// there is no sense using sophisticated algorithm like ode_maffine
	// because the result of the problem has no dependence on initial value
	r = kv::odelong(g, ix, itv(1.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
