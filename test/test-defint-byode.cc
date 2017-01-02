#include <iostream>
#include <limits>

#include <kv/ode.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> T operator() (T x) {
		return 1./x;
	}
};

class Kahaner10 {
	public:
	template <class T> T operator() (T x) {
		return 1./(1. + x);
	}
};

template <class F> class DefintByODE {
	public:
	F f;
	DefintByODE(F f_v) : f(f_v) {}
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t) {
		ub::vector<T> y(1);

		y(0) = f(t);

		return y;
	}
};

int main()
{
	int i;
	ub::vector<itvd> ix;
	bool r;

	itvd end;

	ix.resize(1);
	ix(0) = 0.;

	std::cout.precision(17);

	Func f;
	DefintByODE<Func> g(f);

	// maffineなどの高度な解法は途中経過が初期値に依存しないので意味が無い
	// r = odelong(g, ix, itvd(1.), itvd(3), 12);
	end = 3.;
	r = kv::odelong(DefintByODE<Func>(Func()), ix, itvd(1.), end);
	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << ix << "\n";
		std::cout << end << "\n";
	}
}
