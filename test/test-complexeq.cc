/*
 * Example for solving complex equation by converting to real equation
 */

#include <iostream>
#include "kv/complex.hpp"
#include "kv/kraw-approx.hpp"
#include "kv/allsol.hpp"

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;

struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
		ub::vector<T> y(2);

		y(0) = pow(x(0), 5) * pow(x(1), 9) + pow(1 - x(1), 7);
		y(1) = pow(x(0), 2) * pow(x(1), 2) - (1 - x(0)) * (1 - x(1));

		return y;
	}
};

// convert n-n complex function to 2n-2n real function

template <class F> struct ComplexReal {
	F f;
	ComplexReal(F f) : f(f) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		int n = x.size();
		int m = n / 2;
		int i;
		ub::vector<T> y(n);
		ub::vector< kv::complex<T> > xc(m);
		ub::vector< kv::complex<T> > yc(m);

		for (i=0; i<m; i++) {
			xc(i) = kv::complex<T>(x(i*2), x(i*2+1));
		}

		yc = f(xc);

		for (i=0; i<m; i++) {
			y(i*2) = yc(i).real();
			y(i*2+1) = yc(i).imag();
		}

		return y;
	}
};

int main()
{
	Func f;
	ComplexReal<Func> g(f);
	ub::vector<itv> I(4);
	ub::vector<double> x(4);

	std::cout.precision(17);

	I(0) = itv(-10., 10.);
	I(1) = itv(1e-5, 10.);
	I(2) = itv(-10., 10.);
	I(3) = itv(1e-5, 10.);

	kv::allsol(g, I, 2);

	x(0) = 0.13;
	x(1) = 0.37;
	x(2) = 4.64;
	x(3) = 1.69;

	if (krawczyk_approx(g, x, I, 2, 1)) {
		std::cout << "solution found\n";
	};
}
