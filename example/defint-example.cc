#include <iostream>
#include <kv/defint.hpp>

typedef kv::interval<double> itv;


// simple constant problem

struct Constant_Example {
	template <class T> T operator() (T x) {
		return (T)1.;
	}
};

/*
 *  Knut Petras: "Principles of verified numerical integration".
 * QUADPACK returns wrong solution.
 */

struct Petras {
	template <class T> T operator() (T x) {
		return 5. * sin(x) + (9.*x-4.)*(9*x-8.)*(3*x-4.)*(9.*x-10.)*(kv::constants<T>::pi() - 2.*x)/(1.+(90.*x-110.)*(90.*x-110.)*(90.*x-110.)*(90.*x-110.));
	}
};

/*
 *  http://www.ti3.tuhh.de/intlab/demos/html/dtaylor.html
 */

struct DTaylor {
	template <class T> T operator() (T x) {
		static T pi(kv::constants<itv>::pi());
		return sin(pi * x) - sin(x);
	}
};

struct Sqrt {
	template <class T> T operator() (T x) {
		return sqrt(x);
	}
};

/*
 *  \int_{-1}^1 f(x)dx = pi - 2
 */

struct Func1 {
	template <class T> T operator() (T x) {
		return (1. - x * x) / (1 + x * x);
	}
};

/* 
 * S.M. Rump: Verification methods: Rigorous results using floating-point arithmetic, Acta Numerica, 19, pp. 287-449, 2010
 * http://www.ti3.tu-harburg.de/paper/rump/Ru10.pdf
 * p.86
 */

struct Rump1 {
	template <class T> T operator() (T x) {
		return sin(x + exp(x));
	}
};

struct Nanbu {
	template <class T> T operator() (T x) {
		return log(cos(sqrt(6. * pow(x, 4) + 2. * pow(x, 3) + 1)) + x) / exp(sin(sqrt(pow(x, 2) + 4. * x)) + pow(x, 2));
	}
};


int main() {
	std::cout.precision(17);

	std::cout << kv::defint_autostep(Constant_Example(), (itv)0., (itv)1., 12) << "\n";
	std::cout << kv::defint_autostep(Petras(), (itv)0., kv::constants<itv>::pi(), 12) << "\n";
	std::cout << kv::defint_autostep(DTaylor(), (itv)0., (itv)20., 12) << "\n";
	std::cout << kv::defint_autostep(Sqrt(), (itv)"0.0001", (itv)2., 12) << "\n";
	std::cout << kv::defint_autostep(Func1(), (itv)(-1.), (itv)1., 12) << "\n";
	std::cout << kv::defint_autostep(Rump1(), (itv)(0.), (itv)8., 12) << "\n";
	std::cout << kv::defint_autostep(Nanbu(), (itv)(0.5), (itv)4., 10) << "\n";
}
