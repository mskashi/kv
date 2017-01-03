#include <iostream>
#include <kv/doubleintegral.hpp>

typedef kv::interval<double> itv;

struct Func {
	template <class T> T operator() (const T& x, const T& y) {
		return 1. / (x * x + 2 * y * y + 1.);
	}
};

template <class F, class TT> struct Duffy {
	F f;
	TT l1, h1, l2, h2;
	Duffy(F f, TT l1, TT h1, TT l2, TT h2):
	f(f), l1(l1), h1(h1), l2(l2), h2(h2) {}

	template <class T> T operator() (const T& v1, const T& v2) {
		return f(l1 + (h2-v2)*(v1-l1)/(h2-l2), v2) * (h2-v2)/(h2-l2);
	}
};

int main()
{
	std::cout.precision(17);

	std::cout << kv::doubleintegral(Func(), (itv)(-1.), (itv)1., (itv)(-1.), (itv)1., 8, 20) << "\n";
	std::cout << kv::doubleintegral(Func(), (itv)(-1.), (itv)1., (itv)(-1.), (itv)1., 8, 20, true) << "\n";

	Duffy<Func,itv> g(Func(), (itv)(-1), (itv)1., (itv)(-1.), (itv)1.);

	std::cout << kv::doubleintegral(g, (itv)(-1.), (itv)1., (itv)(-1.), (itv)1., 8, 20) << "\n";

	std::cout << kv::doubleintegral_triangle(Func(), (itv)(-1.), (itv)(-1.), (itv)(-1.), (itv)1., (itv)1., (itv)(-1.), 8, 20) << "\n";
}
