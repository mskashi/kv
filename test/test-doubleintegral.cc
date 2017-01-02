#include <iostream>
#include <kv/doubleintegral.hpp>

typedef kv::interval<double> itvd;

class Func {
	public:
	template <class T> T operator() (T x, T y) {
		return 1. / (x * x + 2 * y * y + 1.);
	}
};

template <class F, class TT> class Duffy {
	public:
	F f;
	TT l1, h1, l2, h2;
	Duffy(F f_v, TT l1_v, TT h1_v, TT l2_v, TT h2_v):
	f(f_v), l1(l1_v), h1(h1_v), l2(l2_v), h2(h2_v) {}

	template <class T> T operator() (T v1, T v2) {
		return f(l1 + (h2-v2)*(v1-l1)/(h2-l2), v2) * (h2-v2)/(h2-l2);
	}
};

int main()
{
	std::cout.precision(17);

	std::cout << kv::doubleintegral(Func(), (itvd)(-1.), (itvd)1., (itvd)(-1.), (itvd)1., 8, 20) << "\n";
	std::cout << kv::doubleintegral(Func(), (itvd)(-1.), (itvd)1., (itvd)(-1.), (itvd)1., 8, 20, true) << "\n";

	Duffy<Func,itvd> g(Func(), (itvd)(-1), (itvd)1., (itvd)(-1.), (itvd)1.);

	std::cout << kv::doubleintegral(g, (itvd)(-1.), (itvd)1., (itvd)(-1.), (itvd)1., 8, 20) << "\n";

	std::cout << kv::doubleintegral_triangle(Func(), (itvd)(-1.), (itvd)(-1.), (itvd)(-1.), (itvd)1., (itvd)1., (itvd)(-1.), 8, 20) << "\n";
}
