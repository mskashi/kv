#include <kv/highderiv.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace ub = boost::numeric::ublas;

struct Func {
	template <class T> T operator()(const T& x) {
		return 1 / (1 + x * x);
	}
};

struct Func2 {
	template <class T> T operator()(const ub::vector<T>& x) {
		return 1 / (1 + x(0) * x(0) + 2 * x(1) * x(1));
	}
};

typedef kv::interval<double> itv;

int main()
{
	// example of function highderiv

	std::cout << kv::highderiv(Func(), 5., 10) << "\n";
	std::cout << kv::highderiv(Func(), itv(5.,5.5), 10) << "\n";

	// example of HighDeriv

	kv::HighDeriv<Func> f(Func(), 10);

	std::cout << f(itv(5., 5.5)) << "\n";

	// example of PartialDeriv

	// \frac{\partial^10}{\partial x0^10}
	kv::PartialDeriv<Func2> g0(Func2(), 0, 10);
	// \frac{\partial^10}{\partial x1^10}
	kv::PartialDeriv<Func2> g1(Func2(), 1, 10);

	ub::vector<itv> v(2);
	v(0) = itv(2, 2.5);
	v(1) = itv(2.5, 3);
	std::cout << g0(v) << "\n";
	std::cout << g1(v) << "\n";
}
