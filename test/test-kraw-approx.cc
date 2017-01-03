#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;

struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x) {
		ub::vector<T> y(2);
		y(0) = x(0) * x(0) - x(1) - 1.;
		y(1) = (x(0) - 2.) * (x(0) - 2.) - x(1) - 1.;
		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix;
	bool b;

	std::cout.precision(17);

	x.resize(2);
	ix.resize(2);

	x(0) = 1.01;
	x(1) = 0.01;

	b = kv::krawczyk_approx(Func(), x, ix, 5, 1);
	if (b == false) {
		std::cout << "fail\n";
	}
}
