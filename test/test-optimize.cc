#include <boost/timer.hpp>
#include <kv/optimize.hpp>

namespace ub = boost::numeric::ublas;


class Func {
	public:
	template <class T> T operator() (ub::vector<T> x){
		T tmp, tmp2;
		tmp = x(0) - 1.;
		tmp2 = x(1) - 2.;
		return tmp * tmp + tmp2 * tmp2;
	}
};


// http://www.msi.co.jp/nuopt/products/derivation/global/index.html
// にあった例
class Func2 {
	public:
	template <class T> T operator() (ub::vector<T> in){
		T x, y, z, w;

		x = in(0);
		y = in(1);
		z = in(2);
		w = in(3);

		return 100. * (y - x*x) * (y - x*x) + (1. - x ) * (1. - x)
		    + 90. * (z - w * w) * (z - w * w) + (1. - w) * (1. - w)
		    + 10.1 * (y - 1.) * (y - 1.) + 10.1 * (z - 1.) * (z - 1.)
		    + 19.8 * (y - 1.) * (z - 1.);
	}
};

class Rosenbrock {
	public:
	template <class T> T operator() (ub::vector<T> x){
		T tmp, tmp2;
		tmp = 1. - x(0);
		tmp2 = x(1) - x(0) * x(0);
		return tmp * tmp + 100. * tmp2 * tmp2;
	}
};


int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	std::cout.precision(17);

	I.resize(4);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-10., 10.);

	t.restart();
	optimize(I, Func2(), 1e-12);
	std::cout << t.elapsed() << " sec\n";
}
