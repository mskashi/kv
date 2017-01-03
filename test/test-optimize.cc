#include <boost/timer.hpp>
#include <kv/optimize.hpp>

namespace ub = boost::numeric::ublas;


struct Func {
	template <class T> T operator() (const ub::vector<T>& x){
		T tmp, tmp2;
		tmp = x(0) - 1.;
		tmp2 = x(1) - 2.;
		return tmp * tmp + tmp2 * tmp2;
	}
};


// example found in
// http://www.msi.co.jp/nuopt/products/derivation/global/index.html
struct Func2 {
	template <class T> T operator() (const ub::vector<T>& in){
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

struct Rosenbrock {
	template <class T> T operator() (const ub::vector<T>& x){
		T tmp, tmp2;
		tmp = 1. - x(0);
		tmp2 = x(1) - x(0) * x(0);
		return tmp * tmp + 100. * tmp2 * tmp2;
	}
};

// taken from: http://ci.nii.ac.jp/naid/110002373498

struct Levy {
	template <class T> T operator() (const ub::vector<T>& x){
		T tmp, tmp2;
		int i;

		tmp = 0;
		for (i=1; i<=5; i++) {
			tmp += i * cos((i-1)*x(0) + i);
		}
		tmp2 = 0;
		for (i=1; i<=5; i++) {
			tmp2 += i * cos((i+1)*x(1) + i);
		}
		
		return tmp * tmp2;
	}
};


typedef kv::interval<double> itv;

int main()
{
	int i;
	boost::timer t;
	ub::vector<itv> I;
	std::list< ub::vector<itv> > result;
	typename std::list< ub::vector<itv> >::iterator p;

	std::cout.precision(17);

	I.resize(4);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	result = optimize(I, Func2(), 1e-5);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}
	std::cout << t.elapsed() << " sec\n";

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = itv(0., 10.);
	result = optimize(I, Levy(), 1e-5);
	p = result.begin();
	while (p != result.end()) {
		std::cout << *(p++) << "\n";
	}
}
