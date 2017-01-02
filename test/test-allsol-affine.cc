#include <boost/timer.hpp>
#include <kv/allsol-affine.hpp>

namespace ub = boost::numeric::ublas;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}
};


struct Yamamura {
	template <class T> T g(T x){
		return 2.5 * x*x*x - 10.5 * x*x + 11.8 * x;
	}
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		int n = x.size();
		ub::vector<T> y(n);
		int i;
		T s;

		s = 0.;
		for (i=0; i<n; i++) s += x(i);

		for (i=0; i<n; i++) {
			y(i) = g(x(i)) + s - (i+1.);
		}

		return y;
	}
};


int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	std::cout.precision(17);

	I.resize(7);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-10., 10.);

	t.restart();
	kv::allsol_affine(Yamamura(), I, 2);
	std::cout << t.elapsed() << " sec\n";
}
