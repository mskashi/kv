#include <boost/timer.hpp>
#include "allsol-affine.hpp"

namespace ub = boost::numeric::ublas;


// 円と直線の交点

class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}
};


// エサキダイオードの怪しい回路。
// ほんとは11.8をちゃんと書く必要あり。

class Func2 {
	template <class T> T g(T x){
		return 2.5 * x*x*x - 10.5 * x*x + 11.8 * x;
	}
	public:
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


class Func3 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		int n = x.size();
		ub::vector<T> y(n);
		int i;

		T z = x(0) * x(0) * x(0);
		for (i=1; i<n; i++) z = z + x(i) * x(i) * x(i);

		for (i=0; i<n; i++) {
			y(i) = x(i) - (z + (i+1.)) * (1./(2.*n));
		}

		return y;
	}
};

int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	I.resize(7);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-10., 10.);

	t.restart();
	kv::allsol_affine(I, Func2(), 2);
	std::cout << t.elapsed() << " sec\n";
}
