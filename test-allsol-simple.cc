#include "allsol-simple.hpp"
#include "interval.hpp"
#include "rdouble2.hpp"
#include <boost/timer.hpp>

using namespace std;

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


class Yamamura {
	template <class T> T g(T x){
		return 2.5 * x*x*x - 10.5 * x*x + 11.8 * x;
	}
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		int n = x.size();
		ub::vector<T> y(n);
		int i;

		T s = 0.;
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

	Func f1;

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-10., 10.);

	t.restart();
	allsol(I, f1, 1);
	cout << t.elapsed() << " sec\n";


	Yamamura f2;

	I.resize(5);
	for (i=0; i<I.size(); i++) I(i) = kv::interval<double>(-10., 10.);

	t.restart();
	allsol(I, f2, 1);
	cout << t.elapsed() << " sec\n";
}
