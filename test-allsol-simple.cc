#include "allsol-simple.hpp"
#include "interval.hpp"
#include "rdouble.hpp"
#include "dd.hpp"
#include "rdd.hpp"
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
		static const T v118 = "11.8";
		return 2.5 * x*x*x - 10.5 * x*x + v118 * x;
	}
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x){
		int n = x.size();
		ub::vector<T> y(n);
		T s;
		int i;

		s = 0.;
		for (i=0; i<n; i++) s += x(i);

		for (i=0; i<n; i++) {
			y(i) = g(x(i)) + s - (i+1.);
		}

		return y;
	}
};

typedef kv::interval<kv::dd> itv;
// typedef kv::interval<double> itv;


int main()
{
	std::cout.precision(32);

	int i;
	boost::timer t;
	ub::vector< itv > I;

	Func f1;

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	allsol(I, f1, 1);
	cout << t.elapsed() << " sec\n";


	Yamamura f2;

	I.resize(5);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	allsol(I, f2, 1);
	cout << t.elapsed() << " sec\n";
}
