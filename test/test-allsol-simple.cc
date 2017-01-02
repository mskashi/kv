#include <kv/allsol-simple.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <boost/timer.hpp>
#ifdef TEST_DD
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#endif

using namespace std;

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
	template <class T> T g(const T& x){
		return 2.5 * x*x*x - 10.5 * x*x + 11.8 * x;
	}
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
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

#ifdef TEST_DD
typedef kv::interval<kv::dd> itv;
#else
typedef kv::interval<double> itv;
#endif


int main()
{
#ifdef TEST_DD
	std::cout.precision(33);
#else
	std::cout.precision(17);
#endif

	int i;
	boost::timer t;
	ub::vector< itv > I;

	Func f1;

	I.resize(2);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	kv::allsol_simple(f1, I, 1);
	cout << t.elapsed() << " sec\n";


	Yamamura f2;

	I.resize(5);
	for (i=0; i<I.size(); i++) I(i) = itv(-10., 10.);

	t.restart();
	kv::allsol_simple(f2, I, 1);
	cout << t.elapsed() << " sec\n";
}
