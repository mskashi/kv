#include <boost/timer.hpp>
#include <kv/allsol.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;

struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * x(0) + x(1) * x(1) - 1.;
		y(1) = x(0) - x(1);

		return y;
	}
};

int main()
{
	int i;
	boost::timer t;
	ub::vector<itv> I;

	std::cout.precision(17);

	I.resize(2);

	for (i=0; i<I.size(); i++) I(i) = itv(-10, 10);
	t.restart();
	kv::allsol(Func(), I, 2);
	std::cout << t.elapsed() << " sec\n";

	for (i=0; i<I.size(); i++) I(i) = itv::whole();
	t.restart();
	kv::allsol(Func(), I, 2);
	std::cout << t.elapsed() << " sec\n";
}
