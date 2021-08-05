#include <kv/allsol.hpp>
#include <chrono>

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
	ub::vector<itv> I;
	std::chrono::system_clock::time_point t;

	std::cout.precision(17);

	I.resize(2);

	for (i=0; i<I.size(); i++) I(i) = itv(-10, 10);
	t = std::chrono::system_clock::now();
	kv::allsol(Func(), I, 2);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	for (i=0; i<I.size(); i++) I(i) = itv::whole();
	t = std::chrono::system_clock::now();
	kv::allsol(Func(), I, 2);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
}
