#include <kv/highderiv.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

struct Func {
	template <class T> T operator()(const T& x) {
		return 1 / (1 + x * x);
	}
};

typedef kv::interval<double> itv;

int main()
{
	std::cout << kv::highderiv(Func(), 5., 10) << "\n";
	std::cout << kv::highderiv(Func(), itv(5.,5.5), 10) << "\n";
}
