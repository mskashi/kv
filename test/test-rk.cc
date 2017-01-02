#include <boost/numeric/ublas/io.hpp>
#include <kv/rk.hpp>

namespace ub = boost::numeric::ublas;

struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};

struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

int main()
{
	ub::vector<double> x;

	// x.resize(2);
	// x(0) = 1.; x(1) = 0.;

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;

	std::cout.precision(17);

	// kv::rk(Func(), x, 0., 1.0);
	kv::rk(Lorenz(), x, 0., 0.0675);

	std::cout << x << "\n";
}
