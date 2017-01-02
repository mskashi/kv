#include "complex.hpp"
#include "interval.hpp"
#include "rdouble.hpp"
#include "dd.hpp"
#include "rdd.hpp"

// typedef kv::complex<double> cp;
// typedef kv::complex<kv::dd> cp;
// typedef kv::complex< kv::interval<double> > cp;
typedef kv::complex< kv::interval<kv::dd> > cp;

int main()
{
	// std::cout.precision(17);
	std::cout.precision(33);

	cp x, y, z;

	x = cp(1., 2.);
	y = cp(3., 4.);
	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";
	std::cout << abs(x) << "\n";
	std::cout << arg(x) << "\n";
	std::cout << sqrt(x) << "\n";
	std::cout << exp(x) << "\n";
	std::cout << log(x) << "\n";
	std::cout << sin(x) << "\n";
	std::cout << cos(x) << "\n";
	std::cout << tan(x) << "\n";
	std::cout << asin(x) << "\n";
	std::cout << acos(x) << "\n";
	std::cout << atan(x) << "\n";
	std::cout << sinh(x) << "\n";
	std::cout << cosh(x) << "\n";
	std::cout << tanh(x) << "\n";
	std::cout << asinh(x) << "\n";
	std::cout << acosh(x) << "\n";
	std::cout << atanh(x) << "\n";
}
