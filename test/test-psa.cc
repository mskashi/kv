#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>

typedef kv::interval<double> itvd;

int main()
{
	kv::psa<itvd> a, b;

	// a = 1  (+ 0t + 0t^2 + ...)
	a = 1.;
	std::cout << a << "\n";

	// a = 1 + 3t + 2t^2
	a.v.resize(3);
	a.v(0) = 1.;
	a.v(1) = 3.;
	a.v(2) = 2.;
	std::cout << a << "\n";

	// b = 2 - t + t^2
	b.v.resize(3);
	b.v(0) = 2.;
	b.v(1) = -1.;
	b.v(2) = 1.;
	std::cout << b << "\n";

	// for basic operations
	std::cout << a + b << "\n";
	std::cout << a - b << "\n";
	std::cout << a * b << "\n";
	std::cout << a / b << "\n";

	// mathematical functions
	std::cout << sqrt(a) << "\n";
	std::cout << exp(a) << "\n";
	std::cout << log(a) << "\n";
	std::cout << sin(a) << "\n";
	std::cout << cos(a) << "\n";

	// integrate
	std::cout << integrate(a) << "\n";

	// set Type-II PSA with domain [0,1]
	kv::psa<itvd>::mode() = 2;
	kv::psa<itvd>::domain() = itvd(0., 1.);

	// for basic operations
	std::cout << a + b << "\n";
	std::cout << a - b << "\n";
	std::cout << a * b << "\n";
	std::cout << a / b << "\n";

	// mathematical functions
	std::cout << sqrt(a) << "\n";
	std::cout << exp(a) << "\n";
	std::cout << log(a) << "\n";
	std::cout << sin(a) << "\n";
	std::cout << cos(a) << "\n";

	// integrate
	std::cout << integrate(a) << "\n";
}
