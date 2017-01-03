#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>

typedef kv::interval<double> itv;

int main()
{
	kv::psa<itv> a, b;

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

	// Type-I PSA (default)
	std::cout << "Type-I PSA\n";

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
	std::cout << sinh(a) << "\n";
	std::cout << cosh(a) << "\n";
	std::cout << pow(a, 2) << "\n";
	std::cout << pow(a, -2) << "\n";
	std::cout << pow(a, a) << "\n";
	std::cout << tan(a / 100) << "\n";
	std::cout << tanh(a / 100) << "\n";
	std::cout << asin(a / 100) << "\n";
	std::cout << acos(a / 100) << "\n";
	std::cout << atan(a / 100) << "\n";
	std::cout << asinh(a / 100) << "\n";
	std::cout << acosh(a + 2) << "\n";
	std::cout << atanh(a / 100) << "\n";

	// integrate
	std::cout << integrate(a) << "\n";

	// evaluate on some point
	std::cout << eval(a, (itv)0.5) << "\n";

	// set order
	std::cout << setorder(a, 1) << "\n";
	std::cout << setorder(a, 3) << "\n";

	// Type-II PSA
	std::cout << "Type-II PSA\n";

	// set Type-II PSA with domain [0,1]
	kv::psa<itv>::mode() = 2;
	kv::psa<itv>::domain() = itv(0., 1.);

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
	std::cout << sinh(a) << "\n";
	std::cout << cosh(a) << "\n";
	std::cout << pow(a, 2) << "\n";
	std::cout << pow(a, -2) << "\n";
	std::cout << pow(a, a) << "\n";
	std::cout << tan(a / 100) << "\n";
	std::cout << tanh(a / 100) << "\n";
	std::cout << asin(a / 100) << "\n";
	std::cout << acos(a / 100) << "\n";
	std::cout << atan(a / 100) << "\n";
	std::cout << asinh(a / 100) << "\n";
	std::cout << acosh(a + 2) << "\n";
	std::cout << atanh(a / 100) << "\n";

	// integrate
	std::cout << integrate(a) << "\n";

	// evaluate on some point
	std::cout << eval(a, (itv)0.5) << "\n";

	// evaluate on whole range
	std::cout << evalrange(a) << "\n";

	// set order
	std::cout << setorder(a, 1) << "\n";
	std::cout << setorder(a, 3) << "\n";
}
