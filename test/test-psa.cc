#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>

typedef kv::interval<double> itvd;

int main()
{
	kv::psa<double> a;
	kv::psa<itvd> b;

	a = 1.;
	b = itvd(1., 2.);

	std::cout << a << "\n";
	std::cout << b << "\n";

	a.v.resize(3);
	a.v(0) = 1.;
	a.v(1) = 3.;
	a.v(2) = 2.;
	std::cout << a << "\n";
	std::cout << a * a << "\n";
	std::cout << a / a << "\n";

	kv::psa<itvd>::mode() = 2;
	kv::psa<itvd>::domain() = itvd(0., 1.);

	b.v.resize(3);
	b.v(0) = 1.;
	b.v(1) = 3.;
	b.v(2) = 2.;
	std::cout << b << "\n";
	std::cout << b * b << "\n";
	std::cout << b / b << "\n";
}
