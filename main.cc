#include "interval.hpp"

typedef kv::interval<double> itvd;

int main()
{
	std::cout.precision(17);
	itvd x, y, z;

	x = 1.;
	y = 10.;
	std::cout << x / y << "\n";

	x = itvd(1., 2.);
	y = itvd(3., 4.);

	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";

	z = sqrt(y);
	std::cout << z << "\n";
}
