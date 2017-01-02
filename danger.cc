#include "interval.hpp"
#include "rdouble2.hpp"

typedef kv::interval<double> itvd;

int main()
{
	std::cout.precision(17);
	itvd x, y;

	x = 1.;
	y = 10.;
	std::cout << x / y << "\n";
}
