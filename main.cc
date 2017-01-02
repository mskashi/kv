#include "interval.hpp"
#include "rdouble2.hpp"

typedef kv::interval<double> itvd;

int main()
{
	std::cout.precision(17);
	std::cout << std::boolalpha;
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

	std::cout << width(z) << "\n";
	std::cout << median(z) << "\n";
	std::cout << norm(z) << "\n";

	std::cout << z.lower() << "\n";
	std::cout << z.upper() << "\n";

	z.lower() = 1.5;
	std::cout << z << "\n";

	const itvd c(1.5);
	std::cout << c.lower() << "\n";

	std::cout << itvd::whole() << "\n";
	std::cout << itvd::hull(1., 2.) << "\n";
	std::cout << itvd::hull(2., 1.) << "\n";
	std::cout << in(1.7, z) << "\n";
	std::cout << in(2.1, z) << "\n";
	std::cout << subset(x, y) << "\n";
	std::cout << subset(z, x) << "\n";
	std::cout << proper_subset(z, x) << "\n";
	std::cout << overlap(x, y) << "\n";
	std::cout << overlap(z, x) << "\n";
	std::cout << intersect(z, x) << "\n";
	std::cout << hull(z, y) << "\n";
	std::cout << hull(z, 1.) << "\n";
	std::cout << hull(1., z) << "\n";
	std::cout << (x < y) << "\n";
	std::cout << (x < 1.) << "\n";
	std::cout << (1. < x) << "\n";
	std::cout << abs(itvd(2., 3.)) << "\n";
	std::cout << abs(itvd(-2., 3.)) << "\n";
	std::cout << abs(itvd(-2., -1.)) << "\n";
	std::cout << mig(itvd(-2., 1.)) << "\n";
	std::cout << mig(itvd(-2., -1.)) << "\n";
	std::cout << mig(itvd(2., 3.)) << "\n";
	bool parted;
	std::cout << division_part1(itvd(1., 2.), itvd(-3., 4.), parted) << "\n";
	if (parted) std::cout << division_part2(itvd(1., 2.), itvd(-3., 4.)) << "\n";
}
