#include "interval.hpp"
#include "rdouble.hpp"

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
	std::cout << itvd::hull(z, y) << "\n";
	std::cout << itvd::hull(z, 1.) << "\n";
	std::cout << itvd::hull(1., z) << "\n";
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
	std::cout << pow(itvd(2., 3.), 2) << "\n";
	std::cout << pow(itvd(2., 3.), -2) << "\n";
	std::cout << pow(itvd(-2., 3.), 2) << "\n";
	std::cout << pow(itvd(-2., 3.), 3) << "\n";
	std::cout << exp(itvd(2.5, 3.5)) << "\n";
	std::cout << expm1(itvd(-0.25, 0.25)) << "\n";
	std::cout << log(itvd(0.75, 1.25)) << "\n";
	std::cout << log1p(itvd(-0.25, 0.25)) << "\n";
	std::cout << sin(itvd(-0.25, 0.25)) << "\n";
	std::cout << cos(itvd(-0.25, 0.25)) << "\n";
	std::cout << tan(itvd(-0.25, 0.25)) << "\n";
	std::cout << atan(itvd(-0.25, 0.25)) << "\n";
	std::cout << asin(itvd(-0.25, 0.25)) << "\n";
	std::cout << acos(itvd(-0.25, 0.25)) << "\n";
	std::cout << atan2(itvd(1.), itvd(1.)) << "\n";
	std::cout << sinh(itvd(-0.25, 0.25)) << "\n";
	std::cout << cosh(itvd(-0.25, 0.25)) << "\n";
	std::cout << tanh(itvd(-0.25, 0.25)) << "\n";
	std::cout << asinh(itvd(-0.25, 0.25)) << "\n";
	std::cout << acosh(itvd(1.5, 2.)) << "\n";
	std::cout << atanh(itvd(-0.25, 0.25)) << "\n";

	x = "0.1";
	std::cout << x << "\n";
	std::cout << "0.1" * x << "\n";
}
