// sample program for "interval.hpp"

#include <kv/interval.hpp>

// define rounding operations for double
// if "rdouble.hpp" is not included, result of computation is not "verified" 
#include <kv/rdouble.hpp>

typedef kv::interval<double> itvd;

int main()
{
	itvd x;
	// itvd y = 1.; // error. constructor is "explicit".
	itvd y = (itvd)1.;
	itvd z(1.);

	x = 1.; // assignment opetator is overloaded.
	y = 10.;
	z = x / y;

	// printed number is rounded to "outward" when "rdouble.hpp" is included.
	std::cout << z << "\n";
	std::cout.precision(17);
	std::cout << z << "\n";

	// copy
	x = itvd(1., 2.);
	y = itvd(3., 4.);
	// assign
	y.assign(3., 4.);

	// basic four operations
	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";

	// operation with constant
	std::cout << x + 1 << "\n";
	std::cout << x + 1. << "\n";

	// compound assignment operator
	z += x;
	z += 1.;

	z = itvd(3., 4.);
	// access to endpoints
	std::cout << z.lower() << "\n";
	std::cout << z.upper() << "\n";
	z.lower() = 3.5;
	std::cout << z << "\n";

	// static functions
	// whole and hull
	std::cout << itvd::whole() << "\n";
	std::cout << itvd::hull(1., 2.) << "\n";
	std::cout << itvd::hull(2., 1.) << "\n";
	std::cout << itvd::hull(1., z) << "\n";
	std::cout << itvd::hull(z, 1.) << "\n";
	std::cout << itvd::hull(x, y) << "\n";

	// friend functions
	std::cout << width(z) << "\n";
	std::cout << median(z) << "\n";
	std::cout << norm(z) << "\n";
	std::cout << std::boolalpha;
	std::cout << in(3.9, z) << "\n";
	std::cout << in(4.1, z) << "\n";
	std::cout << subset(x, y) << "\n";
	std::cout << subset(y, z) << "\n";
	std::cout << proper_subset(y, z) << "\n";
	std::cout << overlap(x, y) << "\n";
	std::cout << overlap(y, z) << "\n";
	std::cout << intersect(y, z) << "\n";
	std::cout << abs(itvd(2., 3.)) << "\n";
	std::cout << abs(itvd(-2., 3.)) << "\n";
	std::cout << abs(itvd(-2., -1.)) << "\n";
	std::cout << mig(itvd(-2., 1.)) << "\n";
	std::cout << mig(itvd(-2., -1.)) << "\n";
	std::cout << mig(itvd(2., 3.)) << "\n";

	// comparison operators
	std::cout << (x < y) << "\n";
	std::cout << (x < 1.) << "\n";
	std::cout << (1. < x) << "\n";
	std::cout << (itvd(1.) == 1.) << "\n";
	std::cout << (itvd(1., 2.) != 3.) << "\n";

	// division_part1, division_part2
	// calculate X / (Y \ 0). the result may be divided into two part.
	bool parted;
	std::cout << division_part1(itvd(1., 2.), itvd(-3., 4.), parted) << "\n";
	// if the result has division_part2, parted is set to true.
	if (parted) std::cout << division_part2(itvd(1., 2.), itvd(-3., 4.)) << "\n";

	// integer power
	std::cout << pow(itvd(2., 3.), 2) << "\n";
	std::cout << pow(itvd(2., 3.), -2) << "\n";
	std::cout << pow(itvd(-2., 3.), 2) << "\n";
	std::cout << pow(itvd(-2., 3.), 3) << "\n";
	// general power
	std::cout << pow(itvd(2., 3.), itvd(2., 3)) << "\n";
	// mathematical functions
	std::cout << sqrt(itvd(2.5, 3.5)) << "\n";
	std::cout << exp(itvd(2.5, 3.5)) << "\n";
	std::cout << exp(itvd(-std::numeric_limits<double>::infinity(), 0.)) << "\n";
	std::cout << exp(itvd(0., std::numeric_limits<double>::infinity())) << "\n";
	std::cout << expm1(itvd(-0.25, 0.25)) << "\n";
	std::cout << log(itvd(0.75, 1.25)) << "\n";
	std::cout << log(itvd(1., std::numeric_limits<double>::infinity())) << "\n";
	std::cout << log(itvd(0., 1.)) << "\n";
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

	// string is converted to interval which includes the number represented
	// by the string when "rdouble.hpp" is included.

	// initialize by string. 
	x = "0.1";
	std::cout << x << "\n";
	// operation with string
	std::cout << "0.1" * x << "\n";
	x += "0.1";
	// comparison with string
	std::cout << ("0.20001" > x) << "\n";

	// numeric constants (interval)
	std::cout << kv::constants<itvd>::pi() << "\n";
	std::cout << kv::constants<itvd>::e() << "\n";
	std::cout << kv::constants<itvd>::ln2() << "\n";
	std::cout << kv::constants<itvd>::str("0.1") << "\n";
	std::cout << kv::constants<itvd>::str("0.1", "0.2") << "\n";

	// numeric constants (double)
	std::cout << kv::constants<double>::pi() << "\n";
	std::cout << kv::constants<double>::e() << "\n";
	std::cout << kv::constants<double>::ln2() << "\n";
	std::cout << kv::constants<double>::str("0.1") << "\n";
}
