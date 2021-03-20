// sample program for "interval.hpp"

#include <kv/interval.hpp>

// define rounding operations for double
// if "rdouble.hpp" is not included, result of computation is not "verified" 
#include <kv/rdouble.hpp>

typedef kv::interval<double> itv;

int main()
{
	itv x;
	// itv y = 1.; // error. constructor is "explicit".
	itv y = (itv)1.;
	itv z(1.);

	x = 1.; // assignment opetator is overloaded.
	y = 10.;
	z = x / y;

	// printed number is rounded to "outward" when "rdouble.hpp" is included.
	std::cout << z << "\n";
	std::cout.precision(17);
	std::cout << z << "\n";

	// copy
	x = itv(1., 2.);
	y = itv(3., 4.);
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

	z = itv(3., 4.);
	// access to endpoints
	std::cout << z.lower() << "\n";
	std::cout << z.upper() << "\n";
	z.lower() = 3.5;
	std::cout << z << "\n";

	// static functions
	// whole and hull
	std::cout << itv::whole() << "\n";
	std::cout << itv::hull(1., 2.) << "\n";
	std::cout << itv::hull(2., 1.) << "\n";
	std::cout << itv::hull(1., z) << "\n";
	std::cout << itv::hull(z, 1.) << "\n";
	std::cout << itv::hull(x, y) << "\n";

	// friend functions
	std::cout << width(z) << "\n";
	std::cout << rad(z) << "\n";
	std::cout << median(z) << "\n";
	std::cout << mid(z) << "\n";
	std::cout << norm(z) << "\n";
	std::cout << mag(z) << "\n";
	std::cout << std::boolalpha;
	std::cout << in(3.9, z) << "\n";
	std::cout << in(4.1, z) << "\n";
	std::cout << subset(x, y) << "\n";
	std::cout << subset(y, z) << "\n";
	std::cout << proper_subset(y, z) << "\n";
	std::cout << overlap(x, y) << "\n";
	std::cout << overlap(y, z) << "\n";
	std::cout << intersect(y, z) << "\n";
	std::cout << abs(itv(2., 3.)) << "\n";
	std::cout << abs(itv(-2., 3.)) << "\n";
	std::cout << abs(itv(-2., -1.)) << "\n";
	std::cout << mig(itv(-2., 1.)) << "\n";
	std::cout << mig(itv(-2., -1.)) << "\n";
	std::cout << mig(itv(2., 3.)) << "\n";
	std::cout << max(itv(2., 4.), itv(3., 5.)) << "\n";
	std::cout << min(itv(2., 4.), itv(3., 5.)) << "\n";
	std::cout << floor(itv(2.5, 3.5)) << "\n";
	std::cout << ceil(itv(2.5, 3.5)) << "\n";

	// test for midrad
	std::cout << "test for midrad\n";
	itv::base_type m, r; // double
	z = itv(1., 1 + std::numeric_limits<itv::base_type>::epsilon() * 3);
	std::cout << mid(z) << "\n";
	std::cout << rad(z) << "\n";
	midrad(z, m, r);
	std::cout << m << "\n";
	std::cout << r << "\n";

	// comparison operators
	std::cout << (x < y) << "\n";
	std::cout << (x < 1.) << "\n";
	std::cout << (1. < x) << "\n";
	std::cout << (itv(1.) == 1.) << "\n";
	std::cout << (itv(1., 2.) != 3.) << "\n";

	// division_part1, division_part2
	// calculate X / (Y \ 0). the result may be divided into two part.
	bool parted;
	std::cout << division_part1(itv(1., 2.), itv(-3., 4.), parted) << "\n";
	// if the result has division_part2, parted is set to true.
	if (parted) std::cout << division_part2(itv(1., 2.), itv(-3., 4.)) << "\n";

	// integer power
	std::cout << pow(itv(2., 3.), 2) << "\n";
	std::cout << pow(itv(2., 3.), -2) << "\n";
	std::cout << pow(itv(-2., 3.), 2) << "\n";
	std::cout << pow(itv(-2., 3.), 3) << "\n";
	// general power
	std::cout << pow(itv(2., 3.), itv(2., 3)) << "\n";
	// mathematical functions
	std::cout << sqrt(itv(2.5, 3.5)) << "\n";
	std::cout << exp(itv(2.5, 3.5)) << "\n";
	// itv::base_type == double
	std::cout << exp(itv(-std::numeric_limits<itv::base_type>::infinity(), 0.)) << "\n";
	std::cout << exp(itv(0., std::numeric_limits<itv::base_type>::infinity())) << "\n";
	std::cout << expm1(itv(-0.25, 0.25)) << "\n";
	std::cout << log(itv(0.75, 1.25)) << "\n";
	std::cout << log(itv(1., std::numeric_limits<itv::base_type>::infinity())) << "\n";
	std::cout << log(itv(0., 1.)) << "\n";
	std::cout << log1p(itv(-0.25, 0.25)) << "\n";
	std::cout << sin(itv(-0.25, 0.25)) << "\n";
	std::cout << cos(itv(-0.25, 0.25)) << "\n";
	std::cout << tan(itv(-0.25, 0.25)) << "\n";
	std::cout << atan(itv(-0.25, 0.25)) << "\n";
	std::cout << asin(itv(-0.25, 0.25)) << "\n";
	std::cout << acos(itv(-0.25, 0.25)) << "\n";
	std::cout << atan2(itv(1.), itv(1.)) << "\n";
	std::cout << sinh(itv(-0.25, 0.25)) << "\n";
	std::cout << cosh(itv(-0.25, 0.25)) << "\n";
	std::cout << tanh(itv(-0.25, 0.25)) << "\n";
	std::cout << asinh(itv(-0.25, 0.25)) << "\n";
	std::cout << acosh(itv(1.5, 2.)) << "\n";
	std::cout << atanh(itv(-0.25, 0.25)) << "\n";

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
	std::cout << kv::constants<itv>::pi() << "\n";
	std::cout << kv::constants<itv>::e() << "\n";
	std::cout << kv::constants<itv>::ln2() << "\n";
	std::cout << kv::constants<itv>::str("0.1") << "\n";
	std::cout << kv::constants<itv>::str("0.1", "0.2") << "\n";

	// numeric constants (double)
	std::cout << kv::constants<itv::base_type>::pi() << "\n";
	std::cout << kv::constants<itv::base_type>::e() << "\n";
	std::cout << kv::constants<itv::base_type>::ln2() << "\n";
	std::cout << kv::constants<itv::base_type>::str("0.1") << "\n";
}
