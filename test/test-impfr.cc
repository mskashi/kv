// sample program for mpfr-interval
// all test is same as test-interval.cc

#include <kv/interval.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

typedef kv::interval< kv::mpfr<106> > itv;

int main()
{
	itv x;
	itv y = (itv)1.;
	itv z(1.);

	x = 1.;
	y = 10.;
	z = x / y;

	std::cout << z << "\n";
	std::cout.precision(33);
	std::cout << z << "\n";

	x = itv(1., 2.);
	y = itv(3., 4.);
	y.assign(3., 4.);

	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";

	std::cout << x + 1 << "\n";
	std::cout << x + 1. << "\n";

	z += x;
	z += 1.;

	z = itv(3., 4.);
	std::cout << z.lower() << "\n";
	std::cout << z.upper() << "\n";
	z.lower() = 3.5;
	std::cout << z << "\n";

	std::cout << itv::whole() << "\n";
	std::cout << itv::hull(1., 2.) << "\n";
	std::cout << itv::hull(2., 1.) << "\n";
	std::cout << itv::hull(1., z) << "\n";
	std::cout << itv::hull(z, 1.) << "\n";
	std::cout << itv::hull(x, y) << "\n";

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
	itv::base_type m, r;
	z = itv(1., 1 + std::numeric_limits<itv::base_type>::epsilon() * 3);
	std::cout << mid(z) << "\n";
	std::cout << rad(z) << "\n";
	midrad(z, m, r);
	std::cout << m << "\n";
	std::cout << r << "\n";

	std::cout << (x < y) << "\n";
	std::cout << (x < 1.) << "\n";
	std::cout << (1. < x) << "\n";

	bool parted;
	std::cout << division_part1(itv(1., 2.), itv(-3., 4.), parted) << "\n";
	if (parted) std::cout << division_part2(itv(1., 2.), itv(-3., 4.)) << "\n";

	std::cout << pow(itv(2., 3.), 2) << "\n";
	std::cout << pow(itv(2., 3.), -2) << "\n";
	std::cout << pow(itv(-2., 3.), 2) << "\n";
	std::cout << pow(itv(-2., 3.), 3) << "\n";
	std::cout << sqrt(itv(2.5, 3.5)) << "\n";
	std::cout << exp(itv(2.5, 3.5)) << "\n";
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

	x = "0.1";
	std::cout << x << "\n";
	std::cout << "0.1" * x << "\n";
	x += "0.1";
	std::cout << ("0.20001" > x) << "\n";

	std::cout << kv::constants<itv>::pi() << "\n";
	std::cout << kv::constants<itv>::e() << "\n";
	std::cout << kv::constants<itv>::ln2() << "\n";
	std::cout << kv::constants<itv>::str("0.1") << "\n";
	std::cout << kv::constants<itv>::str("0.1", "0.2") << "\n";

	std::cout << kv::constants<itv::base_type>::pi() << "\n";
	std::cout << kv::constants<itv::base_type>::e() << "\n";
	std::cout << kv::constants<itv::base_type>::ln2() << "\n";
	std::cout << kv::constants<itv::base_type>::str("0.1") << "\n";
}
