#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <kv/dka.hpp>

namespace ub = boost::numeric::ublas;

int main()
{
	ub::vector< kv::complex<double> > a, r;
	ub::vector< kv::complex< kv::interval<double> > > va, vr;

	std::cout.precision(17);

	// -6 + 11t -6t^2 + t^3
	a.resize(4);
	a(0) = -6.;
	a(1) = 11.;
	a(2) = -6.;
	a(3) = 1.;

	kv::dka(a, r);
	std::cout << r << "\n";

	va = a;
	kv::vdka(va, vr);
	std::cout << vr << "\n";

	// -100 +150t -104t^2 + 43t^3 -10t^4 + t^5
	a.resize(6);
	a(0) = -100.;
	a(1) = 150.;
	a(2) = -104.;
	a(3) = 43.;
	a(4) = -10.;
	a(5) = 1.;

	kv::dka(a, r);
	std::cout << r << "\n";

	va = a;
	kv::vdka(va, vr);
	std::cout << vr << "\n";

	// 2 -7t +9t^2 -5t^3 +t^4  (t-1)^3*(t-2)
	a.resize(5);
	a(0) = 2.;
	a(1) = -7.;
	a(2) = 9.;
	a(3) = -5.;
	a(4) = 1.;

	kv::dka(a, r);
	std::cout << r << "\n";

	va = a;
	kv::vdka(va, vr);
	std::cout << vr << "\n";

	// -1 +3t -3t^2 + t^3  (t-1)^3
	a.resize(4);
	a(0) = -1.;
	a(1) = 3.;
	a(2) = -3.;
	a(3) = 1.;

	kv::dka(a, r);
	std::cout << r << "\n";

	va = a;
	kv::vdka(va, vr);
	std::cout << vr << "\n";
}
