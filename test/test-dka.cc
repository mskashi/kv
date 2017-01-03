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
