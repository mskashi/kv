#include <iostream>
#include <kv/dka.hpp>

int main()
{
	kv::psa< kv::complex<double> > a;

	// -6 + 11t -6t^2 + t^3
	a.v.resize(4);
	a.v(0) = -6.;
	a.v(1) = 11.;
	a.v(2) = -6.;
	a.v(3) = 1.;

	std::cout << dka(a) << "\n";
}
