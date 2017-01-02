#include <iostream> 
#include "interval.hpp"
#include "rdouble.hpp"
#ifdef TEST_DD
#include "dd.hpp"
#include "rdd.hpp"
#endif
#include "affine.hpp"

#ifdef TEST_DD
typedef kv::interval<kv::dd> itvd;
typedef kv::affine<kv::dd> afd;
#else
typedef kv::interval<double> itvd;
typedef kv::affine<double> afd;
#endif

int main()
{
	afd a, b, c, d;
	int i;

	// initialize from interval
	a = itvd(2., 3.);
	std::cout << a << "\n";
	// initialize from string
	b = "0.1";
	std::cout << b << "\n";
	// initialize from constant (don't use new dummy variable)
	c = 1.;
	std::cout << c << "\n";

	// basic four operations
	d = a + b;
	std::cout << d << "\n";
	d = a - b;
	std::cout << d << "\n";
	d = a * b;
	std::cout << d << "\n";
	d = a / b;
	std::cout << d << "\n";

	// standard functions
	d = sqrt(a);
	std::cout << d << "\n";
	d = square(a);
	std::cout << d << "\n";
	d = exp(a);
	std::cout << d << "\n";
	d = log(a);
	std::cout << d << "\n";
	d = abs(a);
	std::cout << d << "\n";

	// maximin number of dummy variables
	std::cout << afd::maxnum() << "\n";

	// change to interval
	std::cout << to_interval(d) << "\n";
	// get radius
	std::cout << rad(d) << "\n";
	// get midpoint
	std::cout << d.get_mid() << "\n";
	// get special error term
	std::cout << d.get_err() << "\n";

	// get coefficients of dummy variables
	for (i=0; i<=afd::maxnum(); i++) {
		std::cout << i << ":" << d.get_coef(i) << "\n";
	}
}
