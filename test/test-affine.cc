#include <iostream> 
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#ifdef TEST_DD
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#endif
#include <kv/affine.hpp>

#ifdef TEST_DD
typedef kv::interval<kv::dd> itv;
typedef kv::affine<kv::dd> afd;
#else
typedef kv::interval<double> itv;
typedef kv::affine<double> afd;
#endif

int main()
{
	afd a, b, c, d;
	int i;

	#ifdef TEST_DD
	std::cout.precision(34);
	#else
	std::cout.precision(17);
	#endif

	// initialize from interval
	a = itv(2., 3.);
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

	// basic four operations with constant
	d = a + 2;
	std::cout << d << "\n";
	d = a - 2;
	std::cout << d << "\n";
	d = a * 2;
	std::cout << d << "\n";
	d = a / 2;
	std::cout << d << "\n";
	d = 2 + a;
	std::cout << d << "\n";
	d = 2 - a;
	std::cout << d << "\n";
	d = 2 * a;
	std::cout << d << "\n";
	d = 2 / a;
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
	d = sin(a); // lazy implementation
	std::cout << d << "\n";
	d = cos(a); // lazy implementation
	std::cout << d << "\n";
	d = tan(a); // lazy implementation
	std::cout << d << "\n";
	d = asin(a/10.); // lazy implementation
	std::cout << d << "\n";
	d = acos(a/10.); // lazy implementation
	std::cout << d << "\n";
	d = atan(a); // lazy implementation
	std::cout << d << "\n";
	d = sinh(a); // lazy implementation
	std::cout << d << "\n";
	d = cosh(a); // lazy implementation
	std::cout << d << "\n";
	d = tanh(a); // lazy implementation
	std::cout << d << "\n";
	d = asinh(a); // lazy implementation
	std::cout << d << "\n";
	d = acosh(a); // lazy implementation
	std::cout << d << "\n";
	d = tanh(a/10); // lazy implementation
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
