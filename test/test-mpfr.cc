#include <iostream> 
#include <kv/mpfr.hpp>

int main()
{
	kv::mpfr<106> x;
	kv::mpfr<106> y(1.);;
	kv::mpfr<106> z("3.14");;
	std::cout.precision(33);

	x = 1.;
	y = "3.14";
	z = x + y;

	std::cout << x << "\n";

	x = 0.;
	std::cout << x << "\n";
	x = -0.1;
	std::cout << x << "\n";
	x = x + 5;
	std::cout << x << "\n";
	x = x + "5";
	std::cout << x << "\n";
	x = 4 + x;
	std::cout << x << "\n";
	x = 4 - x;
	std::cout << x << "\n";
	x = -x;
	std::cout << x << "\n";
	x = x * 0.5;
	std::cout << x << "\n";
	x = sqrt(x);
	std::cout << x << "\n";
	x = floor(x);
	std::cout << x << "\n";
	x = ceil(x);
	std::cout << x << "\n";

	int e;

	x = frexp(x, &e);
	std::cout << x << "\n";
	std::cout << e << "\n";
	std::cout << ldexp(x, 5) << "\n";;

	std::cout << (double)x << "\n";
	std::cout << (int)x << "\n";

	std::cout << std::numeric_limits< kv::mpfr<106> >::epsilon() << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::infinity() << "\n";
	std::cout << - std::numeric_limits< kv::mpfr<106> >::infinity() << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::min() << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::max() << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::digits << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::digits10 << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::min_exponent << "\n";
	std::cout << std::numeric_limits< kv::mpfr<106> >::max_exponent << "\n";
}
