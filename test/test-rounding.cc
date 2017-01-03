#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

typedef kv::interval<double> itv;

int main()
{
	itv x, y, z;

	x = pow(2. ,54);
	y = 1.;
	z = x + y;
	std::cout << "plus: " << (z.lower() == z.upper() ? "error" : "OK") << "\n";

	x = pow(2. ,54);
	y = -1.;
	z = x - y;
	std::cout << "minus: " << (z.lower() == z.upper() ? "error" : "OK") << "\n";

	x = 1. / 3.;
	y = x;
	z = x * y;
	std::cout << "mult: " << (z.lower() == z.upper() ? "error" : "OK") << "\n";

	x = 1.;
	y = 3.;
	z = x / y;
	std::cout << "div: " << (z.lower() == z.upper() ? "error" : "OK") << "\n";

	x = 2.;
	z = sqrt(x);
	std::cout << "sqrt: " << (z.lower() == z.upper() ? "error" : "OK") << "\n";
}
