#include <kv/ddx.hpp>

typedef kv::ddx ddx;

int main()
{
	std::cout.precision(38);
	ddx x, y, z;

	x = 1.;
	y = 3.;
	z = "0.1";
	std::cout << z << "\n";
	std::cout << z * "11.8" << "\n";
	std::cout << x + y << "\n";
	std::cout << x - y << "\n";
	std::cout << x * y << "\n";
	std::cout << x / y << "\n";
	std::cout << sqrt(z) << "\n";
	std::cout << abs(z) << "\n";
	std::cout << pow(x, y) << "\n";
	std::cout << pow(2, y) << "\n";
	std::cout << pow(x, 2) << "\n";
	std::cout << pow(x, 2.) << "\n";
	std::cout << exp(x) << "\n";
	std::cout << log(x) << "\n";
	std::cout << sin(x) << "\n";
	std::cout << cos(x) << "\n";
	std::cout << tan(x) << "\n";
	std::cout << asin(x) << "\n";
	std::cout << acos(x) << "\n";
	std::cout << atan(x) << "\n";
	std::cout << atan2(x, y) << "\n";
	std::cout << sinh(x) << "\n";
	std::cout << cosh(x) << "\n";
	std::cout << tanh(x) << "\n";
	std::cout << asinh(x) << "\n";
	std::cout << acosh(y) << "\n";
	std::cout << atanh(x/y) << "\n";
	std::cout << floor(z + 3) << "\n";
	std::cout << ceil(z + 3) << "\n";
	int i;
	std::cout << frexp(z / 7., &i) << "\n";
	std::cout << i << "\n";
	std::cout << ldexp(z, 5) << "\n";
	std::cout << std::boolalpha;
	std::cout << (z < 0.2) << "\n";
	std::cout << (z < "0.100001") << "\n";

	std::cout << kv::constants<ddx>::pi() << "\n";
	std::cout << kv::constants<ddx>::e() << "\n";
	std::cout << kv::constants<ddx>::ln2() << "\n";
	std::cout << kv::constants<ddx>::str("0.1") << "\n";

	std::cout << std::numeric_limits<ddx>::epsilon() << "\n";
	std::cout << std::numeric_limits<ddx>::infinity() << "\n";
	std::cout << std::numeric_limits<ddx>::max() << "\n";
	std::cout << std::numeric_limits<ddx>::min() << "\n";

	std::cout << std::numeric_limits<ddx>::digits << "\n";
	std::cout << std::numeric_limits<ddx>::digits10 << "\n";
	std::cout << std::numeric_limits<ddx>::radix << "\n";
	std::cout << std::numeric_limits<ddx>::min_exponent << "\n";
	std::cout << std::numeric_limits<ddx>::min_exponent10 << "\n";
	std::cout << std::numeric_limits<ddx>::max_exponent << "\n";
	std::cout << std::numeric_limits<ddx>::max_exponent10 << "\n";
}
