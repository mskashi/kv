#include <kv/dd.hpp>

typedef kv::dd dd;

int main()
{
	std::cout.precision(32);
	dd x, y, z;

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
	std::cout << floor(z + 3) << "\n";
	int i;
	std::cout << frexp(z / 7., &i) << "\n";
	std::cout << i << "\n";
	std::cout << std::boolalpha;
	std::cout << (z < 0.2) << "\n";
	std::cout << (z < "0.100001") << "\n";

	std::cout << kv::constants<dd>::pi() << "\n";
	std::cout << kv::constants<dd>::e() << "\n";
	std::cout << kv::constants<dd>::ln2() << "\n";
	std::cout << kv::constants<dd>::str("0.1") << "\n";

	std::cout << std::numeric_limits<kv::dd>::epsilon() << "\n";
	std::cout << std::numeric_limits<kv::dd>::infinity() << "\n";
	std::cout << std::numeric_limits<kv::dd>::max() << "\n";
	std::cout << std::numeric_limits<kv::dd>::min() << "\n";

	std::cout << std::numeric_limits<kv::dd>::digits << "\n";
	std::cout << std::numeric_limits<kv::dd>::digits10 << "\n";
	std::cout << std::numeric_limits<kv::dd>::radix << "\n";
	std::cout << std::numeric_limits<kv::dd>::min_exponent << "\n";
	std::cout << std::numeric_limits<kv::dd>::min_exponent10 << "\n";
	std::cout << std::numeric_limits<kv::dd>::max_exponent << "\n";
	std::cout << std::numeric_limits<kv::dd>::max_exponent10 << "\n";
}
