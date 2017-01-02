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
}
