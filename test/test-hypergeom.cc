#include <kv/hypergeom.hpp>

typedef kv::interval<double> itv;

int main()
{
	std::cout.precision(17);
	std::cout << kv::hypergeom(itv(2), itv(3), itv(4), itv(0.5)) << "\n";
	std::cout << kv::hypergeom(itv(2), itv(3), itv(4), itv(-0.5)) << "\n";
}
