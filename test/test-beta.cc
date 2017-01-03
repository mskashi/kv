#include <kv/beta.hpp>

typedef kv::interval<double> itv;

int main()
{
	std::cout.precision(17);
	std::cout << kv::beta(itv(0.5), itv(1.5)) << "\n";
	std::cout << kv::beta(itv(50.), itv(50.)) << "\n";
}
