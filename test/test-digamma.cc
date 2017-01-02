#include <kv/digamma.hpp>

int main()
{
	std::cout.precision(17);

	std::cout << kv::digamma(kv::interval<double>(1.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(10.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(-3.5)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(0.5, 3.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(-0.9, -0.1)) << "\n";
}
