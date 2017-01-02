#include <iostream>
#include <kv/gamma.hpp>

int main() {
	int i;
	double x;
	std::cout.precision(17);

	std::cout << kv::gamma(kv::interval<double>(0.4)) << "\n";
	std::cout << kv::gamma(kv::interval<double>(5.7)) << "\n";
	std::cout << kv::gamma(kv::interval<double>(3.2, 4.4)) << "\n";
	std::cout << kv::gamma(kv::interval<double>(-0.5)) << "\n";
	std::cout << kv::gamma(kv::interval<double>(-0.9, -0.1)) << "\n";
}
