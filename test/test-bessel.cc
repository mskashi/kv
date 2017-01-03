#include <iostream>
#include <kv/bessel.hpp>

typedef kv::interval<double> itv;

int main() {
	std::cout.precision(17);

	std::cout << kv::Bessel(1, (itv)3.) << "\n";
	std::cout << kv::Bessel(-1, (itv)3.) << "\n";
	std::cout << kv::Bessel(2, (itv)3.) << "\n";
	std::cout << kv::Bessel(-2, (itv)3.) << "\n";
}
