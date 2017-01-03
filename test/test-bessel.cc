#include <iostream>
#include <kv/bessel.hpp>

typedef kv::interval<double> itv;

int main() {
	std::cout.precision(17);

	std::cout << kv::besselj(1, (itv)3.) << "\n";
	std::cout << kv::besselj(-1, (itv)3.) << "\n";
	std::cout << kv::besselj(2, (itv)3.) << "\n";
	std::cout << kv::besselj(-2, (itv)3.) << "\n";

	std::cout << kv::besselj((itv)0.5, (itv)50.) << "\n";
}
