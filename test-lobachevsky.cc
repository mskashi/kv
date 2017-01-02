#include <iostream>
#include "lobachevsky.hpp"

int main() {
	int i;
	double x;
	std::cout.precision(17);

	for (i=0; i<=100; i++) {
		x = mid(kv::constants<double>::pi() / 100 * i);
		std::cout << x << " ";
		std::cout << kv::loba(kv::interval<double>(x)) << "\n";
	}
}
