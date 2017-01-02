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

	std::cout << kv::digamma(kv::interval<double>(1.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(10.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(-3.5)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(0.5, 3.)) << "\n";
	std::cout << kv::digamma(kv::interval<double>(-0.9, -0.1)) << "\n";

	std::cout << kv::trigamma(kv::interval<double>(1.5)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(2.5)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(4.5)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(-4.5)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(1., 2.5)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(-0.6, -0.4)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(0., 1.)) << "\n";
	std::cout << kv::trigamma(kv::interval<double>(-1., -0.5)) << "\n";

	std::cout << kv::digamma_zero(1.5) << "\n";
	for (i=1; i<=300; i++) {
		std::cout << kv::digamma_zero(-(double)i+0.5) << "\n";
	}
}
