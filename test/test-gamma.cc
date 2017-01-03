#include <iostream>
#include <kv/gamma.hpp>

typedef kv::interval<double> itv;

int main() {
	int i;
	std::cout.precision(17);

	std::cout << kv::gamma(itv(0.4)) << "\n";
	std::cout << kv::gamma(itv(5.7)) << "\n";
	std::cout << kv::gamma(itv(3.2, 4.4)) << "\n";
	std::cout << kv::gamma(itv(-0.5)) << "\n";
	std::cout << kv::gamma(itv(-0.9, -0.1)) << "\n";

	std::cout << kv::lgamma(itv(0.4)) << "\n";
	std::cout << kv::lgamma(itv(5.7)) << "\n";
	std::cout << kv::lgamma(itv(3.2, 4.4)) << "\n";
	std::cout << kv::lgamma(itv(-0.5)) << "\n";
	std::cout << kv::lgamma(itv(-0.9, -0.1)) << "\n";
	std::cout << kv::lgamma(itv(-1.5)) << "\n";
	std::cout << kv::lgamma(itv(10000.)) << "\n";
	std::cout << kv::lgamma(itv(-4.1, 2.)) << "\n";
	std::cout << kv::lgamma(itv(-2.1, 2.)) << "\n";

	std::cout << kv::digamma(itv(1.)) << "\n";
	std::cout << kv::digamma(itv(10.)) << "\n";
	std::cout << kv::digamma(itv(-3.5)) << "\n";
	std::cout << kv::digamma(itv(0.5, 3.)) << "\n";
	std::cout << kv::digamma(itv(-0.9, -0.1)) << "\n";

	std::cout << kv::trigamma(itv(1.5)) << "\n";
	std::cout << kv::trigamma(itv(2.5)) << "\n";
	std::cout << kv::trigamma(itv(4.5)) << "\n";
	std::cout << kv::trigamma(itv(-4.5)) << "\n";
	std::cout << kv::trigamma(itv(1., 2.5)) << "\n";
	std::cout << kv::trigamma(itv(-0.6, -0.4)) << "\n";
	std::cout << kv::trigamma(itv(0., 1.)) << "\n";
	std::cout << kv::trigamma(itv(-1., -0.5)) << "\n";

	std::cout << kv::digamma_zero(1.5) << "\n";
	for (i=1; i<=300; i++) {
		std::cout << kv::digamma_zero(-(double)i+0.5) << "\n";
	}
}
