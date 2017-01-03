#include <iostream>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/geoseries.hpp>

typedef kv::interval<double> itv;

int main()
{
	std::cout.precision(17);

	std::cout << kv::geoseries(itv("0.4", "0.6")) << "\n";
	std::cout << kv::geoseries(itv("-0.6", "-0.4")) << "\n";
	std::cout << kv::geoseries(itv("-0.8", "0.4")) << "\n";
}
