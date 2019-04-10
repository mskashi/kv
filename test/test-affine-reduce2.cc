#include <kv/affine.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::affine<double> aff;
typedef kv::interval<double> itv;

int main()
{
	ub::vector<aff> v;
	int n;

	v.resize(2);

	v(0) = itv(1., 1.5);
	v(1) = itv(2., 2.5);

	v(0) = v(0) * v(1);
	v(1) = v(0) + v(1);

	std::cout << v << "\n";
	std::cout << "maxnum: " << aff::maxnum() << "\n";
	
	n = aff::maxnum();

	v(0) = itv(1,1.25) * v(0) - v(1);
	v(1) = itv(1,1.25) * v(0) + v(1);

	std::cout << v << "\n";
	std::cout << "maxnum: " << aff::maxnum() << "\n";

	epsilon_reduce2(v, n);

	std::cout << v << "\n";
	std::cout << "maxnum: " << aff::maxnum() << "\n";
}
