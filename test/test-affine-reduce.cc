#include <kv/affine.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::affine<double> afd;

int main()
{
	ub::vector<afd> v;

	v.resize(2);

	v(0).a.resize(6);
	v(1).a.resize(6);
	v(0).a(0) = 1.; v(1).a(0) = 2.;
	v(0).a(1) = 5.; v(1).a(1) = 0.;
	v(0).a(2) = 1.; v(1).a(2) = 1.;
	v(0).a(3) = 1.; v(1).a(3) = -1.;
	v(0).a(4) = 0.1; v(1).a(4) = 0.1;
	v(0).a(5) = 0.1; v(1).a(5) = 1.;
	#if AFFINE_SIMPLE >= 1
	v(0).er = 0.; v(1).er = 0.;
	#endif

	afd::maxnum() = 5;

	std::cout << v << "\n";

	epsilon_reduce(v, 4);

	std::cout << v << "\n";
}
