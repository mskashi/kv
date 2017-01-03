#include <iostream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <kv/vleq.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;

int main()
{
	ub::matrix<itv> a(2, 2);
	ub::matrix<itv> b(2, 1);
	ub::matrix<itv> x(2, 1);

	a(0, 0) = 1.; a(0, 1) = 2.;
	a(1, 0) = 3.; a(1, 1) = 4.;

	b(0, 0) = 5.; b(1, 0) = 6.;

	kv::vleq(a, b, x);

	std::cout.precision(20);
	std::cout << x << "\n";
	std::cout << prod(a, x) << "\n";
}
