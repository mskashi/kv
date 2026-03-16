#include <kv/eig.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>


namespace ub = boost::numeric::ublas;

typedef kv::interval< kv::dd > itvdd;
typedef kv::interval<double> itv;

int main()
{
	ub::matrix<double> a;
	ub::matrix< kv::complex<double> > v, d;

	std::cout.precision(17);

	a.resize(2,2);

	a(0,0) = 1.;
	a(0,1) = 2.;
	a(1,0) = 3.;
	a(1,1) = 4.;

	// approximation
	eig(a, v, d);
	std::cout << v << "\n";
	std::cout << d << "\n";

	// verified
	ub::vector< kv::complex<itv> > l;
	veig(a, l);
	std::cout << l << "\n";

	// verified (interval matrix)
	ub::matrix<itv> ia;
	ia = a;
	veig(ia, l);
	std::cout << l << "\n";

	// verified (interval dd matrix)
	ub::matrix<itvdd> iadd;
	ub::vector< kv::complex<itvdd> > ldd;
	std::cout.precision(34);
	iadd = a;
	veig(iadd, ldd);
	std::cout << ldd << "\n";

	// a 3x3 test case where eig tends to fail
	a.resize(3,3);
	a(0,0) = 0.; a(0,1) = 1.; a(0,2) = 0.;
	a(1,0) = 1.; a(1,1) = 0.; a(1,2) = 1.;
	a(2,0) = 0.; a(2,1) = 1.; a(2,2) = 0.;

	std::cout.precision(17);
	eig(a, v, d);
	std::cout << d << "\n";

	veig(a, l);
	std::cout << l << "\n";
}
