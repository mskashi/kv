#include <kv/cardano-ferrari.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST_DD
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
typedef kv::complex< kv::interval<kv::dd> > cp;
#else
typedef kv::complex< kv::interval<double> > cp;
#endif

namespace ub = boost::numeric::ublas;

int main()
{
	ub::vector<cp> in, r;
	int i;

	std::cout.precision(17);

	in.resize(4);
	in(3) = 1.; in(2) = -6.; in(1) = 11.; in(0) = -6.;
	// in(3) = 1.; in(2) = 0.; in(1) = -4.; in(0) = 15.;

	kv::cardano(in, r);

	for (i=0; i<3; i++) {
		std::cout << "solution" << i << ": " << r(i) << "\n";
	}

	in.resize(5);
	in(4) = 1.; in(3) = -10.; in(2) = 35.; in(1) = -50.; in(0) = 24.;
	// in(4) = 1.; in(3) = 2.; in(2) = -30.; in(1) = -4.; in(0) = 15.;
	// in(4) = 1.; in(3) = -4.; in(2) = 6.; in(1) = -4.; in(0) = 1.;
	// in(4) = 1.; in(3) = -2.; in(2) = 0.; in(1) = 2.; in(0) = -1.;
	// in(4) = 1.; in(3) = 0.; in(2) = -5.; in(1) = 0.; in(0) = 4.;
	// in(4) = 1.; in(3) = 0.; in(2) = 2.; in(1) = 0.; in(0) = 1.;
	// in(4) = 1.; in(3) = 2.; in(2) = -2.; in(1) = -6.; in(0) = -3.;
	// in(4) = 1.; in(3) = 0.; in(2) = 3.; in(1) = 0.; in(0) = 2.;

	kv::ferrari(in, r);

	for (i=0; i<4; i++) {
		std::cout << "solution" << i << ": " << r(i) << "\n";
	}
}
