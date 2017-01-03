#include <boost/timer.hpp>
#include <kv/allsol.hpp>

#include "burkardt-non.hpp"

using namespace std;

namespace ub = boost::numeric::ublas;


int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	GenRosen().range(10, I);
	std::cout << "GenRosen\n";
	t.restart();
	kv::allsol(GenRosen(), I);
	cout << t.elapsed() << " sec\n";

	// has singular solution
	Powell().range(I);
	std::cout << "Powell\n";
	t.restart();
	kv::allsol(Powell(), I, 1, 1e-10);
	cout << t.elapsed() << " sec\n";

	Wood().range(I);
	std::cout << "Wood\n";
	t.restart();
	kv::allsol(Wood(), I);
	cout << t.elapsed() << " sec\n";

#if 0
	// too heavy. takes about 400 sec on SX1. has one solution.
	Watson().range(4, I);
	std::cout << "Watson\n";
	t.restart();
	kv::allsol(Watson(), I);
	cout << t.elapsed() << " sec\n";
#endif

	// solutions are different from those in the original page
	Chebyquad().range(4, I);
	std::cout << "Chebyquad\n";
	t.restart();
	kv::allsol(Chebyquad(), I);
	cout << t.elapsed() << " sec\n";

	Brown().range(6, I);
	std::cout << "Brown\n";
	t.restart();
	kv::allsol(Brown(), I);
	cout << t.elapsed() << " sec\n";

	DBVP().range(20, I);
	std::cout << "DBVP\n";
	t.restart();
	kv::allsol(DBVP(), I);
	cout << t.elapsed() << " sec\n";

	DIntEq().range(8, I);
	std::cout << "DIntEq\n";
	t.restart();
	kv::allsol(DIntEq(), I);
	cout << t.elapsed() << " sec\n";

	VDim().range(3, I);
	std::cout << "VDim\n";
	t.restart();
	kv::allsol(VDim(), I);
	cout << t.elapsed() << " sec\n";

	Broyden().range(15, I);
	std::cout << "Broyden\n";
	t.restart();
	kv::allsol(Broyden(), I);
	cout << t.elapsed() << " sec\n";

	BroydenBand().range(10, I);
	std::cout << "BroydenBand\n";
	t.restart();
	kv::allsol(BroydenBand(), I);
	cout << t.elapsed() << " sec\n";

	Hammarling2x2().range(I);
	std::cout << "Hammarling2x2\n";
	t.restart();
	kv::allsol(Hammarling2x2(), I);
	cout << t.elapsed() << " sec\n";

#if 0
	// too heavy
	Hammarling3x3().range(I);
	std::cout << "Hammarling3x3\n";
	t.restart();
	kv::allsol(Hammarling3x3(), I);
	cout << t.elapsed() << " sec\n";
#endif

	P17().range(I);
	std::cout << "P17\n";
	t.restart();
	kv::allsol(P17(), I);
	cout << t.elapsed() << " sec\n";

#if 0
	// has singular solution
	P19().range(I);
	std::cout << "P19\n";
	t.restart();
	kv::allsol(P19(), I, 1e-10);
	cout << t.elapsed() << " sec\n";
#endif

#if 0
	// has singular solution
	P20().range(I);
	std::cout << "P20\n";
	t.restart();
	kv::allsol(P20(), I, 1e-10);
	cout << t.elapsed() << " sec\n";
#endif

	Chandrasekhar().range(20, I);
	std::cout << "Chandrasekhar\n";
	t.restart();
	kv::allsol(Chandrasekhar(), I);
	cout << t.elapsed() << " sec\n";
}
