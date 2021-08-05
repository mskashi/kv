#include <kv/allsol.hpp>
#include "burkardt-non.hpp"
#include <chrono>

namespace ub = boost::numeric::ublas;

int main()
{
	int i;
	ub::vector< kv::interval<double> > I;
	std::chrono::system_clock::time_point t;

	GenRosen().range(10, I);
	std::cout << "GenRosen\n";
	t = std::chrono::system_clock::now();
	kv::allsol(GenRosen(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	// has singular solution
	Powell().range(I);
	std::cout << "Powell\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Powell(), I, 1, 1e-10);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Wood().range(I);
	std::cout << "Wood\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Wood(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

#if 0
	// too heavy. takes about 400 sec on SX1. has one solution.
	Watson().range(4, I);
	std::cout << "Watson\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Watson(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
#endif

	// solutions are different from those in the original page
	Chebyquad().range(4, I);
	std::cout << "Chebyquad\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Chebyquad(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Brown().range(6, I);
	std::cout << "Brown\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Brown(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	DBVP().range(20, I);
	std::cout << "DBVP\n";
	t = std::chrono::system_clock::now();
	kv::allsol(DBVP(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	DIntEq().range(8, I);
	std::cout << "DIntEq\n";
	t = std::chrono::system_clock::now();
	kv::allsol(DIntEq(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	VDim().range(3, I);
	std::cout << "VDim\n";
	t = std::chrono::system_clock::now();
	kv::allsol(VDim(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Broyden().range(15, I);
	std::cout << "Broyden\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Broyden(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	BroydenBand().range(10, I);
	std::cout << "BroydenBand\n";
	t = std::chrono::system_clock::now();
	kv::allsol(BroydenBand(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Hammarling2x2().range(I);
	std::cout << "Hammarling2x2\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Hammarling2x2(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

#if 0
	// too heavy
	Hammarling3x3().range(I);
	std::cout << "Hammarling3x3\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Hammarling3x3(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
#endif

	P17().range(I);
	std::cout << "P17\n";
	t = std::chrono::system_clock::now();
	kv::allsol(P17(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

#if 0
	// has singular solution
	P19().range(I);
	std::cout << "P19\n";
	t = std::chrono::system_clock::now();
	kv::allsol(P19(), I, 1e-10);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
#endif

#if 0
	// has singular solution
	P20().range(I);
	std::cout << "P20\n";
	t = std::chrono::system_clock::now();
	kv::allsol(P20(), I, 1e-10);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
#endif

	Chandrasekhar().range(20, I);
	std::cout << "Chandrasekhar\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Chandrasekhar(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
}
