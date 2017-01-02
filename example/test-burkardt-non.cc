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
	allsol(GenRosen(), I);
	cout << t.elapsed() << " sec\n";

	// singular function
	Powell().range(I);
	std::cout << "Powell\n";
	t.restart();
	allsol(Powell(), I, 1, 1e-10);
	cout << t.elapsed() << " sec\n";

	Wood().range(I);
	std::cout << "Wood\n";
	t.restart();
	allsol(Wood(), I);
	cout << t.elapsed() << " sec\n";

	// 重い。SX1でおよそ400秒。解は1つ。
#if 0
	Watson().range(4, I);
	std::cout << "Watson\n";
	t.restart();
	allsol(Watson(), I);
	cout << t.elapsed() << " sec\n";
#endif

	// 解は出るが文献と異なる
	Chebyquad().range(4, I);
	std::cout << "Chebyquad\n";
	t.restart();
	allsol(Chebyquad(), I);
	cout << t.elapsed() << " sec\n";

	Brown().range(6, I);
	std::cout << "Brown\n";
	t.restart();
	allsol(Brown(), I);
	cout << t.elapsed() << " sec\n";

	DBVP().range(20, I);
	std::cout << "DBVP\n";
	t.restart();
	allsol(DBVP(), I);
	cout << t.elapsed() << " sec\n";

	DIntEq().range(8, I);
	std::cout << "DIntEq\n";
	t.restart();
	allsol(DIntEq(), I);
	cout << t.elapsed() << " sec\n";

	VDim().range(3, I);
	std::cout << "VDim\n";
	t.restart();
	allsol(VDim(), I);
	cout << t.elapsed() << " sec\n";

	Broyden().range(15, I);
	std::cout << "Broyden\n";
	t.restart();
	allsol(Broyden(), I);
	cout << t.elapsed() << " sec\n";

	BroydenBand().range(10, I);
	std::cout << "BroydenBand\n";
	t.restart();
	allsol(BroydenBand(), I);
	cout << t.elapsed() << " sec\n";

	Hammarling2x2().range(I);
	std::cout << "Hammarling2x2\n";
	t.restart();
	allsol(Hammarling2x2(), I);
	cout << t.elapsed() << " sec\n";

#if 0
	// 計算時間長すぎ。
	Hammarling3x3().range(I);
	std::cout << "Hammarling3x3\n";
	t.restart();
	allsol(Hammarling3x3(), I);
	cout << t.elapsed() << " sec\n";
#endif

	P17().range(I);
	std::cout << "P17\n";
	t.restart();
	allsol(P17(), I);
	cout << t.elapsed() << " sec\n";

#if 0
	// singular function
	P19().range(I);
	std::cout << "P19\n";
	t.restart();
	allsol(P19(), I, 1e-10);
	cout << t.elapsed() << " sec\n";
#endif

#if 0
	// singular function
	P20().range(I);
	std::cout << "P20\n";
	t.restart();
	allsol(P20(), I);
	cout << t.elapsed() << " sec\n";
#endif

	Chandrasekhar().range(20, I);
	std::cout << "Chandrasekhar\n";
	t.restart();
	allsol(Chandrasekhar(), I);
	cout << t.elapsed() << " sec\n";
}
