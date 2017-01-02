#include <boost/timer.hpp>
#include "allsol.hpp"
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
	allsol(I, GenRosen());
	cout << t.elapsed() << " sec\n";

	// singular function
	Powell().range(I);
	std::cout << "Powell\n";
	t.restart();
	allsol(I, Powell(), 1, 1e-10);
	cout << t.elapsed() << " sec\n";

	Wood().range(I);
	std::cout << "Wood\n";
	t.restart();
	allsol(I, Wood());
	cout << t.elapsed() << " sec\n";

	// 重い。SX1でおよそ400秒。解は1つ。
#if 0
	Watson().range(4, I);
	std::cout << "Watson\n";
	t.restart();
	allsol(I, Watson());
	cout << t.elapsed() << " sec\n";
#endif

	// 解は出るが文献と異なる
	Chebyquad().range(4, I);
	std::cout << "Chebyquad\n";
	t.restart();
	allsol(I, Chebyquad());
	cout << t.elapsed() << " sec\n";

	Brown().range(6, I);
	std::cout << "Brown\n";
	t.restart();
	allsol(I, Brown());
	cout << t.elapsed() << " sec\n";

	DBVP().range(20, I);
	std::cout << "DBVP\n";
	t.restart();
	allsol(I, DBVP());
	cout << t.elapsed() << " sec\n";

	DIntEq().range(8, I);
	std::cout << "DIntEq\n";
	t.restart();
	allsol(I, DIntEq());
	cout << t.elapsed() << " sec\n";

	VDim().range(3, I);
	std::cout << "VDim\n";
	t.restart();
	allsol(I, VDim());
	cout << t.elapsed() << " sec\n";

	Broyden().range(15, I);
	std::cout << "Broyden\n";
	t.restart();
	allsol(I, Broyden());
	cout << t.elapsed() << " sec\n";

	BroydenBand().range(10, I);
	std::cout << "BroydenBand\n";
	t.restart();
	allsol(I, BroydenBand());
	cout << t.elapsed() << " sec\n";

	Hammarling2x2().range(I);
	std::cout << "Hammarling2x2\n";
	t.restart();
	allsol(I, Hammarling2x2());
	cout << t.elapsed() << " sec\n";

#if 0
	// 計算時間長すぎ。
	Hammarling3x3().range(I);
	std::cout << "Hammarling3x3\n";
	t.restart();
	allsol(I, Hammarling3x3());
	cout << t.elapsed() << " sec\n";
#endif

	P17().range(I);
	std::cout << "P17\n";
	t.restart();
	allsol(I, P17());
	cout << t.elapsed() << " sec\n";

#if 0
	// singular function
	P19().range(I);
	std::cout << "P19\n";
	t.restart();
	allsol(I, P19());
	cout << t.elapsed() << " sec\n";
#endif

#if 0
	// singular function
	P20().range(I);
	std::cout << "P20\n";
	t.restart();
	allsol(I, P20());
	cout << t.elapsed() << " sec\n";
#endif

	Chandrasekhar().range(20, I);
	std::cout << "Chandrasekhar\n";
	t.restart();
	allsol(I, Chandrasekhar());
	cout << t.elapsed() << " sec\n";
}
