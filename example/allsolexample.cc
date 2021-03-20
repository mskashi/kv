#include <boost/timer.hpp>
#include <kv/allsol.hpp>
#include "allsolexample.hpp"

using namespace std;

namespace ub = boost::numeric::ublas;


int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	Matsu1().range(I);
	std::cout << "Matsu1\n";
	t.restart();
	allsol(Matsu1(), I);
	cout << t.elapsed() << " sec\n";

	Matsu2().range(I);
	std::cout << "Matsu2\n";
	t.restart();
	allsol(Matsu2(), I);
	cout << t.elapsed() << " sec\n";

	NoSol().range(I);
	std::cout << "NoSol\n";
	t.restart();
	allsol(NoSol(), I);
	cout << t.elapsed() << " sec\n";

	BadCond().range(I);
	std::cout << "BadCond\n";
	t.restart();
	allsol(BadCond(), I);
	cout << t.elapsed() << " sec\n";

	Hansen1().range(I);
	std::cout << "Hansen1\n";
	t.restart();
	allsol(Hansen1(), I);
	cout << t.elapsed() << " sec\n";

	GE1().range(I);
	std::cout << "GE1\n";
	t.restart();
	allsol(GE1(), I);
	cout << t.elapsed() << " sec\n";

	Shinohara1().range(I);
	std::cout << "Shinohara1\n";
	t.restart();
	allsol(Shinohara1(), I);
	cout << t.elapsed() << " sec\n";

	// has multiple root (2,0)
	Shinohara2().range(I);
	std::cout << "Shinohara2\n";
	t.restart();
	allsol(Shinohara2(), I, 1, 1e-8);
	cout << t.elapsed() << " sec\n";

	Shinohara3().range(I);
	std::cout << "Shinohara3\n";
	t.restart();
	allsol(Shinohara3(), I);
	cout << t.elapsed() << " sec\n";

	ModifiedHimmelblau().range(I);
	std::cout << "ModifiedHimmelblau\n";
	t.restart();
	allsol(ModifiedHimmelblau(), I);
	cout << t.elapsed() << " sec\n";

	Heihachiro().range(I);
	std::cout << "Heihachiro\n";
	t.restart();
	allsol(Heihachiro(), I);
	cout << t.elapsed() << " sec\n";

	Yamamura2().range(7, I);
	std::cout << "Yamamura2(7)\n";
	t.restart();
	allsol(Yamamura2(), I);
	cout << t.elapsed() << " sec\n";

	HydroCarbon().range(I);
	std::cout << "HydroCarbon\n";
	t.restart();
	allsol(HydroCarbon(), I);
	cout << t.elapsed() << " sec\n";
}
