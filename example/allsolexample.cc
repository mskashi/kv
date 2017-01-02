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
	allsol(I, Matsu1());
	cout << t.elapsed() << " sec\n";

	Matsu2().range(I);
	std::cout << "Matsu2\n";
	t.restart();
	allsol(I, Matsu2());
	cout << t.elapsed() << " sec\n";

	NoSol().range(I);
	std::cout << "NoSol\n";
	t.restart();
	allsol(I, NoSol());
	cout << t.elapsed() << " sec\n";

	BadCond().range(I);
	std::cout << "BadCond\n";
	t.restart();
	allsol(I, BadCond());
	cout << t.elapsed() << " sec\n";

	Hansen1().range(I);
	std::cout << "Hansen1\n";
	t.restart();
	allsol(I, Hansen1());
	cout << t.elapsed() << " sec\n";

	GE1().range(I);
	std::cout << "GE1\n";
	t.restart();
	allsol(I, GE1());
	cout << t.elapsed() << " sec\n";

	Shinohara1().range(I);
	std::cout << "Shinohara1\n";
	t.restart();
	allsol(I, Shinohara1());
	cout << t.elapsed() << " sec\n";

	// (2,0)が重解らしい?
	Shinohara2().range(I);
	std::cout << "Shinohara2\n";
	t.restart();
	allsol(I, Shinohara2(), 1, 1e-8);
	cout << t.elapsed() << " sec\n";

	Shinohara3().range(I);
	std::cout << "Shinohara3\n";
	t.restart();
	allsol(I, Shinohara3());
	cout << t.elapsed() << " sec\n";

	ModifiedHimmelblau().range(I);
	std::cout << "ModifiedHimmelblau\n";
	t.restart();
	allsol(I, ModifiedHimmelblau());
	cout << t.elapsed() << " sec\n";

	Heihachiro().range(I);
	std::cout << "Heihachiro\n";
	t.restart();
	allsol(I, Heihachiro());
	cout << t.elapsed() << " sec\n";
}
