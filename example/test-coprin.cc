#include <boost/timer.hpp>
#include <kv/allsol.hpp>

#include "coprin.hpp"

using namespace std;

namespace ub = boost::numeric::ublas;


int main()
{
	int i;
	boost::timer t;
	ub::vector< kv::interval<double> > I;

	std::cout.precision(17);

	Bellido().range(I);
	std::cout << "Bellido\n";
	t.restart();
	allsol(I, Bellido());
	cout << t.elapsed() << " sec\n";

	Bronstein().range(I);
	std::cout << "Bronstein\n";
	allsol(I, Bronstein());
	t.restart();
	for (i=0; i<1000; i++) allsol(I, Bronstein(), 0);
	cout << t.elapsed()/1000. << " sec\n";

	Caprasse().range(I);
	std::cout << "Caprasse\n";
	allsol(I, Caprasse());
	t.restart();
	for (i=0; i<5; i++) allsol(I, Caprasse(), 0);
	cout << t.elapsed()/5. << " sec\n";

	Cyclo().range(I);
	std::cout << "Cyclo\n";
	allsol(I, Cyclo());
	t.restart();
	for (i=0; i<200; i++) allsol(I, Cyclo(), 0);
	cout << t.elapsed()/200. << " sec\n";

	Celestial().range(I);
	std::cout << "Celestial\n";
	allsol(I, Celestial(), 1);
	t.restart();
	for (i=0; i<10; i++) allsol(I, Celestial(), 0);
	cout << t.elapsed()/10. << " sec\n";

	Eco9().range(I);
	std::cout << "Eco9\n";
	t.restart();
	allsol(I, Eco9(), 1);
	cout << t.elapsed() << " sec\n";

	Freudenstein().range(I);
	std::cout << "Freudenstein\n";
	allsol(I, Freudenstein());
	t.restart();
	for (i=0; i<10; i++) allsol(I, Freudenstein(), 0);
	cout << t.elapsed()/10. << " sec\n";

	Geneig().range(I);
	for (i=0; i<I.size(); i++) {
		I(i) = kv::interval<double>(-10., 10.);
	}
	std::cout << "Geneig\n";
	t.restart();
	allsol(I, Geneig());
	cout << t.elapsed() << " sec\n";

	Himmelblau().range(I);
	std::cout << "Himmelblau\n";
	allsol(I, Himmelblau());
	t.restart();
	for (i=0; i<2000; i++) allsol(I, Himmelblau(), 0);
	cout << t.elapsed()/2000. << " sec\n";

	Kincox().range(I);
	std::cout << "Kincox\n";
	allsol(I, Kincox());
	t.restart();
	for (i=0; i<1000; i++) allsol(I, Kincox(), 0);
	cout << t.elapsed()/1000. << " sec\n";

	Redeco8().range(I);
	std::cout << "Redeco8\n";
	t.restart();
	allsol(I, Redeco8());
	cout << t.elapsed() << " sec\n";

	Stenger().range(I);
	std::cout << "Stenger\n";
	allsol(I, Stenger());
	t.restart();
	for (i=0; i<5000; i++) allsol(I, Stenger(), 0);
	cout << t.elapsed()/5000. << " sec\n";

	Yamamura1().range(7, I);
	std::cout << "Yamamura1(7)\n";
	t.restart();
	allsol(I, Yamamura1());
	cout << t.elapsed() << " sec\n";

	Math_Maple().range(I);
	std::cout << "Math_Maple\n";
	t.restart();
	allsol(I, Math_Maple());
	cout << t.elapsed() << " sec\n";

	Collins().range(I);
	std::cout << "Collins\n";
	t.restart();
	allsol(I, Collins());
	cout << t.elapsed() << " sec\n";

	Math_Num1().range(I);
	std::cout << "Math_Num1\n";
	t.restart();
	allsol(I, Math_Num1(), 1, 1e-8);
	cout << t.elapsed() << " sec\n";

	Xu().range(I);
	std::cout << "Xu\n";
	t.restart();
	allsol(I, Xu());
	cout << t.elapsed() << " sec\n";

	Box3().range(I);
	std::cout << "Box3\n";
	t.restart();
	allsol(I, Box3());
	cout << t.elapsed() << " sec\n";

	Rump_univariate().range(I);
	std::cout << "Rump_univariate\n";
	t.restart();
	allsol(I, Rump_univariate());
	cout << t.elapsed() << " sec\n";

	#if 0
	// too difficult for our system
	AOL_cosh1().range(I);
	std::cout << "AOL_cosh1\n";
	t.restart();
	allsol(I, AOL_cosh1(), 2);
	cout << t.elapsed() << " sec\n";
	#endif

	AOL_log1().range(I);
	std::cout << "AOL_log1\n";
	t.restart();
	allsol(I, AOL_log1());
	cout << t.elapsed() << " sec\n";

	DiGregorio().range(I);
	std::cout << "DiGregorio\n";
	t.restart();
	allsol(I, DiGregorio());
	cout << t.elapsed() << " sec\n";

	SMNA90897().range(I);
	std::cout << "SMNA90897\n";
	t.restart();
	allsol(I, SMNA90897());
	cout << t.elapsed() << " sec\n";

	#if 0
	// too difficult for our system
	SMNA92191().range(I);
	std::cout << "SMNA92191\n";
	t.restart();
	allsol(I, SMNA92191(), 2);
	cout << t.elapsed() << " sec\n";
	#endif

	// Bellido().range(I);
	// allsol(I, Bellido());
	// allsol_list(pre_allsol(I, Bellido(), 10000000), Bellido());
}
