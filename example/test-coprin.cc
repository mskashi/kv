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
	allsol(Bellido(), I);
	cout << t.elapsed() << " sec\n";

	Bronstein().range(I);
	std::cout << "Bronstein\n";
	allsol(Bronstein(), I);
	t.restart();
	for (i=0; i<1000; i++) allsol(Bronstein(), I, 0);
	cout << t.elapsed()/1000. << " sec\n";

	Caprasse().range(I);
	std::cout << "Caprasse\n";
	allsol(Caprasse(), I);
	t.restart();
	for (i=0; i<5; i++) allsol(Caprasse(), I, 0);
	cout << t.elapsed()/5. << " sec\n";

	Cyclo().range(I);
	std::cout << "Cyclo\n";
	allsol(Cyclo(), I);
	t.restart();
	for (i=0; i<200; i++) allsol(Cyclo(), I, 0);
	cout << t.elapsed()/200. << " sec\n";

	Celestial().range(I);
	std::cout << "Celestial\n";
	allsol(Celestial(), I);
	t.restart();
	for (i=0; i<10; i++) allsol(Celestial(), I, 0);
	cout << t.elapsed()/10. << " sec\n";

	Eco9().range(I);
	std::cout << "Eco9\n";
	t.restart();
	allsol(Eco9(), I);
	cout << t.elapsed() << " sec\n";

	Freudenstein().range(I);
	std::cout << "Freudenstein\n";
	allsol(Freudenstein(), I);
	t.restart();
	for (i=0; i<10; i++) allsol(Freudenstein(), I, 0);
	cout << t.elapsed()/10. << " sec\n";

	Geneig().range(I);
	for (i=0; i<I.size(); i++) {
		I(i) = kv::interval<double>(-10., 10.);
	}
	std::cout << "Geneig\n";
	t.restart();
	allsol(Geneig(), I);
	cout << t.elapsed() << " sec\n";

	Himmelblau().range(I);
	std::cout << "Himmelblau\n";
	allsol(Himmelblau(), I);
	t.restart();
	for (i=0; i<2000; i++) allsol(Himmelblau(), I, 0);
	cout << t.elapsed()/2000. << " sec\n";

	Kincox().range(I);
	std::cout << "Kincox\n";
	allsol(Kincox(), I);
	t.restart();
	for (i=0; i<1000; i++) allsol(Kincox(), I, 0);
	cout << t.elapsed()/1000. << " sec\n";

	Redeco8().range(I);
	std::cout << "Redeco8\n";
	t.restart();
	allsol(Redeco8(), I);
	cout << t.elapsed() << " sec\n";

	Stenger().range(I);
	std::cout << "Stenger\n";
	allsol(Stenger(), I);
	t.restart();
	for (i=0; i<5000; i++) allsol(Stenger(), I, 0);
	cout << t.elapsed()/5000. << " sec\n";

	Yamamura1().range(7, I);
	std::cout << "Yamamura1(7)\n";
	t.restart();
	allsol(Yamamura1(), I);
	cout << t.elapsed() << " sec\n";

	Math_Maple().range(I);
	std::cout << "Math_Maple\n";
	t.restart();
	allsol(Math_Maple(), I);
	cout << t.elapsed() << " sec\n";

	Collins().range(I);
	std::cout << "Collins\n";
	t.restart();
	allsol(Collins(), I);
	cout << t.elapsed() << " sec\n";

	Math_Num1().range(I);
	std::cout << "Math_Num1\n";
	t.restart();
	allsol(Math_Num1(), I, 1, 1e-8);
	cout << t.elapsed() << " sec\n";

	Xu().range(I);
	std::cout << "Xu\n";
	t.restart();
	allsol(Xu(), I);
	cout << t.elapsed() << " sec\n";

	Box3().range(I);
	std::cout << "Box3\n";
	t.restart();
	allsol(Box3(), I);
	cout << t.elapsed() << " sec\n";

	Rump_univariate().range(I);
	std::cout << "Rump_univariate\n";
	t.restart();
	allsol(Rump_univariate(), I);
	cout << t.elapsed() << " sec\n";

	#if 0
	// too difficult for our system
	AOL_cosh1().range(I);
	std::cout << "AOL_cosh1\n";
	t.restart();
	allsol(AOL_cosh1(), I, 2);
	cout << t.elapsed() << " sec\n";
	#endif

	AOL_log1().range(I);
	std::cout << "AOL_log1\n";
	t.restart();
	allsol(AOL_log1(), I);
	cout << t.elapsed() << " sec\n";

	DiGregorio().range(I);
	std::cout << "DiGregorio\n";
	t.restart();
	allsol(DiGregorio(), I);
	cout << t.elapsed() << " sec\n";

	SMNA90897().range(I);
	std::cout << "SMNA90897\n";
	t.restart();
	allsol(SMNA90897(), I);
	cout << t.elapsed() << " sec\n";

	#if 0
	// too difficult for our system
	SMNA92191().range(I);
	std::cout << "SMNA92191\n";
	t.restart();
	allsol(SMNA92191(), I, 2);
	cout << t.elapsed() << " sec\n";
	#endif

	Sinxx().range(I);
	std::cout << "Sinxx\n";
	t.restart();
	allsol(Sinxx(), I);
	cout << t.elapsed() << " sec\n";

	Sinxx1().range(I);
	std::cout << "Sin1xx\n";
	t.restart();
	allsol(Sinxx1(), I);
	cout << t.elapsed() << " sec\n";

	#if 0
	// very difficult; has 8 solutions
	Sjirk_Boon().range(I);
	std::cout << "Sjirk_Boon\n";
	t.restart();
	allsol(Sjirk_Boon(), I);
	cout << t.elapsed() << " sec\n";
	#endif

	// Bellido().range(I);
	// allsol(Bellido(), I);
	// allsol_list(Bellodo(), pre_allsol(Bellido(), I, 10000000));
}
