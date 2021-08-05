#include <kv/allsol.hpp>
#include "coprin.hpp"
#include <chrono>

namespace ub = boost::numeric::ublas;

int main()
{
	int i;
	ub::vector< kv::interval<double> > I;
	std::chrono::system_clock::time_point t;

	std::cout.precision(17);

	Bellido().range(I);
	std::cout << "Bellido\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Bellido(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Bronstein().range(I);
	std::cout << "Bronstein\n";
	kv::allsol(Bronstein(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<1000; i++) kv::allsol(Bronstein(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 1000 << " sec\n";

	Caprasse().range(I);
	std::cout << "Caprasse\n";
	kv::allsol(Caprasse(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<5; i++) kv::allsol(Caprasse(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 5 << " sec\n";

	Cyclo().range(I);
	std::cout << "Cyclo\n";
	kv::allsol(Cyclo(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<200; i++) kv::allsol(Cyclo(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 200 << " sec\n";

	Celestial().range(I);
	std::cout << "Celestial\n";
	kv::allsol(Celestial(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<10; i++) kv::allsol(Celestial(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 10 << " sec\n";

	Eco9().range(I);
	std::cout << "Eco9\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Eco9(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Freudenstein().range(I);
	std::cout << "Freudenstein\n";
	kv::allsol(Freudenstein(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<10; i++) kv::allsol(Freudenstein(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 10 << " sec\n";

	Geneig().range(I);
	for (i=0; i<I.size(); i++) {
		I(i) = kv::interval<double>(-10., 10.);
	}
	std::cout << "Geneig\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Geneig(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Himmelblau().range(I);
	std::cout << "Himmelblau\n";
	kv::allsol(Himmelblau(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<2000; i++) kv::allsol(Himmelblau(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 2000 << " sec\n";

	Kincox().range(I);
	std::cout << "Kincox\n";
	kv::allsol(Kincox(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<1000; i++) kv::allsol(Kincox(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 1000 << " sec\n";

	Redeco8().range(I);
	std::cout << "Redeco8\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Redeco8(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Stenger().range(I);
	std::cout << "Stenger\n";
	kv::allsol(Stenger(), I);
	t = std::chrono::system_clock::now();
	for (i=0; i<5000; i++) kv::allsol(Stenger(), I, 0);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 / 5000 << " sec\n";

	Yamamura1().range(7, I);
	std::cout << "Yamamura1(7)\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Yamamura1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Math_Maple().range(I);
	std::cout << "Math_Maple\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Math_Maple(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Collins().range(I);
	std::cout << "Collins\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Collins(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Math_Num1().range(I);
	std::cout << "Math_Num1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Math_Num1(), I, 1, 1e-8);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Xu().range(I);
	std::cout << "Xu\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Xu(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Box3().range(I);
	std::cout << "Box3\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Box3(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Rump_univariate().range(I);
	std::cout << "Rump_univariate\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Rump_univariate(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	#if 0
	// too difficult for our system
	AOL_cosh1().range(I);
	std::cout << "AOL_cosh1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(AOL_cosh1(), I, 2);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
	#endif

	AOL_log1().range(I);
	std::cout << "AOL_log1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(AOL_log1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	DiGregorio().range(I);
	std::cout << "DiGregorio\n";
	t = std::chrono::system_clock::now();
	kv::allsol(DiGregorio(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	SMNA90897().range(I);
	std::cout << "SMNA90897\n";
	t = std::chrono::system_clock::now();
	kv::allsol(SMNA90897(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	#if 0
	// too difficult for our system
	SMNA92191().range(I);
	std::cout << "SMNA92191\n";
	t = std::chrono::system_clock::now();
	kv::allsol(SMNA92191(), I, 2);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
	#endif

	Sinxx().range(I);
	std::cout << "Sinxx\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Sinxx(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Sinxx1().range(I);
	std::cout << "Sin1xx\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Sinxx1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	#if 0
	// very difficult; has 8 solutions
	Sjirk_Boon().range(I);
	std::cout << "Sjirk_Boon\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Sjirk_Boon(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
	#endif
}
