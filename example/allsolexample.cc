#include <kv/allsol.hpp>
#include "allsolexample.hpp"
#include <chrono>

namespace ub = boost::numeric::ublas;

int main()
{
	int i;
	ub::vector< kv::interval<double> > I;
	std::chrono::system_clock::time_point t;

	Matsu1().range(I);
	std::cout << "Matsu1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Matsu1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Matsu2().range(I);
	std::cout << "Matsu2\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Matsu2(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	NoSol().range(I);
	std::cout << "NoSol\n";
	t = std::chrono::system_clock::now();
	kv::allsol(NoSol(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	BadCond().range(I);
	std::cout << "BadCond\n";
	t = std::chrono::system_clock::now();
	kv::allsol(BadCond(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Hansen1().range(I);
	std::cout << "Hansen1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Hansen1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	GE1().range(I);
	std::cout << "GE1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(GE1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Shinohara1().range(I);
	std::cout << "Shinohara1\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Shinohara1(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	// has multiple root (2,0)
	Shinohara2().range(I);
	std::cout << "Shinohara2\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Shinohara2(), I, 1, 1e-8);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Shinohara3().range(I);
	std::cout << "Shinohara3\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Shinohara3(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	ModifiedHimmelblau().range(I);
	std::cout << "ModifiedHimmelblau\n";
	t = std::chrono::system_clock::now();
	kv::allsol(ModifiedHimmelblau(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Heihachiro().range(I);
	std::cout << "Heihachiro\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Heihachiro(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	Yamamura2().range(7, I);
	std::cout << "Yamamura2(7)\n";
	t = std::chrono::system_clock::now();
	kv::allsol(Yamamura2(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";

	HydroCarbon().range(I);
	std::cout << "HydroCarbon\n";
	t = std::chrono::system_clock::now();
	kv::allsol(HydroCarbon(), I);
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
}
