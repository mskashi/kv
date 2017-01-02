#include <kv/strobomap.hpp>
#include <kv/kraw-approx.hpp>
#include <kv/allsol.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


struct Ueda {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = -0.1 * x(1) - x(0)*x(0)*x(0)*x(0)*x(0) + 0.2 * cos(t) + 0.05;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix, result;
	bool r;

	std::cout.precision(17);

	x.resize(2);

	Ueda f;

	kv::StroboMap<Ueda,double> g1(f, (itvd)0., kv::constants<itvd>::pi() * 2.);

	kv::FixedPoint< kv::StroboMap<Ueda,double> > h1(g1);

	kv::StroboMap<Ueda,double> g2(f, (itvd)0., kv::constants<itvd>::pi() * 4.);

	kv::FixedPoint< kv::StroboMap<Ueda,double> > h2(g2);

	kv::StroboMap<Ueda,double> g4(f, (itvd)0., kv::constants<itvd>::pi() * 8.);

	kv::FixedPoint< kv::StroboMap<Ueda,double> > h4(g4);


#if 0
	ix.resize(2);
	ix(0) = itvd(-1., 1.);
	ix(1) = itvd(-1., 1.);

	// allsol(h1, ix, 2);
	// allsol(h2, ix, 2);
	// allsol(h4, ix, 2);
#endif

	// D2
	std::cout << "2-periodic-1\n";

	x(0) = 0.074937;
	x(1) = 0.254889;

	r = kv::krawczyk_approx(h2, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// I2
	std::cout << "2-periodic-2\n";

	x(0) = -0.760642;
	x(1) = 0.087372;

	r = kv::krawczyk_approx(h2, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// I2
	std::cout << "2-periodic-3\n";

	x(0) = 0.014128;
	x(1) = 0.190107;

	r = kv::krawczyk_approx(h2, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// I4
	std::cout << "4-periodic-1\n";

	x(0) = 0.172614;
	x(1) = 0.195868;

	r = kv::krawczyk_approx(h4, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// I4
	std::cout << "4-periodic-2\n";

	x(0) = 0.703858;
	x(1) = 0.048399;

	r = kv::krawczyk_approx(h4, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// S
	std::cout << "periodic-1\n";

	x(0) = 0.859130;
	x(1) = 0.942437;

	r = kv::krawczyk_approx(h1, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// D 
	std::cout << "periodic-2\n";

	x(0) = -0.849089;
	x(1) = 0.676844;

	r = kv::krawczyk_approx(h1, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}

	// I
	std::cout << "periodic-3\n";

	x(0) = 0.160008;
	x(1) = 0.033718;

	r = kv::krawczyk_approx(h1, x, result);
	if (r) {
		std::cout << result << "\n";
	} else {
		std::cout << "No Solution Found.\n";
	}
}
