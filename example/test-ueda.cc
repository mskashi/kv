#include <kv/strobomap.hpp>
#include <kv/kraw-approx.hpp>
#include <kv/allsol.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Ueda {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
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

	kv::StroboMap<Ueda,itvd> g1(f, (itvd)0., kv::constants<double>::pi() * 2., 12);

	kv::FixedPoint< kv::StroboMap<Ueda,itvd> > h1(g1);

	kv::StroboMap<Ueda,itvd> g2(f, (itvd)0., kv::constants<double>::pi() * 4., 12);

	kv::FixedPoint< kv::StroboMap<Ueda,itvd> > h2(g2);

	kv::StroboMap<Ueda,itvd> g4(f, (itvd)0., kv::constants<double>::pi() * 8., 12);

	kv::FixedPoint< kv::StroboMap<Ueda,itvd> > h4(g4);


#if 0
	ix.resize(2);
	ix(0) = itvd(-1., 1.);
	ix(1) = itvd(-1., 1.);

	// allsol(ix, h1, 2);
	// allsol(ix, h2, 2);
	// allsol(ix, h4, 2);
#endif

#if 0
	x(0) = 0.074937;
	x(1) = 0.254889;
	ix = x;
	std::cout << h2(ix) << "\n";
	dx = init_dif(ix);
	std::cout << h2(dx) << "\n";
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
