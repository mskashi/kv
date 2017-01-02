#include <iostream>
#include <kv/strobomap.hpp>
#include <kv/allsol.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = - x(0)*x(0)*x(0) - 3.;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix, ix2;
	itvd end;
	bool r;

	std::cout.precision(17);

	Func f;

	kv::StroboMap<Func,double> g(f, (itvd)0., (itvd)1.);

	kv::Shooting_TPBVP< kv::StroboMap<Func,double>, double> h(g, 0., 0., 0, 0);

	ix.resize(1);
	ix(0) = itvd(-700., 700.);
	kv::allsol(ix, h, 2);

	// x.resize(1);
	// kv::rand_newton(x, h, -100., 100.);
	// std::cout << x << "\n";

	// x.resize(1);
	// x(0) = -9.68;
	// x(0) = 0.;
	// x(0) = 972.;
	// x(0) = 622.;
	// x(0) = 787.;
	// r = kv::krawczyk_approx(h, x, ix, 5, 1);
	// if (r) std::cout << ix(0) << "\n";

	#if 0
	ix2.resize(2);
	ix2(0) = 0.;
	ix2(1) = ix(0);
	end = 1.;
	kv::odelong_maffine(f, ix2, itvd(0.), end);
	#endif
}
