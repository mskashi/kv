#include <iostream>
#include <kv/strobomap.hpp>
#include <kv/allsol.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = - x(0)*x(0)*x(0) - 3.;

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix, ix2;
	itv end;
	int r;

	std::cout.precision(17);

	Func f;

	kv::StroboMap<Func,double> g(f, (itv)0., (itv)1.);

	kv::Shooting_TPBVP< kv::StroboMap<Func,double>, double> h(g, 0., 0., 0, 0);

	ix.resize(1);
	ix(0) = itv(-700., 700.);
	kv::allsol(h, ix, 2);

	// x.resize(1);
	// kv::newton_random(h, x);
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
	kv::odelong_maffine(f, ix2, itv(0.), end);
	#endif
}
