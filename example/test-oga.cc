#include <iostream>
#include <kv/strobomap.hpp>
#include <kv/allsol.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


struct Oga {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);
		T tmp;

		y(0) = x(1);
		tmp = 1. + x(1) * x(1);
		y(1) = tmp * sqrt(tmp) * x(0);

		return y;
	}
};

int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix;

	std::cout.precision(17);

	Oga f;

	kv::StroboMap<Oga,double> g(f, (itvd)0., (itvd)10.);

	kv::Shooting_TPBVP< kv::StroboMap<Oga,double>, itvd> h(g, (itvd)1., (itvd)0., 1, 0);

	x.resize(1);
	// [-0.76536687071596288,-0.76536685441430041]
	x(0) = -0.76536686;
	kv::newton(h, x);
	std::cout << x << "\n";

	ix.resize(1);
	ix(0) = itvd(-0.77, -0.76);
	kv::allsol(h, ix, 2);
}
