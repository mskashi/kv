#include <kv/strobomap.hpp>
#include <kv/kraw-approx.hpp>
#include <kv/allsol.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


template <class TT> class Duffing {
	public:
	TT B;

	Duffing(TT B_v): B(B_v) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = -0.1 * x(1) - x(0)*x(0)*x(0) + B * cos(t);

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

	Duffing<double> f(3.5);

	kv::StroboMap<Duffing<double>,double> g(f, (itvd)0., kv::constants<double>::pi() * 2.);

	kv::FixedPoint< kv::StroboMap<Duffing<double>,double> > h(g);

	ix.resize(2);
	ix(0) = itvd(-5., 5.);
	ix(1) = itvd(-5., 5.);

	kv::allsol(ix, h, 2);
}
