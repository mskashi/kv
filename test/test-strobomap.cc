#include <iostream>
#include <kv/strobomap.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itv;
typedef kv::affine<double> aff;


struct Func {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};


int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix;
	ub::vector<aff> ax;
	ub::vector< kv::autodif<double> > dx;
	ub::vector< kv::autodif<itv> > dix;

	std::cout.precision(17);

	x.resize(2);

	Func f;

	kv::StroboMap<Func,double> g(f, (itv)0., (itv)1.);
	kv::FixedPoint< kv::StroboMap<Func,double> > h(g);

	x(0) = 1.;
	x(1) = 1.;
	std::cout << h(x) << "\n";
	dx = kv::autodif<double>::init(x);
	std::cout << h(dx) << "\n";
	ix = x;
	std::cout << h(ix) << "\n";
	ax = x;
	std::cout << h(ax) << "\n";
	dix = kv::autodif<itv>::init(ix);
	std::cout << h(dx) << "\n";
}
