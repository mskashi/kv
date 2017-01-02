#include <iostream>
#include <kv/strobomap.hpp>

namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;
typedef kv::affine<double> afd;


class Func {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1); y(1) = - x(0);

		return y;
	}
};


int main()
{
	ub::vector<double> x;
	ub::vector<itvd> ix;
	ub::vector<afd> ax;
	ub::vector< kv::autodif<double> > dx;
	ub::vector< kv::autodif<itvd> > dix;
	bool r;

	std::cout.precision(17);

	x.resize(2);

	Func f;

	kv::StroboMap<Func,double> g(f, (itvd)0., (itvd)1.);
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
	dix = kv::autodif<itvd>::init(ix);
	std::cout << h(dx) << "\n";
}