/*
 * test program for using rest interval of allsol
 *  (rest interval is important for the problem which has multiple root)
 *  -DUNIFY_REST=1 : unify rest interval (default)
 *  -DUNIFY_REST=0 : do not unify rest interval
 */

#include <kv/allsol.hpp>
#include <boost/timer.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;

// from Yoshitane Shinohara: "suuchi kaiseki no kiso", problem 3.8.3
struct Shinohara2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);
		T p, q;

		p = x(0);
		q = x(1);

		y(0) = p*p*p*p*p - 10.*p*p*p*q*q + 5.*p*q*q*q*q
		       - 3.*p*p*p*p + 18.*p*p*q*q - 3.*q*q*q*q
		       - 2.*p*p*p + 6.*p*q*q + 3.*p*p*q - q*q*q
		       + 12.*p*p - 12.*q*q - 10.*p*q - 8.*p + 8.*q;
		y(1) = 5.*p*p*p*p*q - 10.*p*p*q*q*q + q*q*q*q*q
		       - 12.*p*p*p*q + 12.*p*q*q*q - p*p*p + 3.*p*q*q
		       - 6.*p*p*q + 2.*q*q*q + 5.*p*p - 5.*q*q
		       + 24.*p*q - 8.*p - 8.*q + 4.;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i).assign(-3., 3.);
		}
	}
};

int main()
{
	boost::timer t;
	ub::vector<itv> I;
	std::list< ub::vector<itv> > rest;
	std::list< ub::vector<itv> >::iterator p;

	std::cout.precision(17);

	// (2,0) seems to be multiple root
	Shinohara2().range(I);
	std::cout << "Shinohara2\n";
	t.restart();
	kv::allsol(Shinohara2(), I, 1, 1e-8, &rest);
	std::cout << t.elapsed() << " sec\n";

	p = rest.begin();
	while (p != rest.end()) {
		std::cout << *p << "\n";
		p++;
	}
}
