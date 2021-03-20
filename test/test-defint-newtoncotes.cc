#include <iostream>
#include <kv/defint-newtoncotes.hpp>
#ifdef BENCH
#include <boost/timer.hpp>
#endif

typedef kv::interval<double> itv;


struct Func {
	template <class T> T operator() (const T& x) {
		return 1. / (1. + 10 * x * x);
	}
};

int main() {
	std::cout.precision(17);

	int i;
	itv r;
	#ifdef BENCH
	boost::timer t;
	#endif

	for (i=1; i<=7; i++) {
		std::cout << "[i=" << i << "]\n";
		#ifdef BENCH
		t.restart();
		#endif
		r = kv::defint_newtoncotes(Func(), itv(-1), itv(1), i, 0);
		#ifdef BENCH
		std::cout << "time: " << t.elapsed() << "\n";
		#endif
		std::cout << "width: " << width(r) << "\n";
		std::cout << r << "\n";
	}
}
