#include <iostream>
#include <kv/double-newtoncotes.hpp>
#ifdef BENCH
#include <chrono>
#endif

typedef kv::interval<double> itv;


struct Func {
	template <class T> T operator() (const T& x, const T& y) {
		return 1. / (x * x + 2 * y * y + 1.);
	}
};

int main() {
	std::cout.precision(17);

	int i;
	itv r;
	#ifdef BENCH
	std::chrono::system_clock::time_point t;
	#endif

	for (i=1; i<=7; i++) {
		std::cout << "[i=" << i << "]\n";
		#ifdef BENCH
		t = std::chrono::system_clock::now();
		#endif
		r = kv::double_newtoncotes(Func(), itv(-1), itv(1), itv(-1), itv(1), i, 0, 0);
		#ifdef BENCH
		std::cout << "time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - t).count() / 1e9 << " sec\n";
		#endif
		std::cout << "width: " << width(r) << "\n";
		std::cout << r << "\n";
	}
}
