#include <iostream>
#include <kv/doubleint-curvededge.hpp>

struct Func {
	template <class T> T operator() (const T& x, const T& y) {
		return 1 / (x * x + 2 * y * y + 1);
	}
};

struct Edge_l {
	template <class T> T operator() (const T& x) {
		return -1 + sin(3*x) / 16;
	}
};

struct Edge_h {
	template <class T> T operator() (const T& x) {
		return 1 + sin(5*x) / 16;
	}
};

typedef kv::interval<double> itv;

int main() {
	std::cout.precision(17);

	std::cout << kv::doubleint_curvededge(Func(), Edge_l(), Edge_h(), itv(-1), itv(1), 6, 0, 0) << "\n";
	std::cout << kv::doubleint_curvededge_s(Func(), Edge_l(), Edge_h(), itv(-1), itv(1), 6, 0, 0) << "\n";
}
