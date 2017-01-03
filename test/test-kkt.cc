#include <kv/kkt.hpp>
#include <kv/newton.hpp>
#include <kv/kraw-approx.hpp>

namespace ub = boost::numeric::ublas;

// objective function

struct Func_obj {
	template <class T> T operator()(const ub::vector<T>& x) {
		T x1 = x(0);
		T x2 = x(1);
		T x3 = x(2);

		return pow(x1 - 2, 2) + pow(x2 - 3, 2) + pow(x3 - 4, 2);
	}
	int size() {return 3;}
};

// inequalities

struct Func_inequality {
	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		ub::vector<T> y(1);
		T x1 = x(0);
		T x2 = x(1);
		T x3 = x(2);

		y(0) = x1*x1 + x2*x2 + x3*x3 - 1;

		return y;
	}
};

// equalities

struct Func_equality {
	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		ub::vector<T> y(1);
		T x1 = x(0);
		T x2 = x(1);
		T x3 = x(2);

		y(0) = 4 * x1 + x2 + 2 * x3 - 2;

		return y;
	}
};

// KKT equation

struct Func {
	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		ub::vector<T> y(5);
		T x1 = x(0);
		T x2 = x(1);
		T x3 = x(2);
		T b1 = x(3);
		T mu1 = x(4);

		y(0) = 2*(x1-2) + alpha_plus(b1) * 2 * x1 + mu1 * 4;
		y(1) = 2*(x2-3) + alpha_plus(b1) * 2 * x2 + mu1;
		y(2) = 2*(x3-4) + alpha_plus(b1) * 2 * x3 + mu1 * 2;
		y(3) = alpha_minus(b1) + x1*x1 + x2*x2 + x3*x3 - 1;
		y(4) = 4 * x1 + x2 + 2 * x3 - 2;

		return y;
	}
};


typedef kv::interval<double> itv;

int main()
{
	ub::vector<double> x;
	ub::vector<itv> ix;
	bool b;
	int i;

	std::cout.precision(17);

	x.resize(5);

	// verification using KKT equation generated manually

	Func f1;

	kv::newton_random(f1, x);
	std::cout << x << "\n";
        b = kv::krawczyk_approx(f1, x, ix);
	if (b) {
		std::cout << ix << "\n";
	}

	// verification using KKT equation generated automatically

	Func_obj f;
	Func_inequality g;
	Func_equality h;
	kv::KKT_equation<Func_obj,Func_inequality,Func_equality> f2(f, g, h);

	kv::newton_random(f2, x);
	std::cout << x << "\n";
        b = kv::krawczyk_approx(f2, x, ix);
	if (b) {
		std::cout << ix << "\n";
	}
}
