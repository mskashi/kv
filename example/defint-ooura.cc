/*
 * examples taken from Prof. Ooura's DE integration package FAQ
 *  http://www.kurims.kyoto-u.ac.jp/~ooura/intdefaq-j.html
 */

#include <iostream>
#include <kv/defint-singular.hpp>

typedef kv::interval<double> itv;

/*
 * \int_0^1 1/sqrt(1-x^2) dx
 * = pi/2 = 1.5708...
 */

struct Func1 {
	template <class T> T operator()(const T& x) {
		return pow(1 - x * x, -0.5);
	}
};

struct Func1_s {
	template <class T> T operator()(const T& x) {
		return 1 - x * x;
	}
};

/*
 * \int_0^1 sqrt(x) / (exp(x) - 1) dx
 * = 1.6997...
 */

struct Func2 {
	template <class T> T operator()(const T& x) {
		return sqrt(x) / (exp(x) - 1);
	}
};

struct Func2_s {
	template <class T> T operator()(const T& x) {
		return div_reduce(x, exp(x) - 1, 1);
	}
};

/*
 * \int_1^2 1/sqrt(log(x)) dx
 * = 2.14504...
 *   (impossible now)
 */

/*
 * \int_0^1 log(x)/(x-1) dx
 * = 1.64493...
 */

struct Func4 {
	template <class T> T operator()(const T& x) {
		return log(x) / (x - 1);
	}
};

struct Func4_s0 {
	template <class T> T operator()(const T& x) {
		return 1 / (x - 1);
	}
};

struct Func4_s1 {
	template <class T> T operator()(const T& x) {
		return div_reduce(log(x), x - 1, 1);
	}
};

/*
 * \int_0^1 sqrt(abs(0.7-x)) dx
 * = 0.499986...
 */

struct Func6p {
	template <class T> T operator()(const T& x) {
		static const T c("0.7");
		return sqrt(c - x);
	}
};

struct Func6p_s {
	template <class T> T operator()(const T& x) {
		static const T c("0.7");
		return c - x;
	}
};

struct Func6n {
	template <class T> T operator()(const T& x) {
		static const T c("0.7");
		return sqrt(x - c);
	}
};

struct Func6n_s {
	template <class T> T operator()(const T& x) {
		static const T c("0.7");
		return x - c;
	}
};

/*
 * \int_0^1 1/(1e-10+(0.7-x)^2) dx
 * = 314154....
 */

struct Func7 {
	template <class T> T operator()(const T& x) {
		static const T c1("1e-10");
		static const T c2("0.7");
		return 1 / (c1 + pow(c2 - x, 2));
	}
};

/*
 * \int_0^10000 sin(x)log(x) dx
 * = ....
 */

struct Func8 {
	template <class T> T operator()(const T& x) {
		return sin(x) * log(x);
	}
};

struct Func8_s {
	template <class T> T operator()(const T& x) {
		return sin(x);
	}
};

int main() {
	std::cout.precision(17);

	std::cout << kv::defint_power2_autostep_r(Func1(), Func1_s(), (itv)0., (itv)1., 12, (itv)(-0.5), 1) << "\n";
	std::cout << kv::defint_power_autostep(Func2(), Func2_s(), (itv)0., (itv)1., 12, (itv)(-0.5)) << "\n";
	std::cout << kv::defint_log_autostep(Func4(), Func4_s0(), (itv)0, (itv)0.5, 12) +  kv::defint_singular_autostep_r(Func4(), Func4_s1(), (itv)0.5, (itv)1., 12) << "\n";
	std::cout << kv::defint_power2_autostep_r(Func6p(), Func6p_s(), (itv)0., itv("0.7"), 12, (itv)0.5, 1) + kv::defint_power2_autostep(Func6n(), Func6n_s(), itv("0.7"), (itv)1., 12, (itv)0.5, 1) << "\n";
	std::cout << kv::defint_autostep(Func7(), (itv)0., (itv)1., 12) << "\n";
	std::cout << kv::defint_log_autostep(Func8(), Func8_s(), (itv)0., (itv)10000., 12) << "\n";
}
