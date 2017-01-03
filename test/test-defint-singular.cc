#include <iostream>
#include <kv/defint-singular.hpp>

typedef kv::interval<double> itv;

/*
 * \int_0^1 \sin(x) / x dx
 * = 0.946083...
 */

struct Sinc {
	template <class T> T operator()(const T& x) {
		return sin(x) / x;
	}
};

struct Sinc_s {
	template <class T> T operator()(const T& x) {
		return div_tn(sin(x), 1);
		// return div_reduce(sin(x), x, 1);
	}
};

#define SINC_TH 0.125


/*
 * \int_0^1 x^2 / (1 - \cos(x)) / x dx
 * = 2.05727...
 */

struct RemSing {
	template <class T> T operator()(const T& x) {
		return x * x / (1 - cos(x));
	}
};

struct RemSing_s {
	template <class T> T operator()(const T& x) {
		return div_reduce(x * x, 1 - cos(x), 2);
	}
};

#define REMSING_TH 0.125


/*
 * \int_0^1 \sqrt(x) cos(x) dx
 * = 0.531203...
 */

struct Func1 {
	template <class T> T operator() (const T& x) {
		return sqrt(x) * cos(x);
	}
};

struct Func1_s {
	template <class T> T operator() (const T& x) {
		return cos(x);
	}
};

#define FUNC1_TH 0.125


/*
 * \int_0^1 \sqrt{\sin(x)} dx
 * = 0.642978...
 */

struct SqrtSin {
	template <class T> T operator()(const T& x) {
		return sqrt(sin(x));
	}
};

struct SqrtSin_s {
	template <class T> T operator()(const T& x) {
		return sin(x);
	}
};

#define SQRTSIN_TH 0.125


/*
 * \int_0^1 \sqrt{1 - \cos(x)} dx
 * = 0.34625...
 */

struct Func2 {
	template <class T> T operator()(const T& x) {
		return sqrt(1 - cos(x));
	}
};

struct Func2_s {
	template <class T> T operator()(const T& x) {
		return 1 - cos(x);
	}
};

#define FUNC2_TH 0.125


/*
 * \int_0^0.5 (10 sin(pi x) + sin(3 pi x))^(3/2) sin(pi x) dx
 *  = 7.13876...
 */

struct Tanaka_f {
	template <class T> T operator()(const T& x) {
		static const T p = kv::constants<T>::pi();
		return 10 * sin(p * x) + sin(3 * p * x);
	}
};

struct Tanaka_g {
	template <class T> T operator()(const T& x) {
		static const T p = kv::constants<T>::pi();
		return sin(p * x);
	}
};

struct Tanaka {
	template <class T> T operator()(const T& x) {
		static const T p = kv::constants<T>::pi();
		return pow(10 * sin(p * x) + sin(3 * p * x), 1.5) * sin(p * x);
	}
};

#define TANAKA_TH 0.125


/*
 * \int_0^1 sqrt(sin(x))cos(x) dx
 *  = 0.51460...
 */

struct Func3_f {
	template <class T> T operator()(const T& x) {
		return sin(x);
	}
};

struct Func3_g {
	template <class T> T operator()(const T& x) {
		return cos(x);
	}
};

struct Func3 {
	template <class T> T operator()(const T& x) {
		return sqrt(sin(x)) * cos(x);
	}
};

#define FUNC3_TH 0.125


/*
 * \int_0^1 \log(x) / (x+1) dx
 * = -0.822467033424...
 */

struct LogF {
	template <class T> T operator()(const T& x) {
		return log(x) / (x + 1);
	}
};

struct LogF_s {
	template <class T> T operator()(const T& x) {
		return 1 / (x + 1);
	}
};

#define LOGF_TH 0.125


/*
 * \int_0^1 \log{\sin(x)} dx
 * = -1.05672...
 */

struct LogSin {
	template <class T> T operator()(const T& x) {
		return log(sin(x));
	}
};

struct LogSin_s {
	template <class T> T operator()(const T& x) {
		return sin(x);
	}
};

#define LOGSIN_TH 0.125

/*
 * \int_0^1 \log{1-\cos(x)} dx
 * = -2.72107...
 */

struct LogCos {
	template <class T> T operator()(const T& x) {
		return log(1-cos(x));
	}
};

struct LogCos_s {
	template <class T> T operator()(const T& x) {
		return 1-cos(x);
	}
};

#define LOGCOS_TH 0.125


/*
 * \int_0^1 \log{\sin(x)}\cos(x) dx
 * = -0.98671...
 */

struct LogSinCos {
	template <class T> T operator()(const T& x) {
		return log(sin(x)) * cos(x);
	}
};

struct LogSinCos_f {
	template <class T> T operator()(const T& x) {
		return sin(x);
	}
};

struct LogSinCos_g {
	template <class T> T operator()(const T& x) {
		return cos(x);
	}
};

#define LOGSINCOS_TH 0.125


/*
 * \int_0^1 \sqrt{1 - x^2} dx
 * = pi/4 = 0.7853981633974483096...
 */

struct Quadrant {
	template <class T> T operator() (const T& x) {
		return sqrt(1 - x * x);
	}
};

struct Quadrant_s {
	template <class T> T operator() (const T& x) {
		return 1 - x * x;
	}
};

#define QUADRANT_TH 0.875


int main() {
	std::cout.precision(17);
	itv r1, r2;

	std::cout << "int_0^1 sin(x) / x dx = 0.946083...\n" ;

	r1 = kv::defint_singular(Sinc_s(), (itv)0., (itv)SINC_TH, 12);
	r2 = kv::defint_autostep(Sinc(), (itv)SINC_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout<< kv::defint_singular_autostep(Sinc(), Sinc_s(), (itv)0., (itv)1., 12) << "\n";


	std::cout << "int_0^1 x^2 / (1 - cos(x)) / x dx = 2.05727...\n";

	r1 = kv::defint_singular(RemSing_s(), (itv)0., (itv)REMSING_TH, 12);
	r2 = kv::defint_autostep(RemSing(), (itv)REMSING_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout<< kv::defint_singular_autostep(RemSing(), RemSing_s(), (itv)0., (itv)1., 12) << "\n";


	std::cout << "int_0^1 sqrt(x) cos(x) dx = 0.531203...\n";

	r1 = kv::defint_power(Func1_s(), (itv)0., (itv)FUNC1_TH, 12, (itv)0.5);
	r2 = kv::defint_autostep(Func1(), (itv)FUNC1_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power_autostep(Func1(), Func1_s(), (itv)0., (itv)1., 12, (itv)0.5) << "\n";


	std::cout << "int_0^1 sqrt{sin(x)} dx = 0.642978...\n";

	r1 = kv::defint_power2(SqrtSin_s(), (itv)0., (itv)SQRTSIN_TH, 12, (itv)0.5);
	r2 = kv::defint_autostep(SqrtSin(), (itv)SQRTSIN_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power2_autostep(SqrtSin(), SqrtSin_s(), (itv)0., (itv)1., 12, (itv)0.5) << "\n";


	std::cout << "int_0^1 sqrt{1 - cos(x)} dx = 0.34625...\n";

	r1 = kv::defint_power2(Func2_s(), (itv)0., (itv)SQRTSIN_TH, 12, (itv)0.5, 2);
	r2 = kv::defint_autostep(Func2(), (itv)SQRTSIN_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power2_autostep(Func2(), Func2_s(), (itv)0., (itv)1., 12, (itv)0.5, 2) << "\n";


	std::cout << "int_0^0.5 (10 sin(pi x) + sin(3 pi x))^(3/2) sin(pi x) dx = 7.13876...\n";

	r1 = kv::defint_power3(Tanaka_f(), Tanaka_g(), (itv)0., (itv)TANAKA_TH, 12, (itv)1.5, 1);
	r2 = kv::defint_autostep(Tanaka(), (itv)TANAKA_TH, (itv)0.5, 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power3_autostep(Tanaka(), Tanaka_f(), Tanaka_g(), (itv)0., (itv)0.5, 12, (itv)1.5, 1) << "\n";


	std::cout << "int_0^1 sqrt(sin(x))cos(x) dx = 0.51460...\n";

	r1 = kv::defint_power3(Func3_f(), Func3_g(), (itv)0., (itv)FUNC3_TH, 12, (itv)0.5);
	r2 = kv::defint_autostep(Func3(), (itv)FUNC3_TH, (itv)1, 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power3_autostep(Func3(), Func3_f(), Func3_g(), (itv)0., (itv)1, 12, (itv)0.5) << "\n";


	std::cout << "int_0^1 log(x) / (x+1) dx = -0.822467033424...\n";

	r1 = kv::defint_log(LogF_s(), (itv)0., (itv)LOGF_TH, 12);
	r2 = kv::defint_autostep(LogF(), (itv)LOGF_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log_autostep(LogF(), LogF_s(), (itv)0., (itv)1., 12) << "\n";


	std::cout << "int_0^1 log(sin(x)) dx = -1.05672...\n";

	r1 = kv::defint_log2(LogSin_s(), (itv)0., (itv)LOGSIN_TH, 12);
	r2 = kv::defint_autostep(LogSin(), (itv)LOGSIN_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log2_autostep(LogSin(), LogSin_s(), (itv)0., (itv)1., 12) << "\n";


	std::cout << "int_0^1 log(1-cos(x)) dx = -2.72107...\n";

	r1 = kv::defint_log2(LogCos_s(), (itv)0., (itv)LOGCOS_TH, 12, 2);
	r2 = kv::defint_autostep(LogCos(), (itv)LOGCOS_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log2_autostep(LogCos(), LogCos_s(), (itv)0., (itv)1., 12, 2) << "\n";


	std::cout << "int_0^1 log(sin(x))cos(x) dx = -0.98671...\n";

	r1 = kv::defint_log3(LogSinCos_f(), LogSinCos_g(), (itv)0., (itv)LOGSINCOS_TH, 12);
	r2 = kv::defint_autostep(LogSinCos(), (itv)LOGSINCOS_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log3_autostep(LogSinCos(), LogSinCos_f(), LogSinCos_g(), (itv)0., (itv)1., 12) << "\n";


	std::cout << "int_0^1 sqrt{1 - x^2} dx = pi/4 = 0.7853981633974483096...\n";

	r1 = kv::defint_power2_r(Quadrant_s(), (itv)QUADRANT_TH, (itv)1., 12, (itv)0.5);
	r2 = kv::defint_autostep(Quadrant(), (itv)0., (itv)QUADRANT_TH, 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power2_autostep_r(Quadrant(), Quadrant_s(), (itv)0., (itv)1., 12, (itv)0.5) << "\n";
}
