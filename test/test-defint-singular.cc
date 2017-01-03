#include <iostream>
#include <kv/defint-singular.hpp>

typedef kv::interval<double> itv;

/*
 * \int_0^1 \sqrt{1 - (x-1)^2} dx
 * \int_0^1 \sqrt{x (2 - x)} dx
 * = pi/4 = 0.7853981633974483096...
 */

struct Func1 {
	template <class T> T operator() (T x) {
		return sqrt(x * (2 - x));
	}
};

struct Func1_s {
	template <class T> T operator() (T x) {
		return sqrt(2 - x);
	}
};

#define FUNC1_TH 0.125


/*
 * \int_0^{pi/2} \sqrt{\sin(x)} dx
 * = 1.1981402347...
 */

struct SqrtSin {
	template <class T> T operator()(T& x) {
		return sqrt(sin(x));
	}
};

// \int_0^{pi/2} \sqrt{\sin(x)} dx
// = \int_0^{pi/2} \sqrt{x}\sqrt{\sin(x)/x} dx

struct SqrtSin_s {
	template <class T> T operator()(T& x) {
		return sqrt(div_tn(sin(x), 1));
	}
};

#define SQRTSIN_TH 0.125

/*
 * \int_0^1 \log{\sin(x)} dx
 * = -1.05672...
 */

struct LogSin {
	template <class T> T operator()(T& x) {
		return log(sin(x));
	}
};

struct LogSin_s {
	template <class T> T operator()(T& x) {
		return sin(x);
	}
};

#define LOGSIN_TH 0.125

/*
 * \int_0^1 \log{1-\cos^2(x)} dx
 * = -2.72107...
 */

struct LogCos {
	template <class T> T operator()(T& x) {
		return log(1-cos(x));
	}
};

struct LogCos_s {
	template <class T> T operator()(T& x) {
		return 1-cos(x);
	}
};

#define LOGCOS_TH 0.125


/*
 * \int_0^1 \log(x) / (x+1) dx
 * = -0.822467033424...
 */

struct LogF {
	template <class T> T operator()(T& x) {
		return log(x) / (x + 1);
	}
};

struct LogF_s {
	template <class T> T operator()(T& x) {
		return 1 / (x + 1);
	}
};

#define LOGF_TH 0.125


/*
 * \int_0^1 \sin(x) / x dx
 * = 0.946083...
 */

struct Sinc {
	template <class T> T operator()(T& x) {
		return sin(x) / x;
	}
};

struct Sinc_s {
	template <class T> T operator()(T& x) {
		return div_tn(sin(x), 1);
		// return div_reduce(sin(x), x, 1);
	}
};

#define SINC_TH 0.125


int main() {
	std::cout.precision(17);
	itv r1, r2;

	r1 = kv::defint_power(Func1_s(), (itv)0., (itv)FUNC1_TH, 12, (itv)0.5);
	r2 = kv::defint_autostep(Func1(), (itv)FUNC1_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power_autostep(Func1(), Func1_s(), (itv)0., (itv)1., 12, (itv)0.5) << "\n";


	r1 = kv::defint_power(SqrtSin_s(), (itv)0., (itv)SQRTSIN_TH, 12, (itv)0.5);
	r2 = kv::defint_autostep(SqrtSin(), (itv)SQRTSIN_TH, kv::constants<itv>::pi()/2., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_power_autostep(SqrtSin(), SqrtSin_s(), (itv)0., kv::constants<itv>::pi()/2., 12, (itv)0.5) << "\n";


	r1 = kv::defint_log_singular(LogSin_s(), (itv)0., (itv)LOGSIN_TH, 12, 1);
	r2 = kv::defint_autostep(LogSin(), (itv)LOGSIN_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log_singular_autostep(LogSin(), LogSin_s(), (itv)0., (itv)1., 12, 1) << "\n";

	r1 = kv::defint_log_singular(LogCos_s(), (itv)0., (itv)LOGCOS_TH, 12, 2);
	r2 = kv::defint_autostep(LogCos(), (itv)LOGCOS_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	r1 = kv::defint_log(LogF_s(), (itv)0., (itv)LOGF_TH, 12);
	r2 = kv::defint_autostep(LogF(), (itv)LOGF_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout << kv::defint_log_autostep(LogF(), LogF_s(), (itv)0., (itv)1., 12) << "\n";

	r1 = kv::defint_singular(Sinc_s(), (itv)0., (itv)SINC_TH, 12);
	r2 = kv::defint_autostep(Sinc(), (itv)SINC_TH, (itv)1., 12);

	std::cout << r1 << "\n";
	std::cout << r2 << "\n";
	std::cout << r1 + r2 << "\n";

	std::cout<< kv::defint_singular_autostep(Sinc(), Sinc_s(), (itv)0., (itv)1., 12) << "\n";
}
