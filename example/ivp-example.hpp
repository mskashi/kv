#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace ub = boost::numeric::ublas;

struct SimpleHarmonic {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = - x(0);

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(2);
		x(0) = 0.;
		x(1) = 1.;
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = 100.;
	}
};

struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(3);
		x(0) = 15.;
		x(1) = 15.;
		x(2) = 36.;
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = 20.;
	}
};

struct VdP {
	double mu;

	VdP(double mu = 0.25) : mu(mu) {}

	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = mu * (1. - x(0)*x(0))*x(1) - x(0);

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(2);
		x(0) = 1.;
		x(1) = 1.;
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = 20.;
	}
};

struct Nobi {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = x(0) - x(0)*x(0)*x(0);

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(2);
		x(0) = 0.;
		x(1) = 4.;
		// x(0) = kv::interval<T>(-0.05, 0.05);
		// x(1) = kv::interval<T>(3.95, 4.05);
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = "3.3";
	}
};

struct QuadTest1 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(1);

		y(0) = t * t;

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(1);
		x(0) = 0.;
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = 20.;
	}
};

struct QuadTest2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(1);

		y(0) = 0.;

		return y;
	}

	template <class T>
	void initial_value(ub::vector< kv::interval<T> >&  x) {
		x.resize(1);
		x(0) = 0.;
	}

	template <class T>
	void start_time(kv::interval<T>& x) {
		x = 0.;
	}

	template <class T>
	void stop_time(kv::interval<T>& x) {
		x = 20.;
	}
};
