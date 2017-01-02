#ifndef BURKARDT_INT_HPP
#define BURKARDT_INT_HPP

// Quadrature test problem set from
//   https://people.sc.fsu.edu/~burkardt/f_src/test_int/test_int.html
//   https://people.sc.fsu.edu/~burkardt/f_src/test_int/test_int.f90

#include "interval.hpp"
#include "rdouble.hpp"


class P01 {
	public:
	template <class T> T operator() (T x){
		return exp(x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P03 {
	public:
	template <class T> T operator() (T x){
		return (2./3.) * x * sqrt(x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P04 {
	public:
	template <class T> T operator() (T x){
		return 1. / (x*x*x*x + x*x + 0.9);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = -1.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P08 {
	public:
	template <class T> T operator() (T x){
		return 1. / (1. + x*x*x*x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P10 {
	public:
	template <class T> T operator() (T x){
		return 1. / (1. + x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P11 {
	public:
	template <class T> T operator() (T x){
		return 1. / (1. + exp(x));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P12 {
	public:
	template <class T> T operator() (T x){
		return x / (exp(x) - 1.);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P13 {
	public:
	template <class T> T operator() (T x){
		return sin(x) / x;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 10.;
	}
};

class P14 {
	public:
	template <class T> T operator() (T x){
		return sqrt(50.) * exp(-50. * 3.1415926535897932 * x * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 10.;
	}
};

class P15 {
	public:
	template <class T> T operator() (T x){
		return 25. * exp(-25. * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 10.;
	}
};

class P16 {
	public:
	template <class T> T operator() (T x){
		return 50. / (3.1415926535897932 * (2500. * x * x + 1.));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P17 {
	public:
	template <class T> T operator() (T x){
		T tmp;
		tmp = sin(50. * 3.1415926535897932 * x);
		return tmp * tmp;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P18 {
	public:
	template <class T> T operator() (T x){
		return x / (exp(x) + 1.);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P20 {
	public:
	template <class T> T operator() (T x){
		return 1. / (x * x + 1.005);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = -1.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P22 {
	public:
	template <class T> T operator() (T x){
		return 1. / (x*x*x*x + x*x + 1.);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P39 {
	public:
	template <class T> T operator() (T x){
		return exp(cos(x));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 2. * 3.1415926535897932;
	}
};

class P41 {
	public:
	template <class T> T operator() (T x){
		return exp(-x) * sin(50. * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 2. * 3.1415926535897932;
	}
};

class P44 {
	public:
	template <class T> T operator() (T x){
		return 1. / (1. + x * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = -4.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 4.;
	}
};

class P45 {
	public:
	template <class T> T operator() (T x){
		return exp(x) * cos(x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 3.1415926535897932;
	}
};

class P47 {
	public:
	template <class T> T operator() (T x){
		return (10. * x - 1) * (10. * x - 1.1) * (10. * x - 1.2) * (10. * x - 1.3);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P51 {
	public:
	template <class T> T operator() (T x){
		return cos(8. * sin(x));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 3.1415926535897932;
	}
};

class P57 {
	public:
	template <class T> T operator() (T x){
		return exp(20. * (x - 1)) * sin(4. * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P63 {
	public:
	template <class T> T operator() (T x){
		T alpha;
		alpha = sin(3.1415926535897932 / 12.);
		return (1. - alpha * alpha) / (1. - alpha * cos(x));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 2. * 3.1415926535897932;
	}
};

class P73 {
	public:
	template <class T> T operator() (T x){
		return 2. / (2. + sin(10. * 3.1415926535897932 * x));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P74 {
	public:
	template <class T> T operator() (T x){
		const double c = 3.;
		const double x0 = 0.75;
		return exp( - c * c * (x - x0) * (x - x0));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P75 {
	public:
	template <class T> T operator() (T x){
		return 1. / (x*x*x*x*x*x + 0.9);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = -1.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 1.;
	}
};

class P77 {
	public:
	template <class T> T operator() (T x){
		return cos(x) + 5. * cos(1.6 * x) - 2. * cos(2. * x) + 5. * cos(4.5 * x) + 7. * cos(9. * x);
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 2.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 7.;
	}
};

#endif // BURKARDT_INT_HPP
