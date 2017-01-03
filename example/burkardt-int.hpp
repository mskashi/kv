#ifndef BURKARDT_INT_HPP
#define BURKARDT_INT_HPP

// Quadrature test problem set from
//   https://people.sc.fsu.edu/~burkardt/f_src/test_int/test_int.html
//   https://people.sc.fsu.edu/~burkardt/f_src/test_int/test_int.f90

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>


struct P01 {
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

struct P03 {
	template <class T> T operator() (T x){
		return sqrt(x);
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

struct P04 {
	template <class T> T operator() (T x){
		return 0.92 * cosh(x) - cos(x);
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

struct P05 {
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

struct P07 {
	template <class T> T operator() (T x){
		return 1. / sqrt(x);
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

struct P08 {
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

struct P09 {
	template <class T> T operator() (T x){
		static T pi = kv::constants<T>::pi();

		return 2. / (2. + sin(10 * pi * x));
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

struct P10 {
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

struct P11 {
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

struct P12 {
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

struct P13 {
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

struct P14 {
	template <class T> T operator() (T x){
		static T pi = kv::constants<T>::pi();

		return sqrt((T)50.) * exp(-50. * pi * x * x);
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

struct P15 {
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

struct P16 {
	template <class T> T operator() (T x){
		static T pi = kv::constants<T>::pi();

		return 50. / (pi * (2500. * x * x + 1.));
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

struct P17 {
	template <class T> T operator() (T x){
		T tmp;
		static T pi = kv::constants<T>::pi();

		tmp = sin(50. * pi * x);
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

struct P18 {
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

struct P19 {
	template <class T> T operator() (T x){
		return log(x);
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

struct P20 {
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

struct P22 {
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

struct P39 {
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

struct P41 {
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

struct P44 {
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

struct P45 {
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

struct P47 {
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

struct P51 {
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

struct P57 {
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

struct P63 {
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

struct P73 {
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

struct P74 {
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

struct P75 {
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

struct P77 {
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
