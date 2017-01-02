#ifndef BURKARDT_ODE_HPP
#define BURKARDT_ODE_HPP

// ODE test problem set from
//   https://people.sc.fsu.edu/~burkardt/f_src/test_ode/test_ode.html
//   https://people.sc.fsu.edu/~burkardt/f_src/test_ode/test_ode.f90

#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace ub = boost::numeric::ublas;

class P01 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = - x(0);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P02 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = - x(0)* x(0)*x(0) / 2.;

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P03 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = - cos(t) * x(0);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P04 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = x(0) * (20. - x(0)) / 80.;

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P05 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = (x(0) - t) / (x(0) + t);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 4.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P06 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = 2. * x(0) * (1. - x(1));
		y(1) = - x(1) * (1. - x(0));

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 1.;
		x(1) = 3.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P07 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = -x(0) + x(1);
		y(1) = x(0) - 2. * x(1) + x(2);
		y(2) = x(1) - x(2);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(3);
		x(0) = 2.;
		x(1) = 0.;
		x(2) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P08 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = -x(0);
		y(1) = x(0) - x(1) * x(1);
		y(2) = x(1) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(3);
		x(0) = 1.;
		x(1) = 0.;
		x(2) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P09 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		T norm;

		norm = sqrt(x(0)*x(0) + x(1)*x(1));

		y(0) = - x(1) - x(0) * x(2) / norm;
		y(1) = x(0) - x(1) * x(2) / norm;
		y(2) = x(0) / norm;

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(3);
		x(0) = 3.;
		x(1) = 0.;
		x(2) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P10 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);

		y(0) = x(1) * x(2);
		y(1) = - x(0) * x(2);
		y(2) = -0.51 * x(0) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(3);
		x(0) = 0.;
		x(1) = 1.;
		x(2) = 1.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P11 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(10);
		int i;

		y(0) = - x(0);
		for (i=1; i<9; i++) {
			y(i) = x(i-1) - x(i);
		}
		y(9) = - x(8);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(10);
		x(0) = 1.;
		for (i=1; i<10; i++) {
			x(i) = 0.;
		}
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P12 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(10);
		int i;

		y(0) = - x(0);
		for (i=1; i<9; i++) {
			y(i) = (double)i * x(i-1) - (double)(i+1) * x(i);
		}
		y(9) = - 9. * x(8);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(10);
		x(0) = 1.;
		for (i=1; i<10; i++) {
			x(i) = 0.;
		}
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P13 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(10);
		int i;

		y(0) = - 2. * x(0) + x(1);
		for (i=1; i<9; i++) {
			y(i) = x(i-1) - 2. * x(i) + x(i+1);
		}
		y(9) = x(8) - 2. * x(9);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(10);
		x(0) = 1.;
		for (i=1; i<10; i++) {
			x(i) = 0.;
		}
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P14 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(51);
		int i;

		y(0) = - 2. * x(0) + x(1);
		for (i=1; i<50; i++) {
			y(i) = x(i-1) - 2. * x(i) + x(i+1);
		}
		y(50) = x(49) - 2. * x(50);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(51);
		x(0) = 1.;
		for (i=1; i<51; i++) {
			x(i) = 0.;
		}
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P15 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(30);
		const double k2 = 2.95912208286;
		const double m[6] = {
			1.00000597682,
			0.954786104043e-03,
			0.285583733151e-03,
			0.437273164546e-04,
			0.517759138449e-04,
			0.277777777778e-05
		};
		int i, i3, j, l, ll, mm;
		T p;
		T q[5][5];
		T r[5];

		i = 0;
		for (l=3; l<=15; l+=3) {
			i++;
			p = x(l-3)*x(l-3) + x(l-2)*x(l-2) + x(l-1)*x(l-1);
			r[i-1] = 1. / (p * sqrt(p));
			j = 0;
			for (ll=3; ll<=15; ll+=3) {
				j++;
				if (ll != l) {
					p = (x(l-3) - x(ll-3)) * (x(l-3) - x(ll-3))
					+ (x(l-2) - x(ll-2)) * (x(l-2) - x(ll-2))
					+ (x(l-1) - x(ll-1)) * (x(l-1) - x(ll-1));
					q[i-1][j-1] = 1. / (p * sqrt(p));
					q[j-1][i-1] = q[i-1][j-1];
				}
			}
		}

		i3 = 0;
		for (i=1; i<=5; i++) {
			i3 += 3;
			for (ll = i3-2; ll <= i3; ll++) {
				mm = ll - i3;
				y(ll-1) = x(ll+14);
				p = 0.;
				for (j=1; j<=5; j++) {
					mm += 3;
					if (j != i) {
						p += m[j-1] * (x(mm-1) * (q[i-1][j-1] - r[j-1]) - x(ll-1) * q[i-1][j-1]);
					}
				}
				y(ll+14) = k2 * ( - (m[0] + m[i]) * x(ll-1) * r[i-1] + p);
			}
		}

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		x.resize(30);

		x(0) = 3.42947415189;
		x(1) = 3.35386959711;
		x(2) = 1.35494901715;
		x(3) = 6.64145542550;
		x(4) = 5.97156957878;
		x(5) = 2.18231499728;
		x(6) = 11.2630437207;
		x(7) = 14.6952576794;
		x(8) = 6.27960525067;
		x(9) = -30.1552268759;
		x(10) = 1.65699966404;
		x(11) = 1.43785752721;
		x(12) = -21.1238353380;
		x(13) = 28.4465098142;
		x(14) = 15.3882659679;
		x(15) = -0.557160570446;
		x(16) = 0.505696783289;
		x(17) = 0.230578543901;
		x(18) = -0.415570776342;
		x(19) = 0.365682722812;
		x(20) = 0.169143213293;
		x(21) = -0.325325669158;
		x(22) = 0.189706021964;
		x(23) = 0.0877265322780;
		x(24) = -0.0240476254170;
		x(25) = -0.287659532608;
		x(26) = -0.117219543175;
		x(27) = -0.176860753121;
		x(28) = -0.216393453025;
		x(29) = -0.0148647893090;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};


class P16 {
	public:
	double delta;

	P16(double _delta = 0.1) : delta(_delta) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(4);
		T tmp;

		tmp = sqrt(x(0) * x(0) + x(1) * x(1));
		y(0) = x(2);
		y(1) = x(3);
		y(2) = - x(0) / (tmp * tmp * tmp);
		y(3) = - x(1) / (tmp * tmp * tmp);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {
		x.resize(4);
		x(0) = 1. - delta;
		x(1) = 0.;
		x(2) = 0.;
		x(3) = sqrt((1. + delta)/(1. - delta));
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// P17-P20はP16の初期値違い。
// P16, P17は解けるが、P18-P20は厳しい。

class P17 : public P16 {
	public:
	P17(double _delta = 0.3) {delta = _delta;}
};

class P18 : public P16 {
	public:
	P18(double _delta = 0.5) {delta = _delta;}
};

class P19 : public P16 {
	public:
	P19(double _delta = 0.7) {delta = _delta;}
};

class P20 : public P16 {
	public:
	P20(double _delta = 0.9) {delta = _delta;}
};

class P21 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = - ( 1. - 0.25 / ((t + 1.) * (t + 1.))) * x(0) - 1. / (t + 1.) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 0.6713967071418030;
		x(1) = 0.09540051444747446;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P22 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = (1. - x(0) * x(0)) * x(1) - x(0);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 2.;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P23 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = x(0) * x(0) * x(0) / 6. - x(0) + 2. * sin(2.78535 * t);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 0.;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P24 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = 0.032 - 0.4 * x(1) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 30.;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

class P25 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = sqrt(1. + x(1) * x(1)) / (25. - t);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 0.;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// Lotka-Volterra Predator-Prey Equations
class P31 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);
		const double a = 5.;
		const double b = 1.;
		const double c = 0.5;
		const double d = 2.;

		y(0) = (a - b * x(1)) * x(0);
		y(1) = (c * x(0) - d) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 2.;
		x(1) = 2.;
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

// The Lorenz System
class P32 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(3);
		const double beta = 8./3.;
		const double rho = 28.;
		const double sigma = 10.;

		y(0) = sigma * (x(1) - x(0));
		y(1) = rho * x(0) - x(1) - x(0) * x(2);
		y(2) = - beta * x(2) + x(0) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(3);
		x(0) = 2.;
		x(1) = 2.;
		x(2) = 21.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// The Van der Pol equation
class P33 {
	public:

	double delta;

	P33(double _delta = 1.) : delta(_delta) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = delta * (1. - x(0) * x(0)) * x(1) - x(0);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 2.;
		x(1) = 2.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// The Linearized Damped Pendulum
class P34 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);
		const double g = 32.;
		const double k = 1.;
		const double l = 1.;
		const double m = 1.;

		y(0) = x(1);
		y(1) = - (g / l) * x(0) - (k / m) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 2.;
		x(1) = 2.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// The Nonlinear Damped Pendulum
class P35 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);
		const double g = 32.;
		const double k = 1.;
		const double l = 1.;
		const double m = 1.;

		y(0) = x(1);
		y(1) = - (g / l) * sin(x(0)) - (k / m) * x(1);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 2.;
		x(1) = 2.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 20.;
	}
};

// Duffing's Equation
class P36 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);

		y(0) = x(1);
		y(1) = x(0) * (1. - x(0) * x(0));

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 0.5;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 100.;
	}
};

// Duffing's Equation with Damping and Forcing
class P37 {
	public:
	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(2);
		const double a = 0.3;
		const double k = 0.2;
		const double w = 1.;

		y(0) = x(1);
		y(1) = x(0) * (1. - x(0) * x(0)) - k * x(1) + a * cos(w * t);

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(2);
		x(0) = 0.5;
		x(1) = 0.;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 100.;
	}
};

// Shampine's Ball of Flame
class P38 {
	public:
	double delta;

	P38(double _delta = 0.01) : delta(_delta) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = x(0) * x(0) * (1. - x(0));

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = delta;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = 2./delta;
	}
};

// Polking's first order ODE
class P39 {
	public:

	double a;
	double b;

	P39(double _a = 1., double _b = 0.) : a(_a), b(_b) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = x(0) * x(0) - a * t + b;

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = 0.5;
	}

	template<class T>
	void start_time(kv::interval<T>& x) {

		x = 0.;
	}

	template<class T>
	void stop_time(kv::interval<T>& x) {

		x = (b + 9.) / a;
	}
};

// the Kee problem
class P40 {
	public:

	double eps;

	P40(double _eps = 0.01) : eps(_eps) {}

	template <class T> ub::vector<T> operator() (ub::vector<T> x, T t){
		ub::vector<T> y(1);

		y(0) = x(0) * (x(0) - t) / eps;

		return y;
	}

	template<class T>
	void initial_value(ub::vector< kv::interval<T> >& x) {

		x.resize(1);
		x(0) = -1.;
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

#endif // BURKARDT_ODE_HPP
