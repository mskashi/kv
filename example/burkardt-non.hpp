#ifndef BURKARDT_NON_HPP
#define BURKARDT_NON_HPP

/*
 * nonlinear equation example set from
 *  http://people.sc.fsu.edu/~jburkardt/f_src/test_nonlin/test_nonlin.html
 *  http://people.sc.fsu.edu/~jburkardt/f_src/test_nonlin/test_nonlin.f90
 */

#include <boost/numeric/ublas/vector.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

namespace ub = boost::numeric::ublas;

struct GenRosen {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int i;

		y(0) = 1. - x(0);
		for (i=1; i<s; i++) {
			y(i) = 10. * (x(i) - x(i-1) * x(i-1));
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Powell {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(4);

		y(0) = x(0) + 10. * x(1);
		y(1) = sqrt(5.) * (x(2) - x(3));
		y(2) = (x(1) - 2. * x(2)) * (x(1) - 2. * x(2));
		y(3) = sqrt(10.) * (x(0) - x(3)) * (x(0) - x(3));

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(4);
		for (i=0; i<4; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Wood {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(4);
		T tmp1, tmp2;

		tmp1 = x(1) - x(0) * x(0);
		tmp2 = x(3) - x(2) * x(2);

		y(0) = -200. * x(0) * tmp1 - (1. - x(0));
		y(1) = 200. * tmp1 + 20.2 * (x(1) - 1.) + 19.8 * (x(3) - 1.);
		y(2) = -180. * x(2) * tmp2 - (1. - x(2));
		y(3) = 180. * tmp2 + 20.2 * (x(3) - 1.) + 19.8 * (x(1) - 1.);

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(4);
		for (i=0; i<4; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Watson {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int i, j, k;
		T sum1, sum2, tmp, ti;

		for (i=0; i<s; i++) {
			y(i) = 0;
		}

		for (i=1; i<=29; i++) {
			ti = i / 29.;
			sum1 = 0.;
			tmp = 1.;
			for (j=1; j<s; j++) {
				sum1 += (T)j * tmp * x(j);
				tmp *= ti;
			}
			sum2 = 0;
			tmp = 1.;
			for (j=0; j<s; j++) {
				sum2 += tmp * x(j);
				tmp *= ti;
			}
			tmp = (sum1 - sum2 * sum2 - 1.) / ti;
			for (k=0; k<s; k++) {
				y(k) += tmp * ((T)k - 2. * ti * sum2);
				tmp *= ti;
			}
		}

		y(0) += 3. * x(0) - 2. * x(0) + x(1) + 2. * x(0)*x(0)*x(0);
		y(1) += x(1) - x(0) * x(0) - 1.;

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Chebyquad {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int i, j;
		T t1, t2, t3;

		for (i=0; i<s; i++) {
			y(i) = 0.;
		}

		for (j=0; j<s; j++) {
			t1 = 1.;
			t2 = x(j);
			for (i=0; i<s; i++) {
				y(i) += t2;
				t3 = 2. * x(j) * t2 - t1;
				t1 = t2;
				t2 = t3;
			}
		}

		for (i=0; i<s; i++) {
			y(i) /= (T)s;
		}

		for (i=0; i<s; i++) {
			if ( (i+1) % 2 == 0) {
				y(i) += 1. / (T)((i+1)*(i+1)-1);
			}
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Brown {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int i, j;
		T sum, prod;

		sum = 0.;
		for (i=0; i<s; i++) {
			sum += x(i);
		}

		for (i=0; i<s-1; i++) {
			y(i) = x(i) + sum - (T)(s+1);
		}

		prod = 1.;
		for (i=0; i<s; i++) {
			prod *= x(i);
		}

		y(s-1) = prod - 1.;

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct DBVP {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int k;
		T h, tmp;

		h = 1. / (s + 1);

		for (k=0; k<s; k++) {
			tmp = x(k) + (T)(k+1) * h + 1.;
			y(k) = 2. * x(k) + 0.5 * h * h * tmp*tmp*tmp;
			if (k > 0) {
				y(k) -= x(k-1);
			}
			if (k < s-1) {
				y(k) -= x(k+1);
			}
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct DIntEq {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int j, k;
		T h, tk, tj, sum1, sum2, tmp;

		h = 1. / (s + 1);

		for (k=0; k<s; k++) {
			tk = (k+1.) / (s+1.);
			sum1 = 0.;
			for (j=0; j<k+1; j++) {
				tj = (j+1.) * h;
				tmp = x(j) + tj + 1.;
				sum1 += tj * tmp*tmp*tmp;
			}
			sum2 = 0.;
			for (j=k+1; j<s; j++) {
				tj = (j+1.) * h;
				tmp = x(j) + tj + 1.;
				sum2 += (1. - tj) * tmp*tmp*tmp;
			}
			y(k) = x(k) + h * ( (1. - tk) * sum1 + tk * sum2) / 2.;
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct VDim {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int j;
		T sum1, tmp;

		sum1 = 0.;

		for (j=0; j<s; j++) {
			sum1 += (j+1.) * (x(j) - 1.);
		}

		tmp = sum1 * (1. + 2. * sum1 * sum1);

		for (j=0; j<s; j++) {
			y(j) = x(j) - 1. + (j+1.) * tmp;
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Broyden {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int k;

		for (k=0; k<s; k++) {
			y(k) = (3. - 2. * x(k)) * x(k) + 1.;
			if (k > 0) {
				y(k) -= x(k-1);
			}
			if (k < s-1) {
				y(k) -= 2. * x(k+1);
			}
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct BroydenBand {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);
		int k, k1, k2, j;
		T tmp;

		for (k=0; k<s; k++) {
			k1 = k - 5;
			if (k1 < 0) k1 = 0;
			k2 = k + 1;
			if (k2 > s-1) k2 = s-1;

			tmp = 0.;
			for (j=k1; j<=k2; j++) {
				if (j != k) {
					tmp += x(j) * (1. + x(j));
				}
			}

			y(k) = x(k) + (2. + 5. * x(k)*x(k)) + 1. - tmp;
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Hammarling2x2 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(4);

		y(0) = (x(0) * x(0) + x(1) * x(2)) - 0.0001;
		y(1) = (x(0) * x(1) + x(1) * x(3)) - 1.;
		y(2) = (x(2) * x(0) + x(3) * x(2)) - 0.;
		y(3) = (x(2) * x(1) + x(3) * x(3)) - 0.0001;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(4);
		for (i=0; i<4; i++) {
			x(i).assign(-100., 100.);
		}
	}
};

struct Hammarling3x3 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(9);

		y(0) = (x(0) * x(0) + x(1) * x(3) + x(2) * x(6)) - 0.0001;
		y(1) = (x(0) * x(1) + x(1) * x(4) + x(2) * x(7)) - 1.;
		y(2) = (x(0) * x(2) + x(1) * x(5) + x(2) * x(8));

		y(3) = (x(3) * x(0) + x(4) * x(3) + x(5) * x(6));
		y(4) = (x(3) * x(1) + x(4) * x(4) + x(5) * x(7)) - 0.0001;
		y(5) = (x(3) * x(2) + x(4) * x(5) + x(5) * x(8));

		y(6) = (x(6) * x(0) + x(7) * x(3) + x(8) * x(6));
		y(7) = (x(6) * x(1) + x(7) * x(4) + x(8) * x(7));
		y(8) = (x(6) * x(2) + x(7) * x(5) + x(8) * x(8)) - 0.0001;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(9);
		for (i=0; i<9; i++) {
			x(i).assign(-100., 100.);
		}
	}
};

struct P17 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) + x(1) - 3.;
		y(1) = x(0) * x(0) + x(1) * x(1) - 9.;;

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct P19 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(2);

		y(0) = x(0) * (x(0) * x(0) + x(1) * x(1));
		y(1) = x(1) * (x(0) * x(0) + x(1) * x(1));

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(2);
		for (i=0; i<2; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct P20 {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		ub::vector<T> y(1);

		y(0) = x(0) * (x(0) - 5.) * (x(0) - 5.);

		return y;
	}

	template<class T>
	void range(ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(1);
		for (i=0; i<1; i++) {
			x(i).assign(-10., 10.);
		}
	}
};

struct Chandrasekhar {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x){
		int s = x.size();
		ub::vector<T> y(s);

		ub::vector<T> mu(s);
		T sum, term;
		int i, j;

		for (i=0; i<s; i++) {
			y(i) = x(i);
		}
		for (i=0; i<s; i++) {
			mu(i) = (2.*(i+1.)-1.) / (2. * s);
		}

		for (i=0; i<s; i++) {
			sum = 0.;
			for (j=0; j<s; j++) {
				sum += mu(i) * x(j) / (mu(i) + mu(j));
			}
			term = 1. - 0.9 * sum / (2. * s);
			y(i) -= 1. / term;
		}

		return y;
	}

	template<class T>
	void range(int n, ub::vector< kv::interval<T> >& x) {
		int i;

		x.resize(n);
		for (i=0; i<n; i++) {
			x(i).assign(-3., 3.);
		}
	}
};

#endif // BURKARDT_NON_HPP
