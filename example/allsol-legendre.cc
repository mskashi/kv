#include <kv/allsol.hpp>

namespace ub = boost::numeric::ublas;
typedef kv::interval<double> itv;

struct Legendre {
	int n;

	Legendre(int n) : n(n){
	}

	template <class T> ub::vector<T> operator()(const ub::vector<T>& x) {
		int i;
		T tmp1, tmp2, tmp3;
		ub::vector<T> y(1);

		if (n == 0) {
			y(0) = 1.;
			return y;
		}
		if (n == 1) {
			y(0) = x(0);
			return y;
		}
		tmp1 = 1.;
		tmp2 = x(0);
		for (i=2; i<=n; i++) {
			tmp3 = ((2*i-1) * x(0) * tmp2 - (i-1) * tmp1) / i;
			tmp1 = tmp2;
			tmp2 = tmp3;
		}
		y(0) = tmp3;
		return y;
	}
};

int main()
{
	ub::vector<itv> I(1);
	std::cout.precision(17);

	I(0) = itv(-1, 1);

	allsol(Legendre(20), I, 2);
}
