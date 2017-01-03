#include <kv/interval-conv.hpp>

int main()
{
	kv::mpfr<106> m1;
	kv::dd dd1;
	double d1;
	std::cout.precision(34);

	m1 = 1;
	m1 /= 3;

	std::cout << m1 << "\n";

	kv::mpfrtodd(m1, dd1, 0);
	std::cout << dd1 << "\n";
	kv::mpfrtodd(m1, dd1, -1);
	std::cout << dd1 << "\n";
	kv::mpfrtodd(m1, dd1, 1);
	std::cout << dd1 << "\n";

	kv::mpfrtodouble(m1, d1, 0);
	std::cout << d1 << "\n";
	kv::mpfrtodouble(m1, d1, -1);
	std::cout << d1 << "\n";
	kv::mpfrtodouble(m1, d1, 1);
	std::cout << d1 << "\n";

	dd1 = 1;
	dd1 /= 3;

	std::cout << dd1 << "\n";

	kv::ddtompfr(dd1, m1, 0);
	std::cout << m1 << "\n";
	kv::ddtompfr(dd1, m1, -1);
	std::cout << m1 << "\n";
	kv::ddtompfr(dd1, m1, 1);
	std::cout << m1 << "\n";

	kv::ddtodouble(dd1, d1, 0);
	std::cout << d1 << "\n";
	kv::ddtodouble(dd1, d1, -1);
	std::cout << d1 << "\n";
	kv::ddtodouble(dd1, d1, 1);
	std::cout << d1 << "\n";


	kv::interval< kv::mpfr<106> > m2;
	kv::interval<kv::dd> dd2;
	kv::interval<double> d2;

	d2 = 1.;
	d2 /= 3.;
	std::cout << d2 << "\n";
	kv::idoubletoidd(d2, dd2);
	std::cout << dd2 << "\n";
	kv::idoubletoimpfr(d2, m2);
	std::cout << m2 << "\n";

	dd2 = 1.;
	dd2 /= 3.;
	std::cout << dd2 << "\n";
	kv::iddtoidouble(dd2, d2);
	std::cout << d2 << "\n";
	kv::iddtoimpfr(dd2, m2);
	std::cout << m2 << "\n";

	m2 = 1.;
	m2 /= 3.;
	std::cout << m2 << "\n";
	kv::impfrtoidouble(m2, d2);
	std::cout << d2 << "\n";
	kv::impfrtoidd(m2, dd2);
	std::cout << dd2 << "\n";
}
