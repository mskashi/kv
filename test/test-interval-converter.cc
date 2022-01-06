#include <kv/interval.hpp>
#include <kv/interval-converter.hpp>

int main()
{
	kv::mpfr<140> m1;
	kv::mpfr<53> n1;
	kv::dd dd1;
	double d1;
	std::cout.precision(42);

	std::cout << "[testing rounded_converter]\n";

	std::cout << "dd to double\n";
	dd1 = 1; dd1 /= 3;
	std::cout << dd1 << "\n";
	kv::rounded_converter(dd1, d1, 0);
	std::cout << d1 << "\n";
	kv::rounded_converter(dd1, d1, -1);
	std::cout << d1 << "\n";
	kv::rounded_converter(dd1, d1, 1);
	std::cout << d1 << "\n";

	std::cout << "double to dd\n";
	d1 = 1; d1 /= 3;
	std::cout << d1 << "\n";
	kv::rounded_converter(d1, dd1, 0);
	std::cout << dd1 << "\n";
	kv::rounded_converter(d1, dd1, -1);
	std::cout << dd1 << "\n";
	kv::rounded_converter(d1, dd1, 1);
	std::cout << dd1 << "\n";

	std::cout << "mpfr to double\n";
	m1 = 1; m1 /= 3;
	std::cout << m1 << "\n";
	kv::rounded_converter(m1, d1, 0);
	std::cout << d1 << "\n";
	kv::rounded_converter(m1, d1, -1);
	std::cout << d1 << "\n";
	kv::rounded_converter(m1, d1, 1);
	std::cout << d1 << "\n";

	std::cout << "mpfr to dd\n";
	m1 = 1; m1 /= 3;
	std::cout << m1 << "\n";
	kv::rounded_converter(m1, dd1, 0);
	std::cout << dd1 << "\n";
	kv::rounded_converter(m1, dd1, -1);
	std::cout << dd1 << "\n";
	kv::rounded_converter(m1, dd1, 1);
	std::cout << dd1 << "\n";

	std::cout << "double to mpfr\n";
	d1 = 1; d1 /= 3;
	std::cout << d1 << "\n";
	kv::rounded_converter(d1, m1, 0);
	std::cout << m1 << "\n";
	kv::rounded_converter(d1, m1, -1);
	std::cout << m1 << "\n";
	kv::rounded_converter(d1, m1, 1);
	std::cout << m1 << "\n";

	std::cout << "dd to mpfr\n";
	dd1 = 1; dd1 /= 3;
	std::cout << dd1 << "\n";
	kv::rounded_converter(dd1, m1, 0);
	std::cout << m1 << "\n";
	kv::rounded_converter(dd1, m1, -1);
	std::cout << m1 << "\n";
	kv::rounded_converter(dd1, m1, 1);
	std::cout << m1 << "\n";

	std::cout << "mpfr to mpfr\n";
	m1 = 1; m1 /= 3;
	std::cout << m1 << "\n";
	kv::rounded_converter(m1, n1, 0);
	std::cout << n1 << "\n";
	kv::rounded_converter(m1, n1, -1);
	std::cout << n1 << "\n";
	kv::rounded_converter(m1, n1, 1);
	std::cout << n1 << "\n";

#ifdef __HAVE_FLOAT64X

	_Float64x dx1;
	kv::ddx ddx1;

	std::cout << "float64x to double\n";
	dx1 = 1; dx1 /= 3;
	std::cout << dx1 << "\n";
	kv::rounded_converter(dx1, d1, 0);
	std::cout << d1 << "\n";
	kv::rounded_converter(dx1, d1, -1);
	std::cout << d1 << "\n";
	kv::rounded_converter(dx1, d1, 1);
	std::cout << d1 << "\n";

	std::cout << "float64x to dd\n";
	dx1 = 1; dx1 /= 3;
	std::cout << dx1 << "\n";
	kv::rounded_converter(dx1, dd1, 0);
	std::cout << dd1 << "\n";
	kv::rounded_converter(dx1, dd1, -1);
	std::cout << dd1 << "\n";
	kv::rounded_converter(dx1, dd1, 1);
	std::cout << dd1 << "\n";

	std::cout << "float64x to mpfr\n";
	dx1 = 1; dx1 /= 3;
	std::cout << dx1 << "\n";
	kv::rounded_converter(dx1, m1, 0);
	std::cout << m1 << "\n";
	kv::rounded_converter(dx1, m1, -1);
	std::cout << m1 << "\n";
	kv::rounded_converter(dx1, m1, 1);
	std::cout << m1 << "\n";

	std::cout << "double to float64x\n";
	d1 = 1; d1 /= 3;
	std::cout << dd1 << "\n";
	kv::rounded_converter(d1, dx1, 0);
	std::cout << dx1 << "\n";
	kv::rounded_converter(d1, dx1, -1);
	std::cout << dx1 << "\n";
	kv::rounded_converter(d1, dx1, 1);
	std::cout << dx1 << "\n";

	std::cout << "dd to float64x\n";
	dd1 = 1; dd1 /= 3;
	std::cout << dd1 << "\n";
	kv::rounded_converter(dd1, dx1, 0);
	std::cout << dx1 << "\n";
	kv::rounded_converter(dd1, dx1, -1);
	std::cout << dx1 << "\n";
	kv::rounded_converter(dd1, dx1, 1);
	std::cout << dx1 << "\n";

	std::cout << "mpfr to float64x\n";
	m1 = 1; m1 /= 3;
	std::cout << m1 << "\n";
	kv::rounded_converter(m1, dx1, 0);
	std::cout << dx1 << "\n";
	kv::rounded_converter(m1, dx1, -1);
	std::cout << dx1 << "\n";
	kv::rounded_converter(m1, dx1, 1);
	std::cout << dx1 << "\n";

	std::cout << "ddx to double\n";
	ddx1 = 1; ddx1 /= 3;
	std::cout << ddx1 << "\n";
	kv::rounded_converter(ddx1, d1, 0);
	std::cout << d1 << "\n";
	kv::rounded_converter(ddx1, d1, -1);
	std::cout << d1 << "\n";
	kv::rounded_converter(ddx1, d1, 1);
	std::cout << d1 << "\n";

	std::cout << "ddx to dd\n";
	ddx1 = 1; ddx1 /= 3;
	std::cout << ddx1 << "\n";
	kv::rounded_converter(ddx1, dd1, 0);
	std::cout << dd1 << "\n";
	kv::rounded_converter(ddx1, dd1, -1);
	std::cout << dd1 << "\n";
	kv::rounded_converter(ddx1, dd1, 1);
	std::cout << dd1 << "\n";

	std::cout << "ddx to mpfr\n";
	ddx1 = 1; ddx1 /= 3;
	std::cout << ddx1 << "\n";
	kv::rounded_converter(ddx1, m1, 0);
	std::cout << m1 << "\n";
	kv::rounded_converter(ddx1, m1, -1);
	std::cout << m1 << "\n";
	kv::rounded_converter(ddx1, m1, 1);
	std::cout << m1 << "\n";

	std::cout << "ddx to float64x\n";
	ddx1 = 1; ddx1 /= 3;
	std::cout << ddx1 << "\n";
	kv::rounded_converter(ddx1, dx1, 0);
	std::cout << dx1 << "\n";
	kv::rounded_converter(ddx1, dx1, -1);
	std::cout << dx1 << "\n";
	kv::rounded_converter(ddx1, dx1, 1);
	std::cout << dx1 << "\n";

	std::cout << "double to ddx\n";
	d1 = 1; d1 /= 3;
	std::cout << d1 << "\n";
	kv::rounded_converter(d1, ddx1, 0);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(d1, ddx1, -1);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(d1, ddx1, 1);
	std::cout << ddx1 << "\n";

	std::cout << "dd to ddx\n";
	dd1 = 1; dd1 /= 3;
	std::cout << dd1 << "\n";
	kv::rounded_converter(dd1, ddx1, 0);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(dd1, ddx1, -1);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(dd1, ddx1, 1);
	std::cout << ddx1 << "\n";

	std::cout << "mpfr to ddx\n";
	m1 = 1; m1 /= 3;
	std::cout << m1 << "\n";
	kv::rounded_converter(m1, ddx1, 0);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(m1, ddx1, -1);
	std::cout << ddx1 << "\n";
	kv::rounded_converter(m1, ddx1, 1);
	std::cout << ddx1 << "\n";

#endif // __HAVE_FLOAT64X

	kv::interval< kv::mpfr<140> > m2;
	kv::interval< kv::mpfr<53> > n2;
	kv::interval<kv::dd> dd2;
	kv::interval<double> d2;

	std::cout << "[testing interval_converter]\n";

	std::cout << "interval<dd> to interval<double>\n";
	dd2 = 1.; dd2 /= 3.;
	std::cout << dd2 << "\n";
	d2 = dd2;
	std::cout << d2 << "\n";

	std::cout << "interval<double> to interval<dd>\n";
	d2 = 1.; d2 /= 3.;
	std::cout << d2 << "\n";
	dd2 = d2;
	std::cout << dd2 << "\n";

	std::cout << "interval<mpfr> to interval<double>\n";
	m2 = 1.; m2 /= 3.;
	std::cout << m2 << "\n";
	d2 = m2;
	std::cout << d2 << "\n";

	std::cout << "interval<mpfr> to interval<dd>\n";
	m2 = 1.; m2 /= 3.;
	std::cout << m2 << "\n";
	dd2 = m2;
	std::cout << dd2 << "\n";

	std::cout << "interval<double> to interval<mpfr>\n";
	d2 = 1.; d2 /= 3.;
	std::cout << d2 << "\n";
	m2 = d2;
	std::cout << m2 << "\n";

	std::cout << "interval<dd> to interval<mpfr>\n";
	dd2 = 1.; dd2 /= 3.;
	std::cout << dd2 << "\n";
	m2 = dd2;
	std::cout << m2 << "\n";

	std::cout << "interval<mpfr> to interval<mpfr>\n";
	m2 = 1.; m2 /= 3.;
	std::cout << m2 << "\n";
	m2 = m2;
	std::cout << n2 << "\n";

#ifdef __HAVE_FLOAT64X

	kv::interval<_Float64x> dx2;
	kv::interval<kv::ddx> ddx2;

	std::cout << "interval<float64x> to interval<double>\n";
	dx2 = 1; dx2 /= 3;
	std::cout << dx2 << "\n";
	d2 = dx2;
	std::cout << d2 << "\n";

	std::cout << "interval<float64x> to interval<dd>\n";
	dx2 = 1; dx2 /= 3;
	std::cout << dx2 << "\n";
	dd2 = dx2;
	std::cout << dd2 << "\n";

	std::cout << "interval<float64x> to interval<mpfr>\n";
	dx2 = 1; dx2 /= 3;
	std::cout << dx2 << "\n";
	m2 = dx2;
	std::cout << m2 << "\n";

	std::cout << "interval<double> to interval<float64x>\n";
	d2 = 1.; d2 /= 3.;
	std::cout << d2 << "\n";
	dx2 = d2;
	std::cout << dx2 << "\n";

	std::cout << "interval<dd> to interval<float64x>\n";
	dd2 = 1.; dd2 /= 3.;
	std::cout << dd2 << "\n";
	dx2 = dd2;
	std::cout << dx2 << "\n";

	std::cout << "interval<mpfr> to interval<float64x>\n";
	m2 = 1.; m2 /= 3.;
	std::cout << m2 << "\n";
	dx2 = m2;
	std::cout << dx2 << "\n";

	std::cout << "interval<ddx> to interval<double>\n";
	ddx2 = 1; ddx2 /= 3;
	std::cout << ddx2 << "\n";
	d2 = ddx2;
	std::cout << d2 << "\n";

	std::cout << "interval<ddx> to interval<dd>\n";
	ddx2 = 1; ddx2 /= 3;
	std::cout << ddx2 << "\n";
	dd2 = ddx2;
	std::cout << dd2 << "\n";

	std::cout << "interval<ddx> to interval<mpfr>\n";
	ddx2 = 1; ddx2 /= 3;
	std::cout << ddx2 << "\n";
	m2 = ddx2;
	std::cout << m2 << "\n";

	std::cout << "interval<ddx> to interval<float64x>\n";
	ddx2 = 1; ddx2 /= 3;
	std::cout << ddx2 << "\n";
	dx2 = ddx2;
	std::cout << dx2 << "\n";

	std::cout << "interval<double> to interval<ddx>\n";
	d2 = 1.; d2 /= 3.;
	std::cout << d2 << "\n";
	ddx2 = d2;
	std::cout << ddx2 << "\n";

	std::cout << "interval<dd> to interval<ddx>\n";
	dd2 = 1.; dd2 /= 3.;
	std::cout << dd2 << "\n";
	ddx2 = dd2;
	std::cout << ddx2 << "\n";

	std::cout << "interval<mpfr> to interval<ddx>\n";
	m2 = 1.; m2 /= 3.;
	std::cout << m2 << "\n";
	ddx2 = m2;
	std::cout << ddx2 << "\n";

	std::cout << "interval<float64x> to interval<ddx>\n";
	dx2 = 1.; dx2 /= 3.;
	std::cout << dx2 << "\n";
	ddx2 = dx2;
	std::cout << ddx2 << "\n";

#endif // __HAVE_FLOAT64X

}
