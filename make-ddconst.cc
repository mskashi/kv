#include "interval.hpp"
#include "rdouble.hpp"
#include "dd.hpp"
#include "rdd.hpp"

typedef kv::dd dd;
typedef kv::interval<dd> idd;

void rational_view(double x)
{
	int d;

	frexp(x, &d);
	d = 52 - (d - 1);

	if (d > 53) {
		std::cout << x * ldexp(1., d) << "./" << ldexp(1., 53) << "./" << ldexp(1., d - 53) << ".\n";
	} else {
		std::cout << x * ldexp(1., d) << "./" << ldexp(1., d) << ".\n";
	}
}


int main()
{
	std::cout.precision(33);

	// printf("%.100f\n", pi(1e-1000)) on calc
	idd pi("3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679");
	std::cout << pi << "\n";
	rational_view(pi.upper().a1);
	rational_view(pi.upper().a2);
	rational_view(pi.lower().a1);
	rational_view(pi.lower().a2);

	// printf("%.100f\n", exp(1, 1e-1000)) on calc
	idd e("2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274");
	std::cout << e << "\n";
	rational_view(e.upper().a1);
	rational_view(e.upper().a2);
	rational_view(e.lower().a1);
	rational_view(e.lower().a2);

	// printf("%.100f\n", ln(2, 1e-1000)) on calc
	idd ln2("0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875");
	std::cout << ln2 << "\n";
	rational_view(ln2.upper().a1);
	rational_view(ln2.upper().a2);
	rational_view(ln2.lower().a1);
	rational_view(ln2.lower().a2);

	std::cout << kv::constants<double>::pi() << "\n";
	std::cout << kv::constants<double>::e() << "\n";
	std::cout << kv::constants<double>::ln2() << "\n";

	std::cout << kv::constants<dd>::pi() << "\n";
	std::cout << kv::constants<dd>::e() << "\n";
	std::cout << kv::constants<dd>::ln2() << "\n";
}
