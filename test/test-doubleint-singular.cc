#include <kv/doubleintegral.hpp>
#include <kv/doubleint-singular.hpp>

/*
 * \int_0^{0.125} \int_0^{0.125} \sqrt{x y} \cos{xy} dx dy
 *  = 0.000868036...
 */

struct Func1 {
	template <class T> T operator() (const T& x, const T& y) {
		return sqrt(x * y) * (cos(x * y));
	}
};

struct Func1_f {
	template <class T> T operator() (const T& x, const T& y) {
		return x * y;
	}
};

struct Func1_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};

/*
 * \int_0^{0.125} \int_0^{0.125} \sqrt{x \cos{y}} \cos{xy} dx dy
 *  = 0.00367799...
 */

struct Func2 {
	template <class T> T operator() (const T& x, const T& y) {
		return sqrt(x * cos(y)) * (cos(x * y));
	}
};

struct Func2_f {
	template <class T> T operator() (const T& x, const T& y) {
		return x * cos(y);
	}
};

struct Func2_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};

/*
 * \int_0^{0.125} \int_0^{0.125} \sqrt{(\cos{x}) y} \cos{xy} dx dy
 *  = 0.00367799...
 */

struct Func3 {
	template <class T> T operator() (const T& x, const T& y) {
		return sqrt(cos(x) * y) * (cos(x * y));
	}
};

struct Func3_f {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x) * y;
	}
};

struct Func3_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};

/*
 * \int_0^{0.125} \int_0^{0.125} ((1-\cos{x}) \cos{y})^{1/3} \cos{xy} dx dy
 *  = 0.00185822...
 */

struct Func4 {
	template <class T> T operator() (const T& x, const T& y) {
		return pow((1-cos(x)) * cos(y), 1/T(3)) * (cos(x * y));
	}
};

struct Func4_f {
	template <class T> T operator() (const T& x, const T& y) {
		return (1-cos(x)) * cos(y);
	}
};

struct Func4_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};

/*
 * \int_0^{0.125} \int_0^{0.125} \sqrt{\sin{0.125-x}\cos{y}} \cos{xy} dx dy
 *  = 0.00367596...
 */

struct Func5 {
	template <class T> T operator() (const T& x, const T& y) {
		return sqrt(sin(0.125-x) * cos(y)) * (cos(x * y));
	}
};

struct Func5_f {
	template <class T> T operator() (const T& x, const T& y) {
		return sin(0.125-x) * cos(y);
	}
};

struct Func5_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};

/*
 *  \int_0^{0.1} \int_0^x \sqrt{x + y} dy dx = 0.0154187...
 */

struct Func6 {
	template <class T> T operator() (const T& x, const T& y) {
		return sqrt(x + y);
	}
};

struct Func6_f {
	template <class T> T operator() (const T& x, const T& y) {
		return x + y;
	}
};

struct Func6_g {
	template <class T> T operator() (const T& x, const T& y) {
		return T(1);
	}
};

/*
 * \int_{-1}^1 \int_{-1}^1 (x^2+y^2)^(1/4) * \cos{xy} dx dy
 *  = 3.2003...
 */

struct Func7 {
	template <class T> T operator() (const T& x, const T& y) {
		return pow(x * x + y * y, 0.25) * (cos(x * y));
	}
};

struct Func7_f {
	template <class T> T operator() (const T& x, const T& y) {
		return x * x + y * y;
	}
};

struct Func7_g {
	template <class T> T operator() (const T& x, const T& y) {
		return cos(x * y);
	}
};


typedef kv::interval<double> itv;

int main()
{
	std::cout.precision(17);

	// singularity with edge x=0 and edge y=0

 	std::cout << "\\int_0^{0.125} \\int_0^{0.125} \\sqrt{x y} \\cos{xy} dx dy = 0.000868036...\n";

	std::cout << kv::doubleintegral_power3(Func1_f(), Func1_g(), itv(0.), itv(0.125), itv(0.), itv(0.125), 12, itv(0.5), 1, 1) << "\n";

	// singularity with edge x=0

	std::cout << "\\int_0^{0.125} \\int_0^{0.125} \\sqrt{x \\cos{y}} \\cos{xy} dx dy = 0.00367799...\n";

	std::cout << kv::doubleintegral_power3(Func2_f(), Func2_g(), itv(0.), itv(0.125), itv(0.), itv(0.125), 12, itv(0.5), 1, 0) << "\n";

	// singularity with edge y=0

	std::cout << "\\int_0^{0.125} \\int_0^{0.125} \\sqrt{(\\cos{x}) y} \\cos{xy} dx dy = 0.00367799...\n";

	std::cout << kv::doubleintegral_power3(Func3_f(), Func3_g(), itv(0.), itv(0.125), itv(0.), itv(0.125), 12, itv(0.5), 0, 1) << "\n";

	// singularity with edge x=0 (multiplicity of f_x = 2)

	std::cout << "\\int_0^{0.125} \\int_0^{0.125} ((1-\\cos{x}) \\cos{y})^{1/3} \\cos{xy} dx dy = 0.00185822...\n";

	std::cout << kv::doubleintegral_power3(Func4_f(), Func4_g(), itv(0.), itv(0.125), itv(0.), itv(0.125), 12, 1./itv(3.), 2, 0) << "\n";

	// singularity with right edge x=0.125

	std::cout << "\\int_0^{0.125} \\int_0^{0.125} \\sqrt{\\sin{0.125-x}\\cos{y}} \\cos{xy} dx dy = 0.00367596...\n";

	std::cout << kv::doubleintegral_power3_r1(Func5_f(), Func5_g(), itv(0.), itv(0.125), itv(0.), itv(0.125), 12, itv(0.5), 1, 0) << "\n";

	// singularity with point (0,0)
	std::cout << "\\int_0^{0.1} \\int_0^x \\sqrt{x + y} dy dx = 0.0154187...\n";

	std::cout << kv::doubleintegral_singular_point(Func6_f(), Func6_g(), itv(0.), itv(0.), itv("0.1"), itv(0.), itv("0.1"), itv("0.1"), 12, itv(0.5), 1) << "\n";;


	// singularity with center (0,0)

	std::cout << "\\int_{-1}^1 \\int_{-1}^1 (x^2+y^2)^(1/4) * \\cos{xy} dx dy = 3.2003...\n";

	int i, j;
	itv r, s;
	int div = 4;
	s = 2. / div;
	r = 0.;

	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(s), itv(0.), itv(s), itv(s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(s), itv(s), itv(0.), itv(s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(0.), itv(s), itv(-s), itv(s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(-s), itv(s), itv(-s), itv(0.), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(-s), itv(0.), itv(-s), itv(-s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(-s), itv(-s), itv(0.), itv(-s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(0.), itv(-s), itv(s), itv(-s), 15, itv(0.25), 2);
	r +=  kv::doubleintegral_singular_point(Func7_f(), Func7_g(), itv(0.), itv(0.), itv(s), itv(-s), itv(s), itv(0.), 15, itv(0.25), 2);

	for (i=0; i<div; i++) {
		for (j=0; j<div; j++) {
			if ((i == div/2-1 || i == div/2) && (j == div/2-1 || j == div/2)) continue;
			r += kv::doubleintegral(Func7(), -1. + itv(2.) / div * i, -1. + itv(2.) / div * (i+1), -1. + itv(2.) / div * j, -1. + itv(2.) / div * (j+1), 15, 2);
		}
	}

	std::cout << r << "\n";
}
