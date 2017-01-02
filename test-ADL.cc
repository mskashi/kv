// test program for ADL (Argument Dependent name Lookup) ?
#include <cmath>

namespace hoge {
	class fuga {
		public:
		friend void sqrt(const fuga& x) {}
		friend void func(const fuga& x) {
			double y = 0.;
			// can't compile on VC++ unless below line is uncommented.
			// using std::sqrt;
			sqrt(y);
		}
	};
}

int main()
{
}
