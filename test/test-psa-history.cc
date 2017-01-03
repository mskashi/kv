#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/psa.hpp>

typedef kv::interval<double> itv;

int main()
{
	kv::psa<itv> x;
	int i;

	// approximate solution of dx/dt = -x^2 - sin(x), x(0) = 1

	x = 1.;
	for (i=0; i<5; i++) {
		x = 1. + integrate(-x*x - sin(x));
		std::cout << x << "\n";
	}

	/* calculate same approximate solution using history.
	 *  will be faster and yield same result
	 *         record use
	 *  (i=0)  o      x
	 *  (i=1)  o      o
	 *  (i=2)  o      o
	 *  (i=3)  o      o
	 *  (i=4)  x      o
	 */

	x = 1.;
	kv::psa<itv>::record_history() = true;
	for (i=0; i<5; i++) {
		if (i == 1) {
			kv::psa<itv>::use_history() = true;
		}
		if (i == 4) {
			kv::psa<itv>::record_history() = false;
		}
		x = 1. + integrate(-x*x - sin(x));
		std::cout << x << "\n";
	}
	kv::psa<itv>::use_history() = false;

	// size of history buffer: should be 0
	std::cout << kv::psa<itv>::history().size() << "\n";;
}
