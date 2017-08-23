#include <iostream>
#include <cmath>

int main()
{
	volatile double u, r;

	u = pow(2.0, -52);
	r = fma(1.0 + u, 1.0 - u, -1.0);

	if (r == - pow(2.0, -104)) {
		std::cout << "fma works!\n";
	} else {
		std::cout << "fma is broken!\n";
	}
}
