#include <iostream>
#include <kv/defint.hpp>

#include "burkardt-int.hpp"


typedef kv::interval<double> itv;

int main()
{
	itv start, end;

	std::cout.precision(17);

	P01 p01;

	p01.start_time(start);
	p01.stop_time(end);

	std::cout << kv::defint_autostep(p01, start, end, 12);
	std::cout << "\n";

#if 0
	// not differentiable at x = 0
	P03 p03;

	p03.start_time(start);
	p03.stop_time(end);

	std::cout << kv::defint_autostep(p03, start, end, 12);
	std::cout << "\n";
#endif

	P04 p04;

	p04.start_time(start);
	p04.stop_time(end);

	std::cout << kv::defint_autostep(p04, start, end, 12);
	std::cout << "\n";

	P05 p05;

	p05.start_time(start);
	p05.stop_time(end);

	std::cout << kv::defint_autostep(p05, start, end, 12);
	std::cout << "\n";

#if 0
	// not differentiable at x = 0
	P07 p07;

	p07.start_time(start);
	p07.stop_time(end);

	std::cout << kv::defint_autostep(p07, start, end, 12);
	std::cout << "\n";
#endif

	P08 p08;

	p08.start_time(start);
	p08.stop_time(end);

	std::cout << kv::defint_autostep(p08, start, end, 12);
	std::cout << "\n";

	P09 p09;

	p09.start_time(start);
	p09.stop_time(end);

	std::cout << kv::defint_autostep(p09, start, end, 12);
	std::cout << "\n";

	P10 p10;

	p10.start_time(start);
	p10.stop_time(end);

	std::cout << kv::defint_autostep(p10, start, end, 12);
	std::cout << "\n";

	P11 p11;

	p11.start_time(start);
	p11.stop_time(end);

	std::cout << kv::defint_autostep(p11, start, end, 12);
	std::cout << "\n";
#if 0
	// has a removable singularity at x = 0
	P12 p12;

	p12.start_time(start);
	p12.stop_time(end);

	std::cout << kv::defint_autostep(p12, start, end, 12);
	std::cout << "\n";
#endif

#if 0
	// has a removable singularity at x = 0
	P13 p13;

	p13.start_time(start);
	p13.stop_time(end);

	std::cout << kv::defint_autostep(p13, start, end, 12);
	std::cout << "\n";
#endif

	P14 p14;

	p14.start_time(start);
	p14.stop_time(end);

	std::cout << kv::defint_autostep(p14, start, end, 12);
	std::cout << "\n";

	P15 p15;

	p15.start_time(start);
	p15.stop_time(end);

	std::cout << kv::defint_autostep(p15, start, end, 12);
	std::cout << "\n";

	P16 p16;

	p16.start_time(start);
	p16.stop_time(end);

	std::cout << kv::defint_autostep(p16, start, end, 12);
	std::cout << "\n";

	P17 p17;

	p17.start_time(start);
	p17.stop_time(end);

	std::cout << kv::defint_autostep(p17, start, end, 12);
	std::cout << "\n";

	P18 p18;

	p18.start_time(start);
	p18.stop_time(end);

	std::cout << kv::defint_autostep(p18, start, end, 12);
	std::cout << "\n";

#if 0
	// not differentiable at x = 0
	P19 p19;

	p19.start_time(start);
	p19.stop_time(end);

	std::cout << kv::defint_autostep(p19, start, end, 12);
	std::cout << "\n";
#endif

	P20 p20;

	p20.start_time(start);
	p20.stop_time(end);

	std::cout << kv::defint_autostep(p20, start, end, 12);
	std::cout << "\n";

	P22 p22;

	p22.start_time(start);
	p22.stop_time(end);

	std::cout << kv::defint_autostep(p22, start, end, 12);
	std::cout << "\n";

	P39 p39;

	p39.start_time(start);
	p39.stop_time(end);

	std::cout << kv::defint_autostep(p39, start, end, 12);
	std::cout << "\n";

	P41 p41;

	p41.start_time(start);
	p41.stop_time(end);

	std::cout << kv::defint_autostep(p41, start, end, 12);
	std::cout << "\n";

	P44 p44;

	p44.start_time(start);
	p44.stop_time(end);

	std::cout << kv::defint_autostep(p44, start, end, 12);
	std::cout << "\n";

	P45 p45;

	p45.start_time(start);
	p45.stop_time(end);

	std::cout << kv::defint_autostep(p45, start, end, 12);
	std::cout << "\n";

	P47 p47;

	p47.start_time(start);
	p47.stop_time(end);

	std::cout << kv::defint_autostep(p47, start, end, 12);
	std::cout << "\n";

	P51 p51;

	p51.start_time(start);
	p51.stop_time(end);

	std::cout << kv::defint_autostep(p51, start, end, 12);
	std::cout << "\n";

	P57 p57;

	p57.start_time(start);
	p57.stop_time(end);

	std::cout << kv::defint_autostep(p57, start, end, 12);
	std::cout << "\n";

	P63 p63;

	p63.start_time(start);
	p63.stop_time(end);

	std::cout << kv::defint_autostep(p63, start, end, 12);
	std::cout << "\n";

	P73 p73;

	p73.start_time(start);
	p73.stop_time(end);

	std::cout << kv::defint_autostep(p73, start, end, 12);
	std::cout << "\n";

	P74 p74;

	p74.start_time(start);
	p74.stop_time(end);

	std::cout << kv::defint_autostep(p74, start, end, 12);
	std::cout << "\n";

	P75 p75;

	p75.start_time(start);
	p75.stop_time(end);

	std::cout << kv::defint_autostep(p75, start, end, 12);
	std::cout << "\n";

	P77 p77;

	p77.start_time(start);
	p77.stop_time(end);

	std::cout << kv::defint_autostep(p77, start, end, 12);
	std::cout << "\n";
}
