#include <iostream>
#include <limits>
#include <kv/ode-maffine.hpp>

#include "burkardt-ode.hpp"


namespace ub = boost::numeric::ublas;

typedef kv::interval<double> itvd;


int main()
{
	ub::vector<itvd> x;
	itvd start, end;
	kv::ode_param<double> p;
	bool r;

	std::cout.precision(17);


	P01 p01;

	p01.initial_value(x);
	p01.start_time(start);
	p01.stop_time(end);

	std::cout << "P01\n";
	r = odelong_maffine(p01, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P02 p02;

	p02.initial_value(x);
	p02.start_time(start);
	p02.stop_time(end);

	std::cout << "P02\n";
	r = odelong_maffine(p02, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P03 p03;

	p03.initial_value(x);
	p03.start_time(start);
	p03.stop_time(end);

	std::cout << "P03\n";
	r = odelong_maffine(p03, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P04 p04;

	p04.initial_value(x);
	p04.start_time(start);
	p04.stop_time(end);

	std::cout << "P04\n";
	r = odelong_maffine(p04, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P05 p05;

	p05.initial_value(x);
	p05.start_time(start);
	p05.stop_time(end);

	std::cout << "P05\n";
	r = odelong_maffine(p05, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P06 p06;

	p06.initial_value(x);
	p06.start_time(start);
	p06.stop_time(end);

	std::cout << "P06\n";
	r = odelong_maffine(p06, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P07 p07;

	p07.initial_value(x);
	p07.start_time(start);
	p07.stop_time(end);

	std::cout << "P07\n";
	r = odelong_maffine(p07, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P08 p08;

	p08.initial_value(x);
	p08.start_time(start);
	p08.stop_time(end);

	std::cout << "P08\n";
	r = odelong_maffine(p08, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P09 p09;

	p09.initial_value(x);
	p09.start_time(start);
	p09.stop_time(end);
	// too high order (ex. order=24) causes zero-division.
	p = kv::ode_param<double>().set_order(12);

	std::cout << "P09\n";
	r = odelong_maffine(p09, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P10 p10;

	p10.initial_value(x);
	p10.start_time(start);
	p10.stop_time(end);

	std::cout << "P10\n";
	r = odelong_maffine(p10, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P11 p11;

	p11.initial_value(x);
	p11.start_time(start);
	p11.stop_time(end);
	p = kv::ode_param<double>().set_order(12);

	std::cout << "P11\n";
	r = odelong_maffine(p11, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P12 p12;

	p12.initial_value(x);
	p12.start_time(start);
	p12.stop_time(end);
	// maffineの順序入れ換えで20まで伸びた。
	// end = 8.;

	std::cout << "P12\n";
	r = odelong_maffine(p12, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P13 p13;

	p13.initial_value(x);
	p13.start_time(start);
	p13.stop_time(end);

	std::cout << "P13\n";
	r = odelong_maffine(p13, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P14 p14;

	p14.initial_value(x);
	p14.start_time(start);
	p14.stop_time(end);
	p = kv::ode_param<double>().set_order(18).set_epsilon(1e-10);

	std::cout << "P14\n";
	r = odelong_maffine(p14, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P15 p15;

	p15.initial_value(x);
	p15.start_time(start);
	p15.stop_time(end);
	end = 5.; // 20 is difficult
	p = kv::ode_param<double>().set_order(12).set_restart_max(10);

	std::cout << "P15\n";
	r = odelong_maffine(p15, x, start, end, p);
	// r = odelong(p15, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// Kepler問題

	P16 p16;

	p16.initial_value(x);
	p16.start_time(start);
	p16.stop_time(end);

	std::cout << "P16\n";
	r = odelong_maffine(p16, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	// P17-P20はP16の初期値違い。(離心率の違い)
	// P16, P17は解けるが、P18-P20は厳しい。
	// P18, 19, 20はそれぞれrestart回数2, 3, 6くらいで解ける。

	P17 p17;

	p17.initial_value(x);
	p17.start_time(start);
	p17.stop_time(end);

	std::cout << "P17\n";
	r = odelong_maffine(p17, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P18 p18;

	p18.initial_value(x);
	p18.start_time(start);
	p18.stop_time(end);
	p = kv::ode_param<double>().set_restart_max(2);

	std::cout << "P18\n";
	r = odelong_maffine(p18, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P19 p19;

	p19.initial_value(x);
	p19.start_time(start);
	p19.stop_time(end);
	p = kv::ode_param<double>().set_restart_max(3);

	std::cout << "P19\n";
	r = odelong_maffine(p19, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P20 p20;

	p20.initial_value(x);
	p20.start_time(start);
	p20.stop_time(end);
	p = kv::ode_param<double>().set_order(10).set_restart_max(6);

	std::cout << "P20\n";
	r = odelong_maffine(p20, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P21 p21;

	p21.initial_value(x);
	p21.start_time(start);
	p21.stop_time(end);

	std::cout << "P21\n";
	r = odelong_maffine(p21, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P22 p22;

	p22.initial_value(x);
	p22.start_time(start);
	p22.stop_time(end);

	std::cout << "P22\n";
	r = odelong_maffine(p22, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P23 p23;

	p23.initial_value(x);
	p23.start_time(start);
	p23.stop_time(end);

	std::cout << "P23\n";
	r = odelong_maffine(p23, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P24 p24;

	p24.initial_value(x);
	p24.start_time(start);
	p24.stop_time(end);

	std::cout << "P24\n";
	r = odelong_maffine(p24, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P25 p25;

	p25.initial_value(x);
	p25.start_time(start);
	p25.stop_time(end);

	std::cout << "P25\n";
	r = odelong_maffine(p25, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// P26, P27はstep functionを含む
	// P28は絶対値を含む
	// P29はifを含む
	// P30は特殊な^(1/3)を含む

	P31 p31;

	p31.initial_value(x);
	p31.start_time(start);
	p31.stop_time(end);

	std::cout << "P31\n";
	r = odelong_maffine(p31, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P32 p32;

	p32.initial_value(x);
	p32.start_time(start);
	p32.stop_time(end);

	std::cout << "P32\n";
	r = odelong_maffine(p32, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P33 p33;

	p33.initial_value(x);
	p33.start_time(start);
	p33.stop_time(end);

	std::cout << "P33\n";
	r = odelong_maffine(p33, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P34 p34;

	p34.initial_value(x);
	p34.start_time(start);
	p34.stop_time(end);

	std::cout << "P34\n";
	r = odelong_maffine(p34, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P35 p35;

	p35.initial_value(x);
	p35.start_time(start);
	p35.stop_time(end);

	std::cout << "P35\n";
	r = odelong_maffine(p35, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	// Duffing。なぜか答えが全く違うので要確認--
	P36 p36;

	p36.initial_value(x);
	p36.start_time(start);
	p36.stop_time(end);

	std::cout << "P36\n";
	r = odelong_maffine(p36, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	// Duffing。こっちも全然違う。--
	P37 p37;

	p37.initial_value(x);
	p37.start_time(start);
	p37.stop_time(end);

	std::cout << "P37\n";
	r = odelong_maffine(p37, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

	// 全く進めない。厳しい。--
	// restartするようにしたら、restart一回で簡単に進めた。--
	P38 p38;

	p38.initial_value(x);
	p38.start_time(start);
	p38.stop_time(end);
	// end = 3.001;

	std::cout << "P38\n";
	r = odelong_maffine(p38, x, start, end);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P39 p39;

	p39.initial_value(x);
	p39.start_time(start);
	p39.stop_time(end);
	p = kv::ode_param<double>().set_restart_max(10);

	std::cout << "P39\n";
	r = odelong_maffine(p39, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}


	P40 p40;

	p40.initial_value(x);
	p40.start_time(start);
	p40.stop_time(end);
	p = kv::ode_param<double>().set_restart_max(10);

	std::cout << "P40\n";
	r = odelong_maffine(p40, x, start, end, p);

	if (!r) {
		std::cout << "can't calculate verified solution\n";
	} else {
		std::cout << x << "\n";
		std::cout << end << "\n";
	}

}
