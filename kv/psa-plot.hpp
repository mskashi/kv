/*
 * Copyright (c) 2013-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef PSA_PLOT_HPP
#define PSA_PLOT_HPP

#include <boost/numeric/ublas/vector.hpp>

#include <kv/psa.hpp>
#include <kv/matplotlib.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template<class T1, class T2>
void psa_plot(const psa<T1>& x, const T2& offset, const matplotlib& g, int div = 50, const char * color = "blue") 
{
	double s, e, off, p1, p2, l1, u1, l2, u2;
	int i;

	s = psa<T1>::domain().lower();
	e = psa<T1>::domain().upper();
	// off = mid(offset);
	off = offset;

	for (i=0; i<div; i++) {
		p1 = s + (e - s) * i / div;
		p2 = s + (e - s) * (i+1) / div;
		l1 = eval(x, (T1)p1).lower();
		u1 = eval(x, (T1)p1).upper();
		l2 = eval(x, (T1)p2).lower();
		u2 = eval(x, (T1)p2).upper();
		g.line(off + p1, l1, off + p2, l2, color);
		g.line(off + p1, u1, off + p2, u2, color);
	}
}

} // namespace kv

#endif // PSA_PLOT_HPP
