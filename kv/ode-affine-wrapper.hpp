/*
 * Copyright (c) 2013-2014 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_AFFINE_W_HPP
#define ODE_AFFINE_W_HPP

// Wrapper for ODE Affine

#include <kv/ode-affine.hpp>


namespace kv {

namespace ub = boost::numeric::ublas;

template <class T, class F>
int
ode_wrapper(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>()) {
	int n = init.size();
	int i, j;
	ub::vector<T> init_rad;
	ub::vector< affine<T> > new_init;

	ub::vector< affine<T> > result;
	ub::matrix<T> result_m;
	ub::vector<T> result_c;;
	ub::vector< interval<T> > result_i;

	ub::vector< affine<T> > result2;

	int r;
	interval<T> end2 = end;


	int maxnum_save;

	affine<T> s1, s2;
	interval<T> s2i;
        
	ub::vector< affine<T> > s1_save;
	ub::vector< interval<T> > s2i_save;
        
	init_rad.resize(n);
	for (i=0; i<n; i++) {
		init_rad(i) = rad(init(i));
	}

	maxnum_save = affine<T>::maxnum();
	affine<T>::maxnum() = 0;

	new_init.resize(n);
	for (i=0; i<n; i++) {
		new_init(i) = init(i).a(0) + (affine<T>)(init_rad(i) * interval<T>(-1.,1));
	}

	// ep_reduce must be 0
	ode_param<T> p2 = p;
	p2.ep_reduce = 0;
	r = ode_affine(f, new_init, start, end2, p2);
	if (r == 0) {
		affine<T>::maxnum() = maxnum_save;
		return 0;
	}

	result = new_init;
	
	result_m.resize(n,n);
	result_c.resize(n);
	result_i.resize(n);

	for (i=0; i<n; i++) {
		split(result(i), n, s1, s2);
		result_c(i) = s1.a(0);
		for (j=0; j<n; j++) result_m(i, j) = s1.a(j+1);
		result_i(i) = to_interval(s2);
	}

	affine<T>::maxnum() = maxnum_save;

	result2 = init;
	for (i=0; i<n; i++) {
		result2(i).a(0) = 0.;
		if (init_rad(i) != 0.) result2(i) /= init_rad(i);
	}
	result2 = prod(result_m, result2) + result_c;
	for (i=0; i<n; i++) {
		result2(i) += (affine<T>)result_i(i);
	}

	if (p.ep_reduce == 0) {
		s1_save.resize(n);
		s2i_save.resize(n);
		for (i=0; i<n; i++) {
			split(result2(i), maxnum_save, s1, s2);
			s2i = to_interval(s2);
			s1_save(i) = s1;
			s2i_save(i) = s2i;
		}
		affine<T>::maxnum() = maxnum_save;
		for (i=0; i<n; i++) {
			s1_save(i).resize();
			result2(i) = append(s1_save(i), (affine<T>)s2i_save(i));
		}
	} else {
		epsilon_reduce(result, p.ep_reduce, p.ep_reduce_limit);
	}

	init = result2;
	if (r == 1) end = end2;
	return r;
}

template <class T, class F>
int
odelong_wrapper(F f, ub::vector< affine<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>())
{
	int s = init.size();
	ub::vector< affine<T> > x;
	interval<T> t, t1;
	int r;
	int ret_val = 0;

	x = init;
	t = start;
	p.set_autostep(true);
	while (1) {
		t1 = end;

		r = ode_wrapper(f, x, t, t1, p);
		if (r == 0) {
			if (ret_val == 1) {
				init = x;
				end = t;
			}
			return ret_val;
		}
		ret_val = 1;
		if (p.verbose == 1) {
			std::cout << "t: " << t1 << "\n";
			std::cout << to_interval(x) << "\n";
		}
		if (r == 2) {
			init = x;
			return 2;
		}
		t = t1;
	}
}

template <class T, class F>
int
odelong_wrapper(F f, ub::vector< interval<T> >& init, const interval<T>& start, interval<T>& end, ode_param<T> p = ode_param<T>())
{
	int s = init.size();
	int i;
	ub::vector< affine<T> > x;
	int r;
	interval<T> end2 = end;

	affine<T>::maxnum() = 0;
	x = init;

	r = odelong_wrapper(f, x, start, end2, p);
	if (r == 0) return 0;

	for (i=0; i<s; i++) init(i) = to_interval(x(i));
	if (r == 1) end = end2;

	return r;
}

} // namespace kv

#endif // ODE_AFFINE_W_HPP
