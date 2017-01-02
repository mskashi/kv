/*
 * Copyright (c) 2013 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_CALLBACK
#define ODE_CALLBACK

namespace kv {

template <class T> struct ode_callback {
	virtual void operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
	}
};

} // namespace kv

#endif // ODE_CALLBACK
