#ifndef ODE_PARAM_HPP
#define ODE_PARAM_HPP

#include <cmath>
#include <limits>

namespace kv {


template <class T> struct ode_param {
	int order_v;
	bool autostep_v;
	T tolerance_v;
	int iteration_v;
	int verbose_v;
	int ep_reduce_v;
	int ep_reduce_limit_v;
	int restart_max_v;

	ode_param() :
		order_v(24),
		autostep_v(true),
		tolerance_v(std::numeric_limits<T>::epsilon()),
		iteration_v(2),
		verbose_v(0),
		ep_reduce_v(0),
		ep_reduce_limit_v(0),
		restart_max_v(1)
	{}

	ode_param& set_order(int x) {
		order_v = x;
		return *this;
	}
	ode_param& set_autostep(bool x) {
		autostep_v = x;
		return *this;
	}
	ode_param& set_tolerance(T x) {
		tolerance_v = x;
		return *this;
	}
	ode_param& set_iteration(int x) {
		iteration_v = x;
		return *this;
	}
	ode_param& set_verbose(int x) {
		verbose_v = x;
		return *this;
	}
	ode_param& set_ep_reduce(int x) {
		ep_reduce_v = x;
		return *this;
	}
	ode_param& set_ep_reduce_limit(int x) {
		ep_reduce_limit_v = x;
		return *this;
	}
	ode_param& set_restart_max(int x) {
		restart_max_v = x;
		return *this;
	}
};

} // namespace kv

#endif // ODE_PARAM_HPP
