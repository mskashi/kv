#ifndef ODE_PARAM_HPP
#define ODE_PARAM_HPP

#include <limits>

namespace kv {


template <class T> struct ode_param {
	int order;
	bool autostep;
	T epsilon;
	int iteration;
	int verbose;
	int ep_reduce;
	int ep_reduce_limit;
	int restart_max;

	ode_param() :
		order(24),
		autostep(true),
		epsilon(std::numeric_limits<T>::epsilon()),
		iteration(2),
		verbose(0),
		ep_reduce(0),
		ep_reduce_limit(0),
		restart_max(2)
	{}

	ode_param& set_order(int x) {
		order = x;
		return *this;
	}
	ode_param& set_autostep(bool x) {
		autostep = x;
		return *this;
	}
	ode_param& set_epsilon(T x) {
		epsilon = x;
		return *this;
	}
	ode_param& set_iteration(int x) {
		iteration = x;
		return *this;
	}
	ode_param& set_verbose(int x) {
		verbose = x;
		return *this;
	}
	ode_param& set_ep_reduce(int x) {
		ep_reduce = x;
		return *this;
	}
	ode_param& set_ep_reduce_limit(int x) {
		ep_reduce_limit = x;
		return *this;
	}
	ode_param& set_restart_max(int x) {
		restart_max = x;
		return *this;
	}
};

} // namespace kv

#endif // ODE_PARAM_HPP
