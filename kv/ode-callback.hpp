/*
 * Copyright (c) 2013-2017 Masahide Kashiwagi (kashi@waseda.jp)
 */

#ifndef ODE_CALLBACK
#define ODE_CALLBACK

namespace kv {

// base class of ode_callback

template <class T> struct ode_callback {
	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		return true;
	}
};


// example1: callback function for dense print

template <class T> struct ode_callback_dense_print : ode_callback<T> {
	interval<T> start_g;
	interval<T> step;

	ode_callback_dense_print(interval<T> start, interval<T> step) : start_g(start), step(step) {}

	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		int i, j;
		interval<T> t;
		ub::vector< interval<T> > y;
		psa< interval <T> > tmp;
		int s = result.size();
		y.resize(s);

		for (i = (int)ceil(((start - start_g) / step).lower()); i<=(int)floor(((end - start_g) / step).upper()); i++) {
			t = start_g + step * i - start;
			for (j=0; j<s; j++) {
				tmp = result(j);
				y(j) = eval(tmp, t);
			}
			std::cout << "t: " << start + t << "\n";
			std::cout << y << "\n";
		}
		return true;
	}
};

// example1: callback function for getting all information by list

template <class T> struct ode_callback_list : ode_callback<T> {
	std::list< interval<T> > &start_list;
	std::list< interval<T> > &end_list;
	std::list< ub::vector< psa< interval<T> > > > &psa_list;

	ode_callback_list(std::list< interval<T> >& start_list, std::list< interval<T> >& end_list, std::list< ub::vector< psa< interval<T> > > >& psa_list) : start_list(start_list), end_list(end_list), psa_list(psa_list) {} 

	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		start_list.push_back(start);
		end_list.push_back(end);
		psa_list.push_back(result);
		return true;
	}
};

// example1: callback function for getting dense output by list

template <class T> struct ode_callback_dense_list : ode_callback<T> {
	interval<T> start_g;
	interval<T> step;
	std::list< interval<T> > &time_list;
	std::list< ub::vector< interval<T> > > &value_list;

	ode_callback_dense_list(interval<T> start, interval<T> step, std::list< interval<T> >& time_list, std::list< ub::vector< interval<T> > >& value_list) : start_g(start), step(step), time_list(time_list), value_list(value_list) {} 

	virtual bool operator()(const interval<T>& start, const interval<T>& end, const ub::vector< interval<T> >& x_s, const ub::vector< interval<T> >& x_e, const ub::vector< psa< interval<T> > >& result) const {
		int i, j;
		interval<T> t;
		ub::vector< interval<T> > y;
		psa< interval <T> > tmp;
		int s = result.size();
		y.resize(s);

		for (i = (int)ceil(((start - start_g) / step).lower()); i<=(int)floor(((end - start_g) / step).upper()); i++) {
			t = start_g + step * i - start;
			for (j=0; j<s; j++) {
				tmp = result(j);
				y(j) = eval(tmp, t);
			}
			time_list.push_back(start + t);
			value_list.push_back(y);
		}
		return true;
	}
};

} // namespace kv

#endif // ODE_CALLBACK
