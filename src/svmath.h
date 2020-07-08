#ifndef SVMATH_H
#define SVMATH_H
#include <tuple>
#include <cmath>
#include <algorithm>

// ------- elementary state vector calculus

template <class StateVector>
StateVector operator + (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] += rhs[i];
    return res;
}

template <class StateVector>
StateVector operator - (StateVector lhs,  StateVector rhs)
{
    StateVector res(lhs);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    res[i] -= rhs[i];
    return res;
}


// mult is a component-wise multiplication which is used e.g. in the
// Metropolis algorithm 

inline double mult(const double& a, const double& b)
{
    return a*b;
}

inline double mult(const double& a, const int& b)
{
    return a*double(b);
}

template <class VecType, class StateVector>
inline StateVector mult(const VecType& a, const StateVector& b)
{
    StateVector retval(b);
    for(int i = 0; i < std::tuple_size<StateVector>::value; ++i)
    retval[i] *= a[i];
    return retval;
}

template <class StateVector>
inline StateVector mult(const int& a, const StateVector& b)
{
	StateVector retval(b);
	for(int i = 0; i < std::tuple_size<StateVector>::value; ++i) retval[i] *= a;
	return retval;
}

template <class StateVector>
inline StateVector mult(const double& a, const StateVector& b)
{
	StateVector retval(b);
	for(int i = 0; i < std::tuple_size<StateVector>::value; ++i) retval[i] *= a;
	return retval;
}




inline double dot(const double& a, const double& b)
{
    return a*b;
}

template<class VecType>
inline typename VecType::value_type dot(const VecType& a, const VecType& b)
{
//    typedef typename VecType::value_type FPType;
    return std::inner_product(begin(a), end(a), begin(b), 0.0);
}


template <class StateVector>
inline void reflect(StateVector& vec, const StateVector mirror)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	const double dotp = dot(vec,mirror);

	for (int i=0; i<SymD; i++) vec[i] -= 2*dotp*mirror[i];
}	

template <class Container>
inline void normalize(Container& a)
{
	typename Container::value_type tmp_abs=std::sqrt(dot(a, a));

	for (decltype(a.size()) i = 0; i < a.size(); ++i) a[i] /= tmp_abs;
}



template <class StateVector>
inline void coutsv(StateVector& vec)
{
	const int SymD = std::tuple_size<StateVector>::value;
	
	for (int i=0; i<SymD; i++) cout << vec[i] << "\t";
	cout << endl;
}

#endif
