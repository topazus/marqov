#ifndef HELPERS_H
#define HELPERS_H

#include <iomanip> 
#include <vector>
#include <cmath>
#include "registry.h"

std::vector<int> arange(int lower, int upper)
{
	std::vector<int> retval;
	for (int i=lower; i<upper; i++)
	{	
		retval.push_back(i);
	}
	return retval;
}



template <class T1, class T2, class T3>
class Triple
{
	public:
		Triple(T1 t1, T2 t2, T3 t3) : first(t1), second(t2), third(t3) {}
		T1 first;
		T2 second;
		T3 third;
};

template< class T1, class T2, class T3 >
constexpr auto make_triple( T1&& t, T2&& u, T3&& v) 
{
	return Triple<typename std::decay<T1>::type, typename std::decay<T2>::type, typename std::decay<T3>::type>(t,u,v);
}




void write_logfile(RegistryDB& reg, std::vector<double> loopvar)
{
	std::string logdir  = reg.Get<std::string>("mc.ini", "IO", "logdir" );
	std::string logfile = reg.Get<std::string>("mc.ini", "IO", "logfile" );
	std::ofstream os(logdir+"/"+logfile);
	os << std::setprecision(7);
	for (std::size_t i=0; i<loopvar.size(); i++) os << loopvar[i] << endl;
	os.close();
}

// //C++17 make_from_tuple from cppreference adapted for emplace.
// template <class Cont, class Latt, class Tuple, std::size_t... I>
// constexpr auto emplace_from_tuple_impl(Cont&& cont, Latt&& latt, MARQOV::MARQOVConfig&& mc, Tuple&& t, std::index_sequence<I...> )
// {
//   return cont.emplace_back(std::forward<Latt>(latt), std::forward<decltype(mc)>(mc), std::get<I>(std::forward<Tuple>(t))...) ;
// }
// 
// /** A function to construct an object in a container directly from a tuple
//  * @param cont the container where we append to.
//  * @param t the tuple containing the arguments.
//  */
// template <class Cont, class T, class Tuple, class Latt>
// constexpr auto emplace_from_tuple(Cont&& cont, Latt&& latt, T&& mc, Tuple&& t )
// {
//     return emplace_from_tuple_impl(cont, std::forward<Latt>(latt), std::forward<T>(mc), std::forward<Tuple>(t),
//         std::make_index_sequence<std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
// }




bool startswith(std::string longword, std::string shortword)
{
	if (longword.rfind(shortword, 0) == 0) return true;
	else return false;
}



std::vector<double> create_range(double rangestart, double rangefinal, int steps, std::string type="linear", int endpoint=0)
{
	std::vector<double> range(steps);

	if (type == "linear")
	{
		double rangestep = (rangefinal-rangestart)/double(steps-endpoint);

		for (int i=0; i<steps; i++) range[i] = rangestart + i*rangestep;
	}
	else
	{
		// implement me
		std::cout << "not implemented!" << std::endl;
	}

	return range;
}



// Transforming a one-dimensional, “flattened” index into the N-dimensional vector index of an N-dimensional array

// from: stackoverflow.com/questions/18932339/transforming-a-one-dimensional-flattened-index-into-the-n-dimensional-vector
// explicit formulas for 2D and 3D: stackoverflow.com/questions/18932339/transforming-a-one-dimensional-flattened-index-into-the-n-dimensional-vector

std::vector<int> IndexOf(int k, int nDim, int nBin)
{
	std::vector<int> indices;

	for (int i=0; i<nDim; i++)
	{
		double index = std::fmod( double(k)/std::pow(nBin,i), nBin);
		indices.push_back(int(index));
		k -= index * pow(nBin,i);
	}

	return indices;
}

#endif
