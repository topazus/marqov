#ifndef HELPERS_H
#define HELPERS_H
#include <algorithm>
#include <vector>
#include <cmath>
#include <type_traits>
#include <tuple>
#include "registry.h"

/**
 * A helper to create the cartesian product of a set of containers.
 * @tparam T the head type of the template parameter pack
 * @tparam Ts the remaining parameters.
 */
template <typename T, typename ... Ts>
struct cartprodhelper
{
    typedef decltype(std::tuple_cat(std::make_tuple(std::declval<typename T::value_type>()),
                                    std::declval<typename cartprodhelper<Ts...>::RetType>())) RetType;///< A tuple consisting of our element and recursively collected elements.

    /** Recursively build up a vector of all possible combinations.
     * 
     * @param t A container of elements.
     * @param ts The remaining containers
     * @return a vector of the cartesian product of t with all elements in ts
     */
    static auto call(const T& t, Ts ... ts)
    {
    /* Plan: First we recurse into the remaining template parameters and obtain the resulting vector.
     *       Now we can generate new tuples by taking an entry from the other tuple and combining them with all of our elements.
     *       Then we return the resulting vector in the hope that it might become part of another recursion.
     */
        //get other vector through recursion
        auto othertuplevec = cartprodhelper<Ts...>::call(ts...);
        typedef decltype(std::tuple_cat(std::make_tuple(std::declval<typename T::value_type>()), std::declval<typename decltype(othertuplevec)::value_type>())) NewTupleType;
        std::vector<NewTupleType> retval;
        
        for(auto othertuple : othertuplevec) // for each other vector
        {
            for(auto myelem : t) //create an element that has an element from my vecto attached in front.
                retval.push_back(std::tuple_cat(std::make_tuple(myelem), othertuple) );
        }
        return retval;
    }
};

/** End of recursion: Just a single list.
 * 
 * @tparam T the container of the last argument in cart_prod.
 */
template <typename T>
struct cartprodhelper<T>
{
    typedef decltype(std::make_tuple(std::declval<typename T::value_type>())) RetType; ///< the type of the elements in the last container.
    
    /** The end of the recursion.
     * 
     * Builds a vector of one element tuples.
     * @param inp A container of elements.
     * @return a vector of one-element tuples.
     */
    static auto call(const T& inp)
    {
        typedef std::tuple<typename T::value_type> VecElemType;
        std::vector<VecElemType> retval;
        std::transform(inp.begin(), inp.end(), std::back_inserter(retval),
                       [](typename T::value_type elem) -> VecElemType {return std::make_tuple(elem);});
        return retval;
    }
};

/** This creates a cartesian product of an arbitrary number of input containers. 
 * 
 *  They need to conform to the basic STL Interface. They should expose value_type and provide iterators.
 *  @tparam Ts The template pack for the arbitry containers.
 * 
 *  @param vals The pack of containers.
 *  @return A vector of tuples. Each tuple contains an entry from the cartesian product.
 * FIXME: unclear what happens if one of the parameters themselves is supposed to be a tuple.
 */
template <typename ... Ts>
auto cart_prod(Ts ... vals)
{
    return cartprodhelper<Ts...>::call(vals...);
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

bool startswith(const std::string longword, const std::string shortword)
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
		k -= index * std::pow(nBin,i);
	}

	return indices;
}

#endif
