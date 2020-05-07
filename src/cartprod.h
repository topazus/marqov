#ifndef CARTPROD_H
#define CARTPROD_H
#include <tuple>
#include <vector>
#include <algorithm>

template <typename T, typename ... Ts>
struct cartprodhelper
{
    typedef decltype(std::tuple_cat(std::make_tuple(std::declval<typename T::value_type>()),
                                    std::declval<typename cartprodhelper<Ts...>::RetType>())) RetType;
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

/** End of recursion: Just a single list. */
template <typename T>
struct cartprodhelper<T>
{
    typedef decltype(std::make_tuple(std::declval<typename T::value_type>())) RetType;
    static auto call(const T& inp)
    {
        typedef std::tuple<typename T::value_type> VecElemType;
        std::vector<VecElemType> retval;
        std::transform(inp.begin(), inp.end(), std::back_inserter(retval),
                       [](typename T::value_type elem) -> VecElemType {return std::make_tuple(elem);});
        return retval;
    }
};

/** The creates a cartesian product of an arbitrary number of input containers. 
 *  They need to conform to the basic STL Interface. they should expose value_type and provide iterators.
 *  @return a vector of tuples. Each tuple contains an entry from the cartesian product.
 * FIXME: unclear what happens if one of the parameters themselves is supposed to be a tuple.
 */
template <typename ... Ts>
auto cart_prod(Ts ... vals)
{
    return cartprodhelper<Ts...>::call(vals...);
}
#endif
