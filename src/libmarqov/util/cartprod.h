#ifndef CARTPROD_H
#define CARTPROD_H
#include <algorithm>
#include <vector>
#include <type_traits>
#include <tuple>

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


#endif
