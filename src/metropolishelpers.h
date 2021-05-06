/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef METROPOLISHELPERS_H
#define METROPOLISHELPERS_H
#include <utility>
#include <type_traits>

namespace MARQOV
{
    /**
     * has_bonds utility function.
     * 
     * This function decays to a bool_type to denote whether a lattice provides
     * the getbond function.
     * @tparam L The lattice that we are querying.
     */
	template<class L, class=void> 
	struct has_bonds : std::false_type {};
	
	template<class Lattice> 
	struct has_bonds<Lattice, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Lattice>().bnds(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};
    
    
    
	// what are the input parameters:
	// a: bond family (not yet implemented)
	// i: site to be looked at
	// j: j-th neighbour of i (encodes the site to which the bond goes)
	
	
	// there are three possibilities
	
	//    nbr  |   cpl                      return type
	// -----------------                 ----------------
	//  scalar   scalar                     -> scalar
	//  vector   vector  (same length!)	-> vector
	//  vector   scalar 				-> vector
	

	// note also that "nbr" is always a "statvector" (std::array), 
	// whereas "cpl" can be anything (typically though: int, double, std::vector)
    
    
	template <class Lattice, class NbrType>
	inline auto callbonds_helper(const Lattice& grid, int a, int i, int j, NbrType nbr, std::true_type)
	{
		auto cpl = grid.bnds(a,i)[j];  
		// improve me: if only one particular bond is needed, 
		// why load all in the first place

		return mult(cpl, nbr);
	}
	
	template <class Lattice, class NbrType>
	inline auto callbonds_helper(const Lattice& grid, int a, int i, int j, NbrType nbr, std::false_type)
	{
		return nbr;
	}
	
	template <class Lattice, class NbrType>
	inline auto callbonds(const Lattice& grid, int a, int i, int j, NbrType nbr)
	{
		return callbonds_helper(grid, a, i, j, nbr, typename has_bonds<Lattice>::type());
	}
	
	/**
     * Promote_Array utility class.
     * 
     * This function takes two types and tries to figure out,which one is wider
     * @tparam A the first type.
     * @tparam B the other type.
     */
	template<class A, class B>
	struct Promote_Array
	{
		//Get Element types
		typedef decltype(dot(std::declval<A>(), std::declval<A>())) AElemType;
		typedef decltype(dot(std::declval<B>(), std::declval<B>())) BElemType;
		// get result type of addition
		typedef typename std::common_type<AElemType, BElemType>::type CommonType;
		// return A if the common type is A, else B
		typedef typename std::conditional<std::is_same<AElemType, CommonType>::value, A, B>::type CommonArray;
	};

	// A helper to check whether the lattice provides a termselector method
	/**
     * has_terms utility class.
     * This class tries to figure out whether there is a has_terms function.
     * @tparam Grid the Grid which we check.
     */
	template<class Grid , class = void> 
	struct has_terms : std::false_type {};
	
	template<class Grid>
	struct has_terms<Grid, MARQOV::detail::type_sink_t<
		decltype(std::declval<Grid>().termselector(std::declval<int>()))
		>> : std::true_type {};
	// in C++17 (which we don't use), this can be solved with void_t
	
	template <class Grid>
	std::vector<int> get_terms_helper(const Grid& grid, int idx, std::false_type) {return {-1};}
	
	template <class Grid>
	std::vector<int> get_terms_helper(const Grid& grid, int idx, std::true_type) {return grid.termselector(idx);}
	
	template <class Grid>
	std::vector<int> get_terms(const Grid& grid, int idx) {	return get_terms_helper<Grid>(grid, idx, has_terms<Grid>{}); }

	/** has_nbrs utility struct
     * 
     * This struct checks whether the lattice provides a getnbrs method.
     * @tparam Grid the lattice that we query.
     */
	template<class Grid, class = void> 
	struct has_nbrs : std::false_type {};
	
	template<class Grid>
	struct has_nbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().nbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto getnbrs_helper(const Grid& grid, int fam, int idx, std::false_type)
	{
		cout << "nbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto getnbrs_helper(const Grid& grid, int fam, int idx, std::true_type)
	{
		return grid.nbrs(fam,idx); 
	}
	
	template <class Grid>
	auto getnbrs(const Grid& grid, int fam, int idx)
	{
		return getnbrs_helper<Grid>(grid, fam, idx, has_nbrs<Grid>{}); 
	}

    /** has_flexnbrs utility struct.
     * 
     * This utility class checks whether the lattice provides a getflexnbrs
     * method.
     * @tparam Grid The grid that we query.
     */
	template<class Grid, class = void> 
	struct has_flexnbrs : std::false_type {};
	
	template<class Grid>
	struct has_flexnbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().getflexnbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto getflexnbrs_helper(const Grid& grid, int fam, int idx, std::false_type) 
	{
		cout << "flexnbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto getflexnbrs_helper(const Grid& grid, int fam, int idx, std::true_type)  
	{
		return grid.getflexnbrs(fam,idx); 
	}
	
	template <class Grid>
	auto getflexnbrs(const Grid& grid, int fam, int idx) 
	{
		return getflexnbrs_helper<Grid>(grid, fam, idx, has_flexnbrs<Grid>{}); 
	}
    
    /**
     * Is_Container utility struct.
     * 
     * Helpers to determine if the interactions are container-like.
     * By that we mean whether it has a .size() method and an array access
     * operator.
     * See metropolis.h for an example.
     * @see metropolis.h
     * @tparam Cont The Container that we query.
     */
    template <class Cont, class = void, class = void>
    struct Is_Container
    {
        static constexpr bool value = false;
    };
    
    template <class Cont>
    struct Is_Container<Cont,
    MARQOV::detail::type_sink_t<decltype(std::declval<Cont>().size())>,
    MARQOV::detail::type_sink_t<decltype(std::declval<Cont>().operator[](std::declval<std::size_t>()))>
    >
    {
        static constexpr bool value = true;
    };
};


#endif
