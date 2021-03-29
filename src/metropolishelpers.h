#ifndef METROPOLISHELPERS_H
#define METROPOLISHELPERS_H

namespace MARQOV
{
	//A helper to decide in the Metropolis code whether a lattice provides the getbond function
	template<class L, class=void> 
	struct has_bonds : std::false_type {};
	
	template<class Lattice> 
	struct has_bonds<Lattice, MARQOV::type_sink_t< 
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
	auto callbonds_helper(const Lattice& grid, int a, int i, int j, NbrType nbr, std::true_type)
	{
		auto cpl = grid.bnds(a,i)[j];  
		// improve me: if only one particular bond is needed, 
		// why load all in the first place

		return mult(cpl, nbr);
	}
	
	template <class Lattice, class NbrType>
	auto callbonds_helper(const Lattice& grid, int a, int i, int j, NbrType nbr, std::false_type)
	{
		return nbr;
	}
	
	template <class Lattice, class NbrType>
	auto callbonds(const Lattice& grid, int a, int i, int j, NbrType nbr)
	{
		return callbonds_helper(grid, a, i, j, nbr, typename has_bonds<Lattice>::type());
	}
	
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
	
	template<class, class = void> 
	struct has_terms : std::false_type {};
	
	template<class Grid>
	struct has_terms<Grid, MARQOV::detail::type_sink_t<
		decltype(std::declval<Grid>().termselector(std::declval<int>()))
		>> : std::true_type {};
	// in C++17 (which we don't use), this can be solved with void_t
	
	template <class Grid>
	std::vector<int> get_terms_helper(Grid& grid, int idx, std::false_type) {return {-1};}
	
	template <class Grid>
	std::vector<int> get_terms_helper(Grid& grid, int idx, std::true_type) {return grid.termselector(idx);}
	
	template <class Grid>
	std::vector<int> get_terms(Grid& grid, int idx) {	return get_terms_helper<Grid>(grid, idx, has_terms<Grid>{}); }
	
	


	// A helper to check whether the lattice provides a getnbrs method
	
	template<class, class = void> 
	struct has_nbrs : std::false_type {};
	
	template<class Grid>
	struct has_nbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().nbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto getnbrs_helper(Grid& grid, int fam, int idx, std::false_type) 
	{
		cout << "nbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto getnbrs_helper(Grid& grid, int fam, int idx, std::true_type)  
	{
		return grid.nbrs(fam,idx); 
	}
	
	template <class Grid>
	auto getnbrs(Grid& grid, int fam, int idx) 
	{
		return getnbrs_helper<Grid>(grid, fam, idx, has_nbrs<Grid>{}); 
	}




	// A helper to check whether the lattice provides a getflexnbrs method
	
	template<class, class = void> 
	struct has_flexnbrs : std::false_type {};
	
	template<class Grid>
	struct has_flexnbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().getflexnbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto getflexnbrs_helper(Grid& grid, int fam, int idx, std::false_type) 
	{
		cout << "flexnbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto getflexnbrs_helper(Grid& grid, int fam, int idx, std::true_type)  
	{
		return grid.getflexnbrs(fam,idx); 
	}
	
	template <class Grid>
	auto getflexnbrs(Grid& grid, int fam, int idx) 
	{
		return getflexnbrs_helper<Grid>(grid, fam, idx, has_flexnbrs<Grid>{}); 
	}


};


#endif
