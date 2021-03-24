#ifndef METROPOLIS_H
#define METROPOLIS_H
#include "rngcache.h"

namespace MARQOV
{
    //A helper to decide in the Metropolis code whether a lattice provides the getbond function
    template<class L, class=void> 
    struct has_bonds : std::false_type {};
    
    template<class Lattice> 
    struct has_bonds<Lattice, MARQOV::type_sink_t< decltype( std::declval<Lattice>().getbnds(std::declval<int>(), std::declval<int>()) ) > > : std::true_type {};
    
    
    
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
    
    // note also that "nbr" is always a "statvector" (std::array), whereas "cpl" can be anything (typically though: int, double, std::vector)
    
    
    template <class Lattice, class NbrType>
    auto callbonds_helper(const Lattice& grid, int a, int i, int j, NbrType nbr, std::true_type)
    {
        auto cpl = grid.getbnds(a,i)[j];  // improve me: if only one particular bond is needed, why load all in the first place
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
	struct has_getnbrs : std::false_type {};
	
	template<class Grid>
	struct has_getnbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().getnbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto get_nbrs_helper(Grid& grid, int fam, int idx, std::false_type) 
	{
		cout << "getnbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto get_nbrs_helper(Grid& grid, int fam, int idx, std::true_type)  
	{
		return grid.getnbrs(fam,idx); 
	}
	
	template <class Grid>
	auto get_nbrs(Grid& grid, int fam, int idx) 
	{
		return get_nbrs_helper<Grid>(grid, fam, idx, has_getnbrs<Grid>{}); 
	}




	// A helper to check whether the lattice provides a getflexnbrs method
	
	template<class, class = void> 
	struct has_getflexnbrs : std::false_type {};
	
	template<class Grid>
	struct has_getflexnbrs<Grid, MARQOV::detail::type_sink_t< 
		decltype( std::declval<Grid>().getflexnbrs(std::declval<int>(), std::declval<int>()) ) 
		>> : std::true_type {};

	
	template <class Grid>
	auto get_flexnbrs_helper(Grid& grid, int fam, int idx, std::false_type) 
	{
		cout << "getflexnbrs not implement!" << flush;
		exit(0); // improve me
		return std::vector<int>{};
	}
	
	template <class Grid>
	auto get_flexnbrs_helper(Grid& grid, int fam, int idx, std::true_type)  
	{
		return grid.getflexnbrs(fam,idx); 
	}
	
	template <class Grid>
	auto get_flexnbrs(Grid& grid, int fam, int idx) 
	{
		return get_flexnbrs_helper<Grid>(grid, fam, idx, has_getflexnbrs<Grid>{}); 
	}







// -------------------------- Metropolis Algorithm -----------------------------

    template <class Hamiltonian, class Lattice>
    struct Metropolis
    {
        template <class StateSpace, class M, class RNGType>
        static int move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNGCache<RNGType>& rng, double beta, int rsite);
    };
    
    
    template <class Hamiltonian, class Lattice>
    template <class StateSpace, class M, class RNGType>
    int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNGCache<RNGType>& rng, double beta, int rsite)
    {
        typedef typename Hamiltonian::StateVector StateVector;
        
        // old state vector at rsite
        StateVector& svold = statespace[rsite];
        // propose new configuration
        StateVector svnew = metro.newsv(svold);
        
        // interaction part
        double interactionenergydiff = 0;
        for (typename std::remove_cv<decltype(ham.Nalpha)>::type a=0; a<ham.Nalpha; ++a)
        {
			typedef decltype(ham.interactions[a]->get(statespace[0])) InteractionType;
			typedef decltype(callbonds<Lattice>(grid, a, rsite, 0, ham.interactions[a]->get(statespace[0]))) BondType;
			typedef typename Promote_Array<InteractionType, BondType>::CommonArray CommonArray;
			            
			auto nbrs = get_nbrs<Lattice>(grid, a, rsite);
			            
			CommonArray averagevector = {0};
            
			// sum over neighbours
			for (std::size_t i = 0; i < nbrs.size(); ++i)
			{
				// index of the neighbour
				auto idx = nbrs[i];
				                
				// configuration of the neighbour
				auto nbr = ham.interactions[a]->get(statespace[idx]);
				                
				// add neighbours (also accounting for bond strength if available)
				averagevector = averagevector + MARQOV::callbonds<Lattice>(grid, a, rsite, i, nbr);
			}

			interactionenergydiff += ham.interactions[a]->J * (dot(svnew - svold, averagevector));
		}
        
        
        // onsite energy part
        auto terms = get_terms<Lattice>(grid, rsite);
        if (terms[0] == -1) terms = arange(0, ham.Nbeta);
        
        double onsiteenergydiff = 0;
        for (typename std::remove_cv<decltype(ham.Nbeta)>::type b=0; b<terms.size(); ++b)
        {
            
            // select on-site term
            const int tidx = terms[b]; 
            
            // compute the difference
            auto diff = ham.onsite[tidx]->get(svnew) - ham.onsite[tidx]->get(svold);
            
            // multiply the constant
            onsiteenergydiff += dot(ham.onsite[tidx]->h, diff);
        }
        
        // flex term energy
	   double flexenergydiff = 0;
        for (typename std::remove_cv<decltype(ham.Ngamma)>::type c=0; c<ham.Ngamma; ++c)
        {
			auto nbrs = get_flexnbrs<Lattice>(grid, c, rsite);
			auto diff = ham.multisite[c]->diff(rsite, svold, svnew, nbrs, statespace, grid);
			flexenergydiff += dot(ham.multisite[c]->k, diff);
        }
        

        // collect everything
        double dE 	= interactionenergydiff + onsiteenergydiff + flexenergydiff;
        
	   // evaluate
        int retval = 0;
        if ( dE <= 0 )
        {
            svold = svnew;
            retval = 1;
        }
        else if (rng.real() < exp(-beta*dE))
        {
            svold = svnew;
            retval = 1;
        }
        return retval;
    }
    
    // Single Metropolis update step statevectors on a lattice
    // returns an integer which encodes whether the flip attempt was successful (1) or not (0)
    template <class Grid, class Hamiltonian, template<class> class RefType>
    inline int Core<Grid, Hamiltonian, RefType>::metropolisstep(int rsite)
    {
        return Metropolis<Hamiltonian, Grid>::move(this->ham, this->grid, statespace, this->metro, this->rngcache, this->beta, rsite);
    }
};


#endif
