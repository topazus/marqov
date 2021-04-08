#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <type_traits>
#include <cmath>
#include "rngcache.h"
#include "metropolishelpers.h"

// -------------------------- Metropolis Algorithm -----------------------------

namespace MARQOV
{
/**
 * A class to encapsulate the Metropolis update.
 * Using the power of partial template class specializations it is possible to define moves
 * peculiar to your model.
 */
	template <class Hamiltonian, class Lattice>
	struct Metropolis
	{
		template <class StateSpace, class M, class RNGType>
		static int move(const Hamiltonian& ham, 
		   				const Lattice& grid, 
						StateSpace& statespace, 
						M& metro, 
						RNGCache<RNGType>& rng, 
						double beta, 
						int rsite);
	};

    template <class Hamiltonian, class Lattice>
    template <class StateSpace, class M, class RNGType>
    int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
    										const Lattice& grid, 
										StateSpace& statespace, 
										M& metro, RNGCache<RNGType>& rng, 
										double beta, 
										int rsite)
    {
        static_assert(Is_Container<decltype(std::declval<Hamiltonian>().interactions)>::value, "[MARQOV::Metropolis] COMPILATION FAILED: interactions are not a container.");
		typedef typename Hamiltonian::StateVector StateVector;
        
		// old state vector at rsite
		StateVector& svold = statespace[rsite];
		// propose new configuration
		StateVector svnew = metro.newsv(svold);
		        
		// interaction part
		double interactionenergydiff = 0;
		for (typename std::remove_cv<decltype(ham.interactions.size())>::type a=0; a<ham.interactions.size(); ++a)
		{
			typedef decltype(ham.interactions[a]->get(statespace[0])) InteractionType;
			typedef decltype(callbonds<Lattice>(grid, a, rsite, 0, ham.interactions[a]->get(statespace[0]))) BondType;
			typedef typename Promote_Array<InteractionType, BondType>::CommonArray CommonArray;
			            
			auto nbrs = getnbrs<Lattice>(grid, a, rsite);
			            
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
			auto nbrs = getflexnbrs<Lattice>(grid, c, rsite);
			auto diff = ham.multisite[c]->diff(rsite, svold, svnew, nbrs, statespace, grid);
			flexenergydiff += dot(ham.multisite[c]->k, diff);
		}
        

		// collect everything
		double dE = interactionenergydiff + onsiteenergydiff + flexenergydiff;
		
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
		return Metropolis<Hamiltonian, Grid>::move(this->ham,
											this->grid, 
											statespace, 
											this->metro, 
											this->rngcache, 
											this->beta, 
											rsite);
	}
};


#endif
