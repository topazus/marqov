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

#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <type_traits>
#include <cmath>
#include "rngcache.h"
#include "metropolishelpers.h"

// -------------------------- Metropolis Algorithm -----------------------------

namespace MARQOV
{


/** Helper function for the Metropolis algorithm, sums over local neighbourhood of spin
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * 
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has function grid.nbrs
 * @param std::true_type lattice has function grid.bnds
 *
 * @return sum over neighbours after applying the interaction term and weighted by respective bond strengths
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::true_type, 
											  std::true_type) 
{
	typedef typename Hamiltonian::StateVector StateVector;

	// gather neighbours and bonds
	const auto nbrs = getnbrs<Lattice>(grid, a, rsite);
	const auto bnds = getbnds<Lattice>(grid, a, rsite);

	StateVector neighbourhood = {0};

	// sum over neighbours
	for (std::size_t i=0; i<nbrs.size(); ++i)
	{
		// index of the neighbour
		auto idx = nbrs[i];
		                
		// configuration of the neighbour after applying the interaction
		auto nbr = ham.interactions[a]->get(statespace[idx]);
		                
		// sum up contributions from neighbourbood
		neighbourhood = neighbourhood + mult(bnds[i], nbr);
	}

	return neighbourhood;
}


/**
 * Helper function for the Metropolis algorithm, sums over local neighbourhood of spin
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has function grid.nbrs
 * @param std::false_type lattice has no function grid.bnds
 *
 * @return sum over neighbours after applying the interaction term
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::true_type, 
											  std::false_type) 
{
	typedef typename Hamiltonian::StateVector StateVector;

	// gather neighbours
	const auto nbrs = getnbrs<Lattice>(grid, a, rsite);

	StateVector neighbourhood = {0};

	// sum over neighbours
	for (std::size_t i=0; i<nbrs.size(); ++i)
	{
		// index of the neighbour
		auto idx = nbrs[i];
		                
		// configuration of the neighbour after applying the interaction
		auto nbr = ham.interactions[a]->get(statespace[idx]);
		                
		// sum up contributions from neighbourbood
		neighbourhood = neighbourhood + nbr;
	}

	return neighbourhood;
}


/**
 * Helper function for the Metropolis algorithm, sums over local neighbourhood of spin.
 * this overloading covers the case when the lattice has no function "nbrs".
 * in this case the general Metropolis algorithm cannot be used hence the user is given an error message
 *
 * @tparam Lattice the type of the lattice
 * @tparam Hamiltonian the type of the Hamiltonian
 * @tparam StateSpace the type of the statespace
 * @param grid the lattice
 * @param ham the Hamiltonian
 * @param statespace the statespace
 * @param a the index of the interaction term under consideration
 * @param rsite the site under consideration
 * @param std::true_type lattice has no function grid.nbrs
 * @param std::false_type lattice has no function grid.bnds
 *
 * @return does not apply 
 */
template <class Lattice, class Hamiltonian, class StateSpace>
typename Hamiltonian::StateVector nbrhoodloop(const Lattice& grid, 
											  const Hamiltonian& ham, 
											  const StateSpace& statespace, 
											  int a, 
											  int rsite, 
											  std::false_type, 
											  std::false_type) 
{
	cout << "[MARQOV] Error: The lattice does not provide the following function: nbrs" << endl;
	cout << "In order to use the general Metropolis algorithm, this function must be implemented" << endl;
	cout << "Alternatively, you can write you own specialization of the Metropolis algorithm" << endl;
	exit(0);
}




/**
 * A class to encapsulate the Metropolis update.
 * Using the power of partial template class specializations it is possible to
 * define moves peculiar to your model.
 * @tparam Hamiltonian The Hamiltonian used to generate a Metropolis move
 * @tparam Lattice The lattice used for the neighbourhood relations.
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

	/**
	 * The actual Metropolis move attempt
	 *
	 * @tparam Hamiltonian the type of the Hamiltonian
	 * @tparam Lattice the type of the lattice
	 * @tparam StateSpace the type of the state space
	 * @tparam M the type of the initializer
	 * @tparam RNG the type of the random number generator
	 * @param ham the Hamiltonian
	 * @param grid the lattice
	 * @param statespace the statespace
	 * @param metro initalizer, proposes new state vector configuration
	 * @param rng the random number generator
	 * @param beta inverse temperature
	 * @param rsite site to be considered for an update
	 *
	 * @return integer, encoding whether the update was accepted or rejected
	 */
    template <class Hamiltonian, class Lattice>
    template <class StateSpace, class M, class RNGType>
    int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
    											const Lattice& grid, 
												StateSpace& statespace, 
												M& metro, RNGCache<RNGType>& rng, 
												double beta, 
												int rsite)
    {
        static_assert(Is_Container<decltype(std::declval<Hamiltonian>().interactions)>::value, 
			"[MARQOV::Metropolis] COMPILATION FAILED: interactions are not a container.");
        static_assert(Is_Container<decltype(std::declval<Hamiltonian>().onsite)>::value, 
			"[MARQOV::Metropolis] COMPILATION FAILED: onsite terms are not a container.");
        static_assert(Is_Container<decltype(std::declval<Hamiltonian>().multisite)>::value, 
			"[MARQOV::Metropolis] COMPILATION FAILED: multisite terms are not a container.");

		// definitions
		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename std::remove_cv<decltype(ham.interactions.size())>::type InteractionSizeType;
		typedef typename std::remove_cv<decltype(ham.onsite.size())>::type OnSiteSizeType;
		typedef typename std::remove_cv<decltype(ham.multisite.size())>::type FlexSizeType; 

		// availablity of certain functions in the lattice
		typedef typename has_nbrs<Lattice>::type HasNbrs;
		typedef typename has_bnds<Lattice>::type HasBnds;
		typedef typename has_trms<Lattice>::type HasTrms;
        
		// old state vector at index rsite
		StateVector& svold = statespace[rsite];
		// propose new configuration
		StateVector svnew = metro.newsv(svold);
		        

		// I. interaction part
		double interactionenergydiff = 0;
		for (InteractionSizeType a=0; a<ham.interactions.size(); a++)
		{
			// sum up neighbourhood contributions
			auto nbrhood = nbrhoodloop<Lattice,Hamiltonian,StateSpace>(grid, ham, statespace, a, rsite, HasNbrs(), HasBnds());

			// compute interaction energy difference
			interactionenergydiff += ham.interactions[a]->J * (dot(svnew-svold, nbrhood));
		}
        
        
		// II. onsite energy part  -----
		double onsiteenergydiff = 0;

		// find Hamiltonian terms associated to "rsite"
		auto terms = get_terms<Lattice,Hamiltonian>(grid, ham, rsite, HasTrms());
		
		// loop over these terms
		for (OnSiteSizeType b=0; b<terms.size(); b++)
		{
			// select on-site term
			const int tidx = terms[b]; 
			
			// compute the difference
			auto diff = ham.onsite[tidx]->get(svnew) - ham.onsite[tidx]->get(svold);
			
			// multiply the constant
			onsiteenergydiff += dot(ham.onsite[tidx]->h, diff);
		}
        

		// III. flex term energy
		double flexenergydiff = 0;
		for (FlexSizeType c=0; c<ham.multisite.size(); c++)
		{
			auto nbrs = getflexnbrs<Lattice>(grid, c, rsite);
// 			auto diff = ham.multisite[c]->diff(rsite, svold, svnew, nbrs, statespace, grid);
// 			flexenergydiff += dot(ham.multisite[c]->k, diff);
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
		else if (rng.real() < std::exp(-beta*dE))
		{
			svold = svnew;
			retval = 1;
		}
		return retval;
	}
};


#endif
