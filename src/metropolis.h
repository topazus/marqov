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

		// definitions
		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename HasInteractions<Hamiltonian>::type HasInt;
		typedef typename HasFlexTerms<Hamiltonian>::type HasFlex;
		typedef typename HasOnsite<Hamiltonian>::type HasOns;
		
		// old state vector at index rsite
		StateVector& svold = statespace[rsite];
		// propose new configuration
		StateVector svnew = metro.newsv(svold);
		        

		// I. interaction part
		double interactionenergydiff = compute_interactionenergydiff(grid, ham, statespace, svnew, svold, rsite, HasInt());
        
		// II. onsite energy part
		auto onsiteenergydiff = compute_onsiteenergydiff(grid, ham, svnew, svold, rsite, HasOns());

		// III. flex term energy
		auto flexenergydiff = compute_flexenergydiff(grid, ham, statespace, svnew, svold, rsite, HasFlex());


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
