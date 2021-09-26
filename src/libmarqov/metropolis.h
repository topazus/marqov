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
#include "../hamiltonian/util/initializers.h"
#include <cassert>
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
		template <class StateSpace, class RNGType>
		static int move(const Hamiltonian& ham, 
		   				const Lattice& grid, 
						StateSpace& statespace, 
						RNGCache<RNGType>& rng, 
						double beta, 
						int rsite);
	};
    
//     double P4(double x)
//     {
//         return
//     }
    template <class RNGType>
    inline bool update_accepted(double dE, double beta, RNGCache<RNGType>& rng)
    {
        bool accept = false;
        if ( dE <= 0 )
		{
            accept = true;
		}
		else
        {
            // if not, accept with probability depending on Boltzmann weight
            //try to decide depending on the magnitude:
            double rngnum = rng.real();
            int rndexpo = 0;
            double bintemp = M_LOG2E * beta;
            double action2 = - bintemp*dE;
            auto rdmantissa = std::frexp(rngnum, &rndexpo);

//             std::cout<<rngnum<<" "<<rndexpo<<" "<<action2<<" "<<std::exp(-beta*dE)<<'\n';
//             std::cout<<rngnum<<" "<<std::exp(-beta*dE)<<" | "<<rndexpo<<" "<<std::trunc(action2)<<std::endl;
            if (rndexpo == lround(std::trunc(action2)))// if both numbers are of equal magnitude
            {//check harder
                if ( rngnum < std::exp(-beta*dE)) 
                {
                    accept = true;
                }
            }
            else // decide based on the magnitudes
                accept = (rndexpo < std::trunc(action2));
//         std::cout<<(accept?"true":"false")<<" "<<(rngnum < std::exp(-beta*dE)? "exp true": "exp false")<<std::endl;
//         assert(accept == (rngnum < std::exp(-beta*dE)));
        }
		return accept;
    }
	/**
	 * The actual Metropolis move attempt.
	 *
     * The behaviour of the basic local MARQOV step, how to generate a new state vector
     * from an old one, can be set by the initializers, @see initializers.h
	 * @tparam Hamiltonian the type of the Hamiltonian
	 * @tparam Lattice the type of the lattice
	 * @tparam StateSpace the type of the state space
	 * @tparam M the type of the initializer
	 * @tparam RNG the type of the random number generator
	 *
	 * @param ham the Hamiltonian
	 * @param grid the lattice
	 * @param statespace the statespace
	 * @param rng the random number generator
	 * @param beta inverse temperature
	 * @param rsite site to be considered for an update
	 *
	 * @return integer, encoding whether the update was accepted or rejected
	 */
    template <class Hamiltonian, class Lattice>
    template <class StateSpace/*, class M*/, class RNGType>
    int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
    											const Lattice& grid, 
												StateSpace& statespace, 
												RNGCache<RNGType>& rng, 
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
		StateVector svnew(Initializer<Hamiltonian>::newsv(svold, rng) );

		// I. interaction part
		double interactionenergydiff = compute_interactionenergydiff(grid, ham, statespace, svnew, svold, rsite, HasInt());
        
		// II. onsite energy part
		auto onsiteenergydiff = compute_onsiteenergydiff(grid, ham, svnew, svold, rsite, HasOns());

		// III. flex term energy
		auto flexenergydiff = compute_flexenergydiff(grid, ham, statespace, svnew, svold, rsite, HasFlex());


		// collect energy differences
		double dE = interactionenergydiff + onsiteenergydiff + flexenergydiff;
		int retval = 0;
		if (update_accepted(dE, beta, rng))
		{
			svold = svnew;
			retval = 1;
		}
		return retval;
	}
};


#endif
