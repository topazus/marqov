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

#ifndef WOLFF_H
#define WOLFF_H
#include <vector>
#include <cmath>
#include "embedder.h"

namespace MARQOV
{
    template <class RNGType>
    inline bool wolff_update_accepted(double dE, double beta, RNGCache<RNGType>& rng)
    {
        constexpr double twolog2e =  2.0 * 1.4426950408889634074; /* log_2 e */
        bool accept = true;
        if ( dE <= 0 )
        {
            accept = false;
        }
	else
        {// if not, accept with probability depending on Boltzmann weight
            double rngnum = rng.real();
            double action = -beta * dE;
            
            double action2 = twolog2e * action;// transform e^x to 2^x
            int64_t i64;
            std::memcpy(&i64, &rngnum, sizeof(rngnum));
            if ((( i64 >>52)  )-1022  != detail::trunctoi64(action2))// if both numbers have different magnitudes.
            {// decide based on the magnitudes. This is likely, hence first in the branch
                accept = (((i64>>52)  )-1022 > detail::trunctoi64(action2));
            }
            else
            {//both numbers are of same magnitude -> evaluate the remainder
                i64 = (i64 & ~0xFFF0000000000000) | 0x3FE0000000000000;

                double remainder = action2-std::trunc(action2);
                std::memcpy(&rngnum, &i64, sizeof(rngnum));
                if ( rngnum < detail::exp2app6(remainder) )
                {
                    accept = false;
                }
            }
        }
	return accept;
    }
    
    
    /**
     * This class serves as an entry point for easily defining your own
     * specializations of the Wolff Algorithm of your Hamiltonians.
     * To that end it has the two prototypical template parameters:
     * @tparam Hamiltonian The Hamiltonian that the Wolff algo will use.
     * @tparam Lattice The Lattice, that the Wolff algo should use.
     */


    template <class Hamiltonian, class Lattice>
    struct Wolff
    {
        template <class RNGType, class StateSpace>
        static inline int move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNGCache<RNGType>& rng, double beta, int rsite);
    };
    


    template <class Hamiltonian, class Lattice>
    template <class RNGType, class StateSpace>
    int Wolff<Hamiltonian, Lattice>::move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNGCache<RNGType>& rng, double beta, int rsite)
    {
		// set up embedder
		Embedder<Hamiltonian, Lattice> embd(ham, grid, statespace);
		embd.draw(rng,statespace[rsite]);

        // prepare stack
        typedef typename Hamiltonian::StateVector StateVector;
        std::vector<int> cstack(grid.size(), 0);

        // add cluster seed and flip it
        int q = 0;
        cstack[q] = rsite;
        int clustersize = 1;
		embd.flip(statespace[rsite]);

        // loop over stack as long as non-empty
        while (q>=0)
        {
            // extract last sv in stack
            const int currentidx = cstack[q];
            q--;
            
            // get its neighbours
			for (decltype(ham.interactions.size()) a = 0; a < ham.interactions.size(); a++)
			{
        		const double gcpl = ham.interactions[a]->J;
            	auto nbrs = grid.nbrs(a, currentidx);

            	// loop over neighbours
            	for (std::size_t i = 0; i < nbrs.size(); ++i)
            	{
            	    // extract corresponding sv
            	    const auto currentnbr = nbrs[i];
            	    StateVector& candidate = statespace[currentnbr];
            	    
					const auto lcpl = embd.coupling(currentidx, currentnbr);
					const double cpl = gcpl*lcpl;

					
            	    // test whether site is added to the cluster
                    if (wolff_update_accepted(cpl, beta, rng))
                    {
                        q++;
                        cstack[q] = currentnbr;
                        clustersize++;
                        embd.flip(candidate);
                    }
            	}
            }
        }
        return clustersize;
    }
};



#endif
