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
    /** Determine if the Wolff would accept the move with the proposed energy change.
     * 
     * @param dE The energy change
     * @param beta The inverse temperature
     * @param rng the Random numer generator
     * @return true if we extend the cluster, else false.
     */
    template <class RNGType>
    inline bool wolff_update_accepted(double dE, double beta, RNGCache<RNGType>& rng)
    {
        return !update_accepted<2, RNGType>(dE, beta, rng);
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
