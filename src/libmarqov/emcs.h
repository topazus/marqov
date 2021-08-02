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

#ifndef EMCS_H
#define EMCS_H
#include "rngcache.h"
#include "metropolis.h"
#include "wolff.h"
/** @file emcs.h
 * Defines the Elementary Monte Carlo Step (EMCS).
 * @see core.h , the main file, for docs.
 */

namespace MARQOV
{
    /** Entry point for overwriting the EMCS.
     * 
     * This class serves as an entry point for easily defining your own
     * specializations of the notion of a basic sweep of your Hamiltonians.
     * Further uses are e.g. adding your own moves
     * To that end it has the two prototypical template parameters:
     * @tparam Hamiltonian The Hamiltonian that the Wolff algo will use.
     * @tparam Lattice The Lattice, that the Wolff algo should use.
     */
    template <class Hamiltonian, class Lattice>
    struct EMCS
    {
        template <class RNGType, class StateSpace, class Timer>
        static inline double move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNGCache<RNGType>& rng, double beta, int ncluster, int nmetro, Timer&);
    };
    
    template <class Hamiltonian, class Lattice>
    template <class RNGType, class StateSpace, class Timer>
    double EMCS<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
    											const Lattice& grid, 
												StateSpace& statespace, 
												RNGCache<RNGType>& rngcache, 
												double beta, int ncluster, int nmetro, Timer& mrqvt)
    {
        // cluster updates
        mrqvt.switch_clock("cluster");
        double avgclustersize = 0;
        for (int j=0; j < ncluster; j++)
        {
            const int seed = rngcache.integer(grid.size());

            avgclustersize += Wolff<Hamiltonian, Lattice>::move(ham, grid, statespace, rngcache, beta, seed);
        }

        // Metropolis sweeps
        mrqvt.switch_clock("metrop");
		double avgacceptance;
        for (int j=0; j < nmetro; j++)
        {
            // loop sites
            for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
            {
                const int rsite = rngcache.integer(grid.size());
                int avgacceptance = Metropolis<Hamiltonian, Lattice>::move(ham,
                                                grid, 
                                                statespace, 
                                                rngcache, 
                                                beta, 
                                                rsite);
            }
        }

        return avgacceptance/grid.size()/double(nmetro);
//        return avgclustersize/ncluster;
    }
    
template <class Grid, class Hamiltonian, class MutexType, class RNGType, template<class> class RefType>
double Core<Grid, Hamiltonian, MutexType, RNGType, RefType>::elementaryMCstep()
{
    return EMCS<Hamiltonian, Grid>::move(this->ham,
											this->grid, 
											statespace, 
											this->rngcache, 
											this->beta, mcfg.ncluster, mcfg.nmetro, mrqvt);
	}
};
#endif
