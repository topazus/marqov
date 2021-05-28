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
#include "metropolis.h"
#include "wolff.h"
// Defines the Elementary Monte Carlo Step (EMCS)
// See main file, core.h, for docs.

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


// //     template <class Hamiltonian, class Lattice>
// //     struct EMCS
// //     {
// //         template <class RNG, class M, class StateSpace>
// //         static inline double move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta);
// //     };
// //     
// //     template <class Hamiltonian, class Lattice>
// //     template <class StateSpace, class M, class RNGType>
// //     double EMCS<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
// //     											const Lattice& grid, 
// // 												StateSpace& statespace, 
// // 												M& metro, RNGCache<RNGType>& rngcache, 
// // 												double beta)
// //     {
// //         
// //         
// //         
// //         
// //     // cluster updates
// // 	mrqvt.switch_clock("cluster");
// // 	double avgclustersize = 0;
// // 	for (int j=0; j < mcfg.ncluster; j++)
// // 	{
// // 		const int seed = rngcache.integer(this->grid.size());
// // 
// // 		avgclustersize += 
// // 		Wolff<Hamiltonian, Grid>::move(ham, grid, statespace, rngcache, beta, seed);
// // 		//wolffstep(seed);
// // 	}
// // 
// // 
// // 	// Metropolis sweeps
// // 	mrqvt.switch_clock("metrop");
// // 	for (int j=0; j<mcfg.nmetro; j++)
// // 	{
// // 		// loop sites
// // 		for(decltype(this->grid.size()) i = 0; i < this->grid.size(); ++i)
// // 		{
// // 			const int rsite = rngcache.integer(this->grid.size());
// //             Metropolis<Hamiltonian, Grid>::move(ham,
// // 											grid, 
// // 											statespace, 
// // 											metro, 
// // 											rngcache, 
// // 											beta, 
// // 											rsite);
// //             
// //             
// //             
// // // 			metropolisstep(rsite);
// // 		}
// // 	}
// // 
// // 	return avgclustersize/mcfg.ncluster;
// //         
// //         
// //         
// //     }
    
template <class Grid, class Hamiltonian, template<class> class RefType>
double Core<Grid, Hamiltonian, RefType>::elementaryMCstep()
{
	// cluster updates
	mrqvt.switch_clock("cluster");
	double avgclustersize = 0;
	for (int j=0; j < mcfg.ncluster; j++)
	{
		const int seed = rngcache.integer(this->grid.size());

		avgclustersize += wolffstep(seed);
	}


	// Metropolis sweeps
	mrqvt.switch_clock("metrop");
	for (int j=0; j<mcfg.nmetro; j++)
	{
		// loop sites
		for(decltype(this->grid.size()) i = 0; i < this->grid.size(); ++i)
		{
			const int rsite = rngcache.integer(this->grid.size());
			metropolisstep(rsite);
		}
	}

	return avgclustersize/mcfg.ncluster;
}
};
#endif
