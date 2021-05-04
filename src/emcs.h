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

namespace MARQOV
{
template <class Grid, class Hamiltonian, template<class> class RefType>
double Core<Grid, Hamiltonian, RefType>::elementaryMCstep()
{
	const int SymD = std::tuple_size<StateVector>::value;

	// cluster updates
	marqovtime.switch_clock("cluster");
	double avgclustersize = 0;
	for (int j=0; j < mcfg.ncluster; j++)
	{
		const int rsite = rngcache.integer(this->grid.size());

		// Heisenberg; random direction
		const auto rdir = rnddir<RNGCache<RNGType>, typename StateVector::value_type, SymD>(rngcache);

		avgclustersize += wolffstep(rsite, rdir);
	}

	// Metropolis sweeps
	marqovtime.switch_clock("local");
	for (int j=0; j<mcfg.nsweeps; j++)
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
