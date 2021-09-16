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

#ifndef REPLICATE_H
#define REPLICATE_H
#include <vector>
#include <tuple>
#include "../core.h"

/** A helper to replicate parameter sets (to be used for disordered systems).
*
* @tparam Params the type of a parameter set
* @param params	the parameter set
* @param nreplicas number of replicas
*
* @return returns vector of replicated parameters
*/
template <class Params>
inline std::vector<Params> replicator(const std::vector<Params>& params, int nreplicas)
{
	std::vector<Params> newparams;

	int mcid = 0;
	for(std::size_t i=0; i<params.size(); ++i)
	{
		auto&& l  = std::get<0>(params[i]);
		auto mp   = std::get<1>(params[i]);
		auto&& hp = std::get<2>(params[i]);

		for (int j = 0; j < nreplicas; ++j)
		{
			auto mpr(mp);
			mpr.setrepid(j);
			mpr.setid(mcid++);
			newparams.emplace_back(l, mpr, hp);
		}
	}

	return newparams;
}



/** Set up a vector with the parameters for the simulation.
 * 
 * @tparam L the argument for the lattice information. If a Lattice or a reference to a full lattice is given for lp
 *           the universal reference decays to a plain reference. If a temporary object is given, move semantcs is invoked.
 * @tparam HArgs A template parameter pack with the parameters of the Hamiltonian
 * 
 * @param l_or_lp Depending on context this is usually either a reference to a lattice, or a tuple of arguments for a lattice.
 * @param mp The parameters for MARQOV.
 * @param hp The array with the hamiltonian parameters.
 * 
 * @return a vector of tuples with the full arguments for MARQOV::Core
 */
template <class L, class HArgs>
inline auto finalize_parameter(L&& l_or_lp, const MARQOV::Config& mp, const std::vector<HArgs>& hp)
{
	typedef std::tuple<L, MARQOV::Config, HArgs> RetType;
	std::vector<RetType> params;

	for(std::size_t i=0; i < hp.size(); ++i) 
	{
		params.emplace_back(l_or_lp, mp, hp[i]);
	}

	return params;
}

#endif
