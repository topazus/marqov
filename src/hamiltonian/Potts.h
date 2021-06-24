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

#ifndef POTTS_H
#define POTTS_H
#include <array>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include <array>
#include "util/hamparts.h"
#include "util/observables.h"
#include "util/termcollection.h"


// ------------------------------ OBSERVABLES ---------------------------


// ------------------------------ HAMILTONIAN ---------------------------

/**
 * Ising Hamiltonian
 * This defines the Ising Hamiltonian. It only consists of a single part,
 * namely the interaction.
 * @tparam SpinType the type in which to store the binary magnetization values.
 */
template <int Q>
class Potts
{
	public:
		
		//  ----  Parameters  ----

		static constexpr int q = Q;
		double J;
		static constexpr int SymD = 1;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<int, SymD> StateVector;

		//  ----  Hamiltonian terms  ---- 

        std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};

		Potts(double J) : J(J), name("Ising")
		{}
		~Potts() {delete interactions[0];}


		//  ----  Observables ----

		Magnetization  obs_m;
        decltype(std::make_tuple(obs_m)) observables = {std::make_tuple(obs_m)};


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (decltype(grid.size()) i = 0; i < grid.size(); i++)
			{
				statespace[i][0] = rng.integer(Q);
			}
		}

};


template <int Q>
class Initializer<Potts<Q>>
{
	public:

		template <class RNGCache>
		static typename Potts<Q>::StateVector newsv(typename Potts<Q>::StateVector& svold, RNGCache& rng) 
		{
			typedef typename Potts<Q>::StateVector StateVector;
			StateVector retval(svold);
			retval[0] = rng.integer(Q);
			return retval;
		}
};


#endif
