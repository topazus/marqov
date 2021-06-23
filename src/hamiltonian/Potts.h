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
template <typename SpinType = int>
class Potts
{
	public:
		
		//  ----  Parameters  ----

		int q;
		double J;
		static constexpr int SymD = 1;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;

		//  ----  Hamiltonian terms  ---- 

        std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};

		Potts(int q, double J) : q(q), J(J), name("Ising")
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
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}

};


template <typename SpinType>
class Initializer<Potts<SpinType>>
{
	public:

		template <class RNGCache>
		static typename Potts<SpinType>::StateVector newsv(typename Potts<SpinType>::StateVector& svold, RNGCache& rng) 
		{
			typedef typename Potts<SpinType>::StateVector StateVector;
			StateVector retval(svold);
			retval[0] = rng.integer(Potts<SpinType>::q);
			return retval;
		}
};


#endif
