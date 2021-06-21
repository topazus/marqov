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

#ifndef BLUMECAPELBIPARTITE_H
#define BLUMECAPELBIPARTITE_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include <vector>
#include "util/hamparts.h"


// ------------------------------ OBSERVABLES ---------------------------

#include "util/observables.h"

// ------------------------------ INITIALIZER ---------------------------

#include "util/initializers.h"


// ------------------------------ HAMILTONIAN ---------------------------

template <class StateVector>
class BlumeCapelBipartite_onsite
{
	public:
        const double& h;
		BlumeCapelBipartite_onsite(const double& D) : h(D) {}
		double get (const StateVector& phi) {return dot(phi, phi);};
};


template <typename SpinType = int>
class BlumeCapelBipartite
{
	public:

		//  ----  Parameter  ----

		double J, DA, DB;
		static constexpr int SymD = 1;
		const std::string name;


		//  ----  Definitions  ----

		typedef std::array<SpinType, SymD> StateVector;
//		template <typename RNG>
//		using MetroInitializer = Spin1_Initializer<StateVector, RNG>;

		//  ----  Hamiltonian terms  ----

        std::array<Standard_Interaction<StateVector>*,1> interactions = {new Standard_Interaction<StateVector>(J)};
        std::array<BlumeCapelBipartite_onsite<StateVector>*, 2> onsite = {new BlumeCapelBipartite_onsite<StateVector>(DA), new BlumeCapelBipartite_onsite<StateVector>(DB)};

		BlumeCapelBipartite(double J, double DA, double DB) : J(J), DA(DA), DB(DB), name("BlumeCapelBipartite"), observables(obs_m)
		{
		}
		
		~BlumeCapelBipartite(){delete interactions[0]; delete onsite[0]; delete onsite[1];}


		//  ----  Observables ----

		Magnetization obs_m;
        std::tuple<Magnetization> observables;


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}

};

template <typename SpinType>
class Initializer<BlumeCapelBipartite<SpinType> > : public  Spin1_Initializer<typename BlumeCapelBipartite<SpinType>::StateVector > {};
#endif
