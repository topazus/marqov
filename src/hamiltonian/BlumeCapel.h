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

#ifndef BLUMECAPEL_H
#define BLUMECAPEL_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "../hamparts.h"
#include "termcollection.h"


template <class StateVector, class RNG>
class BlumeCapel_Initializer
{
	public:
		BlumeCapel_Initializer()   {}
		BlumeCapel_Initializer(RNG& rn) : rng(rn) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 

			int state = retval[0];

			if (state == 0)
			{
				if (rng.real() < 0.5)  state = -1;
				else				state = +1;
			}
			else // +1/-1
			{
				if (rng.real() < 0.5) state *= -1;
				else state = 0;
			}

			retval[0] = state;

			return retval;
		};

	private:
		RNG& rng;
};


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class BlumeCapel
{
	public:

		//  ----  Parameters  ----

		double J, D;
		constexpr static int SymD = 1;
		const std::string name;


		
		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = BlumeCapel_Initializer<StateVector, RNG>;



		//  ----  Hamiltonian terms  ----

		std::array<Standard_Interaction<StateVector>*, 1>    interactions = {new Standard_Interaction<StateVector>(J)};
		std::array<Onsite_Quadratic<StateVector>*, 1>        onsite       = {new Onsite_Quadratic<StateVector>(D)};
		std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;
	
		BlumeCapel(double J, double D) : J(J), D(D), name("BlumeCapel"), observables(obs_m) {}



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



		//  ----  Wolff  ----

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0) const
		{
			if (sv1[0] == 0) return 0.0;
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}


		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0) const
		{
			sv[0] *= -static_cast<SpinType>(1.0);
		}

};
#endif
