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

#ifndef BLUMEEMERYGRIFFITHS_H
#define BLUMEEMERYGRIFFITHS_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "../hamparts.h"
#include "termcollection.h"

// ------------------------------ OBSERVABLES ---------------------------

// ...



// ------------------------------ INITIALIZER ---------------------------

template <class StateVector, class RNG>
class BlumeEmeryGriffiths_Initializer
{
	public:
		BlumeEmeryGriffiths_Initializer(RNG& rn) : rng(rn) {}

		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			int state = retval[0];

			if (state == 0)
			{
				if (rng.real() < 0.5) state = -1;
				else state = +1;
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

template <class StateSpace, class StateVector>
class BiquadraticExchangeInteraction
//class BiquadraticExchangeInteraction : public FlexTerm<StateSpace, StateVector>
{

	double k = 1;

	public:

		BiquadraticExchangeInteraction(double k) : k(k) {}


		template <class Lattice>
		double diff (const int rsite,
					const StateVector& svold, 
					const StateVector& svnew, 
					std::vector<int>& nbrs, 
					StateSpace& s,
					Lattice& grid)
		{

			double neighbourhood = 0;

			for (std::size_t i=0; i<nbrs.size(); ++i)
			{
				auto idx = nbrs[i];
				auto nbr = s[idx];
				neighbourhood += mult(nbr,nbr);
			}

			auto svdiff = mult(svnew,svnew)-mult(svold,svold);

			return dot(svdiff, neighbourhood);
		}
};


template <typename SpinType = int>
class BlumeEmeryGriffiths
{
	public:

		//  ----  Parameters  ----

		double J, D, K;
		static constexpr int SymD = 1;
		const std::string name;


		
		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = BlumeEmeryGriffiths_Initializer<StateVector, RNG>;



		//  ----  Hamiltonian terms  ----

		std::array<Standard_Interaction<StateVector>*, 1>    interactions = {new Standard_Interaction<StateVector>(J)};
		std::array<Onsite_Quadratic<StateVector>*, 1>        onsite       = {new Onsite_Quadratic<StateVector>(D)};
		std::array<BiquadraticExchangeInteraction<StateVector*,StateVector>,0> multisite = {};//new BiquadraticExchangeInteraction<StateVector*,StateVector>(K)};
//		std::array<FlexTerm<StateVector*,  StateVector>*, 1> multisite = {new BiquadraticExchangeInteraction<StateVector*,StateVector>(K)};
	
		BlumeEmeryGriffiths(double J, double D, double K) : J(J), D(D), K(K), name("BlumeCapel"), observables(obs_m) {}



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





#endif
