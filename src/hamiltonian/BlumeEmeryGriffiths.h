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
class BiquadraticExchangeInteraction : public FlexTerm<StateSpace, StateVector>
{
	public:

		BiquadraticExchangeInteraction(double k)
		{
			this->k = k;
		}
		~BiquadraticExchangeInteraction() {};


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
				neighbourhood = neighbourhood + dot(nbr,nbr);
			}

			auto svdiff = dot(svnew,svnew)-dot(svold,svold);
			return dot(svdiff, neighbourhood);
		}


		template <class Grid>
		double energy(const StateSpace& statespace, const Grid& grid, int c)
		{
			const int N = grid.size();
			double retval = 0;

			for (int idx=0; idx<N; idx++)
			{
				auto nbrs = grid.flexnbrs(c, idx);
				auto self = statespace[idx];
				
				for (std::size_t i=0; i<nbrs.size(); ++i)
				{
					auto idx = nbrs[i];
					auto nbr = statespace[idx];

					retval += dot(self,self) * dot(nbr,nbr);
				}
			}
			return 0.5*retval; // account for double counting
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

		
		BiquadraticExchangeInteraction<StateVector*,StateVector> biquadratic_exchange_int;

		std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};
		std::array<Onsite_Quadratic<StateVector>*, 1>     onsite       = {new Onsite_Quadratic<StateVector>(D)};
		std::vector<FlexTerm<StateVector*,StateVector>*>  multisite ;
	
		BlumeEmeryGriffiths(double J, double D, double K) : J(J), 
															D(D), 
															K(K), 
															name("BlumeCapel"), 
															observables(obs_m),
															biquadratic_exchange_int(K)
							{
								multisite.push_back(&biquadratic_exchange_int);
							}



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




// ------------------------------ SPECIALIZATIONS ---------------------------

namespace MARQOV
{
	/** Specialization of the Embedding class for the Blume Capel model */

	template <class SpinType, class Lattice>
	class Embedder<BlumeEmeryGriffiths<SpinType>,Lattice>
	{
		typedef BlumeEmeryGriffiths<SpinType> Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;

		private:
			const Hamiltonian& ham;
			const Lattice& lat;
			StateSpace& statespace;

		public:
			/** Constructs a BEG embedding object.
			 *
			 * @param ham The corresponding Hamiltonian
			 * @param lat The corresponding lattice
			 * @param statespace The statespace of the simulation
			 */
			Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			/** Set new embedding variable.
			  * @note For this specific model this is not necessary, as the embedding is fixed
			  */
			template <class RNG>
			void draw(RNG& rng) { } // nothing to draw 


			/** Computes the Wolff coupling when attempting to add a spin to the cluster
			*
			* @param pos1 The position (index) of the current state vector (which is already in the cluster)
			* @param pos2 The position (index) of a neighbour being checked whether it will become part of the cluster as well
			* @return The scalar Wolff coupling (a double)
			*
			* @note For this model, only clusters on the subset of +1/-1 spins are created. Spin-0 site are left untouched. The update will hence not be ergodic!
			*/
			double coupling(int pos1, int pos2) const
			{
				const SpinType s1 = statespace[pos1][0];
				const SpinType s2 = statespace[pos2][0];

				if (s1 == 0) return 0.0; // no cluster if seed state is zero
				if (s1 == s2) return 0.0;
				else return -1.0;
			}


			/** Specifies how a spin flip is performed */
			void flip(StateVector& sv) const {sv[0] *= -1;}

	};
}


#endif

