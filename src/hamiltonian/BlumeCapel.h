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


// ------------------------------ OBSERVABLES ---------------------------

#include "../observables.h"


// ------------------------------ INITIALIZER ---------------------------

#include "initializers.h"


// ------------------------------ HAMILTONIAN ---------------------------

/**
* The Blume-Capel Hamiltonian
*
* @tparam SpinType the type of the state vector
*/
template <typename SpinType = int>
class BlumeCapel
{
	public:

		//  ----  Parameters  ----

		double J, D;
		static constexpr int SymD = 1;
		const std::string name;


		
		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = Spin1_Initializer<StateVector, RNG>;



		//  ----  Hamiltonian terms  ----

		std::array<Standard_Interaction<StateVector>*, 1>    interactions = {new Standard_Interaction<StateVector>(J)};
		std::array<Onsite_Quadratic<StateVector>*, 1>        onsite       = {new Onsite_Quadratic<StateVector>(D)};

		/** Constructor of the Blume-Capel model
		*
		* @param J standard interaction parameter
		* @param D zero-field splitting parameter
		*/
		BlumeCapel(double J, double D) : J(J), D(D), name("BlumeCapel"), observables(obs_m) {}



		//  ----  Observables ----

		Magnetization obs_m;
        std::tuple<Magnetization> observables;

};




// ------------------------------ SPECIALIZATIONS ---------------------------

namespace MARQOV
{

	/** Specialization of the Embedding class for the Blume Capel model
	*
	* @tparam SpinType the type of the spin
	* @tparam Lattice the type of the lattice
	*/
	template <class SpinType, class Lattice>
	class Embedder<BlumeCapel<SpinType>,Lattice>
	{
		typedef BlumeCapel<SpinType> Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef MARQOV::Space<StateVector, Lattice> StateSpace;

		private:
			const Hamiltonian& ham;
			const Lattice& lat;
			StateSpace& statespace;

		public:
			/** Constructs a Blume Capel embedding object.
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
