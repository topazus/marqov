/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2021, The MARQOV Project
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

#ifndef EMBEDDER_H
#define EMBEDDER_H
#include "svmath.h"

namespace MARQOV 
{
/** The default Embedder.
  * projects the model onto one of the SymD Cartesian dimensions
  *
  * @tparam Hamiltonian the type of the Hamiltonian
  * @tparam Lattice the type of the Lattice
  * @tparam StateSpace The type of the Statespace.
  *
  * @note This default Embedder will also work, i.e. lead to a code which compiles
  * and runs, however it might not be physically reasonable for models other than
  * standard O(N) models. Therefore, make sure to implement your own embedder when
  * using more complex models -- or disable cluster updates, in which case it will 
  * not be executed
  */
template <class Hamiltonian, class Lattice>
class Embedder
{
	// definitions
	typedef typename Hamiltonian::StateVector StateVector;
    typedef Space<typename Hamiltonian::StateVector, Lattice> StateSpace;
	constexpr static int SymD = Hamiltonian::SymD;

	private:
		const Hamiltonian& ham;
		const Lattice& lat;
		const StateSpace& statespace;

	public:

		int rdir; // encodes the random direction

		Embedder(const Hamiltonian& ham, const Lattice& lat, const StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};

		/** Set new embedding variable.
		  *
		  * Typically, this function is executed once before every cluster update. The  variable
		  * can be drawn randomly (for which case an RNG is provided), but of course can also follow
		  * some sequental scheme.
		  *
		  * @tparam RNG the type of the random number generator
		  * @param rng the random number generator
		  */
		template <class RNG>
		void draw(RNG& rng)
		{
			rdir = rng.integer(SymD);
		}


		/** Computes the Wolff coupling when attempting to add a spin to the cluster
		*
		* @param pos1 The position (index) of the current state vector (which is already in the cluster)
		* @param pos2 The position (index) of a neighbour being checked whether it will become part of the cluster
		*
		* @return The scalar Wolff coupling (a double)
		*/
		auto coupling(int pos1, int pos2) const
		{
			return statespace[pos1][rdir] * statespace[pos2][rdir];
		}


		/** Specifies how a spin flip in the embedded (reduced) model is performed 
		*
		* @param sv the spin to flipped
		*/
		void flip(StateVector& sv) const
		{
			sv[rdir] = -sv[rdir];
			normalize(sv);
		}
	};
}
#endif
