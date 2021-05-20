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

#ifndef EMBEDDER_H
#define EMBEDDER_H


namespace MARQOV 
{

// the default embedder
// projects the model onto one of the SymD Cartesian dimensions
template <class Hamiltonian, class Lattice>
class Embedder
{
	// definitions
	typedef typename Hamiltonian::StateVector StateVector;
	typedef StateVector* StateSpace;
	constexpr static int SymD = Hamiltonian::SymD;

	private:
		const Hamiltonian& ham;
		const Lattice& lat;
		const StateSpace& statespace;

	public:

		int rdir; // encodes the random direction

		Embedder(const Hamiltonian& ham, const Lattice& lat, const StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};

		template <class RNG>
		void draw(RNG& rng)
		{
			rdir = rng.integer(SymD);
		}


		auto coupling(int pos1, int pos2)
		{
			return statespace[pos1][rdir] * statespace[pos2][rdir];
		}


		void flip(StateVector& sv)
		{
			sv[rdir] = -sv[rdir];
			normalize(sv);
		}
	};
}
#endif
