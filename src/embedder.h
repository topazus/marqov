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

template <class Hamiltonian, class Lattice>
class Embedder
{
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;
		constexpr static int SymD = Hamiltonian::SymD;

	private:
		const Hamiltonian& ham; // why const?
		const Lattice& lat;
		StateSpace& statespace;

	public:

		int randomvar = 0;


		Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};

		template <class RNG>
		void draw(RNG& rng)
		{
			randomvar = rng.integer(SymD);
		}


		double coupling(int pos1, int pos2)
		{
			
			return statespace[pos1][randomvar] * statespace[pos2][randomvar];
		}


		void flip(StateVector& sv)
		{
			sv[randomvar] = -sv[randomvar];
			normalize(sv);
		}
	};
}
#endif
