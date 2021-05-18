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



template <class Hamiltonian>
class Embedder
{
	private:
		const Hamiltonian& ham; // why const?

	public:

		typedef typename Hamiltonian::StateVector StateVector;
		Embedder(const Hamiltonian& ham) : ham(ham) {};

		template <class RNG>
		void draw(RNG& rng)
		{
			
		}


		double coupling(StateVector& currentsv, StateVector& candidate)
		{
			return ham.J;
		}


		void flip(StateVector& sv)
		{
			sv[0] = -sv[0];
		}
};

#endif
