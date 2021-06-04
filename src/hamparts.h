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

#ifndef HAMILTONIANPARTS_H
#define HAMILTONIANPARTS_H

/**
 * A generic interface for a two-body interaction.
 * @tparam StateVector the type of the StateVector that the interaction
 *                     acts upon.
 */
template <class StateVector>
class Interaction
{
	public:
		double J;
		Interaction(double J) : J(J) {}
		virtual ~Interaction() {};
		virtual StateVector get(const StateVector& phi_i) = 0;
};




/**
 * A generic interface for an site interaction.
 * @tparam StateVector the type of the StateVector that the interaction
 *                     acts upon.
 * @tparam CouplingType
 */
template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
		CouplingType h;
		OnSite(CouplingType h) : h(h) {}
		virtual ~OnSite(){};
		virtual CouplingType get(const StateVector& phi) = 0;
};




/**
 * A generic interface for a flexible term.
 * @tparam StateSpace these terms my act upon the full state space, hence we
 *                    need its type.
 * @tparam StateVector the type of the StateVector that the interaction
 *                     acts upon.
 */
template <class StateSpace, class StateVector>
class FlexTerm
{
	public:
		double k;
		FlexTerm(double k) : k(k) {};
		virtual ~FlexTerm() {};

		/** Interface for the energy difference to be used in local update algorithms
		*
		* @param rsite site under consideration
		* @param svold state before the update
		* @param svnew proposed state after the update
		* @param nbrs neighbours of the site
		* @param the statespace
		*
		* @return energy difference (double)
		*/
		virtual double diff (const int rsite,
					const StateVector& svold,
					const StateVector& svnew,
					std::vector<int>& nbrs,
					StateSpace& s) = 0;

		template <class Grid>


		/** Computes the total energy of this Hamiltonian term on the lattice
		*
		* @param s the statespace
		* @param grid the lattice
		* @param c the sublattice / bond family
		*
		* @return energy difference (double)
		*/
		double energy(const StateSpace& s, const Grid& grid, int c) {return 0;}
		// TODO: remove this default implementation and instead check for existence at compile-time
};
#endif
