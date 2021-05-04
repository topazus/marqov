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

template <class StateVector>
class Interaction
{
	public:
		double J;
		virtual StateVector get(const StateVector& phi_i) = 0;
		virtual ~Interaction() {};
};

template <class StateVector, typename CouplingType>
class OnSite
{
	public: 
		CouplingType h;
		virtual CouplingType get(const StateVector& phi) = 0;
		virtual ~OnSite(){};
};


template <class StateSpace, class StateVector>
class FlexTerm
{
	public:
		double k;
		virtual double get(const StateVector& sv, int svpos, StateSpace s) = 0;
		virtual ~FlexTerm() {};
		template <class Lattice>
		double diff (const int rsite,
					const StateVector& svold,
					const StateVector& svnew,
					std::vector<int>& nbrs,
					StateSpace& s,
					Lattice& grid) {return 0;}

};


#endif
