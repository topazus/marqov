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

#ifndef TERMCOLLECTION_H
#define TERMCOLLECTION_H

#include "../hamparts.h"

template <class StateVector, typename CouplingType = double>
class Standard_Interaction : public Interaction<StateVector>
{
	public:
		// what is correct?
		// this:

//		const CouplingType& J;
//		Standard_Interaction(const CouplingType &J) : J(J) {}

		// or this:

		Standard_Interaction(const CouplingType &J) {this->J = J;}



		StateVector get (const StateVector& phi) {return phi;}

		// TODO: rename J to something more generic, like "constant"
		
};

template <class StateVector, typename CouplingType = double>
class Onsite_Quadratic : public OnSite<StateVector, CouplingType> 
{
	public:
		Onsite_Quadratic(CouplingType constant)
		{
			this->h = constant;
		}
		CouplingType get (const StateVector& phi) {return dot(phi,phi);}
};

#endif
