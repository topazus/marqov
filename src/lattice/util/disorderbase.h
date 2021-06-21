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

#ifndef GRID_H
#define GRID_H

#include <vector>
#include <random>
#include "points.h"
#include "distance.h"

/**
 * Disordered Grid Base Class.
 * provides default implementations for the most important interfaces of a disordered lattice
 */
class DisorderType
{
	protected:
		int npoints;
		DisorderType(){}
		DisorderType(int npoints) : npoints(npoints) {}

	public:
		std::vector<std::vector<int>> neighbours;

		std::vector<int> nbrs(const int a, const int i)
		{
			return neighbours[i];
		}

		std::size_t size() const {return npoints;}

		inline int identify(int i) {return 0;};
		inline std::vector<int> termselector(int sublattice) const {return {-1};}
};




#endif
