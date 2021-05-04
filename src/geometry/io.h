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

#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

// export coordinates and neighbour relations

template <class Grid>
void save_geometry(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

	const int dim = grid.dim;

	for (int i=0; i<grid.size(); i++)
	{
		auto nbrs = grid.nbrs(0,i);
		auto crds = grid.crds(i);

		os << nbrs.size() + crds.size();

		for (int j=0; j<crds.size(); j++) os << "\t" << crds[j];
		for (int j=0; j<nbrs.size(); j++) os << "\t" << nbrs[j];

		os << std::endl;
	}
}



// export coordinates and neighbour relations and bond strenghts (only scalars!)

template <class Grid>
void save_geometry_deluxe(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

	const int dim = grid.dim;

	for (int i=0; i<grid.size(); i++)
	{
		auto nbrs = grid.nbrs(0,i);
		auto crds = grid.crds(i);
		auto bnds = grid.getbnds(0,i);

		os << nbrs.size() + crds.size() + bnds.size();

		for (int j=0; j<crds.size(); j++) os << "\t" << crds[j];
		for (int j=0; j<nbrs.size(); j++) os << "\t" << nbrs[j];
		for (int j=0; j<bnds.size(); j++) os << "\t" << bnds[j];

		os << std::endl;
	}
}

#endif
