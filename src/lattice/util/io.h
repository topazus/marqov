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

#ifndef IO_H
#define IO_H

#include <string>
#include <fstream>

/**
* export coordinates and neighbour relations
*
* @tparam Grid the type of the lattice
* @param grid the lattice
* @param path directory in which the geometry is to be stored
*/
template <class Grid>
void save_geometry(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

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



/**
* export coordinates, neighbour relations and bond strengths
*
* @tparam Grid the type of the lattice
* @param grid the lattice
* @param path directory in which the geometry is to be stored
*/

template <class Grid>
void save_geometry_deluxe(Grid& grid, const std::string path)
{
	std::ofstream os;
	os.open(path.c_str());
	os << std::fixed << std::setprecision(16);

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


int import_geometry_ncoords(std::string path)
{
	std::ifstream in(NULL);
	in.open(path.c_str());
	std::string row;
    std::getline(in, row);
	return std::stoi(row);
}


/**
* @brief Import coordinates and neighbour relations from CSV file
*
* @param path filename including full path
* @param grid here the coordinates will be stored
* @param nbrs here the neighbour relations will be stored
* @param ncoords number of coordinates per point
* @param idxoffs set to 1 if vertex indices in the CSV are counted from 1. Set 0 zero if they are counted from 0
*
* @return the number of vertices
*
* @note for format specifications see general documentation-
*/
int import_geometry(const std::string path, std::vector<std::vector<double>>& grid, std::vector<std::vector<int>>& nbrs, int ncoords, int idxoffs = 0)
{

	std::ifstream in(NULL);
	in.open(path.c_str());
	std::string row;

        std::getline(in, row);
        std::getline(in, row);

	// loop over lines
    while (!in.eof()) 
	{
        std::getline(in, row);
        if (in.bad() || in.fail()) 
		{
            break;
        }
		std::istringstream ss(row);
		std::string substr;

		int counter = 0;
		std::vector<double> g;
		std::vector<int> n;

		// loop over "words" in a line
		while(std::getline(ss, substr, '\t'))
		{
			// if coordinate, transform to double
			if (counter < ncoords) g.push_back(std::stod(substr));	
			// if vertex, transform to int
			else n.push_back(std::stoi(substr)-idxoffs);
			counter++;
    	}
		grid.push_back(g);
		nbrs.push_back(n);
	}
	return nbrs.size();
}





#endif
