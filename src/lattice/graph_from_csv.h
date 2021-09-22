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



#ifndef GRAPHFROMCSV_H 
#define GRAPHFROMCSV_H

#include "util/io.h"

#include <vector>
#include <string>

/** Read minimal graph from csv file.
* @note for format specifications see general documentation
*/
class GraphFromCSV
{
    public:
    int npoints, ncoords;

    std::vector<std::vector<int>> neighbours;
    std::vector<std::vector<double>> grid;

    GraphFromCSV(std::string filename) 
    {
		ncoords = import_geometry_ncoords(filename);
        npoints = import_geometry(filename, grid, neighbours, ncoords);
    }
            
    std::vector<int> nbrs(const int a, const int i) const
    {
        return neighbours[i];
    }

    std::vector<double> crds(const int i) const
    {
        return grid[i];
    }

    std::size_t size() const {return npoints;}
};

/** This function provides the necessary overload for MARQOV::Core
 *  such that we can dump useful info into the HDF5 file
 * 
 * @param h5loc the group of the lattice
 * @param l the lattice that we are dumping info about.
 */
void writelat(H5::Group& h5loc, const GraphFromCSV& l)
{
	dumpscalartoH5(h5loc, "name", std::string("GraphFromCSV"));
	dumpscalartoH5(h5loc, "npoints", l.npoints);
}

#endif
