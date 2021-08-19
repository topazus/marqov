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
    int npoints;

    std::vector<std::vector<int>> neighbours;

    GraphFromCSV(std::string filename) 
    {
        npoints = import_geometry(filename, neighbours);
    }
            
    std::vector<int> nbrs(const int a, const int i) const
    {
        return neighbours[i];
    }

    std::size_t size() const {return npoints;}
};


/* Read graph with coordinates from csv file.
* @note for format specifications see general documentation
*/
class GraphFromCSVwithCoords
{
    public:
    int npoints;

    std::vector<std::vector<int>> neighbours;
    std::vector<std::vector<double>> grid;

    GraphFromCSVwithCoords(std::string filename, int ncoords) 
    {
        npoints = import_geometry(filename, grid, neighbours, ncoords, false);
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
#endif
