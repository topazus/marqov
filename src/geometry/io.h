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


void import_geometry( const int N, std::vector<std::vector<double> >& grid, std::vector< std::vector<int> >& n, const std::string path, int dim )
{
		grid.clear();
		n.clear();
		n.resize( N );
		
		std::ifstream is;
		is.open( path.c_str() );
		if (!is) std::cout << "\n#######\nError in import_geom: Failed to fetch grid!\n#######\n";
		else
		{
				double val;
				int lineLength;
				int neighbour;	
				int ct = 0;
				double x,y;
				while( is >> val )
				{
						lineLength = val;
                        std::vector<double> xt(dim);
                        for(int i = 0; i < dim; ++i)
                            is >> xt[i];
                        grid.push_back(xt);
						for (int i=0; i<lineLength-2; i++)
						{
								is >> neighbour;
								n[ct].push_back( neighbour );
						}
						ct++;
						if (ct>=N) break;
				}
		}
}


#include <sstream>

std::vector<std::vector<std::string>> readfile(std::istream &in) {
    std::vector<std::vector<std::string>> table;
    std::string row;
    while (!in.eof()) 
	{
        std::getline(in, row);
        if (in.bad() || in.fail()) 
		{
            break;
        }
		std::istringstream ss(row);
		std::string substr;
		while(std::getline(ss, substr, '\t'))
		{
			cout << substr << endl;
    	}
	}
    return table;
}




#endif
