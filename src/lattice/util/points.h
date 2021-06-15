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


#ifndef POINTS_H
#define POINTS_H

#include <vector>
#include <random>

/**
 * Point cloud base class.
 * Provides a way to store coordinates as well as interaces to access them.
 */
class PointCloud
{
	protected:
		std::vector<std::vector<double>> pccrds; ///< this is where the actual coordinates are stored

	public:
		int dim; ///< dimensionality
		int npoints; ///< number of points
		int len; ///< linear length

		/**
		* Point cloud base class constructor
		*/
		PointCloud(){}
		PointCloud(int npoints, int len, int dim) : dim(dim), npoints(npoints), len(len) {}

		/**
		* Getter function for the coordinates of a site
		*
		* @param i index of site
		*
		* @return std::vector<double> containing coordinates of site with index i
		*/
		const std::vector<double>& crds(const int i) const { return pccrds[i];}

		double getx(const int i) const { return pccrds[i][0]; }
		double gety(const int i) const { return pccrds[i][1]; }

		std::size_t size() const {return npoints; }
};


/**
 * Random Poissonian point cloud
 */
class Poissonian : public PointCloud
{
	public:
		/** Constructor of the Poissonian point cloud
		*
		* @param npoints number of points
		* @param len linear length of the hypercube containing the points
		* @param dim dimensionality of the point cloud
		*/
		Poissonian(int npoints, int len, int dim) : PointCloud(npoints, len, dim)
		{
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, 1.0);
			for (int i=0; i<npoints; i++)
			{
				std::vector<double> site;
				for (int j=0; j<dim; j++) {site.push_back(dis(gen));}
				pccrds.push_back(site);
			}
		}
};



/**
 * Regular square point cloud
 * FIXME: calculate coordinates on demand
 */
class RegularSquare : public PointCloud
{
	public:
		RegularSquare(int len) : PointCloud(len*len, len, 2)
		{
			dim = 2;
			npoints = len*len;
			
			const double delta = 1.0 / (double)(len);
			for (int i=0; i<len; i++)  
			{
				for (int j=0; j<len; j++)
				{
					const std::vector<double> site = {(i+0.5)*delta, (j+0.5)*delta};
					pccrds.push_back(site);
				}
			}
		}
};

#endif
