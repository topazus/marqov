#ifndef CC_LATTICE_H
#define CC_LATTICE_H

#include "../geometry/points.h"	
#include "../geometry/distance.h"	
#include "helpers/constantcoordination2D.h"

/**
 * Constant Coordination Lattice
 * A class that can create a Constant Coordination lattice on a user defined
 * point cloud.
 * @tparam PointCloud The Point Cloud that is used by the user.
 */
template <class PointCloud>
class ConstantCoordinationLattice
{
	private:
		PointCloud cloud;
	
	public:
		int npoints, len, dim;
		std::vector<std::vector<int>> neighbours;
		ConstantCoordinationLattice(const int len, const int dim) : cloud(len*len, len, dim), 
																	npoints(pow(len,dim)), 
																	len(len), 
																	dim(dim)
		{
			if (dim > 2) throw std::invalid_argument("The CC lattice is currently not implemented for d="+std::to_string(dim));
			constant_coordination_lattice(cloud, neighbours);
		}

		// implement nbrs
		std::vector<int> nbrs(const int a, const int i) const
		{
			return neighbours[i];
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return cloud.crds(i);
		}

		std::size_t size() const {return npoints;}
};

#endif
