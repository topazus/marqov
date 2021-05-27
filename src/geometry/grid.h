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
#include "regular_lattice.h"
#include "constantcoordination2D.h"

/**
 * Disordered Grid Base Class
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
		inline std::vector<int> termselector(int sublattice){return {-1};}
};

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
		int counter = 0; // only for test, remove me later

		int npoints, len, dim;
		std::vector<std::vector<int>> neighbours;
		ConstantCoordinationLattice(const int len, const int dim) : cloud(len*len, len, dim), npoints(pow(len,dim)), len(len), dim(dim)
		{
			if (dim > 2) throw std::invalid_argument("The CC lattice is currently not implemented for d="+std::to_string(dim));
			constant_coordination_lattice(cloud, neighbours);
		}

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

/**
 * The Regular Hypercubic lattice class
 * This class provides routines for encapsulating all neighbour
 * and coordinate relations for n-dimensional hypercubic lattices.
 */
class RegularHypercubic
{
	private:
		RegularLattice lattice;

	public:
		int len, dim, npoints;


		RegularHypercubic(int len, int dim) : lattice(len, dim), len(len), dim(dim), npoints(pow(len, dim)) {}

		// override nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		// only for debug, delete me later
		std::vector<double> bnds(const int alpha, const int i) const
		{
			return {1,1,1,1};
		}
		std::size_t size() const {return npoints;}
};

/**
 * Simple Bipartite Lattice
 */
class SimpleBipartite
{
	private:
		RegularLattice lattice;

	public:
		int len, dim, npoints;

		SimpleBipartite(int len, int dim) : lattice(len,dim), len(len), dim(dim), npoints(pow(len,dim))
		{
			if (len%2 != 0) cout << "ERROR: linear lattice size must be even!" << endl;
		}


		inline int identify(int i) // is this correct?
		{
			auto index = IndexOf(i, dim, len);
			
			int quersumme = 0;
			for (decltype(index.size()) j = 0; j < index.size(); j++) quersumme += index[j];

			if (quersumme%2 == 0) return 0;
			else return 1;
		}


		inline std::vector<int> termselector(int rsite)
		{
			return {this->identify(rsite)};
		}


		// override nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		std::size_t size() const {return npoints;}
};

/**
 * Super Chaos Lattice
 * @tparam PointCloud the point cloud that we use
 * @tparam bond_type 
 */
template <class PointCloud, typename bond_type>
class SuperChaos : public DisorderType
{
private:
	int symD;
public:
	std::vector<std::vector<std::vector<bond_type>>> bonds;

	SuperChaos(const PointCloud& cloud) : DisorderType(cloud.size())
	{
		// prepare neighbour array
		const int npoints = cloud.size;
		this->neighbours.resize(npoints);
		this->bonds.resize(npoints);
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> disreal(0.0, 1.0);
        std::uniform_int_distribution<> disint(0, npoints-1);

		// draw bonds
		for (int i=0; i<2*npoints; i++) // make me variable
		{
			const int j = disint(gen);
			const int k = disint(gen);

			auto jcoordinates = cloud.crds(j);
			auto kcoordinates = cloud.crds(k);

			if (i!=k && distancePBSQ_nD(jcoordinates,kcoordinates) < 0.1) // make me 1/L dependent
			// improve me: make sure bond does not already exist!
			{
				this->neighbours[j].push_back(k);
				this->neighbours[k].push_back(j);
			}
		}

		/* under construction
		
		// compute weights
		for (int k=0; k<npoints; k++)
		{
			std::vector<bond_type> temp;
			for (int m=0; m<this->nbrs[k].size(); m++)
			{
				bond_type subtemp;
				for (int n=0; n<sizeof(bond_type); n++)
				{
					bond_type[n] =disreal(gen);
				}
				temp.push_back(subtemp);
			}
			bonds[k] = temp;
		}
				*/
	}

	// override bnds
	std::vector<bond_type> bnds(const int i)
	{
		return bonds[i];
	}
};

/**
 * Erdos-Renyj Graph
 */
class ErdosRenyi : public DisorderType
{
	private:
		int p;
	public:
		ErdosRenyi(int npoints, double p) : DisorderType(npoints), p(p)
		{
            std::random_device rd;  //Will be used to obtain a seed for the random number engine
            std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
            std::uniform_real_distribution<> dis(0.0, 1.0);
			// prepare neighbour array
			this->neighbours.resize(npoints);
	
			// the actual implementation
			// of the Erdos-Renyj Graph
			for (int i=0; i<npoints; i++)
			{
				for (int j=i+1; j<npoints; j++)
				{
					if (dis(gen) < this->p)
					{
						this->neighbours[i].push_back(j);
						this->neighbours[j].push_back(i);
					}
				}
			}
		}
};

/**
 * A helper class for getting random numbers from a gaussian.
 */
class GaussianPDF
{
	private:
		std::mt19937 gen;
		std::normal_distribution<double> d;

	public:
		GaussianPDF() : gen(std::random_device{}()), d(std::normal_distribution<>{0,1}) {} 
		GaussianPDF(double mu, double sigma) : gen(std::random_device{}()), d(std::normal_distribution<double>{mu,sigma}) {} 

		double draw() {return(d(gen));}
};

/**
 * A helper class for getting random numbers from a bimodal distribution.
 */
class BimodalPDF
{
	private:
		std::mt19937 gen;
		std::discrete_distribution<int> d;
	public:
		BimodalPDF() : gen(std::random_device{}()), d(std::discrete_distribution<int>{1,0,1}) {}
		BimodalPDF(double weight) : gen(std::random_device{}()), d(std::discrete_distribution<int>{weight,0,1-weight}) {}

		int draw() {return(d(gen)-1);}
};

/**
 * A lattice with bond disorder
 * @tparam PDFType The type of disorder distribution to use.
 */
template <class PDFType>
class RegularRandomBond :  public DisorderType
{
	private:
		RegularLattice lattice;
		PDFType PDF;
	
	public:
		int len, dim, npoints;

		using bond_type = decltype(PDF.draw());
		std::vector<std::vector<bond_type>> bonds;
		
		RegularRandomBond(int len, int dim) : lattice(len,dim), len(len), dim(dim), npoints(pow(len,dim))
		{
			// construct bonds
			bonds.resize(lattice.size());
			for (decltype(lattice.size()) i=0; i<lattice.size(); i++)
			{
				for (decltype(lattice[i].size()) j = 0; j < lattice[i].size(); j++) // why does lattice[i].size even work?
				{
					bonds[i].push_back(PDF.draw());
				}
			}

			// "symmetrize" bonds
			for (decltype(lattice.size()) i = 0; i < lattice.size(); i++)
			{
				auto lnbrs = lattice.nbrs(1,i);

				for (decltype(lattice[i].size()) j = 0; j < lattice[i].size(); j++)
				{
					// find i in bonds[lnbr] and replace its value by bonds[i][j]
					auto lnbr = lnbrs[j];
					auto nbrs_temp = lattice.nbrs(1, lnbr);

					auto it  = std::find(nbrs_temp.begin(), nbrs_temp.end(), i);
					auto idx = std::distance(nbrs_temp.begin(), it);

					bonds[lnbr][idx] = bonds[i][j];
				}
			}
		}

		// override nbrs
		std::vector<int> nbrs(const int alpha, const int i) const
		{
			return lattice.nbrs(alpha, i);
		}

		// implement bnds
		std::vector<bond_type> bnds(const int alpha, const int i) const
		{
			return bonds[i];
		}

		// implement crds
		std::vector<double> crds(const int i) const
		{
			return lattice.crds(i);
		}

		std::size_t size() const {return npoints;}

};

#endif
