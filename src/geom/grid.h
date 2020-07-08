#ifndef GRID_H
#define GRID_H


#include "points.h"
#include "distance.h"
#include "regular_lattice.h"
#include "constantcoordination2D.h"


// ---------------- Disordered Grid Base Class ----------------

class DisorderType
{
	protected:
		int npoints;
		DisorderType(){}
		DisorderType(int npoints) : npoints(npoints) {}

	public:
		std::vector<std::vector<int>> nbrs;

		std::vector<int> getnbrs(const int a, const int i)
		{
			return nbrs[i];
		}

		std::size_t size() const {return npoints;}

};



// ----------------- Constant Coordination --------------------

template <class PointCloud>
class ConstantCoordinationLattice
{
	private:
		PointCloud cloud;
	
	public:
		int counter = 0; // only for test, remove me later

		int npoints, len, dim;
		std::vector<std::vector<int>> nbrs;
		ConstantCoordinationLattice(const int len, const int dim) : cloud(len*len, len, dim), npoints(pow(len,dim)), len(len), dim(dim)
		{
			if (dim > 2) throw std::invalid_argument("The CC lattice is currently not implemented for d="+std::to_string(dim));
			constant_coordination_lattice(cloud, nbrs);
		}

		std::vector<int> getnbrs(const int a, const int i) const
		{
			return nbrs[i];
		}
		// implement getcrds
		std::vector<double> getcrds(const int i) const
		{
			return cloud.getcrds(i);
		}

		std::size_t size() const {return npoints;}
};




// ----------------- Regular Hypercubic --------------------

class RegularHypercubic
{
	private:
		RegularLattice lattice;

	public:
		int len, dim, npoints;


		RegularHypercubic(int len, int dim) : len(len), dim(dim), npoints(pow(len,dim)), lattice(len,dim) {};


		// override getnbrs
		std::vector<int> getnbrs(const int alpha, const int i) const
		{
			return lattice.getnbrs(alpha, i);
		}

		// implement getcrds
		std::vector<double> getcrds(const int i) const
		{
			return lattice.getcrds(i);
		}

		std::size_t size() const {return npoints;}
};





// ----------------- "Super Chaos" --------------------

template <class PointCloud, typename bond_type>
class SuperChaos : public DisorderType
{
private:
	int symD;
	RND rng;

public:
	std::vector<std::vector<std::vector<bond_type>>> bnds;

	SuperChaos(const PointCloud& cloud) : DisorderType(cloud.size()), rng(0,1)
	{
		// prepare neighbour array
		const int npoints = cloud.size;
		this->nbrs.resize(npoints);
		this->bnds.resize(npoints);

		// prepare random number generator
		rng.seed(time(NULL)+std::random_device{}());	
		rng.set_integer_range(npoints);


		// draw bonds
		for (int i=0; i<2*npoints; i++) // make me variable
		{
			const int j = rng.i();
			const int k = rng.i();

			auto jcrds = cloud.getcrds(j);
			auto kcrds = cloud.getcrds(k);

			if (i!=k && distancePBSQ_nD(jcrds,kcrds) < 0.1) // make me 1/L dependent
			// improve me: make sure bond does not already exist!
			{
				this->nbrs[j].push_back(k);
				this->nbrs[k].push_back(j);
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
					bond_type[n] = rng.d();
				}
				temp.push_back(subtemp);
			}
			bnds[k] = temp;
		}
				*/
	}

	// override getbnds
	std::vector<bond_type> getbnds(const int i)
	{
		return bnds[i];
	}
};





// ----------------- Erdos-Renyj Graph --------------------

class ErdosRenyi : public DisorderType
{
	private:
		int p;
		RND rng;
	
	public:
		ErdosRenyi(int npoints, double p) : DisorderType(npoints), p(p), rng(0,1)
		{
			// prepare neighbour array
			this->nbrs.resize(npoints);
	
			// seed random number generator
			rng.seed(time(NULL)+std::random_device{}());	
	
			// the actual implementation
			// of the Erdos-Renyj Graph
			for (int i=0; i<npoints; i++)
			{
				for (int j=i+1; j<npoints; j++)
				{
					if (rng.d() < p)
					{
						this->nbrs[i].push_back(j);
						this->nbrs[j].push_back(i);
					}
				}
			}
		}
};




// ----------------- Random Bond Disorder --------------------


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


template <class PDFType>
class RegularRandomBond :  public DisorderType
{
	private:
		RegularLattice lattice;
		PDFType PDF;
	
	public:
		int len, dim, npoints;

		using bond_type = decltype(PDF.draw());
		std::vector<std::vector<bond_type>> bnds;
		
		RegularRandomBond(int len, int dim) : dim(dim), len(len), lattice(len,dim), npoints(pow(len,dim))
		{
			// construct bonds
			bnds.resize(lattice.size());
			for (int i=0; i<lattice.size(); i++)
			{
				for (int j=0; j<lattice[i].size(); j++) // why does lattice[i].size even work?
				{
					bnds[i].push_back(PDF.draw());
				}
			}

			// "symmetrize" bonds
			for (int i=0; i<lattice.size(); i++)
			{
				auto lnbrs = lattice.getnbrs(1,i);

				for (int j=0; j<lattice[i].size(); j++)
				{
					// find i in bnds[lnbr] and replace its value by bnds[i][j]
					auto lnbr = lnbrs[j];
					auto nbrs_temp = lattice.getnbrs(1, lnbr);

					auto it  = std::find(nbrs_temp.begin(), nbrs_temp.end(), i);
					auto idx = std::distance(nbrs_temp.begin(), it);

					bnds[lnbr][idx] = bnds[i][j];
				}
			}
		}

		// override getnbrs
		std::vector<int> getnbrs(const int alpha, const int i) const
		{
			return lattice.getnbrs(alpha, i);
		}

		// implement getbnds
		std::vector<bond_type> getbnds(const int alpha, const int i) const
		{
			return bnds[i];
		}

		// implement getcrds
		std::vector<double> getcrds(const int i)
		{
			return lattice.getcrds(i);
		}

		std::size_t size() const {return npoints;}

};

#endif
