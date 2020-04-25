#include "distance.h"

// Point cloud base class
class PointCloud
{
	public:
		int dim, size;
		PointCloud(){}
		PointCloud(int size, int dim) : size(size), dim(dim) {}

		std::vector<std::vector<double>> crds;

		std::vector<double> getcrds(const int i)
		{
			return crds[i];
		}
};


// random Poissonian point cloud
class Poissonian : public PointCloud
{
	private:
		RND rng;

	public:
		Poissonian(int dim, int size) : PointCloud(size, dim), rng(0,1)
		{
			rng.seed(time(NULL)+std::random_device{}());

			for (int i=0; i<size; i++)
			{
				std::vector<double> site;

				for (int j=0; j<dim; j++)
				{
					site.push_back(rng.d());
				}
				crds.push_back(site);
			}
		}
};


// regular square point cloud
// improve me: calculate coordinates on demand
class RegularSquare : public PointCloud
{
	public:
		int len;

		RegularSquare(int len) : len(len), PointCloud(len*len, 2)
		{
			dim = 2;
			size = len*len;
			
			const double delta = 1.0 / (double)(len);
			for (int i=0; i<len; i++)  
			{
				for (int j=0; j<len; j++)
				{
					const std::vector<double> site = {(i+0.5)*delta, (j+0.5)*delta};
					crds.push_back(site);
				}
			}
		}
};



// -----------------------------------------
// ---------- Disorder base class ----------

template <typename bond_type = int>
class DisorderType
{
	protected:
		int npoints;
		DisorderType(){}
		DisorderType(int npoints) : npoints(npoints) {}

	public:
		std::vector<std::vector<int>> nbrs;

		std::vector<int> getnbrs(const int i)
		{
			return nbrs[i];
		}

		std::size_t size() const {return npoints;}

};




// ------------------------------
// Disorder Example 1: SuperChaos
// ------------------------------

template <class PointCloud, typename bond_type>
class SuperChaos : public DisorderType<bond_type>
{
private:
	int symD;
	RND rng;

public:
	std::vector<std::vector<std::vector<bond_type>>> bnds;

	SuperChaos(const PointCloud& cloud) : DisorderType<bond_type>(cloud.size()), rng(0,1)
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



// -------------------------------------
// Disorder Example 2: Erdos-Renyj Graph
// -------------------------------------
class ErdosRenyi : public DisorderType<int>
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



// --------------------------------------------------------------------
// Disorder Example 3: Regular square lattice with random bond disorder
// --------------------------------------------------------------------
template <typename bond_type> 
class RegularRandomBond:  public DisorderType<bond_type>
{
	private:
		RND rng;
		int len, dim;
		RegularLattice lattice;
		RegularSquare cloud; // only for 2D!!!
		double p;
	
	public:
		std::vector<std::vector<std::vector<bond_type>>> bnds;
		
		RegularRandomBond(int dim, int len, double p) : dim(dim), len(len), p(p), rng(0,1), lattice(len,dim), cloud(len)
		{
			// improve me:
			// calculate point coordinates on demand

			// prepare random number generator
			rng.seed(time(NULL)+std::random_device{}());	


			// construct bonds
			bnds.resize(lattice.size());
			for (int i=0; i<lattice.size(); i++)
			{
				for (int j=0; j<lattice[i].size(); j++) // why does lattice[i].size even work?
				{
					bond_type bndval;

					if (rng.d() < p) bndval = 1;
					else             bndval = -1;

					bnds[i].push_back(std::vector<bond_type>{bndval});
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
		std::vector<int> getnbrs(const int alpha, const int i)
		{
			return lattice.getnbrs(alpha, i);
		}

		// override getbnds
		std::vector<bond_type> getbnds(const int alpha, const int i, const int j)
		{
			return bnds[i][j];
		}

		// implement getcrds
		std::vector<double> getcrds(const int i)
		{
			return cloud.getcrds(i);
		}

		std::size_t size() const {return pow(len,dim);}

};
