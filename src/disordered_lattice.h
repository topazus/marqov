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
		
		inline std::vector<bond_type> getbnds(const int i)
		{
			return 1;
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
	std::vector<std::vector<bond_type>> bnds;

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
			{
				this->nbrs[j].push_back(k);
				this->nbrs[k].push_back(j);
			}
		}

		
		// compute weights
		for (int k=0; k<npoints; k++)
		{
			std::vector<bond_type> temp;
			for (int m=0; m<this->nbrs[k].size(); m++)
			{
				bond_type subtemp;
				/* under construction
				for (int n=0; n<sizeof(bond_type); n++)
				{
					bond_type[n] = rng.d();
				}
				*/
				temp.push_back(subtemp);
			}
			bnds[k] = temp;
		}
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
		PointCloud cloud;
		double p;
	
	public:
		std::vector<std::vector<bond_type>> bnds;

		RegularRandomBond(int dim, int len, double p) : dim(dim), len(len), p(p), rng(0,1)
		{
			// improve me
			// 1. get rid of unelegant "move"
			// 2. calculate point coordinates on demand

			RegularLattice templattice(len, dim);
			lattice = std::move(templattice);

			RegularSquare tempcloud(len);
			cloud = std::move(tempcloud);

			// prepare random number generator
			rng.seed(time(NULL)+std::random_device{}());	

			// construct bonds
			for (int i=0; i<lattice.size(); i++)
			{
				std::vector<bond_type> bnd;
				for (int j=0; j<lattice[i].size(); j++)
				{
					if (rng.d() < p) bnd.push_back(1.0);
					else             bnd.push_back(-1.0);
				}
				bnds.push_back(bnd);
			}
		}

		// override getnbrs
		std::vector<int> getnbrs(const int alpha, const int i)
		{
			return lattice.getnbrs(alpha, i);
		}

		// override getbnds
		std::vector<bond_type> getbnds(const int alpha, const int i)
		{
			return bnds[i];
		}

		// implement getcrds
		std::vector<double> getcrds(const int i)
		{
			return cloud.getcrds(i);
		}

		std::size_t size() const {return pow(len,dim);}

};

