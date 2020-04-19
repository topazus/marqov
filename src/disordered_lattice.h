
// base class for point clouds which is nothing
// but a set of coordinate vectors
class PointCloud
{
	public:
		int dim, size;
		PointCloud(){}

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
		Poissonian(int dim, int size) : rng(0,1)
		{
			dim = dim;
			size = size;

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


// coordinates of a regular square lattice
// improve me: calculate coordinates on demand
class RegularSquare : public PointCloud
{
	public:
		int len;

		RegularSquare(int len) : len(len)
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



// Base class encoding the links and their weights
template <typename bond_type = int>
class DisorderType
{
	public:
		std::vector<std::vector<int>> nbrs;

		DisorderType(){}

		std::vector<int> getnbrs(const int i)
		{
			return nbrs[i];
		}
		
		inline bond_type getweights(const int i)
		{
			return 1;
		}

};



// Disorder Class implementation, example 1:
// get a point cloud and draw some random connections 
template <class PointCloud, typename bond_type>
class SomeRandomConnections : public DisorderType<bond_type>
{
	private:
		RND rng;
	
	public:
		SomeRandomConnections(const PointCloud& cloud) : rng(0,1)
		{

			// prepare neighbour array
			const int npoints = cloud.size;
			this->nbrs.resize(npoints);

			// prepare random number generator
			rng.seed(time(NULL)+std::random_device{}());	
			rng.set_integer_range(npoints);

			// the actual implementation of the geometry
			for (int i=0; i<2*npoints; i++)
			{
				const int j = rng.i();
				const int k = rng.i();
	
				if (i!=k)
				{
					this->nbrs[j].push_back(k);
					this->nbrs[k].push_back(j);
				}
			}
		}

		// override getweights
		std::vector<std::vector<std::vector<bond_type>>> bnds;
		bond_type getweights(const int i)
		{
			return bnds[i];
		}
};


// Disorder Class implementation, example 2:
// regular square lattice with random bond disorder
template <typename bond_type> 
class RegularRandomBond:  public DisorderType<bond_type>
{
	private:
		RND rng;
		int len, dim;
		RegularLattice lattice;
		PointCloud cloud;
	
	public:
		std::vector<std::vector<bond_type>> bnds;

		RegularRandomBond(int dim, int len) : dim(dim), len(len), rng(0.9, 1.1)
		{
			// improve me
			// 1. unelegant move
			// 2. calculate point coordinates on demand

			RegularLattice templattice(len, dim);
			lattice = std::move(templattice);

			RegularSquare tempcloud(len);
			cloud = std::move(tempcloud);

			for (int i=0; i<lattice.size(); i++)
			{
				std::vector<bond_type> bnd;
				for (int j=0; j<lattice[i].size(); j++)
				{
					bnd.push_back(rng.d());
				}
				bnds.push_back(bnd);
			}
		}

		// override getnbrs
		std::vector<int> getnbrs(const int i)
		{
			const int alpha = 1;
			return lattice.getnbrs(alpha, i);
		}

		// override getweights
		std::vector<bond_type> getweights(const int i)
		{
			return bnds[i];
		}

		// implement getcrds
		std::vector<double> getcrds(const int i)
		{
			return cloud.getcrds(i);
		}
};

