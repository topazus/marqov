
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
class Regular2D : public PointCloud
{
	public:
		int len;

		Regular2D(int len) : len(len)
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



// Base class encoding the links
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



// some dummy implementation of DisorderType
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


//template <class RegularType, typename bond_type>
template <typename bond_type>
class RegularRandomBond:  public DisorderType<bond_type>
{
	private:
		RND rng;
		int len, dim;
		RegularLattice lattice;
	
	public:
		RegularRandomBond(int dim, int len) : dim(dim), len(len), rng(0,1)
		{
			RegularLattice templattice(len, dim);
			lattice = std::move(templattice);
		}

		// override getnbrs
		std::vector<int> getnbrs(const int i)
		{
			const int alpha = 1;
			return lattice.getnbrs(alpha, i);
		}

		// override getweights
		std::vector<std::vector<std::vector<bond_type>>> bnds;
		bond_type getweights(const int i)
		{
			return bnds[i];
		}
};

// this class is supposed to collect everything that is needed for
// for quenched disordered lattices; let's see if this works out ...
template <class PointCloud, class DisorderType, class BondWeights>
class DisorderedLattice
{
	private:
		int SymD;
		RND rng;
		PointCloud cloud;
		DisorderType disorder;


	public:
		DisorderedLattice(PointCloud cloud, int SymD) : cloud(cloud), symD(symD), rng(0,1) 
		{
			disorder.construct(cloud);

		}

		std::vector<int> getnbrs(const int i)
		{
			return disorder.getnbrs(i);
		}

		int symD;

//		std::vector<int> getnbrs(int a, int i) {return 0;}
//		std::size_t size() const {return 0;}
};
