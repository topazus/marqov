
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
class DisorderType
{
	public:
		std::vector<std::vector<int>> nbrs;

		DisorderType(){}

		std::vector<int> getnbrs(const int i)
		{
			return nbrs[i];
		}
};



// some dummy implementation of DisorderType
template <class PointCloud>
class SomeRandomConnections : public DisorderType
{
	private:
		RND rng;
	
	public:
		SomeRandomConnections(const PointCloud cloud) : rng(0,1)
		{
			rng.seed(time(NULL)+std::random_device{}());	

			const int npoints = cloud.size;
			nbrs.resize(npoints);
			rng.set_integer_range(npoints);
			for (int i=0; i<2*npoints; i++)
			{
				const int j = rng.i();
				const int k = rng.i();
	
				if (i!=k)
				{
					nbrs[j].push_back(k);
					nbrs[k].push_back(j);
				}
			}
		}
};



template <class PointCloud, class DisorderType>
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

//		std::vector<std::vector<int>> nbrs;
//		std::vector<std::vector<double>> crds;
//		std::vector<std::vector<std::vector<double>>> bnds;


		std::vector<int> getnbrs(int a, int i) {return 0;}
		std::size_t size() const {return 0;}
};
