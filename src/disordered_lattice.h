

class Poissonian
{
	public:
		int dim, size;
		RND rng;
		std::vector<std::vector<double>> crds;

		Poissonian(int dim, int size) : dim(dim), size(size), rng(0,1)
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


class Regular2D
{
	public:
		int dim = 2;
		int size;
		int len;
		std::vector<std::vector<double>> crds;

		Regular2D(int len) : len(len)
		{
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


template <class PointCloud>
class SomeRandomConnections
{

private:
	RND rng;

public:

	std::vector<std::vector<int>> nbrs;

	SomeRandomConnections() : rng(0,1) {}

	void construct(const PointCloud cloud)
	{
		const int npoints = cloud.size();
		nbrs.resize(npoints);
		rng.set_integer_range(npoints);
		for (int i=0; i<4*npoints; i++)
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

	std::vector<int> getnbrs(const int i)
	{
		return nbrs[i];
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
