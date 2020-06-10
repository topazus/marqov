#ifndef POINTS_H
#define POINTS_H


// implement simple class for "site" which only carries the coordinate


// Point cloud base class
class PointCloud
{
	public:
		int dim, npoints;
		PointCloud(){}
		PointCloud(int npoints, int dim) : npoints(npoints), dim(dim) {}

		std::vector<std::vector<double>> crds;

		std::vector<double> getcrds(const int i) const { return crds[i];	}

		double getx(const int i) const { return crds[i][0]; }
		double gety(const int i) const { return crds[i][1]; }

		std::size_t size() const {return npoints; }
};


// random Poissonian point cloud
class Poissonian : public PointCloud
{
	private:
		RND rng;

	public:
		Poissonian(int dim, int npoints) : PointCloud(npoints, dim), rng(0,1)
		{
			rng.seed(time(NULL)+std::random_device{}());

			for (int i=0; i<npoints; i++)
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
			npoints = len*len;
			
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

#endif
