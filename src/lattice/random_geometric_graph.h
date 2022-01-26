#ifndef RGG_LATTICE_H
#define RGG_LATTICE_H

#include "util/points.h"	
#include "util/distance.h"	
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>


void construct_RGG_3D(const double R, const PointCloud& cloud, std::vector<std::vector<int>>& neighbours)
{
	const double Rsq = R*R;
	double upper[3]={0,0,0};
	double lower[3]={1,1,1};	
	
	//parametrization of cube surfaces
	double tcds[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};

	//parametrization of cube edges
	double tcde[12][6]={{0,0,1, 0,0,0},{0,0,1, 1,0,0},{0,0,1, 0,1,0},{0,0,1, 1,1,0},
				     {0,1,0, 0,0,0},{0,1,0, 1,0,0},{0,1,0, 0,0,1},{0,1,0, 1,0,1},
				     {1,0,0, 0,0,0},{1,0,0, 0,1,0},{1,0,0, 0,0,1},{1,0,0, 0,1,1}};

	//parametrization of cube corners
	double tcdc[8][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};	

	int number_of_nn;
	int size = cloud.size();
	// arrays to collect the search results
	std::vector<std::vector<int>> pointIdxRadiusSearch;
	std::vector<std::vector<float>> pointRadiusSquaredDistance;


	// construct PCL point cloud
	pcl::PointCloud<pcl::PointXYZ>::Ptr pclcloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	pcl::PointXYZ searchPoint;

	// transform input cloud into PCL format
	pclcloud->width = size;
	pclcloud->height= 1;
	pclcloud->points.resize(pclcloud->width * pclcloud->height);		
	for(int i=0; i<pclcloud->points.size(); i++)
	{
		auto crds = cloud.crds(i);
		pclcloud->points[i].x = crds[0];
		pclcloud->points[i].y = crds[1];
		pclcloud->points[i].z = crds[2];
	}

	// set up k-d tree
	kdtree.setInputCloud(pclcloud);

	// loop through points
	for(int i=0; i< pclcloud->points.size(); i++)
	{
		auto crds = cloud.crds(i);
		searchPoint.x = crds[0];
		searchPoint.y = crds[1];
		searchPoint.z = crds[2];
		
		const auto spx = searchPoint.x;
		const auto spy = searchPoint.y;
		const auto spz = searchPoint.z;

		pointIdxRadiusSearch.push_back(std::vector<int>());
		pointRadiusSquaredDistance.push_back(std::vector<float>());
		int pt=0;
		number_of_nn = kdtree.radiusSearch(searchPoint, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);		
		int surface_cuts=0;
	
		// surfaces
		for(int j=0; j<6; j++)
		{
			bool c1x = (searchPoint.x + tcds[j][0] < upper[0]);
			bool c2x = (searchPoint.x + tcds[j][0] > lower[0]);
			bool c1y = (searchPoint.y + tcds[j][1] < upper[1]);
			bool c2y = (searchPoint.y + tcds[j][1] > lower[1]);
			bool c1z = (searchPoint.z + tcds[j][2] < upper[2]);
			bool c2z = (searchPoint.z + tcds[j][2] > lower[2]);

			if( !( c1x and c2x and c1y and c2y and c1z and c2z) ) 
			{
			 	pointIdxRadiusSearch.push_back(std::vector<int>());
			 	pointRadiusSquaredDistance.push_back(std::vector<float>());
			 	pt++;
			 	pcl::PointXYZ tmp_searchPoint(	spx-tcds[j][0], 
												spy-tcds[j][1], 
												spz-tcds[j][2]);
	
			 	number_of_nn += kdtree.radiusSearch(tmp_searchPoint, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);				
			 	surface_cuts++;
			 }
		}

		// edges
		int edge_cuts=0;
		if(surface_cuts >= 2)
		{
			for(int j=0; j<12; j++)
			{
				pointIdxRadiusSearch.push_back(std::vector<int>());
				pointRadiusSquaredDistance.push_back(std::vector<float>());
				pt++;
				pcl::PointXYZ tmp_searchPoint2(	spx+((tcde[j][0]==1)?(0.0):(1.0))*((tcde[j][3]==1)?(-1.0):(1.0)),
										 		spy+((tcde[j][1]==1)?(0.0):(1.0))*((tcde[j][4]==1)?(-1.0):(1.0)),
										 		spz+((tcde[j][2]==1)?(0.0):(1.0))*((tcde[j][5]==1)?(-1.0):(1.0)));
				number_of_nn+=kdtree.radiusSearch(tmp_searchPoint2, R, pointIdxRadiusSearch[pt],pointRadiusSquaredDistance[pt]);
				edge_cuts++;
			}
		}

		// corners
		if(edge_cuts >= 2)
		{
			for(int j=0; j<8; j++)
			{
				const auto distsq = (spx-tcdc[j][0])*(spx-tcdc[j][0])
									+(spy-tcdc[j][1])*(spy-tcdc[j][1])
									+(spz-tcdc[j][2])*(spz-tcdc[j][2]);

				if(distsq > Rsq)
				{
					pointIdxRadiusSearch.push_back(std::vector<int>());
					pointRadiusSquaredDistance.push_back(std::vector<float>());
					pt++;
					pcl::PointXYZ tmp_searchPoint3(	spx+((tcdc[j][0]==1)?(-1.0):(1.0)),
													spy+((tcdc[j][1]==1)?(-1.0):(1.0)),
													spz+((tcdc[j][2]==1)?(-1.0):(1.0)));
					
					number_of_nn+=kdtree.radiusSearch(tmp_searchPoint3, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);
				}
			}
		}

		for(int j=1; j<=pt; j++)
		{
			pointIdxRadiusSearch[0].insert(pointIdxRadiusSearch[0].end(), pointIdxRadiusSearch[j].begin(), pointIdxRadiusSearch[j].end());
		}
		

		// erase self
		pointIdxRadiusSearch[0].erase(pointIdxRadiusSearch[0].begin());

		// store
		neighbours.push_back(pointIdxRadiusSearch[0]);		
	}

}


// not working
void construct_RGG_3D_brute(const double R, const PointCloud& points, std::vector<std::vector<int>>& neighbours)
{
	std::cout << "Not Working!!!" << std::endl;
	int number_of_nn; 
	int size = points.npoints;
	int len = points.len;
	std::vector<int> pointIdxRadiusSearch;
	std::vector<float> pointRadiusSquaredDistance;

	std::cout << ">> " << size << std::endl;

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZ>);	// point cloud object
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;									// kd-Tree object
	pcl::PointXYZ searchPoint;											// search point (point_xy object)

	// Adjust point cloud size
	cloud->width = 3*3*3*size;
	cloud->height = 1;
	cloud->points.resize (cloud->width * cloud->height);

	// Populate point cloud object with the coordinates
	std::vector<double> offset = {-1, 0, 1};

	int counter = 0;

	for (size_t i=0; i<size; i++)
	{
					auto crds = points.crds(i);
		for (int ix=0; ix<3; ix++)
		{
			for (int iy=0; iy<3; iy++)
			{
				for (int iz=0; iz<3; iz++)
				{
					
//					std::cout << counter << "  " << ix << "  " << iy << "  " << iz << std::endl;
					counter++;
//					std::cout << counter << "  " << crds[0]+offset[ix] << "  " << crds[1] << "  " << crds[2] << std::endl;
	    			cloud->points[i].x = crds[0]+offset[ix];
	    			cloud->points[i].y = crds[1]+offset[iy];
	    			cloud->points[i].z = crds[2]+offset[iz];
//					std::cout << "done" << std::endl;
				}
			}
		}
	}


	// Construct kd-tree
	kdtree.setInputCloud(cloud);

	// Search tree
	for (int j=0; j<size; j++) 
	{
		auto crds = points.crds(j);
		searchPoint.x = crds[0];
		searchPoint.y = crds[1];
		searchPoint.z = crds[2];
		std::cout << counter << "  " << crds[0] << "  " << crds[1] << "  " << crds[2] << std::endl;

		// perform search
		number_of_nn = kdtree.radiusSearch (searchPoint, R, pointIdxRadiusSearch, pointRadiusSquaredDistance);

//		std::cout << R << "  " << pointIdxRadiusSearch.size() << std::endl;

		// erase the first element since it is the search point itself
		pointIdxRadiusSearch.erase(pointIdxRadiusSearch.begin());

		// apply modulo
		for (int k=0; k<pointIdxRadiusSearch.size(); k++) pointIdxRadiusSearch[k] = pointIdxRadiusSearch[k] % size;
		neighbours.push_back(pointIdxRadiusSearch);
		// fixme: für sehr kleine gitter kann der "gleiche" punkte zweimal der nachbar werden
		// wohl kein problem für große gitter
	}
}





/**
 * Random Geometric Graph
 * A class that can create a Random Geometric Graph on a user defined
 * point cloud.
 * @tparam PointCloud The Point Cloud that is used by the user.
 */
template <class PointCloud>
class RandomGeometricGraph
{
	private:
		PointCloud cloud;
	
	public:
		int npoints, len, dim;
		std::vector<std::vector<int>> neighbours;
		RandomGeometricGraph(const int len, const int dim, const double search_radius) : 
																	cloud(pow(len,dim), len, dim), 
																	npoints(pow(len,dim)), 
																	len(len), 
																	dim(dim)
		{
//			construct_RGG_3D_brute(search_rad, cloud, neighbours);
			construct_RGG_3D(search_radius, cloud, neighbours);
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
