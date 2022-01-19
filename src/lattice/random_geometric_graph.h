#ifndef RGG_LATTICE_H
#define RGG_LATTICE_H

#include "util/points.h"	
#include "util/distance.h"	
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>


bool construct_RGG_3D(const double R, const PointCloud& cloud, std::vector<std::vector<int>>& neighbours)
{
	double upper[3]={0,0,0};
	double lower[3]={1,1,1};	
	
	double tcds[6][3]={{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}}; 				//parametrization of cube surfaces

	double tcde[12][6]={{0,0,1, 0,0,0},{0,0,1, 1,0,0},{0,0,1, 0,1,0},{0,0,1, 1,1,0},		//
				     {0,1,0, 0,0,0},{0,1,0, 1,0,0},{0,1,0, 0,0,1},{0,1,0, 1,0,1},		//parametrization of cube edges
				     {1,0,0, 0,0,0},{1,0,0, 0,1,0},{1,0,0, 0,0,1},{1,0,0, 0,1,1}};		//

	double tcdc[8][3]={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1},{1,1,1}};		//parametrization of cube corners

	int number_of_nn;
	int size = cloud.size();
	std::vector< std::vector<int> > pointIdxRadiusSearch;
	std::vector< std::vector<float> > pointRadiusSquaredDistance;

	pcl::PointCloud<pcl::PointXYZ>::Ptr pclcloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	pcl::PointXYZ searchPoint;

	// transform cloud into PCL format
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
		
		pointIdxRadiusSearch.push_back(std::vector<int>());
		pointRadiusSquaredDistance.push_back(std::vector<float>());
		int pt=0;
		number_of_nn = kdtree.radiusSearch(searchPoint, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);		
		int surface_cuts=0;
	
		for(int j=0; j<6; j++)
		{
			if( !( ( (searchPoint.x + tcds[j][0] < upper[0]) and (searchPoint.x + tcds[j][0] > lower[0]) ) and
			       ( (searchPoint.y + tcds[j][1] < upper[1]) and (searchPoint.y + tcds[j][1] > lower[1]) ) and
			       ( (searchPoint.z + tcds[j][2] < upper[2]) and (searchPoint.z + tcds[j][2] > lower[2]) ) ) )
				   {
						pointIdxRadiusSearch.push_back(std::vector<int>());
						pointRadiusSquaredDistance.push_back(std::vector<float>());
						pt++;
						pcl::PointXYZ tmp_searchPoint(searchPoint.x-tcds[j][0] , searchPoint.y-tcds[j][1] , searchPoint.z-tcds[j][2]);
	
						number_of_nn += kdtree.radiusSearch(tmp_searchPoint, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);				
						surface_cuts++;
					}
		}

		int edge_cuts=0;
		if(surface_cuts >= 2)
		{
			for(int j=0; j<12; j++)
			{
				if(true)
				{
		
					pointIdxRadiusSearch.push_back(std::vector<int>());
					pointRadiusSquaredDistance.push_back(std::vector<float>());
					pt++;
					pcl::PointXYZ tmp_searchPoint2(searchPoint.x+((tcde[j][0]==1)?(0.0):(1.0))*((tcde[j][3]==1)?(-1.0):(1.0)),
											 searchPoint.y+((tcde[j][1]==1)?(0.0):(1.0))*((tcde[j][4]==1)?(-1.0):(1.0)),
											 searchPoint.z+((tcde[j][2]==1)?(0.0):(1.0))*((tcde[j][5]==1)?(-1.0):(1.0)));
					number_of_nn+=kdtree.radiusSearch(tmp_searchPoint2, R, pointIdxRadiusSearch[pt],pointRadiusSquaredDistance[pt]);
					edge_cuts++;
				} 
			}
		}

		if(edge_cuts >= 2)
		{
			for(int j=0; j<8; j++)
			{
				if( ((searchPoint.x-tcdc[j][0])*(searchPoint.x-tcdc[j][0])+(searchPoint.y-tcdc[j][1])*(searchPoint.y-tcdc[j][1])+(searchPoint.z-tcdc[j][2])*(searchPoint.z-tcdc[j][2])) > R*R)
				{
										
					
					pointIdxRadiusSearch.push_back(std::vector<int>());
					pointRadiusSquaredDistance.push_back(std::vector<float>());
					pt++;
					pcl::PointXYZ tmp_searchPoint3(searchPoint.x+((tcdc[j][0]==1)?(-1.0):(1.0)),searchPoint.y+((tcdc[j][1]==1)?(-1.0):(1.0)),searchPoint.z+((tcdc[j][2]==1)?(-1.0):(1.0)));
					
					number_of_nn+=kdtree.radiusSearch(tmp_searchPoint3, R, pointIdxRadiusSearch[pt], pointRadiusSquaredDistance[pt]);
				}
			}
		}

		for(int j=1;j<=pt;j++)
		{
			pointIdxRadiusSearch[0].insert(pointIdxRadiusSearch[0].end(), pointIdxRadiusSearch[j].begin(), pointIdxRadiusSearch[j].end());
		}
		
		pointIdxRadiusSearch[0].erase(pointIdxRadiusSearch[0].begin());
		neighbours.push_back(pointIdxRadiusSearch[0]);		
	}
}


// not working
bool construct_RGG_3D_brute(const double R, const PointCloud& points, std::vector<std::vector<int>>& neighbours)
{
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
