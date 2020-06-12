#include "prng.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using std::cout;
using std::endl;


// auxiliary functions
inline std::size_t myrndnbr(LCG& lcg, int16_t len)
{
	if (len < 2) return 0;
	
	std::size_t ret = 0;
	uint mr = lcg.rnd();
	switch(len)
	{
		case 1:
		ret = 0;
		break;
		case 2:
		ret = mr%2;
		case 3:
		ret = mr%3;
		case 4:
		ret = mr%4;
		break;
		case 5:
		ret = mr%5;
		break;
		case 6:
		ret = mr%6;
		break;
		case 7:
		ret = mr%7;
		break;
		case 8:
		ret = mr%8;
		break;
		default:
		ret = mr%len;
		break;
	}
	return ret;
}


// shuffle ordinary std::vector

void shuffle_vector(std::vector<int>& vec, int seed=1)
{
	std::mt19937_64 g( seed+std::random_device{}()+time(NULL));
	std::shuffle(std::begin(vec), std::end(vec), g);
}




inline void updateneighbours(std::vector<std::vector<int>>& neighbours, const int i, const int j, const int i1, const int j1)
{
	for (int q=0; q<neighbours[i].size();  q++) if (neighbours[i][q]  == i1) {neighbours[i][q]  = j; continue;}
	for (int q=0; q<neighbours[i1].size(); q++) if (neighbours[i1][q] == i)  {neighbours[i1][q] = j1; continue;}
	for (int q=0; q<neighbours[j].size();  q++) if (neighbours[j][q]  == j1) {neighbours[j][q]  = i; continue;}
	for (int q=0; q<neighbours[j1].size(); q++) if (neighbours[j1][q] == j)  {neighbours[j1][q] = i1; continue;}
}




bool has_idx(std::vector<int>& n,int16_t len, int j)
{
    for (int i = 0; i < len; ++i)
    if(n[i] == j) return true;
    return false;
}



// use only, if the points are not sorted in some way
// if they are shuffle them before. this is crucial!
bool force_fixed_neighbours_naive_box(const int K, std::vector< std::vector<int> >& neighbours, const std::vector<int>& box, int offset )
{
     int size = box.size();

	// K must be even (Handshaking Lemma)
     if (K%2 == 1) { throw std::invalid_argument("error! K must be even!"); }

	// check if there are enough points
     if (size < K+1+2*offset) { throw std::invalid_argument("error! Found box with too few sites!"); }

     for (int i=0; i<size; i++)
     {
          for (int j=-K/2-offset; j<0-offset; j++)
          {
               neighbours[box[i]].push_back(box[(i+j+size)%size]);
          }
          for (int j=1+offset; j<K/2+1+offset; j++)
          {
               neighbours[box[i]].push_back(box[(i+j)%size]);
          }
     }

	
	bool duplicate;
     for (int i=0; i<size; i++)
	{
		std::vector<int> local_copy = neighbours[box[i]];

		std:sort(local_copy.begin(), local_copy.end());
		auto it = std::unique(local_copy.begin(), local_copy.end());
		duplicate = !(it == local_copy.end());
		if (duplicate) break;
	}



	if (duplicate)
	{
     	for (int i=0; i<size; i++) 
		{
			for (int j=0; j<K; j++) neighbours[box[i]].pop_back();
		}
	}
	if (duplicate) return false;
	else return true;
}


double avgNeighbours(std::vector<std::vector<int>> neighbours)
{
     int sum = 0;
     int number = neighbours.size();
     for (int i =0; i<number; i++) sum += neighbours[i].size();
     return sum / double(number);
}



int geometric_simulated_annealing_box(const PointCloud& cloud, std::vector<std::vector<int>>& neighbours, const std::vector<int>& box, int K, long int Nsteps, int seed=1)
{
	const int nparticles = box.size();

	const bool verbose = false;

	// general procedure: 
	// choose two bonds i--i1 and j--j1 (that is four points, i,i1,j,j1)
	// remove those bonds and reconnect the points i--j and i1--j1 
	// (important remark: choosing i--j1, i1--j as the new bonds does not work in all cases!)
	// if the new setting is better than the old one, accept it; else accept it only with a certain probability
	
	int ifc = 0;
	int elc = 0;
	int llc = 0;
	int conc = 0;
	
	// prepare random number generator for choosing sites and neighbours
	std::mt19937_64 g( seed+std::random_device{}()+time(NULL));
	std::uniform_int_distribution<> dI(0,nparticles-1);
	std::uniform_int_distribution<> dK(0,K-1);

	for (int k=0; k<Nsteps; k++)
	{
		// choose two random site on the point cloud
		int i = box[dI(g)];
		int j = box[dI(g)];
		
		// i and j must neither be the same points nor neighbours of each other
		if ((i==j) || has_idx(neighbours[i], K, j)) {conc++; continue;}
		
		// choose a neighbour of i and j, respectively
		// this determines two bonds i--i1 and j--j1
		int i1 = neighbours[i][dK(g)];
		int j1 = neighbours[j][dK(g)];
		
		// i1 and j1 must neither be the same points nor neighbours of each other
		if ((i1==j1) || has_idx(neighbours[i1], K, j1)) {conc++; continue; }
		
		// no two points of {i, j, i1, j1} must be the same
		// this is fullfilled implicitly
		
		// calculate the (squared!) distances of the bonds i--i1, j--j1, i--j and i1--j1 
		const double d1 = distancePBSQ_nD(cloud.getcrds(i),  cloud.getcrds(i1));
		const double d2 = distancePBSQ_nD(cloud.getcrds(j),  cloud.getcrds(j1));
		const double d3 = distancePBSQ_nD(cloud.getcrds(i),  cloud.getcrds(j));
		const double d4 = distancePBSQ_nD(cloud.getcrds(i1), cloud.getcrds(j1));
		
		if ( (d3+d4) < (d1+d2) )
		{
			ifc++;
			updateneighbours(neighbours, i, j, i1, j1);
		}
		else {llc++;}
	}


	if (verbose)
	{
		double sum = double(ifc+elc+llc);
		cout << endl << ifc/sum <<  "\t" << llc/sum << endl;
		cout         << ifc     <<  "\t" << llc     << "\t" << conc << endl;
	}

	return ifc;
}





// new (third) approach to CC lattice using local boxes
// supposed to scale as O(L^2) while being  scale free and very isotropic

bool constant_coordination_lattice(const PointCloud& cloud, std::vector<std::vector<int>>& neighbours)
{
	neighbours.clear();

	const int K1 = 2;						// coordination number first step
	const int K2 = 2;						// coordination number second step

	const int K = K1 + K2;					// total coordination number

	const int size = cloud.size();				// number of particles

	const int L = int(sqrt(size));			// linear system length

	std::vector<int> nbrsizes;
	nbrsizes.reserve(size);
	for (uint k = 0; k < size; ++k) nbrsizes.push_back(K);
	neighbours.resize(size);

	// set up boxes for initial connection
	const int boxlen = 4;
	const int nboxes1D = L / boxlen;
	const int nboxes2D = nboxes1D*nboxes1D;
	std::vector<std::vector<std::vector<int>>> boxes(nboxes1D, std::vector<std::vector<int>>(nboxes1D));

	// set up boxes for dynamical rewiring
	const int boxlen_r = 4;
	const int nboxes1D_r = L / boxlen_r;
	const int nboxes2D_r = nboxes1D_r*nboxes1D_r;
	std::vector<std::vector<std::vector<int>>> boxes_r(nboxes1D_r, std::vector<std::vector<int>>(nboxes1D_r));

	// dynamical rewiring attempts (per box)

//	double a = 100;	
	double a = 50;	

	const long int Nsteps = int(K*a*pow(boxlen_r,4));
     if (L%(2*boxlen)   != 0) { throw std::invalid_argument("error! L needs to be a multiple of two times the box length!"); }
     if (L%(2*boxlen_r) != 0) { throw std::invalid_argument("error! L needs to be a multiple of two times the box length!"); }
	
	// bin the points into boxes
	for (int i=0; i<L*L; i++)
	{
		int ix = floor(cloud.getx(i) * nboxes1D);
		int iy = floor(cloud.gety(i) * nboxes1D);
		boxes[ix][iy].push_back(i);

		int ix_r = floor(cloud.getx(i) * nboxes1D_r);
		int iy_r = floor(cloud.gety(i) * nboxes1D_r);
		boxes_r[ix_r][iy_r].push_back(i);
	}
	
    std::vector<int> workbox;
	// loop over boxes and rewire
	for (int ix=0; ix<nboxes1D; ix+=2) 
	{
		for (int iy=0; iy<nboxes1D; iy+=2) 
		{
	
			int downleft_idx = ix;
			int downleft_idy = iy;
	
			int upright_idx   = downleft_idx + 1;
			int upright_idy   = downleft_idy + 1;
			int downright_idx = downleft_idx + 1;
			int downright_idy = downleft_idy;
			int upleft_idx    = downleft_idx;
			int upleft_idy    = downleft_idy + 1;
	
			// account for periodic boundaries (nboxes1D+1=0)
			if (upright_idx   >= nboxes1D) upright_idx   = 0;
			if (downright_idx >= nboxes1D) downright_idx = 0;
			if (upleft_idx    >= nboxes1D) upleft_idx    = 0;
			if (upright_idy   >= nboxes1D) upright_idy   = 0;
			if (downright_idy >= nboxes1D) downright_idy = 0;
			if (upleft_idy    >= nboxes1D) upleft_idy    = 0;
            
            	workbox.clear();

			for (int i=0; i<boxes[upleft_idx][upleft_idy].size(); i++)       
				workbox.push_back(boxes[upleft_idx][upleft_idy][i]);
			for (int i=0; i<boxes[downleft_idx][downleft_idy].size(); i++)   
				workbox.push_back(boxes[downleft_idx][downleft_idy][i]);
			for (int i=0; i<boxes[upright_idx][upright_idy].size(); i++)     
				workbox.push_back(boxes[upright_idx][upright_idy][i]);
			for (int i=0; i<boxes[downright_idx][downright_idy].size(); i++) 
				workbox.push_back(boxes[downright_idx][downright_idy][i]);
	
			shuffle_vector(workbox);
			bool success = force_fixed_neighbours_naive_box(K1, neighbours, workbox, 0);
		}
	}

	// save_geometry(cloud, neighbours, geomdir+"/geom-"+to_string(L)+"-"+to_string(1)+"-1.dat");

	for (int ix=1; ix<nboxes1D; ix+=2) 
	{
		for (int iy=1; iy<nboxes1D; iy+=2) 
		{
			int downleft_idx = ix;
			int downleft_idy = iy;
	
			int upright_idx   = downleft_idx + 1;
			int upright_idy   = downleft_idy + 1;
			int downright_idx = downleft_idx + 1;
			int downright_idy = downleft_idy;
			int upleft_idx    = downleft_idx;
			int upleft_idy    = downleft_idy + 1;
	
			// account for periodic boundaries (nboxes1D+1=0)
			if (upright_idx   >= nboxes1D) upright_idx   = 0;
			if (downright_idx >= nboxes1D) downright_idx = 0;
			if (upleft_idx    >= nboxes1D) upleft_idx    = 0;
			if (upright_idy   >= nboxes1D) upright_idy   = 0;
			if (downright_idy >= nboxes1D) downright_idy = 0;
			if (upleft_idy    >= nboxes1D) upleft_idy    = 0;

            	workbox.clear();

			for (int i=0; i<boxes[upleft_idx][upleft_idy].size(); i++)       
				workbox.push_back(boxes[upleft_idx][upleft_idy][i]);
			for (int i=0; i<boxes[downleft_idx][downleft_idy].size(); i++)   
				workbox.push_back(boxes[downleft_idx][downleft_idy][i]);
			for (int i=0; i<boxes[upright_idx][upright_idy].size(); i++)     
				workbox.push_back(boxes[upright_idx][upright_idy][i]);
			for (int i=0; i<boxes[downright_idx][downright_idy].size(); i++) 
				workbox.push_back(boxes[downright_idx][downright_idy][i]);

			bool success = false;
			while (!success)
			{
				shuffle_vector(workbox);
				success = force_fixed_neighbours_naive_box(K2, neighbours, workbox, 0);
			}
		}
	}

	// save_geometry(cloud, neighbours, geomdir+"/geom-"+to_string(L)+"-"+to_string(1)+"-2.dat");

	for (int k=0; k<4; k++)
	{
		int xinitial, yinitial;
		if (k==0) { xinitial=0; yinitial=0; } 
		if (k==1) { xinitial=1; yinitial=1; } 
		if (k==2) { xinitial=0; yinitial=1; } 
		if (k==3) { xinitial=1; yinitial=0; } 
		for (int ix=xinitial; ix<nboxes1D_r; ix+=2) 
		{
			for (int iy=yinitial; iy<nboxes1D_r; iy+=2) 
			{
				int downleft_idx = ix;
				int downleft_idy = iy;
		
				int upright_idx   = downleft_idx + 1;
				int upright_idy   = downleft_idy + 1;
				int downright_idx = downleft_idx + 1;
				int downright_idy = downleft_idy;
				int upleft_idx    = downleft_idx;
				int upleft_idy    = downleft_idy + 1;
		
				// account for periodic boundaries (nboxes1D+1=0)
				if (upright_idx   >= nboxes1D_r) upright_idx   = 0;
				if (downright_idx >= nboxes1D_r) downright_idx = 0;
				if (upleft_idx    >= nboxes1D_r) upleft_idx    = 0;
				if (upright_idy   >= nboxes1D_r) upright_idy   = 0;
				if (downright_idy >= nboxes1D_r) downright_idy = 0;
				if (upleft_idy    >= nboxes1D_r) upleft_idy    = 0;
		
				workbox.clear();
				for (int i=0; i<boxes_r[upleft_idx][upleft_idy].size(); i++)       
					workbox.push_back(boxes_r[upleft_idx][upleft_idy][i]);
				for (int i=0; i<boxes_r[downleft_idx][downleft_idy].size(); i++)   
					workbox.push_back(boxes_r[downleft_idx][downleft_idy][i]);
				for (int i=0; i<boxes_r[upright_idx][upright_idy].size(); i++)     
					workbox.push_back(boxes_r[upright_idx][upright_idy][i]);
				for (int i=0; i<boxes_r[downright_idx][downright_idy].size(); i++) 
					workbox.push_back(boxes_r[downright_idx][downright_idy][i]);
		
				int ww = geometric_simulated_annealing_box(cloud, neighbours, workbox, K, Nsteps);
			}
		}
	}

	return true;
}
