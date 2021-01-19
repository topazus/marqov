#ifndef SSHLATTE_H
#define SSHLATTE_H

// covers the 1+1 dimensional and the 2+1 dimensional case
// note that the dimension is set via the proceessor statement SSH_2D

class SSHLattice
{
	public:
	
		int counter = 0;
	
		typedef std::vector<int> value_type;
		
		SSHLattice() {}
		
		SSHLattice(int l, int ltau, int d) : len(l), lentime(ltau), dim(d) // dim is the spatial dimension!
		{
			// dim from config file will be ignored (improve me!)

			#ifdef SSH_2D
				dim = 2;
			#else
				dim = 1;
			#endif

			vol = dim * pow(len, dim); // spatial volume (number of bonds in a time slice)
			nsites = vol * lentime; // total volume (number of sites)
		}


		//    illustration of the bond arrangement in 2+1 dimensions
		//    horizontal and vertical bonds are arranged in an alternating fashion
		//    shown are the time slices t=0 and t=1 for a lattice of size L=3
		//    X's denote the coordinate grid, numbers represent the bonds
		// 
		//    o       o       o                o       o       o
		//    |       |       |                |       |       |
		//    15      16      17               33      34      35
		//    |       |       |                |       |       |
		//    X---12--X---13--X---14--o        X---30--X---31--X---32--o
		//    |       |       |                |       |       |
		//    9       10      11               27      28      29
		//    |       |       |                |       |       |
		//    X---6---X---7---X---8---o        X---24--X---25--X---26--o
		//    |       |       |                |       |       |
		//    3       4       5                21      22      23
		//    |       |       |                |       |       |
		//    X---0---X---1---X---2---o        X---18--X---19--X---20--o


		
		value_type getnbrs(int a, int i) const 
		{
			value_type retval;
  			int tslice = i / vol; // time slice
  			int offset = i % vol; // position (offset) in time slice

			switch(a) 
			{
				// returns the (two) neighbours along the time direction
				case 0:
  					retval = {((tslice-1+lentime)%lentime)*vol+offset, ((tslice+1)%lentime)*vol+offset};
					break;

				// returns the entire lattice
				case 1:
					retval.reserve(this->nsites);
					for(int j = 0; j < i; ++j)
						retval.push_back(j);
					for(int j = i+1; j < this->nsites; ++j)
						retval.push_back(j);
					break;
			}

			return retval;
		}


		std::vector<double> getcrds(int i) const 
		{
			// returns bond coordinates, which is five numbers in 2+1 dimension
   			// (xstart, xend, ystart, yend, time)
			// and only three numbers in 1+1 dimension (xstart, xend, time)
   			
   			const int tslice = i / vol;  // time slice
   			const int offset = i % vol;  // position (offset) in time slice

			#ifdef SSH_2D
   			
   				const int yslice = offset / len;  // y-row
   				const int xpos   = offset % len;  // position in y-row
   				
   				// since horizontal and vertical bonds are arranged in an alternating fashion (designers choice!)
   				// we require to extract the actual y coordinate
   				int yreal = yslice/2;

				if (yslice%2 == 0)
				{
					const int left = xpos;
					const int right = (xpos+1)%len;  // account for p.b.c in x
					return {left, right, yreal, yreal, tslice};
				}
				else
				{
					const int lower = yreal;
					const int upper = (yreal+1)%len; // account for p.b.c in y
					return {xpos, xpos, lower, upper, tslice};
				}

			#else
				const int left = offset;
				const int right = (offset+1)%len;
				return {left, right, tslice};
			#endif
		 }




		
		 value_type getbnds(int a, int i) const 
		 {
		    value_type retval;
		    switch(a) 
		    {
		        case 0:
		            retval = {1,1};
		            break;
		        case 1:
			   		cout << "This should not have happened!" << endl;
		            retval.reserve(this->nsites);
		            for(int j = 0; j < i; ++j)
		            	retval.push_back(1.0);
		            for(int j = i+1; j < this->nsites; ++j)
		            	retval.push_back(1.0);
		            break;
		    }
		    return retval;
		}
		 
		
		std::size_t size() const {return nsites;}
		int len, lentime, vol;
		int dim;
		std::size_t nsites;
};

#endif
