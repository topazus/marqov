#ifndef SSHLATTE_H
#define SSHLATTE_H

class SSHLattice
{
	public:
	
		int counter = 0;
	
		typedef std::vector<int> value_type;
		
		SSHLattice() {}
		
		SSHLattice(int l, int ltau, int d) : len(l), lentime(ltau), dim(d)
		{
		    vol = pow(len, dim-1); // spatial volume
		    nsites = vol * lentime; // total volume
		}
		
		value_type getnbrs(int a, int i) const 
		{
			value_type retval;

			const int k = i / vol;
			const int offset = i % vol;

			switch(a) 
			{
				case 0:
					retval = {((k-1+lentime)%lentime)*vol+offset, ((k+1)%lentime)*vol+offset};
					break;
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


// spatial ordering for a 2+1 dimensional setting (L=3, Ltime>3)

//    | | |    | | |   | | |
//    | | .    | | .   | | . 
//    | . 24   | . 25  | . 26
//    . 21.    . 22.   . 23.
//    18. |    19. |   20. |
//    . | |    . | |   . | |
//    | | |    | | |   | | |
//    | | .    | | .   | | .
//    | . 15   | . 16  | . 17
//    . 12.    . 13    . 14.
//    9 . |    10. |   11 .|
//    . | |    . | |   . | |
//    | | |    | | |   | | |
//    | | .    | | .   | | .
//    | . 6    | . 7   | . 8
//    . 3      . 4     . 5
//    0        1       2

// e.g. getnbrs(a=0, 10): 1, 19


		
		std::vector<double> getcrds(int k) const 
		{
		 	// transform one-dimensional index to n-d coordinates
			std::vector<int> indices = IndexOfRect(k, dim, len, lentime);
		 	std::vector<double> retval(dim,0);

		 	for (int i=0; i<retval.size(); i++) 
		 	{ 
		 		retval[i] = double(indices[i]); ///len; 
//		 		retval[i] += 0.5/len;
		 	}
		
			return retval;
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
		            retval.reserve(this->nsites);
		            for(int j = 0; j < i; ++j)
		            	retval.push_back(1.0/j);
		            for(int j = i+1; j < this->nsites; ++j)
		            	retval.push_back(1.0/j);
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
