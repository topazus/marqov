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
		    nsites = pow(len, dim-1);
		    nsites = nsites * lentime;
		}
		
		value_type getnbrs(int a, int i) const 
		{
			value_type retval;
			const int j = i / lentime;
			const int offset = j * lentime;
			switch(a) 
			{
				case 0:
					retval = {(i-1+lentime)%lentime+offset, (i+1)%lentime+offset};
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
		
		/*
		std::vector<double> getcrds(int k) const 
		{
		 	// transform one-dimensional index to n-d coordinates
				std::vector<int> indices = IndexOf(k, dim, len);
		 	std::vector<double> retval(dim,0);
		
		 	// transform to double and normalize to unit hypercube
		 	for (int i=0; i<retval.size(); i++) 
		 	{ 
		 		retval[i] = double(indices[i])/len; 
		 		retval[i] += 0.5/len;
		 	}
		
				return retval;
		 }
		 */
		
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
		int len, lentime;
		int dim;
		std::size_t nsites;
};

#endif
