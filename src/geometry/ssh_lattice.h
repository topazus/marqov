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

			const int k = i / len;
			const int offset = i % len;

			switch(a) 
			{
				case 0:
					retval = {((k-1+lentime)%lentime)*len+offset, ((k+1)%lentime)*len+offset};
					break;
				case 1:
					retval.reserve(this->nsites);
					for(int j = 0; j < i; ++j)
						retval.push_back(j);
					for(int j = i+1; j < this->nsites; ++j)
						retval.push_back(j);
					break;
			}

//			cout << i << " ---   " << endl;
//			for (int i=0; i<retval.size(); i++) cout << retval[i] << "\t";
//			cout << endl;
//			cout << endl;

			return retval;
		}





//  3  oooooooo
//  2  oooooooo
//  1  oooooooo
//  0  oooooooo
//     01234556...


		
		std::vector<double> getcrds(int k) const 
		{
		 	// transform one-dimensional index to n-d coordinates
			std::vector<int> indices = IndexOfRect(k, dim, len, lentime);
		 	std::vector<double> retval(dim,0);

//			cout << k << "  " << indices[0] << "  " << indices[1] << endl;
		
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
		int len, lentime;
		int dim;
		std::size_t nsites;
};

#endif
