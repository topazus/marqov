#ifndef SSHLATTE_H
#define SSHLATTE_H

class SSHLattice
{
public:

	int counter = 0;

//    friend class NArray_Iterator;
    typedef std::vector<int> value_type;

    SSHLattice() {}

    SSHLattice(int l, int d) : length(l), dim(d), pows(dim) 
    {
        numberatoms = 1;
        for(int i = 0; i < dim; ++i)
        {
            pows[i] = numberatoms;
            numberatoms *= length;
        }
    }

    value_type getnbrs(int a, int i) const {return this->operator[](i);}

    std::vector<double> getcrds(int k) const 
    {
		// transform one-dimensional index to n-d coordinates
    		std::vector<int> indices = IndexOf(k, dim, length);
		std::vector<double> retval(dim,0);

		// transform to double and normalize to unit hypercube
		for (int i=0; i<retval.size(); i++) 
		{ 
			retval[i] = double(indices[i])/length; 
			retval[i] += 0.5/length;
		}

    		return retval;
	}

		value_type operator[](int i) const
		{
			int j = i / length;

			int offset = j * length;



			return std::vector<int> {(i-1+length)%length+offset, (i+1)%length+offset};

		}

    std::size_t size() const {return numberatoms;}
    int length;
    int dim;
private:
    std::size_t numberatoms;
    std::vector<int> pows;
};

#endif
