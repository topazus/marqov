#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "lattice/util/regular_lattice.h"

using namespace std;

std::vector<int> refIndexOf(int k, int nDim, int nBin)
{
    std::vector<int> indices;

    for (int i=0; i<nDim; i++)
    {
        double index = std::fmod( double(k)/std::pow(nBin,i), nBin);
        indices.push_back(int(index));
        k -= index * std::pow(nBin,i);
    }
    return indices;
}


class RefLattice
{
public:
    typedef std::vector<int> value_type;
    RefLattice(int l, int d) : length(l), dim(d), pows(dim) 
    {
        numberatoms = 1;
        for(int i = 0; i < dim; ++i)
        {
            pows[i] = numberatoms;
            numberatoms *= length;
        }
    }

    std::vector<double> crds(int k) const 
    {
		// transform one-dimensional index to n-d coordinates
    		std::vector<int> indices = refIndexOf(k, dim, length);
		std::vector<double> retval(dim,0);

		// transform to double and normalize to unit hypercube
		for (decltype(retval.size()) i=0; i<retval.size(); i++) 
		{ 
			retval[i] = double(indices[i])/length; 
			retval[i] += 0.5/length;
		}

    		return retval;
	}

    std::size_t size() const {return numberatoms;}
    int length{0};
    int dim{0};
private:
    std::size_t numberatoms{0};
    std::vector<int> pows{0};
};


int main()
{
    int l = 16;
    for(uint ndim = 1; ndim < 5; ++ndim)
    {
        long nsize = pow(l, ndim);
        RefLattice ref(l, ndim);
        RegularLattice rl(l, ndim);
        for(int k = 0; k < nsize; ++k)
            if(ref.crds(k) != rl.crds(k)) {exit(-1);}
    }
}
