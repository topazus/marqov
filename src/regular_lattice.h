#ifndef REGULARLATTE_H
#define REGULARLATTE_H

class RegularLattice
{
public:
//    friend class NArray_Iterator;
    typedef std::vector<int> value_type;

    RegularLattice() {}

    RegularLattice(int l, int d) : length(l), dim(d), pows(dim) 
    {
        numberatoms = 1;
        for(int i = 0; i < dim; ++i)
        {
            pows[i] = numberatoms;
            numberatoms *= length;
        }
    }

    value_type getnbrs(int a, int i) const {return this->operator[](i);}

    std::vector<double> getcrds(int i) const 
    {

    		// Todo: 
		// - test me!
		// - closed formula for the n-d case?
		// - length of retval is known at compile time -> probably not use vectors


    		std::vector<double> retval(dim);

		// transform one-dimensional index to n-d coordinates

		if (dim == 2)
		{
			retval[0] = i % length;
			retval[1] = i / length; // integer division
		}
		else if (dim == 3)
		{
			int j = i % (length*length);

			retval[0] = i / (length*length);
			retval[1] = j / length;;
			retval[2] = j % length;
		}
		else
		{
			std::cout << "error! not implemented!" << std::endl;
		}

    		return retval;
	}

		value_type operator[](int i) const
		{
			std::vector<int> temp(2*dim);
			//calculate neighbours for site i
			for(int j = 0; j < dim; ++j)
			{
				int pl = pows[j]*length;

				//positive additions
				int c = i + pows[j];

				//test positive additions for PBCs
				if (c >= (c/pl)*pl)
				{
					temp[2*j] = (i/pl)*pl + c % pl;
				}
				else
				{
					temp[2*j] = c;
				}
				
				//negative additions
				c = i - pows[j];

				//test negative additions for PBCs
				if (c < (i/pl)*pl)
				{
					temp[2*j+1] = (i/pl)*pl + (c + pl) % pl;
				}
				else
				{
					temp[2*j+1] = c;
				}
			}
			return temp;
		}

    std::size_t size() const {return numberatoms;}
    int length;
    int dim;
private:
    std::size_t numberatoms;
    std::vector<int> pows;
};

#endif
