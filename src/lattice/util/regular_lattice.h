/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2022, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef REGULAR_LATTICE_H
#define REGULAR_LATTICE_H

#include <vector>

class RegularLattice
{
public:
    typedef std::vector<int> value_type;
    RegularLattice(int l, int d) : length(l), dim(d), pows(dim) 
    {
        numberatoms = 1;
        for(int i = 0; i < dim; ++i)
        {
            pows[i] = numberatoms;
            numberatoms *= length;
        }
    }

    value_type nbrs(int a, int i) const {return this->operator[](i);}

    std::vector<double> crds(int k) const 
    {
        std::vector<double> retval(dim);
        double mya = static_cast<double>(1)/length;
        for (int i=0; i < dim; i++)
        {
            int index = (k/pows[i])%length;
            k -= index * pows[i];
            retval[i] = index * mya;
            retval[i] += 0.5 * mya;
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
    int length{0};
    int dim{0};
private:
    std::size_t numberatoms{0};
    std::vector<int> pows{0};
};

#endif
