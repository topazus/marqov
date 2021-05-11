/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
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

#ifndef SSHLATTE_H
#define SSHLATTE_H

#include <vector>
#include <cmath>

class SSHLattice
{
	public:
	
		int counter = 0;
	
		typedef std::vector<int> value_type;
		
		SSHLattice() {}
		
		SSHLattice(int l, int ltau, int d) : len(l), lentime(ltau), dim(d)
		{
		    nsites = std::pow(len, dim-1);
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
//					retval = {(i-1+lentime)%lentime+offset, (i+1)%lentime+offset};
					retval = {((k-1+lentime)%lentime)*len+offset, ((k+1)%lentime)*len+offset};
					break;
				case 1:
					retval.reserve(this->nsites);
					for(int j = 0; j < i; ++j)
						retval.push_back(j);
					for(int j = i+1; j < static_cast<int>(this->nsites); ++j)
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
		            for(int j = i+1; j < static_cast<int>(this->nsites); ++j)
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
