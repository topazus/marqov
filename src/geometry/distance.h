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

#ifndef DISTANCE_H
#define DISTANCE_H

#include <cmath>

// returns squared nD distances, respecting periodic boundaries in the unit hypercube.
template <class Container>
double distancePBSQ_nD(const Container& a, const Container& b)
{
    Container diff(a);
    for(std::size_t i = 0; i < a.size(); ++i)//form the deltas
        diff[i] -= b[i];
    for(std::size_t i = 0; i < a.size(); ++i)//check PBCs
        if(fabs(diff[i]) > 0.5) diff[i] = 1.0 - fabs(diff[i]);
    
    typename Container::value_type retval = 0.0;
    for(decltype(diff.size()) i = 0; i < diff.size(); ++i)
        retval += diff[i] * diff[i];
    return retval;
}

#endif
