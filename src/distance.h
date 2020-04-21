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
