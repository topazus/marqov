#ifndef VECTORHELPERS_H
#define VECTORHELPERS_H
#include <cmath>
#include <array>

// Generate random vector on the SymD-dimensional unit sphere

template <class RND, typename FPType, int SymD> 
auto rnddir(RND& rn) -> std::array<FPType, SymD>
{

	// Spherical coordinates according to:
	// https://sites.math.washington.edu/~morrow/335_12/sphericalCoords.pdf

    std::array<FPType, SymD> retval;
    retval[0] = 1;//beginning of recursion
    if(SymD > 1) //poor mans constexpr if
    {
        //set up that auxiliary data that might be erased later
        FPType angles[SymD-1];
        for(int i = 0; i < SymD - 1; ++i)
            angles[i] = 2.0*rn.d();
        //recursion for the polar(?) angles
        for(int j = 1; j < SymD - 1; ++j)
            retval[j] = retval[j-1] * angles[j-1]*(2.0 - angles[j-1]);
        for(int j = 0; j < SymD - 2; ++j)
            retval[j] = std::sqrt(retval[j]) * (1.0-angles[j]);
        retval[SymD - 2] = std::sqrt(retval[SymD - 2]);
        //fix up the azimuthal angle
        //note the dependency on retval[SymD-2] here. Hence the order is important
        retval[SymD - 1] = retval[SymD - 2] * cos(M_PI*angles[SymD - 2]);
        retval[SymD - 2] = retval[SymD - 2] * sin(M_PI*angles[SymD - 2]);
    }
    return retval;
}

#endif
