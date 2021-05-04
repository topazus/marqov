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

#ifndef VECTORHELPERS_H
#define VECTORHELPERS_H
#include <cmath>
#include <array>
#include <type_traits>

/**
 * A helper class to generate a random vector on an SymD-dimensional unit
 * sphere
 * @tparam RNG The rng that is used for drawing random numbers
 * @tparam valuetype The type of the elements of the vector.
 * @tparam SymD The dimensionality
 * @tparam Enable A helper to distinguish cardinal types from floating point
 *                types and choose a different implementation.
 */
template <class RNG, typename valuetype, int SymD, typename Enable = void> 
struct Rnddir_Helper
{
static auto rnddir(RNG& rn) -> std::array<valuetype, SymD>
{
	// Spherical coordinates according to:
	// https://sites.math.washington.edu/~morrow/335_12/sphericalCoords.pdf
    std::array<valuetype, SymD> retval;
    retval[0] = 1;//beginning of recursion
    if(SymD > 1) //poor mans constexpr if
    {
        //set up that auxiliary data that might be erased later
        valuetype angles[SymD-1];
        for(int i = 0; i < SymD - 1; ++i)
            angles[i] = 2.0*rn.real();
        //recursion for the polar(?) angles
        for(int j = 1; j < SymD - 1; ++j)
            retval[j] = retval[j-1] * angles[j-1]*(2.0 - angles[j-1]);
        for(int j = 0; j < SymD - 2; ++j)
            retval[j] = std::sqrt(retval[j]) * (1.0-angles[j]);
        retval[SymD - 2] = std::sqrt(retval[SymD - 2]);
        //fix up the azimuthal angle
        //note the dependency on retval[SymD-2] here. Hence the order is important
        retval[SymD - 1] = retval[SymD - 2] * std::cos(M_PI*angles[SymD - 2]);
        retval[SymD - 2] = retval[SymD - 2] * std::sin(M_PI*angles[SymD - 2]);
    }
    return retval;
}
};

template <class RND, typename inttype>
struct Rnddir_Helper<RND, inttype, 1, 
typename std::enable_if<std::is_integral<inttype>::value>::type> // enable only for integers
{
static auto rnddir(RND& rn) -> std::array<inttype, 1>
{
    std::array<inttype, 1> retval = {1};
    if (rn.real() < 0.5) retval[0] = -1;
    return retval;
}
};

template <class RNG, typename valuetype, int SymD> 
auto rnddir(RNG& rn) -> std::array<valuetype, SymD>
{
    return Rnddir_Helper<RNG, valuetype, SymD>::rnddir(rn);
}
#endif
