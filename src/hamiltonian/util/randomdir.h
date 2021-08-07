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

/** Helper for random directions on an nD unit sphere.
 * 
 * A helper class to draw a random vector on an SymD-dimensional unit sphere.
 * @tparam RNG The rng that is used for drawing random numbers
 * @tparam valuetype The type of the elements of the vector.
 * @tparam SymD The dimensionality
 * @tparam Enable A helper to distinguish cardinal types from floating point
 *                types and choose a different implementation.
 */
template <class RNG, typename T, int SymD, typename Enable = void> 
struct Rnddir_Helper
{
    /** Draw a random directions on an n-D unit sphere.
     * 
     * Formulas are from this
     * [external PDF](https://sites.math.washington.edu/~morrow/335_12/sphericalCoords.pdf) with some hints.
     * @note For the used rng, we assume it behaves like the RNGCache ...
     * @param rn where to get the random numbers from.
     * @param retval On input: storage to use for the new random vector. On return a random direction uniformly distributed on the unit sphere.
     */
    static void set_to_rnddir(RNG& rn, T& retval) noexcept
    {
        //beginning of recursion already defined outside.
        if(SymD > 1) //poor mans constexpr if
        {
            //set up that auxiliary data that might be erased later
            typename T::value_type angles[SymD-1];
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
    }
};

/** Specialization for arrays of 1 arithmetic type.
 * 
 * If the state space is a scalar, the only directions are +-1.
 * @see Rnddir_Helper
 * @tparam RND The rng that is used for drawing random numbers

 */
template <class RND, typename T>
struct Rnddir_Helper<RND, T, 1, 
typename std::enable_if<std::is_arithmetic<typename T::value_type>::value>::type> // enable only for arithmetic types
{
    /** Draw a random 1D direction.
     * 
     * So effectively fot all arithmetic types, this is either +1 or -1.
     * @note For the used rng, we assume it behaves like the RNGCache ...
     * @param rn where to get the random numbers from.
     * @param t where to write the bit. On return: Either +1 or -1, chosen randomly.
     */
    static void set_to_rnddir(RND& rn, T& t) noexcept
    {
        if (rn.real() < 0.5) t[0] = -t[0];
    }
};

/** Create a random direction.
 * 
 * This creates a random direction in SymD space with proper type.
 * 
 * @see Rnddir_Helper
 * @param rn the RNG to draw the required random numbers from.
 * @param t 
 * @return the random direction in the respective space.
 */
template <class RNG, typename T>
void set_to_rnddir(RNG&& rn, T& t) noexcept
{
    t[0] = typename T::value_type{1};
    Rnddir_Helper<RNG, T, std::tuple_size<T>::value>::set_to_rnddir(std::forward<RNG>(rn), t);
}
#endif
