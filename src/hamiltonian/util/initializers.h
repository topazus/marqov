/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2021-2022, The MARQOV Project
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

#ifndef INITIALIZERS_H
#define INITIALIZERS_H
#include <array>
#include <type_traits>
#include "randomdir.h"

/***********************************************************************
 * @file initializers.h
 * This file contains two possible default initializers.
 **********************************************************************/

template <class SV, class Enable = void> class SVInitializer;

/** Initializer for Ising-like models.
 * draws a random spin site +1/-1.
 *
 * @tparam IntType the integer type used for storage.
 */
template <typename IntType>
class SVInitializer<std::array<IntType, 1>, typename std::enable_if<
std::is_integral<IntType>::value &&
std::is_signed<IntType>::value>::type>
{
    typedef std::array<IntType, 1> StateVector;
public:
    /** Flip the Ising spin.
     * 
     * @tparam RNGCache the Random Number Generator.
     * 
     * @param svold the old state vector.
     * @param rng a reference to the RNG of the Monte Carlo simulation.
     * @return the new state vector with flipped spin
     */
    template <class RNGCache>
    static StateVector newsv(const StateVector& svold, RNGCache& rng)
    {
        StateVector retval(svold);
        retval[0] = -retval[0];
        return retval;
    }
};

template <typename T>
class SVInitializer<T, typename std::enable_if<std::is_floating_point<typename T::value_type>::value>::type>
{
    using StateVector = T;
public:
    /** Flip the Spin.
     * 
     * @tparam RNGCache the Random Number Generator.
     * 
     * @param svold The old state vector.
     * @param rng A reference to the RNG of the Monte Carlo simulation.
     * @return The new state vector with flipped spin
     */
    template <class RNGCache>
    static StateVector newsv(const StateVector& svold, RNGCache& rng)
    {
        StateVector retval;
        set_to_rnddir(rng, retval);
        return retval;
    }
};

/** Initializer for models with integer spin states -1, 0, 1
*
* @tparam StateVector the type of the state vector.
*/
template <class StateVector>
class Spin1_Initializer
{
	public:
        /** Draw a new state vector based on the old one.
         * 
         * @tparam RNG The Random Number Generator.
         * 
         * @param svold The old state vector.
         * @param rng A reference to the RNG of the Monte Carlo simulation.
         * @return The new state vector
         */
        template <class RNGCache>
        static StateVector newsv(const StateVector& svold, RNGCache& rng)
		{
			StateVector retval(svold);
			int state = retval[0];

			if (state == 0)
			{
				if (rng.real() < 0.5) state = -1;
				else state = +1;
			}
			else // +1/-1
			{
				if (rng.real() < 0.5) state *= -1;
				else state = 0;
			}

			retval[0] = state;
			return retval;
		};
};

/** The generic case that can be specified by a user for a Hamiltonian.
 * 
 * @tparam Hamiltonian The Hamiltonian that you use.
 */
template <class Hamiltonian>
class Initializer : public SVInitializer<typename Hamiltonian::StateVector>
{};

#endif
