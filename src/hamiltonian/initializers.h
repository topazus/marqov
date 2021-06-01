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

#ifndef INITIALIZERS_H
#define INITIALIZERS_H
#include <array>
#include <tuple>
#include <string>
#include <functional>

// ------------------------------ INITIALIZER ---------------------------

/** Initializer for O(N)-like models.
 * draws a random vector on the N-dimensional unit shell
 *
 * @tparam StateVector the type of the state vector
 * @tparam RNG the type of the random number generator
 */
template <class StateVector, class RNG>
class NVector_Initializer
{
    public:
        // provide the spin dimension as a compile-time constant expression
        static constexpr int SymD = std::tuple_size<StateVector>::value;

        // constructor
        NVector_Initializer(RNG& rn) : rng(rn) {}

        // generate new statevector
        StateVector newsv(const StateVector&)
        {
            return rnddir<RNG, double, SymD>(rng);
        };

    private:
        RNG& rng;
};


/** Initializer for models with integer spin states -1, 0, 1
*
* @tparam StateVector the type of the state vector
* @tparam RNG the type of the random number generator
*/
template <class StateVector, class RNG>
class Spin1_Initializer
{
	public:
		Spin1_Initializer(RNG& rn) : rng(rn) {}

		StateVector newsv(const StateVector& svold)
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

	private:
		RNG& rng;
};

#endif
