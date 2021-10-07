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

#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <type_traits>
#include <cmath>
#include <cstring>
#include <cinttypes>
#include "rngcache.h"
#include "metropolishelpers.h"
#include "../hamiltonian/util/initializers.h"

// -------------------------- Metropolis Algorithm -----------------------------

namespace MARQOV
{

	/**
	 * A class to encapsulate the Metropolis update.
	 * Using the power of partial template class specializations it is possible to
	 * define moves peculiar to your model.
	 * @tparam Hamiltonian The Hamiltonian used to generate a Metropolis move
	 * @tparam Lattice The lattice used for the neighbourhood relations.
	 */
	template <class Hamiltonian, class Lattice>
	struct Metropolis
	{
		template <class StateSpace, class RNGType>
		static int move(const Hamiltonian& ham, 
		   				const Lattice& grid, 
						StateSpace& statespace, 
						RNGCache<RNGType>& rng, 
						double beta, 
						int rsite);
	};
    namespace detail
    {
        /** calculate 2^(-x) for x in [0,1] with a sixth order polynomial approximation.
        * 
        * Determined with lolremez to minimize the absolut error.
        * 
        * @param x a positive real value in [0,1]
        * @return 2^(-x)
        */
        static inline double exp2app6(double x) noexcept
        {
            double u = 1.0944614513161076e-4;
            u = u * x + 1.2758920456516547e-3;
            u = u * x + 9.5802562122641916e-3;
            u = u * x + 5.5491137793956362e-2;
            u = u * x + 2.4022437948159954e-1;
            u = u * x + 6.9314704935511511e-1;
            return u * x + 9.9999999867786389e-1;
        }

#ifdef __SSE2__
#include <emmintrin.h>
#endif

        /** Truncate a floating value to an integer and truncate in the process.
         * 
         * @param f a double
         * @return an integer representation of the truncated value of f.
         */
        static inline int64_t trunctoi64(double f) noexcept
        {
#ifdef __SSE2__
#if (defined(__PGIC__) && __PGIC__ <= 18)
            return _mm_cvttsd_si64x(_mm_set_sd(f));
#else
            return _mm_cvttsd_si64(_mm_set_sd(f));
#endif
#else
            return (int64_t)std::trunc(f);
#endif
        }
    }

    /** Determine whether the move with the proposed energy change dE is accepted.
     * 
     * @tparam RNGType the type of the underlying RNG.
     * 
     * @param dE the energy change of the proposed move.
     * @param beta the inverse temperature beta
     * @param rng the RNGCache object that we use for drawing random numbers
     * @return true if the move is accepted, else false
    */
    template <class RNGType>
    inline bool update_accepted(double dE, double beta, RNGCache<RNGType>& rng)
    {
        constexpr double log2e =  1.4426950408889634074; /* log_2 e */
        bool accept = false;
        if ( dE <= 0 )
        {
            accept = true;
        }
	else
        {// if not, accept with probability depending on Boltzmann weight
            double rngnum = rng.real();
            double action = -beta * dE;
            
            double action2 = log2e * action;// transform e^x to 2^x
            int64_t i64;
            std::memcpy(&i64, &rngnum, sizeof(rngnum));
            if ((( i64 >>52)  )-1022  != detail::trunctoi64(action2))// if both numbers have different magnitudes.
            {// decide based on the magnitudes. This is likely, hence first in the branch
                accept = (((i64>>52)  )-1022 < detail::trunctoi64(action2));
            }
            else
            {//both numbers are of same magnitude -> evaluate the remainder
                i64 = (i64 & ~0xFFF0000000000000) | 0x3FE0000000000000;

                double remainder = action2-std::trunc(action2);
                std::memcpy(&rngnum, &i64, sizeof(rngnum));
                if ( rngnum < detail::exp2app6(remainder) )
                {
                    accept = true;
                }
            }
        }
	return accept;
    }

	/**
	 * The actual Metropolis move attempt.
	 *
     * The behaviour of the basic local MARQOV step, how to generate a new state vector
     * from an old one, can be set by the initializers, @see initializers.h
	 * @tparam Hamiltonian the type of the Hamiltonian
	 * @tparam Lattice the type of the lattice
	 * @tparam StateSpace the type of the state space
	 * @tparam M the type of the initializer
	 * @tparam RNG the type of the random number generator
	 *
	 * @param ham the Hamiltonian
	 * @param grid the lattice
	 * @param statespace the statespace
	 * @param rng the random number generator
	 * @param beta inverse temperature
	 * @param rsite site to be considered for an update
	 *
	 * @return integer, encoding whether the update was accepted or rejected
	 */
    template <class Hamiltonian, class Lattice>
    template <class StateSpace/*, class M*/, class RNGType>
    int Metropolis<Hamiltonian, Lattice>::move(const Hamiltonian& ham, 
    											const Lattice& grid, 
												StateSpace& statespace, 
												RNGCache<RNGType>& rng, 
												double beta, 
												int rsite)
    {

		// definitions
		typedef typename Hamiltonian::StateVector StateVector;
		typedef typename HasInteractions<Hamiltonian>::type HasInt;
		typedef typename HasFlexTerms<Hamiltonian>::type HasFlex;
		typedef typename HasOnsite<Hamiltonian>::type HasOns;
		
		// old state vector at index rsite
		StateVector& svold = statespace[rsite];
		// propose new configuration
		StateVector svnew(Initializer<Hamiltonian>::newsv(svold, rng) );

		// I. interaction part
		double interactionenergydiff = compute_interactionenergydiff(grid, ham, statespace, svnew, svold, rsite, HasInt());
        
		// II. onsite energy part
		auto onsiteenergydiff = compute_onsiteenergydiff(grid, ham, svnew, svold, rsite, HasOns());

		// III. flex term energy
		auto flexenergydiff = compute_flexenergydiff(grid, ham, statespace, svnew, svold, rsite, HasFlex());


		// collect energy differences
		double dE = interactionenergydiff + onsiteenergydiff + flexenergydiff;
		int retval = 0;
		if (update_accepted(dE, beta, rng))
		{
			svold = svnew;
			retval = 1;
		}
		return retval;
	}
};


#endif
