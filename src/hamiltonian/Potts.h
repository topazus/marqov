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

#ifndef POTTS_H
#define POTTS_H
#include <array>
#include <vector>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include <array>
#include "util/hamparts.h"
#include "util/observables.h"
#include "util/termcollection.h"


// ------------------------------ OBSERVABLES ---------------------------


/** Magnetization of the Q-state Potts model
*
* @tparam Q the number of spin values
*/
template <int Q>
class PottsMagnetization
{
	public:
	std::string name, desc;

	template <class StateSpace, class Grid>
	double measure(const StateSpace& statespace, const Grid& grid)
	{
		const auto N = grid.size();

		// count different states
		std::array<int,Q> magarray;
		for (int i=0; i<Q; i++) magarray[i] = 0;
		for (std::size_t i=0; i<N; i++) magarray[statespace[i][0]]++;

		// find maximum
		int maxelem = *std::max_element(magarray.begin(), magarray.end());

		// definition of the magnetization
		return double(Q*maxelem-N)/double(N*(Q-1));
	}

	PottsMagnetization() : name("m"), desc("Magnetization of the Q-color Potts model") {}

};
			


// ------------------------------ HAMILTONIAN ---------------------------

/** The Potts interaction term
*
* @tparam StateSpace the type of the state space
* @tparam StateVector the type of the state vector
*
* @note This term does not match the required canonical interaction form and his
* hence implented as a Flex Term
*/
template <class StateVector>
class PottsInteraction
{
	public:
		const double k;///< Scalar prefactor
		PottsInteraction(const double k) : k(k) {}

		template <class StateSpace>
		double diff (const int rsite,
					const StateVector& svold,
					const StateVector& svnew,
					const decltype(std::declval<typename StateSpace::Lattice>().nbrs(0,0))& nbrs,
					StateSpace& s)
		{
			double energy_old = 0;
			double energy_new  = 0;

			for (std::size_t i=0; i<nbrs.size(); ++i)
			{
				auto idx = nbrs[i];
				auto nbr = s[idx]; 

				// Dirac delta interaction 
				if (svold[0] == nbr[0]) energy_old++;
				if (svnew[0] == nbr[0]) energy_new++;
			}

			return energy_new-energy_old;
		}
		template <class StateSpace>
		double energy(const StateSpace& s, const typename StateSpace::Lattice& grid, int c) {return 0;}
};


/**
 * Potts Hamiltonian with Q states.
 *
 * @tparam Q defines the number of spin values 0,1,...,Q-1 
 */
template <int Q>
class Potts
{
	public:
		
		//  ----  Parameters  ----

		static constexpr int q = Q; // number of states
		static constexpr int SymD = 1;

		typedef std::array<int, SymD> StateVector;

		//  ----  Hamiltonian terms  ---- 

		// note: since this Hamilitonian has no standard interaction term, a specialization of
		// the Wolff algorithm necessarily has to be provided // FIXME

		// The Potts interaction is not of the standard form and hence categorized as a Flex Term
		PottsInteraction<StateVector> potts_interaction;
		std::array<decltype(potts_interaction)*, 1> multisite{&potts_interaction};

		Potts(double J) : name("Potts"), potts_interaction(J)
		{}
		//  ----  Observables ----

		PottsMagnetization<Q>  obs_m;
		decltype(std::make_tuple(obs_m)) observables = {std::make_tuple(obs_m)};


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (decltype(grid.size()) i = 0; i < grid.size(); i++)
			{
				statespace[i][0] = rng.integer(Q);
			}
		}

};


/** Specialization of the Metropolis initialier for the Q-state Potts model
*
* @tparam Q the number of the spin values
*/
template <int Q>
class Initializer<Potts<Q>>
{
	public:

		template <class RNGCache>
		static typename Potts<Q>::StateVector newsv(typename Potts<Q>::StateVector& svold, RNGCache& rng) 
		{
			typedef typename Potts<Q>::StateVector StateVector;
			StateVector retval(svold);
			retval[0] = rng.integer(Q);
			return retval;
		}
};



namespace MARQOV
{
	/** Specialization of the Embedding class for the Potts model 
	*
	* @tparam Q number of spin values
	* @tparam Lattice the type of the lattice
	*/
	template <int Q, class Lattice>
	class Embedder<Potts<Q>, Lattice>
	{
		typedef Potts<Q> Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef Space<typename Hamiltonian::StateVector, Lattice> StateSpace;
		static constexpr int SymD = Hamiltonian::SymD;

		private:
			const Hamiltonian& ham;
			const Lattice& lat;
			const StateSpace& statespace;

			int newq = -1;
			int oldq;

		public:
			/** Constructs a Potts embedding object
			*
			* @param ham The corresponding Hamiltonian
			* @param lat The corresponding lattice
			* @param statespace The statespace of the simulation
			*/
			Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			/** Set new embedding variable.
			* @tparam RNG the type of the random number generator
			* @param rng reference to the random number generator 
			* @param sv current state of the cluster seed
			*/
			template <class RNG>
			void draw(RNG& rng, StateVector& sv) 
			{
				// store old spin value
				oldq = sv[0];

				// draw new spin state which is different (!) from the old one
				do {newq = rng.integer(Q);}
				while (newq == oldq);
			}

			/** Computes the Wolff coupling when attempting to add a spin to the cluster
			*
			* @param pos1 The position (index) of the current state vector
			* @param pos2 The position (index) of a neighbour being checked whether it will become part of the cluster
			* @return The scalar Wolff coupling (a double)
			*/
			double coupling(int pos1, int pos2) const
			{
				// compare against original q-value of the seed
				if (oldq != statespace[pos2][0]) return 0;
				else return -0.5;
			}

			/** Specifies how a spin flip in the embedded (reduced) model is performed
			*
			* @param sv the spin to flipped
			*/
			void flip(StateVector& sv) const
			{
				sv[0] = newq;
			}
		};


		/** Specialization of the Wolff cluster algorithm for the Q-state Potts model.
		*
		* @tparam Lattice the type of the lattice
		* @tparam Q number of spin values
		*
		* @note this is an example how you can apply a Wolff update to a flex term
		*/
        template <class Lattice, int Q>
        struct Wolff<Potts<Q>, Lattice>
        {
            std::vector<int> cstack = std::vector<int>(4096/sizeof(int), 0);///< the size of the stack is meant to be preserved across different cluster processes.

            template <class RNG, class StateSpace>
            inline int move(const Potts<Q>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
            {
				typedef Potts<Q> Hamiltonian;
				typedef typename Hamiltonian::StateVector StateVector;
				constexpr static int SymD = Hamiltonian::SymD;
                            
				// set up embedder
				Embedder<Hamiltonian, Lattice> embd(ham, grid, statespace);
				embd.draw(rng,statespace[rsite]);
		
		        // add cluster seed and flip it
		        int q = 0;
		        cstack[q] = rsite;
		        int clustersize = 1;
				embd.flip(statespace[rsite]);
		
		        // loop over stack as long as non-empty
		        while (q>=0)
		        {
		            // extract last sv in stack
		            const int currentidx = cstack[q];
		            q--;
		            
		            // get its neighbours
					const int c = 0;
					{
		        		const double gcpl = ham.multisite[c]->k;
		            	auto nbrs = grid.flexnbrs(c, currentidx);
                        if(q + nbrs.size() > cstack.size()) 
                            cstack.resize(2*cstack.size());

		            	// loop over neighbours
		            	for (std::size_t i = 0; i < nbrs.size(); ++i)
		            	{
		            	    // extract corresponding sv
		            	    const auto currentnbr = nbrs[i];
		            	    StateVector& candidate = statespace[currentnbr];
		            	    
							const auto lcpl = embd.coupling(currentidx, currentnbr);
							const double cpl = gcpl*lcpl;
                            
                            if (wolff_update_accepted(cpl, beta, rng))
                            {
                                q++;
                                cstack[q] = currentnbr;
                                clustersize++;
                                embd.flip(candidate);
                            }
		            	}
		            }
		        }
		        return clustersize;
		    }    
                            
        };
}

#endif
