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

#ifndef ISING_H
#define ISING_H
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

class IsingGenericVectorValuedObs
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid>
		std::vector<double> measure(const StateSpace& statespace, const Grid& grid)
		{
			std::vector<double> retval;

			for (int i=0; i<5; i++) retval.push_back(42+0.1*i);

			return retval;
		}
		IsingGenericVectorValuedObs() : name("dummy"), desc("testing vector-valued observables ...") {}
};




// ------------------------------ HAMILTONIAN ---------------------------

/**
 * Ising Hamiltonian
 * This defines the Ising Hamiltonian. It only consists of a single part,
 * namely the interaction.
 * @tparam SpinType the type in which to store the binary magnetization values.
 */
template <typename SpinType = int>
class Ising
{
	public:
		
		//  ----  Parameters  ----

		double J;
		static constexpr int SymD = 1;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;

		//  ----  Hamiltonian terms  ---- 

        std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};

		Ising(double J) : J(J), name("Ising"), obs_e(*this), obs_fx(0), obs_fy(1)
		{}
		~Ising() {delete interactions[0];}


		//  ----  Observables ----

		Magnetization  obs_m;
		Energy<Ising>  obs_e;
		MagFTComp      obs_fx;
		MagFTComp      obs_fy;
		IsingGenericVectorValuedObs dummy;
        decltype(std::make_tuple(obs_m, obs_e, obs_fx, obs_fy, dummy)) observables = {std::make_tuple(obs_m, obs_e, obs_fx, obs_fy, dummy)};


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (decltype(grid.size()) i = 0; i < grid.size(); i++)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}

};




// ------------------------------ SPECIALIZATIONS ---------------------------

namespace MARQOV 
{
	/** Specialization of the Wolff algorithm for the Ising model.
	* The principle structure is identical to the general version of the algorithm
	* See wolff.h for details
	*/
	template <class Lattice>
	struct Wolff<Ising<int>, Lattice>
	{
		template <class RNG, class StateSpace>
		static inline int move(const Ising<int>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
		 {
			// prepare stack and add seed
            std::vector<int> cstack(grid.size(), 0);
            int q = 0;
            auto& seed = statespace[rsite];
			const int val = seed[0];
            cstack[q] = rsite;
			seed[0] *= -1;
            int clustersize = 1;
            int current;
				
			const auto gcpl = ham.interactions[0]->J; 
			const auto prob = -std::expm1(+2.0*beta*gcpl);

			// work through stack
            while (q>=0)
            {
                current = cstack[q];
                q--;

				// loop over neighbours
				const auto nbrs = grid.nbrs(0, current);
				for (std::size_t i = 0; i < nbrs.size(); ++i)
				{
					const auto currentnbr = nbrs[i];
					auto& candidate = statespace[currentnbr];

					// decide whether to add to cluster
					if (candidate[0] == val)
					{
						if (rng.real() < prob)
						{
							q++;
							cstack[q] = currentnbr;
							clustersize++;
							candidate[0] *= -1;
						}
					}
				}
			}
			return clustersize;
		}
	};

}



#endif
