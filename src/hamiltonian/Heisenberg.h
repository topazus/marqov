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

#ifndef HEISENBERG_H
#define HEISENBERG_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include <array>
#include <vector>
#include "../vectorhelpers.h"
#include "../hamparts.h"
#include "termcollection.h"


// ------------------------------ OBSERVABLES ---------------------------

#include "../observables.h"


// ------------------------------ INITIALIZER ---------------------------

#include "initializers.h"


// ------------------------------ HAMILTONIAN ---------------------------


template <class StateVector>
class Heisenberg_interaction
{
	public:
		Heisenberg_interaction(const double myJ) : J(-myJ) {}
		StateVector get (const StateVector& phi) {return phi;};
        const double J;
};


/** Heisenberg Hamiltonian.
* This defines the Heisenberg Hamiltonian. It only consists of a single part,
* namely the standard interaction.
*
* @tparam SpinType the type in which to store the vector-valued magnetization values.
* @tparam CouplingType the type in which the coupling of the on-site term would be stored (if there was one)
*/
template <typename SpinType, typename CouplingType=double>
class Heisenberg
{
	public:

		//  ----  Parameters  ----

		double J;
		static constexpr int SymD = 3;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer =  NVector_Initializer<StateVector, RNG>; 


		//  ----  Hamiltonian terms  ----

        std::vector<Heisenberg_interaction<StateVector>*> interactions;
        std::array<OnSite<StateVector, CouplingType>*, 0> onsite;
        std::array<FlexTerm<Space<StateVector, RegularHypercubic>,  StateVector>*, 0> multisite;

		Heisenberg(double J) : J(J), name("Heisenberg"), observables(obs_m)
        {
            interactions.push_back(new Heisenberg_interaction<StateVector>(J));
        }
//		~Heisenberg() {delete interactions[0];} // fixme
		

		//  ----  Observables ----

		Magnetization obs_m;
        std::tuple<Magnetization> observables;


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
};




// ------------------------------ SPECIALIZATIONS ---------------------------


namespace MARQOV
{

	/** Specialization of the Embedding class for the Heisenberg model
	*
	* @tparam SpinType the type in which to store the magnetization values.
	* @tparam CouplingType the type of the coupling of the on-site term (in case there is one)
	*/

	template <class SpinType, class CouplingType, class Lattice>
	class Embedder<Heisenberg<SpinType, CouplingType>, Lattice>
	{
		typedef Heisenberg<SpinType,CouplingType> Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;
		static constexpr int SymD = Hamiltonian::SymD;

		private:

			const Hamiltonian& ham;
			const Lattice& lat;
			const StateSpace& statespace;

			std::array<SpinType,SymD> rdir;

		public:
			/** Constructs a Heisenberg embedding object.
			*
			* @param ham The corresponding Hamiltonian
			* @param lat The corresponding lattice
			* @param statespace The statespace of the simulation
			*/
			Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			/** Set new embedding variable.
			*
			* Typically, this function is executed once before every cluster update. The variable
			* can be drawn randomly (for which case an RNG is provided), but of course can also follow
			* some sequential scheme.
			*
			* @tparam RNG the type of the random number generator
			* @param rng reference to the random number generator
			*/
			template <class RNG>
			void draw(RNG& rng)	{rdir = rnddir<RNG, double, SymD>(rng);}


			/** Computes the Wolff coupling when attempting to add a spin to the cluster
			*
			* @param pos1 The position (index) of the current state vector (which is already in the cluster)
			* @param pos2 The position (index) of a neighbour being checked whether it will become part of the cluster as well
			*
			* @return The scalar Wolff coupling (a double)
			*/
			double coupling(int pos1, int pos2) const
			{
				return dot(statespace[pos1], rdir) * dot(statespace[pos2], rdir);
			}



			/** Specifies how a spin flip in the embedded (reduced) model is performed
			*
			* @param sv the spin to flipped
			*/
			void flip(StateVector& sv)
			{
				const double dotp = dot(sv, rdir);
				for (int i=0; i<SymD; i++) sv[i] -= 2*dotp*rdir[i];
				normalize(sv);  // necessary?
			}
	};




    template <class Lattice, class SpinType, class CouplingType>
    struct Wolff<Heisenberg<SpinType, CouplingType>, Lattice>
    {
        template <class RNG, class StateSpace>
        static inline int move(const Heisenberg<SpinType, CouplingType>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
        {
			typedef Heisenberg<SpinType, CouplingType> Hamiltonian;
            typedef typename Hamiltonian::StateVector StateVector;
			constexpr static int SymD = Hamiltonian::SymD;

            std::vector<int> cstack(grid.size(), 0);
            int q = 0;
            auto& seed = statespace[rsite];
            cstack[q] = rsite;

			int rdir = rng.integer(SymD); // random Cartesian direction
			seed[rdir] *= -1;
            
            int clustersize = 1;
            int current = 0;
				
			// plain Heisenberg model has only one interaction term
			const auto gcpl = ham.interactions[0]->J; 

            while (q>=0)
            {
                current = cstack[q];
                q--;

				const auto nbrs = grid.nbrs(0, current);
				for (std::size_t i = 0; i < nbrs.size(); ++i)
				{
					const auto currentnbr = nbrs[i];
					StateVector& candidate = statespace[currentnbr];

					const auto lcpl = statespace[current][rdir] * candidate[rdir];
					const auto cpl  = gcpl*lcpl;

					if (cpl > 0)
					{
						if (rng.real() < -std::expm1(-2.0*beta*cpl))
						{
							q++;
							cstack[q] = currentnbr;
							clustersize++;
							candidate[rdir] *= -1;
						}
					}
				}
			}
			return clustersize;
        }
    };
}
#endif
