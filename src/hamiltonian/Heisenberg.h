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
#include "../obsparts.h"


// ----------------------------------------------------------------------

template <class StateVector, class RNG>
class Heisenberg_Initializer
{
	public:
		// provide the spin dimension as a compile-time constant expression
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		Heisenberg_Initializer()   {}
		Heisenberg_Initializer(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector&) 
		{
			return rnddir<RNG, double, SymD>(rng);
		};

	private:
		RNG& rng;
};



template <class StateVector>
class Heisenberg_interaction
{
	public:
		Heisenberg_interaction(const double myJ) : J(-myJ) {}
		StateVector get (const StateVector& phi) {return phi;};
        const double J;
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename CouplingType=double>
class Heisenberg
{
	public:

		double J;
		constexpr static int SymD = 3;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		
        // the next construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)
		template <typename RNG>
		using MetroInitializer =  Heisenberg_Initializer<StateVector, RNG>; 

		// requires pointers
        std::vector<Heisenberg_interaction<StateVector>*> interactions;
        std::array<OnSite<StateVector, CouplingType>*, 0> onsite;
        std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;

		Heisenberg(double J) : J(J), name("Heisenberg"), observables(obs_m)
        {
            interactions.push_back(new Heisenberg_interaction<StateVector>(J));
        }
		~Heisenberg() {delete interactions[0];} // fixme
		

		// instantiate and choose observables
		Magnetization obs_m;
        std::tuple<Magnetization> observables;

		// state space initializer
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
};






namespace MARQOV
{

	/** Specialization of the Embedding class for the Heisenberg model
	*
	* @tparam SpinType the type in which to store the magnetization values.
	* @tparam CouplingType the type of the coupling of the on-site term (in case there is one)
	*/
	template <class SpinType>
	class Embedder<Heisenberg<SpinType>,RegularHypercubic>
	{
			typedef Heisenberg<SpinType> Hamiltonian;
			typedef RegularHypercubic Lattice;
			typedef typename Hamiltonian::StateVector StateVector;
			typedef StateVector* StateSpace;
			constexpr static int SymD = Hamiltonian::SymD;

		private:

			const Hamiltonian& ham; // why const?
			const Lattice& lat;
			const StateSpace& statespace;

			std::array<SpinType,SymD> rdir;

		public:
			/** Constructs a Heisenberg embedding object.
			*
			* @param ham the corresponding Hamiltonian
			*/
			Embedder(const Hamiltonian& ham, const Lattice& lat, const StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			/** Set new embedding variable.
			*
			* Typically, this function is executed once before every cluster update. The random variable
			* can be drawn randomly (for which an RNG is provided), but of course can also follow some
			* sequential scheme.
			*
			* @tparam RNG the type of the random number generator
			* @param rng reference to the random number generator
			*/
			template <class RNG>
			void draw(RNG& rng)
			{
				rdir = rnddir<RNG, double, SymD>(rng);
			}


			/** Compute the Wolff coupling when attempting to add a spin to the cluster
			*
			* @param currentsv the current state vector (which is already in the cluster)
			* @param candidate its neighbour being checked whether it will become part of the cluster as well
			*
			* @return The scalar Wolff coupling (a double)
			*/
			double coupling(StateVector& currentsv, StateVector& candidate)
			{
				dot(currentsv, rdir) * dot(candidate, rdir);
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




			/* update me

    template <class Lattice, class SpinType, class FPType>
    struct Wolff<Heisenberg<SpinType, FPType>, Lattice>
    {
        template <class RNG, class StateSpace>
        static inline int move(const Heisenberg<SpinType, FPType>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
        {
            typedef typename Heisenberg<SpinType, FPType>::StateVector StateVector;
            std::vector<int> cstack(grid.size(), 0);
            int q = 0;
            ham.wolff_flip(statespace[rsite], rdir);
            cstack[q] = rsite;
            
            int clustersize = 1;
            int current = 0;
            
            while (q>=0)
            {
                current = cstack[q];
                q--;

				// plain Heisenberg model has only one interaction term
				auto coupling = ham.interactions[0]->J; 
				const auto proj1 = coupling*dot(statespace[current], rdir);

				const auto nbrs = grid.nbrs(0, current);
				for (std::size_t i = 0; i < nbrs.size(); ++i)
				{
					const auto currentidx = nbrs[i];
					StateVector& candidate = statespace[currentidx];

					const auto proj2 = dot(candidate, rdir);

					if (proj1*proj2 > 0)
					{
						const auto prob = -std::expm1(-2.0*beta*proj1*proj2);

						if (rng.real() < prob)
						{
							q++;
							cstack[q] = currentidx;
							clustersize++;
							ham.wolff_flip(candidate, rdir);
						}
					}
				}
			}
			return clustersize;
        }
    };
*/
}
#endif
