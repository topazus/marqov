/** @file BlumeEmeryGriffiths.h
  *
  * @brief The Blume Emergy Griffiths (BEG) Hamiltonian
  *
  * introduced in M. Blume, V. J. Emery, and R. B. Griftiths, Phys. Rev. A 4, 1071 (1971)
  * phase diagram: L. Wang, F. Lee, and J. D. Kimel, Phys. Rev. B 36, 8945 (1987)
  */

#ifndef BLUMEEMERYGRIFFITHS_H
#define BLUMEEMERYGRIFFITHS_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "util/hamparts.h"
#include "util/termcollection.h"


// ------------------------------ OBSERVABLES ---------------------------

#include "util/observables.h"


// ------------------------------ INITIALIZER ---------------------------

#include "util/initializers.h"


// ------------------------------ HAMILTONIAN ---------------------------

template <class StateVector>
class BiquadraticExchangeInteraction
{
	public:
		double k;
		BiquadraticExchangeInteraction(double k) : k(k) {}

		template <class StateSpace>
		double diff (const int rsite,
					const StateVector& svold, 
					const StateVector& svnew, 
					const decltype(std::declval<typename StateSpace::Lattice>().nbrs(0,0))& nbrs, 
					StateSpace& s)
		{
			double neighbourhood = 0;

			for (std::size_t i=0; i<nbrs.size(); ++i)
			{
				auto idx = nbrs[i];
				auto nbr = s[idx];
				neighbourhood = neighbourhood + dot(nbr,nbr);
			}

			auto svdiff = dot(svnew,svnew)-dot(svold,svold);
			return dot(svdiff, neighbourhood);
		}


		template <class StateSpace>
		double energy(const StateSpace& statespace, const typename StateSpace::Lattice& grid, int c)
		{
			const int N = grid.size();
			double retval = 0;

			for (int idx=0; idx<N; idx++)
			{
				auto nbrs = grid.flexnbrs(c, idx);
				auto self = statespace[idx];
				
				for (std::size_t i=0; i<nbrs.size(); ++i)
				{
					auto idx = nbrs[i];
					auto nbr = statespace[idx];

					retval += dot(self,self) * dot(nbr,nbr);
				}
			}
			return 0.5*retval; // account for double counting
		}
};

template <typename SpinType = int>
class BlumeEmeryGriffiths
{
	public:

		//  ----  Parameters  ----

		double J, D, K;
		static constexpr int SymD = 1;
		const std::string name;


		
		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
//		template <typename RNG>
//		using MetroInitializer = Spin1_Initializer<StateVector, RNG>;



		//  ----  Hamiltonian terms  ----

		
		BiquadraticExchangeInteraction<StateVector> biquadratic_exchange_int;

		std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};
		std::array<Onsite_Quadratic<StateVector>*, 1>     onsite       = {new Onsite_Quadratic<StateVector>(D)};
		std::array<decltype(biquadratic_exchange_int)*, 1> multisite {&biquadratic_exchange_int};
	
		BlumeEmeryGriffiths(double J, double D, double K) : J(J), 
															D(D), 
															K(K), 
															name("BlumeCapel"), 
															observables(obs_m),
															biquadratic_exchange_int(K)
							{
#ifdef __PGI
            	//The following is necessary to make PGI-19.10 happy
        MARQOV::Space<StateVector, RegularHypercubic> dummy(10);
        StateVector dummy1;
        std::vector<int> temp; 
        biquadratic_exchange_int.diff(1, dummy1, dummy1, temp, dummy);
#endif
							}

		~BlumeEmeryGriffiths()
		{
			delete interactions[0];
			delete onsite[0];
		}

		//  ----  Observables ----

		Magnetization obs_m;
        std::tuple<Magnetization> observables;


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}
};

template <typename SpinType>
class Initializer<BlumeEmeryGriffiths<SpinType> > : public  Spin1_Initializer<typename BlumeEmeryGriffiths<SpinType>::StateVector > {};


// ------------------------------ SPECIALIZATIONS ---------------------------

namespace MARQOV
{
	/** Specialization of the Embedding class for the Blume Capel model */

	template <class SpinType, class Lattice>
	class Embedder<BlumeEmeryGriffiths<SpinType>,Lattice>
	{
		typedef BlumeEmeryGriffiths<SpinType> Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef MARQOV::Space<typename Hamiltonian::StateVector, Lattice> StateSpace;

		private:
			const Hamiltonian& ham;
			const Lattice& lat;
			StateSpace& statespace;

		public:
			/** Constructs a BEG embedding object.
			 *
			 * @param ham The corresponding Hamiltonian
			 * @param lat The corresponding lattice
			 * @param statespace The statespace of the simulation
			 */
			Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			/** Set new embedding variable.
			  * @note For this specific model this is not necessary, as the embedding is fixed
			  */
			template <class RNG>
			void draw(RNG& rng, StateVector& sv) { } // nothing to draw 


			/** Computes the Wolff coupling when attempting to add a spin to the cluster
			*
			* @param pos1 The position (index) of the current state vector (which is already in the cluster)
			* @param pos2 The position (index) of a neighbour being checked whether it will become part of the cluster as well
			* @return The scalar Wolff coupling (a double)
			*
			* @note For this model, only clusters on the subset of +1/-1 spins are created. Spin-0 site are left untouched. The update will hence not be ergodic!
			*/
			double coupling(int pos1, int pos2) const
			{
				const SpinType s1 = statespace[pos1][0];
				const SpinType s2 = statespace[pos2][0];

				if (s1 == 0) return 0.0; // no cluster if seed state is zero
				if (s1 == s2) return 0.0;
				else return -1.0;
			}


			/** Specifies how a spin flip is performed.
             * 
             * @param sv The statevector that is flipped.
             */
			void flip(StateVector& sv) const {sv[0] = -sv[0];}

	};
}


#endif

