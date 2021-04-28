#ifndef BLUMECAPEL_H
#define BLUMECAPEL_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "../hamparts.h"
// ----------------------------------------------------------------------

template <class StateVector>
class BlumeCapel_interaction
{
    public:
        const double& J;
		BlumeCapel_interaction(const double& myJ) : J(myJ) {}
		StateVector get (const StateVector& phi) {return phi;}
};


template <class StateVector>
class BlumeCapel_onsite
{
    public:
        const double& h;
		BlumeCapel_onsite(double& D) : h(D) {}
		double get (const StateVector& phi) {return dot(phi,phi);};
};


template <class StateVector, class RNG>
class BlumeCapel_Initializer
{
	public:
		BlumeCapel_Initializer()   {}
		BlumeCapel_Initializer(RNG& rn) : rng(rn) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 

			int state = retval[0];

			if (state == 0)
			{
				if (rng.real() < 0.5)  state = -1;
				else				state = +1;
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


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class BlumeCapel
{
	public:
		double J, D;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = BlumeCapel_Initializer<StateVector, RNG>;

		// instantiate interaction terms (requires pointers)
        std::array<BlumeCapel_interaction<StateVector>*, 1> interactions = {new BlumeCapel_interaction<StateVector>(J)};
        std::array<BlumeCapel_onsite<StateVector>*, 1> onsite = {new BlumeCapel_onsite<StateVector>(D)};
        std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;
	
		BlumeCapel(double J, double D) : J(J), D(D), name("BlumeCapel") {}

		// instantiate and choose observables
		Magnetization obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}


		// state space initializer
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(int i=0; i<grid.size(); ++i)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}



		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0) const
		{
			if (sv1[0] == 0) return 0.0;
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}


		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0) const
		{
			sv[0] *= -static_cast<SpinType>(1.0);
		}

};
#endif
