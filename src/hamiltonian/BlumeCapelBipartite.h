#ifndef BLUMECAPELBIPARTITE_H
#define BLUMECAPELBIPARTITE_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "../hamparts.h"


// ------------------------------ OBSERVABLES ---------------------------

// ...

// ----------------------------------------------------------------------

template <class StateVector>
class BlumeCapelBipartite_interaction : public Interaction<StateVector> 
{
	public:
		BlumeCapelBipartite_interaction(double J)
		{
			this->J = J;
		}
		StateVector get (const StateVector& phi) {return phi;};
};


template <class StateVector>
class BlumeCapelBipartite_onsite : public OnSite<StateVector, double>
{
	public:
		BlumeCapelBipartite_onsite(double D)
		{
			this->h = D;
		}
		double get (const StateVector& phi) {return dot(phi,phi);};
};


template <class StateVector, class RNG>
class BlumeCapelBipartite_Initializer
{
	public:
		BlumeCapelBipartite_Initializer()   {}
		BlumeCapelBipartite_Initializer(RNG& rn) : rng(rn) {}

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
class BlumeCapelBipartite
{
	public:
		double J, DA, DB;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = BlumeCapelBipartite_Initializer<StateVector, RNG>;
        
        std::array<BlumeCapelBipartite_interaction<StateVector>*, 1> interactions = {new BlumeCapelBipartite_interaction<StateVector>(J)};
		static constexpr uint Nbeta = 2;
		static constexpr uint Ngamma = 0;
		
		// instantiate interaction terms (requires pointers)
		OnSite<StateVector, double>* onsite[Nbeta];
		FlexTerm<StateVector*,  StateVector>* multisite[Ngamma];
	
		BlumeCapelBipartite(double J, double DA, double DB) : J(J), DA(DA), DB(DB), name("BlumeCapelBipartite")
		{	
			onsite[0]       = new BlumeCapelBipartite_onsite<StateVector>(DA);		
			onsite[1]       = new BlumeCapelBipartite_onsite<StateVector>(DB);		
		}
		
		~BlumeCapelBipartite(){delete interactions[0]; }

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
