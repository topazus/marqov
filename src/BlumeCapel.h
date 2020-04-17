#ifndef BLUMECAPEL_H
#define BLUMECAPEL_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class BlumeCapelMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			double mag = 0.0;

			for (int i=0; i<N; i++)
			{
					mag += statespace[i][0];
			}

			return std::abs(mag)/double(N);
		}
		BlumeCapelMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class BlumeCapel_interaction : public Interaction<StateVector> 
{
public:
	BlumeCapel_interaction()
	{
		this->J = -1;	// +1 ferro, -1 antiferro
	}
	StateVector operator() (const StateVector& phi) {return phi;};
};


template <class StateVector>
class BlumeCapel_onsite : public OnSite<StateVector, double>
{
	public:
		BlumeCapel_onsite(double D, double beta)
		{
			this->h = D/beta;
		}
		double operator() (const StateVector& phi) {return dot(phi,phi);};
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
				if (rng.d() < 0.5)  state = -1;
				else				state = +1;
			}
			else // +1/-1
			{
				if (rng.d() < 0.5) state *= -1;
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
		const double D = 0.655;
		double beta;
		constexpr static int SymD = 1;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = BlumeCapel_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 1;
		static constexpr uint Ngamma = 0;
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, double>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		BlumeCapel(double mybeta) : beta(mybeta) 
		{	
			interactions[0] = new BlumeCapel_interaction<StateVector>(); 
			onsite[0]       = new BlumeCapel_onsite<StateVector>(D, beta);		
		}
		
	
		// instantiate and choose observables
		BlumeCapelMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}

		// not yet implemented

		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0)
		{
			if (sv1[0] == 0) return 0.0;
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}


		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0)
		{
			sv[0] *= -static_cast<SpinType>(1.0);
		}

};
#endif
