#ifndef ISING_H
#define ISING_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"

// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class IsingMag
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
		IsingMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
	Ising_interaction(double J)
	{
		this->J = J;
	}
	StateVector operator() (const StateVector& phi) {return phi;};
};


template <class StateVector, class RNG>
class Ising_Initializer
{
	public:
		Ising_Initializer()   {}
		Ising_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			retval[0] = -retval[0];
			return retval;
		};
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class Ising
{
	public:
		double J;
		constexpr static int SymD = 1;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = Ising_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		Ising(double J) : J(J) {	interactions[0] = new Ising_interaction<StateVector>(J); }
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		IsingMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}


		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0)
		{
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}

		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0)
		{
			sv[0] *= -1;
		}


		// think about this ... (see update.h)
		template <typename bond_type>
		inline double wolff_scalarize(const std::vector<bond_type>& bond)
		{
			return bond[0];
		}

	
};
#endif
