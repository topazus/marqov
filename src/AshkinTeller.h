#ifndef ASHKINTELLER_H
#define ASHKINTELLER_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"

// 3-color the Ashkin-Teller model (prototype)

// numerically treated as embedded Ising models (compare Zhu et. al 2015)
// make sure you set up a the EMCS properly

// todo: pass Hamiltonian parameters through the constructor



// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class AshkinTellerMag
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
					mag += statespace[i][1];
					mag += statespace[i][2];
			}

			return std::abs(mag)/double(3*N);
		}
		AshkinTellerMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class AshkinTeller_interaction : public Interaction<StateVector> 
{
public:
	AshkinTeller_interaction()
	{
		this->J = -1;
	}
	StateVector operator() (const StateVector& phi) {return phi;};
};


template <class redStateVector, class RNG>
class AshkinTeller_Initializer
{
	public:
		AshkinTeller_Initializer()   {}
		AshkinTeller_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		redStateVector newsv(const redStateVector& svold) 
		{
			redStateVector retval(svold); 
			retval = -retval;			// 
//			retval[0] = -retval[0];		// if redStateVector was array<int>
			return retval;
		};
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class AshkinTeller
{
	public:
		double J = 1;
		double K = 0.5;

		constexpr static int SymD = 3;
		constexpr static int rSymD = 1;
		typedef std::array<SpinType, SymD> StateVector;
		typedef int redStateVector; // reduced StateVector

		template <typename RNG>
		using MetroInitializer = AshkinTeller_Initializer<redStateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		AshkinTeller()  {	interactions[0] = new AshkinTeller_interaction<StateVector>(); }
//		AshkinTeller(double myj) : j(myj) {	interactions[0] = new AshkinTeller_interaction<StateVector>(); }
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		AshkinTellerMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}



		// we need to use the filtered Metropolis method in order to access the 
		// embedded Ising models

		template <class color = int>
		inline double metro_coupling(StateVector& sv1, StateVector& sv2, const color a=0)
		{
			switch (a)
			{
				case 0: return J + K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
				case 1: return J + K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
				case 2: return J + K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
			}
		}

		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class color = int>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const color a=0) const
		{
			if (sv1[a] == sv2[a]) return 0.0;
			else return - metro_coupling(sv1, sv2, a);
		}

		template <class color = int>
		inline void wolff_flip(StateVector& sv, const color a=0) const
		{
			sv[a] *= -1;
		}

};

		// filter functions
		// obvious improvements that should be made:
		// - unify us ...
		// - move into Hamiltonian and use the StateVector typdefs
		// - should return array<int,1> rather than int ...

		int& reduce_ref(std::array<int,3>& sv, int comp) {return sv[comp];}
		int  reduce_cpy(std::array<int,3>  sv, int comp) {return sv[comp];}

#endif
