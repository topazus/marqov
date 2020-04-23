#ifndef ASHKINTELLER_H
#define ASHKINTELLER_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"

// 3-color the Ashkin-Teller model (prototype)

// numerically treated as embedded Ising models (compare Zhu et. al 2015)
// at the moment Metropolis is not working! only Wolff!
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
		this->J = -1;	// +1 ferro, -1 antiferro
	}
	StateVector operator() (const StateVector& phi) {return phi;};
};


template <class rStateVector, class RNG>
class AshkinTeller_Initializer
{
	public:
		AshkinTeller_Initializer()   {}
		AshkinTeller_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		rStateVector newsv(const rStateVector& svold) 
		{
			rStateVector retval(svold); 
			retval = -retval;			// if int
//			retval[0] = -retval[0];		// if array<int>
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
//		typedef std::array<SpinType, rSymD> rStateVector; // reduced StateVector
		typedef int rStateVector; // reduced StateVector

		template <typename RNG>
		using MetroInitializer = AshkinTeller_Initializer<rStateVector, RNG>;

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

		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class color = int>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const color a=0)
		{
//			cout << sv1[0] << "  " << sv1[1] << "  " << sv1[2] << endl;
			if (sv1[a] == sv2[a]) return 0.0;
			else
			{
				if      (a==0) return - (J + K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]));
				else if (a==1) return - (J + K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]));
				else if (a==2) return - (J + K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]));
			}

		}

		template <class color = int>
		inline void wolff_flip(StateVector& sv, const color a=0)
		{
			sv[a] *= -1;
		}

		/* 
		rStateVector reduce(StateVector sv, int comp)
		*/

};
#endif

		// obvious improvements that should be made:
		// - unify us ...
		// - move into Hamiltonian and use the StateVector typdefs
		// - should return array<int,1> rather than int ...

		int& reduce_ref(std::array<int,3>& sv, int comp) {return sv[comp];}
		int  reduce_cpy(std::array<int,3>  sv, int comp) {return sv[comp];}

