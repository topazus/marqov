#ifndef XXZANTIFERRO_H
#define XXZANTIFERRO_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "vectorhelpers.h"


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class XXZAntiferroMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			constexpr static int SymD = 3;	// improve me
			const     int N = grid.size();

			std::vector<double> mag(SymD,0) ;

			for (int i=0; i<N; i++)
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j];
				}
			}
			
			double retval = 0;
			for (int j=0; j<SymD; j++)
			{
				retval += mag[j]*mag[j];
			}
			return sqrt(retval)/double(N);
		}
		XXZAntiferroMag() : name("m") {}
};

// Magnetization Squared - > some dummy for testing two observables
class XXZAntiferroMagSq
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		int measure(const StateSpace& statespace, const Grid& grid)
		{
			constexpr static int SymD = 3;	// improve me
			const     int N = grid.size();

			std::vector<double> mag(SymD,0) ;

			for (int i=0; i<N; i++)
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j];
				}
			}
			
			double retval = 0;
			for (int j=0; j<SymD; j++)
			{
				retval += mag[j]*mag[j];
			}
			return (int) retval;
		}
		XXZAntiferroMagSq() : name("msq") {}
};


// ----------------------------------------------------------------------

template <class StateVector, class RNG>
class XXZAntiferro_Initializer
{
	public:
		// provide the spin dimension as a compile-time constant expression
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		XXZAntiferro_Initializer()   {}
		XXZAntiferro_Initializer(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(StateVector&) 
		{
			return rnddir<RNG, double, SymD>(rng);
		};

	private:
		RNG& rng;
};



template <class StateVector>
class XXZAntiferro_interaction : public Interaction<StateVector> 
{
	public:
		double Delta = 0.8; // uniaxial exchange anisotropy

		XXZAntiferro_interaction()
		{
	 		this->J = 1;
		}
		StateVector operator() (StateVector& phi) 
		{
			StateVector retval;

			retval[0] = Delta*phi[0];
			retval[1] = Delta*phi[1];
			retval[2] = phi[2];

			return retval;
		};
};


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class XXZAntiferro
{
	public:

		double beta;
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		
		template <typename RNG>
		using MetroInitializer =  XXZAntiferro_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0; //1;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		XXZAntiferro(double mybeta) : beta(mybeta) {   interactions[0] = new XXZAntiferro_interaction<StateVector>(); }
		
		typedef std::tuple<XXZAntiferroMag, XXZAntiferroMagSq> ObsTs;
		XXZAntiferroMag obs_m;
        XXZAntiferroMagSq obs_msq;

		auto getobs() { return std::make_tuple(obs_m, obs_msq); }
		

		// using the Wolff cluster algorithm requires to implement 
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A> 
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a)
		{
			return sv1[2]*sv2[2]; // perform the cluster update only in the z-components
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a)
		{
			sv[2] *= -sv[2];
			normalize(sv);  // necessary?
		}

};
#endif
