#ifndef HEISENBERG_H
#define HEISENBERG_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "vectorhelpers.h"


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class HeisenbergMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			constexpr static int SymD = 3;	// improve me
			const     static int N = grid.size();

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
		HeisenbergMag() : name("m") {}
};


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
		StateVector newsv(StateVector&) 
		{
			return rnddir<RNG, double, SymD>(rng);
		};

	private:
		RNG& rng;
};



template <class StateVector>
class Heisenberg_interaction : public Interaction<StateVector> 
{
	public:
		Heisenberg_interaction()
		{
	 		this->J = -1;	// -1 ferro, +1 antiferro
		}
		StateVector operator() (StateVector& phi) {return phi;};
};


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class Heisenberg
{
	public:

		double beta;
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		
		template <typename RNG>
		using MetroInitializer =  Heisenberg_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		Heisenberg(double mybeta) : beta(mybeta) {   interactions[0] = new Heisenberg_interaction<StateVector>(); }

		HeisenbergMag obs_m;

		auto getobs() { return std::make_tuple(obs_m); }
		

		// using the Wolff cluster algorithm requires to implement 
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A> 
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a)
		{
			return dot(sv1, a) * dot(sv2, a);
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a)
		{
			const double dotp = dot(sv, a);
			for (int i=0; i<SymD; i++) sv[i] -= 2*dotp*a[i];

			normalize(sv);  // necessary?
		}

};
#endif
