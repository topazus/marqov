#ifndef XXZANTIFERRO_H
#define XXZANTIFERRO_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "vectorhelpers.h"
#include "hamiltonianparts.h"

// ------------------------------ OBSERVABLES ---------------------------

// Staggered magnetization easy axis (z)
class XXZAntiferroStaggeredMagZ
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			const int L = grid.length;

			if ((L+1) % 2 == 0) cout << "error"<< endl;

			double magA = 0;
			double magB = 0;

			if (grid.dim == 2)
			{
				// implement me
			}

			if (grid.dim == 3)
			{
				for (int i=0; i<L; i++)
				{
					for (int j=0; j<L; j++)
					{
						for (int k=0; k<L; k++)
						{
							const int idx = i*L*L + j*L + k;

							// sublattice A
							if ((i+j+k) % 2 == 0) 	
								magA += statespace[idx][2];
							// sublattice B
							else					
								magB += statespace[idx][2];
						}
					}
				}
			}

			return fabs(magA-magB) / double(0.5*N);
		}

		XXZAntiferroStaggeredMagZ() : name("mstagz") {}
};



// Staggered magnetization perpendicular to easy axis (x-y plane)
class XXZAntiferroStaggeredMagXY
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			const int L = grid.length;

			if ((L+1) % 2 == 0) cout << "error"<< endl;

			double magAx = 0;
			double magBx = 0;
			double magAy = 0;
			double magBy = 0;

			if (grid.dim == 2)
			{
				// implement me
			}

			if (grid.dim == 3)
			{
				for (int i=0; i<L; i++)
				{
					for (int j=0; j<L; j++)
					{
						for (int k=0; k<L; k++)
						{
							const int idx = i*L*L + j*L + k;

							if ((i+j+k) % 2 == 0) // sublattice A
							{
								magAx += statespace[idx][0];
								magAy += statespace[idx][1];
							}
							else // sublattice B
							{
								magBx += statespace[idx][0];
								magBy += statespace[idx][1];
							}
						}
					}
				}
			}

			return std::sqrt(pow(magAx-magBx,2) + pow(magAy-magBy,2)) / double(0.5*N);
		}

		XXZAntiferroStaggeredMagXY() : name("mstagxy") {}
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
		StateVector newsv(const StateVector&) 
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
		StateVector operator() (const StateVector& phi) 
		{
			StateVector retval;

			retval[0] = Delta*phi[0];
			retval[1] = Delta*phi[1];
			retval[2] = phi[2];

			return retval;
		};
};

template <class StateVector>
class XXZAntiferro_extfield : public OnSite<StateVector, double>
{
	public:
		XXZAntiferro_extfield()
		{
			this->h = -4.2;
		}
		double operator() (const StateVector& phi) {return phi[2];};
};


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class XXZAntiferro
{
	public:
        double j;
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		
		template <typename RNG>
		using MetroInitializer =  XXZAntiferro_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta  = 1;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		XXZAntiferro(double myj) : j(myj)
		{
			interactions[0] = new XXZAntiferro_interaction<StateVector>(); 
			onsite[0]       = new XXZAntiferro_extfield<StateVector>();
		}
		
		XXZAntiferroStaggeredMagZ  obs_mstagz;
		XXZAntiferroStaggeredMagXY obs_mstagxy;

		auto getobs() { return std::make_tuple(obs_mstagz, obs_mstagxy); }
		

		// using the Wolff cluster algorithm requires to implement 
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A> 
		inline auto wolff_coupling(StateVector& sv1, StateVector& sv2, const A a)
		{
			return sv1[2]*sv2[2]; // perform the cluster update only in the z-components
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a)
		{
			sv[2] = -sv[2];
			normalize(sv);  // necessary?
		}
};
#endif
