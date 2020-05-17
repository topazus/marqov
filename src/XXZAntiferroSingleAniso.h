#ifndef XXZANTIFERROSINGLEANISO_H
#define XXZANTIFERROSINGLEANISO_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "vectorhelpers.h"
#include "hamiltonianparts.h"

// ------------------------------ OBSERVABLES ---------------------------

// Staggered magnetization easy axis (z)
class XXZAntiferroSingleAnisoStaggeredMagZ
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
				for (int i=0; i<L; i++)
				{
					for (int j=0; j<L; j++)
					{
						const int idx = i*L + j;

						// sublattice A
						if ((i+j) % 2 == 0) 	
							magA += statespace[idx][2];
						// sublattice B
						else					
							magB += statespace[idx][2];
					}
				}
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

		XXZAntiferroSingleAnisoStaggeredMagZ() : name("mstagz") {}
};



// Staggered magnetization perpendicular to easy axis (x-y plane)
class XXZAntiferroSingleAnisoStaggeredMagXY
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
				for (int i=0; i<L; i++)
				{
					for (int j=0; j<L; j++)
					{
						const int idx = i*L + j;

						// sublattice A
						if ((i+j) % 2 == 0) 	
						{
							magAx += statespace[idx][0];
							magAy += statespace[idx][1];
						}
						// sublattice B
						else					
						{
							magBx += statespace[idx][0];
							magBy += statespace[idx][1];
						}
					}
				}
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

		XXZAntiferroSingleAnisoStaggeredMagXY() : name("mstagxy") {}
};

// ----------------------------------------------------------------------

template <class StateVector, class RNG>
class XXZAntiferroSingleAniso_Initializer
{
	public:
		// provide the spin dimension as a compile-time constant expression
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		XXZAntiferroSingleAniso_Initializer()   {}
		XXZAntiferroSingleAniso_Initializer(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector&) 
		{
			return rnddir<RNG, double, SymD>(rng);
		};

	private:
		RNG& rng;
};



template <class StateVector>
class XXZAntiferroSingleAniso_interaction : public Interaction<StateVector> 
{
	public:
		double Delta; // uniaxial exchange anisotropy

		XXZAntiferroSingleAniso_interaction(double myDelta) : Delta(myDelta)
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
class XXZAntiferroSingleAniso_extfield : public OnSite<StateVector, double>
{
	public:

		double H;

		XXZAntiferroSingleAniso_extfield(double myH) : H(myH)
		{
			this->h = H;
		}
		double operator() (const StateVector& phi) {return phi[2];};
};

template <class StateVector>
class XXZAntiferroSingleAniso_onsiteaniso : public OnSite<StateVector, double>
{
	public:

		double D;

		XXZAntiferroSingleAniso_onsiteaniso(double myD) : D(myD)
		{
			this->h = D;
		}
		double operator() (const StateVector& phi) {return pow(phi[2],2);};
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class XXZAntiferroSingleAniso
{
	public:
        	double Delta, H, D, id;
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		
		template <typename RNG>
		using MetroInitializer =  XXZAntiferroSingleAniso_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta  = 2;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		XXZAntiferroSingleAniso(double id, double myH, double myDelta, double myD) : Delta(myDelta), H(myH), D(myD), id(id)
		{
			interactions[0] = new XXZAntiferroSingleAniso_interaction<StateVector>(Delta); 
			onsite[0]       = new XXZAntiferroSingleAniso_extfield<StateVector>(H);
			onsite[1]       = new XXZAntiferroSingleAniso_onsiteaniso<StateVector>(D);
		}
		
		XXZAntiferroSingleAnisoStaggeredMagZ  obs_mstagz;
		XXZAntiferroSingleAnisoStaggeredMagXY obs_mstagxy;

		auto getobs() { return std::make_tuple(obs_mstagz, obs_mstagxy); }
		

		// using the Wolff cluster algorithm requires to implement 
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A> 
		inline auto wolff_coupling(StateVector& sv1, StateVector& sv2, const A a) const
		{
			return sv1[2]*sv2[2]; // perform the cluster update only in the z-components
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a) const
		{
			sv[2] = -sv[2];
			normalize(sv);  // necessary?
		}
};
#endif
