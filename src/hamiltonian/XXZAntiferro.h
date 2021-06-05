#ifndef XXZANTIFERRO_H
#define XXZANTIFERRO_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "../vectorhelpers.h"
#include "../hamparts.h"

// ------------------------------ OBSERVABLES ---------------------------

/** Staggered magnetization along the easy axis.
 * Note that the implementation of staggered magnetizations is always
 * characteristic to the lattice being used, in this case RegularHypercubic
*/
class XXZAntiferroStaggeredMagZ
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			const int L = grid.len;

			if ((L+1) % 2 == 0) cout << "[MARQOV] Error: Lattice size must be a multiple of two!" << endl; // TODO

			double magA = 0;
			double magB = 0;

			if (grid.dim == 2)
			{
				cout << "[MARQOV] Error: This is not implemented!" << endl; // TODO
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



/** Staggered magnetization perpendicular to the easy axis.
 * Note that the implementation of staggered magnetizations is always
 * characteristic to the lattice being used, in this case RegularHypercubic
*/
class XXZAntiferroStaggeredMagXY
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			const int L = grid.len;

			if ((L+1) % 2 == 0) cout << "error"<< endl;

			double magAx = 0;
			double magBx = 0;
			double magAy = 0;
			double magBy = 0;

			if (grid.dim == 2)
			{
				cout << "[MARQOV] Error: This is not implemented!" << endl; // TODO
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

			return std::sqrt(pow(magAx-magBx,2) + std::pow(magAy-magBy,2)) / double(0.5*N);
		}

		XXZAntiferroStaggeredMagXY() : name("mstagxy") {}
};



// ------------------------------ INITIALIZER ---------------------------

#include "initializers.h"



// ------------------------------ HAMILTONIAN ---------------------------

/** Anistropic XXZ interaction.
  * interaction strength in z-direction is different from x,y
  */
template <class StateVector>
class XXZAntiferro_interaction
{
	public:
		const double& Delta; // uniaxial exchange anisotropy
		static constexpr double J = 1; // global coupling (TODO)
		XXZAntiferro_interaction(const double& myDelta) : Delta(myDelta) {}
		StateVector get (const StateVector& phi) 
		{
			StateVector retval;

			retval[0] = Delta*phi[0]; // x
			retval[1] = Delta*phi[1]; // y
			retval[2] = phi[2];       // z

			return retval;
		};
};


/** Anisotropic external field.
 * The field only acts along the z-direction
 *
 * @tparam StateVector the type of the state vectors
 */
template <class StateVector>
class XXZAntiferro_extfield : public OnSite<StateVector, double>
{
	public:
		/** Constructor of anisotropic external field term
		*
		* @param H the coupling constant of the field in z-direction
		*/
		XXZAntiferro_extfield(const double& H) : OnSite<StateVector,double>(H) {}
		double get (const StateVector& phi) {return phi[2];};
};



/** Hamiltonian of the antiferromagnetic XXZ O(3) model
  * 
  * @tparam SpinType the type in which to store the vector-valued magnetization values.
  */
template <typename SpinType>
class XXZAntiferro
{
	public:

		//  ----  Parameters  ----

        double Delta, H;
		static constexpr int SymD = 3;
		const std::string name;


		//  ---- Definitions  -----
		
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = NVector_Initializer<StateVector, RNG>; 


		//  ----  Hamiltonian terms  ----

        std::array<XXZAntiferro_interaction<StateVector>*, 1> interactions = {new XXZAntiferro_interaction<StateVector>(Delta)};
        std::array<XXZAntiferro_extfield<StateVector>*, 1> onsite = {new XXZAntiferro_extfield<StateVector>(H)};
        std::array<FlexTerm<Space<StateVector, RegularHypercubic>,  StateVector>*, 0> multisite;

		XXZAntiferro(double Delta, double H) : Delta(Delta), H(H), name("XXZAntiferro") {}


		//  ----  Observables ----
		
		XXZAntiferroStaggeredMagZ  obs_mstagz;
		XXZAntiferroStaggeredMagXY obs_mstagxy;
        decltype(std::make_tuple(obs_mstagz, obs_mstagxy)) observables = {std::make_tuple(obs_mstagz, obs_mstagxy)};

};
#endif
