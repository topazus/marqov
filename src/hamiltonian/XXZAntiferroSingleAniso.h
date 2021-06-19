#ifndef XXZANTIFERROSINGLEANISO_H
#define XXZANTIFERROSINGLEANISO_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "util/randomdir.h"
#include "util/hamparts.h"


// ------------------------------ OBSERVABLES ---------------------------

#include "XXZAntiferro.h"



// ------------------------------ INITIALIZER ---------------------------

#include "util/initializers.h"



// ------------------------------ HAMILTONIAN ---------------------------

#include "XXZAntiferro.h"


/** Single-site anistropic term for the XXZ model
 * @tparam StateVector the type of the state vectors 
 */
template <class StateVector>
class XXZAntiferro_onsiteaniso : public OnSite<StateVector, double>
{
	public:

		double D; ///< coupling strength of the anisotropic term

		/** Constructor of the single-site anistropic term for the XXZ model
		*
		* @param D coupling strength of the anisotropic term
		*/
		XXZAntiferro_onsiteaniso(double D) : OnSite<StateVector,double>(D) {}
		double get (const StateVector& phi) {return pow(phi[2],2);};
};



/** Hamiltonian of the antiferromagnetic XXZ O(3) model with an additional single-site anisotropy.
  * see XXZAntiferr.h for more documentation
  * @tparam SpinType type in which the spins are stored
  * @param Delta uniaxial exchange anisotropy
  * @param H external field strength
  * @param D single-ion anisotropy
  */
template <typename SpinType>
class XXZAntiferroSingleAniso
{
	public:

		//  ----  Parameters  ----

        double Delta, H, D;
		static constexpr int SymD = 3;
		const std::string name;


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = NVector_Initializer<StateVector, RNG>; 


		//  ----  Hamiltonian terms  ----	

        std::array<XXZAntiferro_interaction<StateVector>*, 1> interactions = {new XXZAntiferro_interaction<StateVector>(Delta)};
        std::vector<OnSite<StateVector, double>*> onsite;

		/** Constructor for the Hamiltonian of the antiferromagnetic XXZ O(3) model with 
		  * an additional single-site anisotropy.
		  *
		  * @param H external field strength 
		  * @param Delta uniaxial exchange anisotropy
		  * @param D single-ion anisotropy
		  */
		XXZAntiferroSingleAniso(double H, double Delta, double D) : Delta(Delta), H(H), D(D), name("XXZAntiferroSingleAniso")
		{
            onsite.push_back(new XXZAntiferro_extfield<StateVector>(H));
            onsite.push_back(new XXZAntiferro_onsiteaniso<StateVector>(D));
		}

		~XXZAntiferroSingleAniso()
        {
            delete onsite[1]; delete onsite[0]; delete interactions[0];
        }


		//  ----  Parameter Names  ----

		std::string paramname(int i) {//A helper function to have nice names for the I/O
            std::string retval;
            switch(i)
            {
                case (0): retval = "H"; break;
                case (1): retval = "Delta"; break;
                case (2): retval = "D"; break;
                default : break;
            }
            return retval;
        }
		
		//  ----  Observables ----

		XXZAntiferroStaggeredMagZ  obs_mstagz;
		XXZAntiferroStaggeredMagXY obs_mstagxy;
        decltype(std::make_tuple(obs_mstagz, obs_mstagxy)) observables = {std::make_tuple(obs_mstagz, obs_mstagxy)};

};
#endif
