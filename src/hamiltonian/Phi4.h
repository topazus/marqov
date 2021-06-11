#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include <vector>
#include "../vectorhelpers.h"
#include "../hamparts.h"
#include "termcollection.h" 


// ------------------------------ OBSERVABLES ---------------------------

#include "../observables.h"


// ------------------------------ INITIALIZER ---------------------------
// here we have two different choices of initializers

template <class StateVector, class RNG>
class Phi4_Initializer_Radial
{
	public:
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		Phi4_Initializer_Radial(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector& osv) 
		{
			double amp = 0.5; // Amplitude (TODO: make me a class parameter)

			auto newdir = rnddir<RNG, double, SymD>(rng);
			auto nsv = osv + mult(amp,newdir);
			return nsv;
		};

	private:
		RNG& rng;
};


template <class StateVector, class RNG>
class Phi4_Initializer_Cartesian
{
	public:
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		Phi4_Initializer_Cartesian(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector& osv) 
		{
			double amp = 0.5; // Amplitude (TODO: make me a class parameter)

			const int comp =	rng.integer(SymD);
			double r = rng.real(-1.0, 1.0);
			double oldval = osv[comp];
			double newval = oldval + amp*r;
			
			auto nsv = osv;
			nsv[comp] = newval;
			return nsv;
		};

	private:
		RNG& rng;
};


// ------------------------------ HAMILTONIAN --------------------------

/** A possible implementation of a fourth-power on-site term.
* This is specific to the Phi4 model
*
* @tparam StateVector the type of the state vector
*/
template <class StateVector>
class Onsite_Fourth_Minus_One : public OnSite<StateVector, double> 
{
	public:
		Onsite_Fourth_Minus_One(double constant) : OnSite<StateVector,double>(constant) {}
		inline double get (const StateVector& phi) {return pow(dot(phi,phi)-1.0, 2);}
};




/** Phi4 Hamiltonian.
 * This defines a O(N) Hamiltonian with additional mass and fourth-order term
 * @tparam SpinType the type in which to store the magnetization values.
 * @tparam CouplingType the type in which to store the coupling constants of the on-site terms
 */
template <typename SpinType, typename CouplingType=double>
class Phi4
{
	public:
		//  ----  Parameters  ----

		double beta, lambda, mass;
		const double J = -1;
		static constexpr int SymD = 1;
		const std::string name = "Phi4";


		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer =  Phi4_Initializer_Cartesian<StateVector, RNG>; 



		//  ----  Hamiltonian terms  ----

		Standard_Interaction<StateVector> phi4interaction;
		Onsite_Quadratic<StateVector> onsite_standard;
		Onsite_Fourth_Minus_One<StateVector> onsite_fourth_minus_one;

		std::array<Standard_Interaction<StateVector>*, 1>      interactions = {new Standard_Interaction<StateVector>(J)};
		std::vector<OnSite<StateVector, CouplingType>*>        onsite; // empty here, to be filled in the constructor!

		Phi4(double beta, double lambda, double mass) : beta(beta), 
												lambda(lambda), 
												mass(mass), 
												name("Phi4"), 
                                                phi4interaction(J),
												onsite_standard(mass/beta),
												onsite_fourth_minus_one(lambda/beta),
												obs_fx(0), 
												obs_fy(1), 
												obs_fz(2)
		{
			onsite.push_back(&onsite_standard);
			onsite.push_back(&onsite_fourth_minus_one);


			#ifdef __PGI
            	//The following three lines are necessary to make PGI-19.10 happy
            	StateVector dummy;
            	onsite_standard.get(dummy);
            	onsite_fourth_minus_one.get(dummy);
			#endif
		}



		//  ----  Observables  ----

		Magnetization obs_m;
		MagFTComp obs_fx;
		MagFTComp obs_fy;
		MagFTComp obs_fz;
        decltype(std::make_tuple(obs_m, obs_fx, obs_fy, obs_fz)) observables = {std::make_tuple(obs_m, obs_fx, obs_fy, obs_fz)};



		//  ---- Parameter Names ----

		/** Allows to give the Hamiltonian parameter names
		*
		* @param i index of the parameter
		*
		* @return the parameter name (string)
		*/
		std::string paramname(int i)
		{
			std::string name;
			switch (i)
			{
				case (0): name = "beta"; break;
				case (1): name = "lambda"; break;
				case (2): name = "mass"; break;
				default: break;
			}
			return name;
		}

		//  ----  Initializer  ----

#ifdef STATIC_BOUNDARY

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				auto nnbrs = grid.nbrs(0,i).size();
				if (nnbrs < 7) 
				{
					auto coords = grid.crds(i);
					if (coords[0] > 0 && coords[1] > 0) statespace[i][0] = 2;
					else if (coords[0] < -0.5 && coords[1] < 0) statespace[i][0] = 2;
					else statespace[i][0] = -2;
				}
				else statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
#endif

};
#endif
