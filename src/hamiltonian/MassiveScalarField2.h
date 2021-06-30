#ifndef MASSIVESCALARFIELD2_H
#define MASSIVESCALARFIELD2_H
#include <array>
#include <cmath>
#include <vector>
#include "util/randomdir.h"
#include "util/hamparts.h"
#include "util/termcollection.h" 


// ------------------------------ OBSERVABLES ---------------------------

#include "util/observables.h"


// ------------------------------ HAMILTONIAN --------------------------


/** Action of a massive scalar field
 */
class MassiveScalarField2
{
	public:
		//  ----  Parameters  ----

		int q = 7; // hardcoded!!
		double beta, k, mass, sqrtg;
		static constexpr int SymD = 1;
		const std::string name = "MassiveScalarField2";


		//  ---- Definitions  -----

		typedef std::array<double, SymD> StateVector;
		typedef MARQOV::Space<StateVector,GraphFromCSV> StateSpace;



		//  ----  Hamiltonian terms  ----

		std::array<Standard_Interaction<StateVector>*, 1> interactions 
			= {new Standard_Interaction<StateVector>(-k)};

		std::array<Onsite_Quadratic<StateVector>*, 1> onsite 
			= {new Onsite_Quadratic<StateVector>(0.5*(q*k+sqrtg*mass))};

		MassiveScalarField2(double k, double sqrtg, double mass) : k(k), 
																  mass(mass), 
																  sqrtg(sqrtg),
																  name("MassiveScalarField")
		{
		}



		//  ----  Observables  ----

		Magnetization obs_m;
        decltype(std::make_tuple(obs_m)) observables = {std::make_tuple(obs_m)};



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
				case (0): name = "k"; break;
				case (1): name = "sqrtg"; break;
				case (2): name = "mass"; break;
				default: break;
			}
			return name;
		}


		//  ----  Initializer  ----

		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
//				auto nnbrs = grid.nbrs(0,i).size();
//				if (nnbrs < 7) statespace[i][0] = 1;
//				else statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
				statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
};


// ------------------------------ INITIALIZER ---------------------------

template <>
class Initializer<MassiveScalarField2>
{
	public:
		Initializer() {}

		template <class RNGCache>
		// auto?
		static typename MassiveScalarField2::StateVector newsv(const typename MassiveScalarField2::StateVector& svold, RNGCache& rng)
		{
			double amp = 0.5; // Amplitude (TODO: make me a class parameter)
			double r = rng.real(-1.0, 1.0);

			double oldval = svold[0];
			double newval = oldval + amp*r;
			auto nsv = svold;
			nsv[0] = newval;
			return nsv;
		}
};



#endif
