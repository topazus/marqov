#ifndef PHI4_H
#define PHI4_H
#include <array>
#include <cmath>
#include <vector>
#include "util/randomdir.h"
#include "util/hamparts.h"
#include "util/termcollection.h" 


auto hyperfilter = [](auto p)
{
	auto& mp = std::get<1>(p);
	auto& hp = std::get<2>(p);

	std::string str_repid = std::to_string(mp.repid);
	std::string str_mass  = "mass"+std::to_string(std::get<3>(hp));
	mp.outname = str_mass + "_" + str_repid;

	return p;
};


// ------------------------------ OBSERVABLES ---------------------------

#include "util/observables.h"

class BulkMagnetization
{
	public:
		std::string name; ///< The name of the observable
        std::string desc; ///< A helpful description that will be used in the HDF5 output files.

        /** Construct a bulk magnetization object
         * 
         */
		BulkMagnetization() : name("mb"), desc("bulk magnetization") {}

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
		 	// q used to distinguish boundary and bulk sites; bulk sites have q neighbours, boundary ones less
			const int q = 7; // hard-coded (fixme)

			const/*expr static*/ int SymD = statespace[0].size();
			std::vector<double> mag(SymD, 0);
			int counter = 0;

			for (decltype(grid.size()) i = 0; i < grid.size(); i++) 
			{
				auto nbrs = grid.nbrs(i,0);
				int nnbrs = nbrs.size();
				if (nnbrs == q)
				{
					counter++;
					for (int j=0; j<SymD; j++)
					{
						mag[j] += statespace[i][j]; 
					}
				}
			}

			double retval = 0;
			for (int j=0; j<SymD; j++) retval += mag[j]*mag[j];

			return sqrt(retval)/double(counter);
		}
};


// ------------------------------ HAMILTONIAN --------------------------


template <class StateSpace, class StateVector> 
class FiniteDifferencesInteraction
{

	public:
		const double k;
		FiniteDifferencesInteraction(double k) : k(k) {}
		~FiniteDifferencesInteraction() {}

        double diff (const int rsite,
						const StateVector& svold,
						const StateVector& svnew,
						std::vector<int>& nbrs,
						StateSpace& s) 

		{
			double energy_old = 0;
			double energy_new  = 0;
			for (std::size_t i=0; i<nbrs.size(); ++i)
			{                                     
				auto idx = nbrs[i];                
				auto nbr = s[idx];                  

				energy_old += pow((svold[0]-nbr[0]),2);
				energy_new += pow((svnew[0]-nbr[0]),2);
			}
			return energy_new-energy_old;
		}

};




/** Action of a massive scalar field
 */
class MassiveScalarField
{
	public:
		//  ----  Parameters  ----

		double beta, k, mass, sqrtg;
		static constexpr int SymD = 1;
		const std::string name = "MassiveScalarField";


		//  ---- Definitions  -----

		typedef std::array<double, SymD> StateVector;
		typedef MARQOV::Space<StateVector,GraphFromCSV> StateSpace;



		//  ----  Hamiltonian terms  ----

		Onsite_Quadratic<StateVector> onsite_standard;
		FiniteDifferencesInteraction<StateSpace, StateVector> maininteraction;

		std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(0)};
		std::vector<OnSite<StateVector, double>*>        onsite; // empty here, to be filled in the constructor!
		std::vector<FiniteDifferencesInteraction<StateSpace,StateVector>*> multisite;

		MassiveScalarField(double k, double sqrtg, double mass) : k(k), 
																  mass(-mass), 
																  sqrtg(sqrtg),
																  name("MassiveScalarField"), 
																  onsite_standard(-0.5*sqrtg*mass), 
																  maininteraction(0.5*k)
		{
			onsite.push_back(&onsite_standard);
			multisite.push_back(&maininteraction);
		}



		//  ----  Observables  ----

		Magnetization obs_m;
		BulkMagnetization obs_mb;
        decltype(std::make_tuple(obs_m, obs_mb)) observables = {std::make_tuple(obs_m, obs_mb)};



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
				auto nnbrs = grid.nbrs(0,i).size();
				if (nnbrs < 7) statespace[i][0] = 1;
				else statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
//				statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
};


// ------------------------------ INITIALIZER ---------------------------

//template <class StateVector>
template <>
class Initializer<MassiveScalarField>
{
	public:
		Initializer() {}

		template <class RNGCache>
//		static StateVector newsv(const StateVector& svold, RNGCache& rng)
		static typename MassiveScalarField::StateVector newsv(const typename MassiveScalarField::StateVector& svold, RNGCache& rng)
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
