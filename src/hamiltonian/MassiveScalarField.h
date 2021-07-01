#ifndef MASSIVESCALARFIELD_H
#define MASSIVESCALARFIELD_H
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


// Finite Differences Interaction
template <class StateSpace, class StateVector> 
class FinDiffInt
{

	public:
		const double k;
		FinDiffInt(double k) : k(k) {}
		~FinDiffInt() {}

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
		typedef MARQOV::Space<StateVector,HyperbolicRegularFromCSV> StateSpace;



		//  ----  Hamiltonian terms  ----

		std::array<OnSite<StateVector, double>*,1> onsite 
			= {new Onsite_Quadratic<StateVector>(0.5*sqrtg*mass)};
		
		std::array<FinDiffInt<StateSpace,StateVector>*,1> multisite 
			= {new FinDiffInt<StateSpace, StateVector>(0.5*k)};

		MassiveScalarField(double k, double sqrtg, double mass) : k(k), 
																  mass(mass), 
																  sqrtg(sqrtg),
																  name("MassiveScalarField") {}



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
			double boundary_value = 5;

			for(decltype(grid.size()) i = 0; i < grid.size(); ++i)
			{
				// random initialization with mean zero is not a good idea!
				// it can introduce frustration effects to the system

				if (grid.is_boundary_site(i)) statespace[i][0] = boundary_value;
				else 
				{
					statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
					statespace[i][0] += boundary_value;
				}
			}
		}
};


// ------------------------------ INITIALIZER ---------------------------

template <>
class Initializer<MassiveScalarField>
{
	public:
		Initializer() {}

		template <class RNGCache>
		static typename MassiveScalarField::StateVector newsv(const typename MassiveScalarField::StateVector& svold, RNGCache& rng)
		{
			double amp = 0.5;
			double r = rng.real(-1.0, 1.0);

			double oldval = svold[0];
			double newval = oldval + amp*r;
			auto nsv = svold;
			nsv[0] = newval;
			return nsv;
		}
};



namespace MARQOV
{

	template <class Lattice>
	struct Metropolis<MassiveScalarField, Lattice>
	{
	    template <class StateSpace, class RNG>
	    static int move(const MassiveScalarField& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
	    {

			// boundary cells are not being updated
			if (grid.is_boundary_site(rsite)) return 0;

			// typedefs
			typedef MassiveScalarField Hamiltonian;
			typedef typename Hamiltonian::StateVector StateVector;

			// old state vector at index rsite
			StateVector& svold = statespace[rsite];
			// propose new configuration
			StateVector svnew(Initializer<Hamiltonian>::newsv(svold, rng) );

			// no interaction term

			// only one on-site term
			const int b = 0;
			auto diffb = ham.onsite[b]->get(svnew) - ham.onsite[b]->get(svold);
			double onsiteenergydiff = dot(ham.onsite[b]->h, diffb);

			// one flex term
			const int c = 0; 
			auto nbrs = getflexnbrs<Lattice>(grid, c, rsite);
			auto diffc = ham.multisite[c]->diff(rsite, svold, svnew, nbrs, statespace);
			double flexenergydiff = dot(ham.multisite[c]->k, diffc);



	    	// sum up energy differences
	    	double dE 	= onsiteenergydiff + flexenergydiff;
	
	    	int retval = 0;
	    	if ( dE <= 0 )
	    	{
	    	    svold = svnew;
	    	    retval = 1;
	    	}
	    	else if (rng.real() < exp(-beta*dE))
	    	{
	    	    svold = svnew;
	    	    retval = 1;
	    	}
	    	return retval;
	    }
	};



	// write specialization for Wolff in order to get rid of "interactions" in Hamiltonian

	template <class Lattice>
	struct Wolff<MassiveScalarField, Lattice>
	{
	    template <class StateSpace, class RNG>
	    static int move(const MassiveScalarField& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite)
	    {

			return 0;
	    }
	};
}


#endif
