#ifndef ISING_H
#define ISING_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include <array>
#include "../hamparts.h"
#include "../obsparts.h"
#include "../metropolis.h"
#include "termcollection.h"


// ------------------------------ OBSERVABLES ---------------------------

class IsingGenericVectorValuedObs
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid>
		std::vector<double> measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();
			std::vector<double> retval;

			for (int i=0; i<5; i++) retval.push_back(42+0.1*i);

			return retval;
		}
		IsingGenericVectorValuedObs() : name("dummy"), desc("testing vector-valued observables ...") {}
};

// ----------------------------------------------------------------------

template <class StateVector>
class Ising_interaction
{
public:
	Ising_interaction(const double& myJ) : J(myJ) {}
	StateVector get (const StateVector& phi) {return phi;};
    const double& J;
};


template <class StateVector, class RNG>
class Ising_Initializer
{
	public:
		Ising_Initializer()   {}
		Ising_Initializer(RNG&) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			retval[0] = -retval[0];
			return retval;
		};
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class Ising
{
	public:
		double J;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = Ising_Initializer<StateVector, RNG>;

        // instantiate interaction terms (requires pointers)
        std::array<standard_interaction<StateVector>*, 1> interactions = {new standard_interaction<StateVector>(J)};
//        std::array<Ising_interaction<StateVector>*, 1> interactions = {new Ising_interaction<StateVector>(J)};
        std::array<OnSite<StateVector, int>*, 0> onsite;
        std::array<FlexTerm<StateVector*,  StateVector>*, 0> multisite;


		Ising(double J) : J(J), name("Ising"), obs_e(*this), obs_fx(0), obs_fy(1)
		{}
		~Ising() {delete interactions[0];}

		// instantiate and choose observables
		ScalarMagnetization  obs_m;
		Energy<Ising>		 obs_e;
		ScalarMagFTComp      obs_fx;
		ScalarMagFTComp      obs_fy;
		IsingGenericVectorValuedObs dummy;

		auto getobs()	{return std::make_tuple(obs_m, obs_e, obs_fx, obs_fy, dummy);}


		// initialize state space
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (int i=0; i<grid.size(); i++)
			{
				if (rng.real() > 0.5) statespace[i][0] = 1;
				else statespace[i][0] = -1;
			}
		}


		// using the Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A = bool>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const A a=0) const 
		{
			if (sv1[0] == sv2[0]) return 0.0;
			else return -1.0;
		}

		template <class A = bool>
		inline void wolff_flip(StateVector& sv, const A a=0) const 
		{
			sv[0] *= -1;
		}


};



/* // comment out for debug


namespace MARQOV {

    //Work around GCC Bug for specializations of things in namespaces:
    //https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56480
    // still occurs on 6.4.0

	template <class Lattice>
	struct Metropolis<Ising<int>, Lattice>
	{
	    template <class StateSpace, class M, class RNG>
	    static int move(const Ising<int>& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta, int rsite)
	    {
			typedef typename Ising<int>::StateVector StateVector;
		    	// old state vector at rsite
			StateVector& svold = statespace[rsite];
			// propose new configuration
			StateVector svnew = metro.newsv(svold);
			
			// interaction part
			double interactionenergydiff = 0;
			for(typename std::remove_cv<decltype(ham.Nalpha)> ::type a = 0; a < ham.Nalpha; ++a)
			{
				auto nbrs = grid.nbrs(a, rsite);
				typedef decltype(ham.interactions[a]->get(statespace[0])) InteractionType;
				typedef decltype(MARQOV::callbonds<Lattice>(grid, a, rsite, 0, ham.interactions[a]->get(statespace[0]))) BondType;
				typename MARQOV::Promote_Array<InteractionType, BondType>::CommonArray averagevector = {0};

				// sum over neighbours
				for (std::size_t i = 0; i < nbrs.size(); ++i)
				{
					auto idx = nbrs[i];
					auto nbr = ham.interactions[a]->get(statespace[idx]);
					averagevector = averagevector + MARQOV::callbonds<Lattice>(grid, a, rsite, i, nbr);
				}
				interactionenergydiff += ham.interactions[a]->J * (dot(svnew - svold, averagevector));
			}

	    		// sum up energy differences
	    		double dE 	= interactionenergydiff;
	
	    		 // improve me: what about models with discrete statevectors where the acceptance probability should be
	    		 // looked up in tables? -> specialized Metropolis routine for this case??
	
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


	template <class Lattice>
	struct Wolff<Ising<int>, Lattice>
	{
		template <class DirType, class RNG, class StateSpace>
		static inline int move(const Ising<int>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite, const DirType&)
		{
			typedef typename Ising<int>::StateVector StateVector;
			// prepare stack
			std::vector<int> cstack(grid.size(), 0);
		
			// add initial site and flip it
			int q = 0;
			cstack[q] = rsite;
			const int val = statespace[rsite][0];
			ham.wolff_flip(statespace[rsite]);
			int clustersize = 1;
		
			// compute 'Wolff probability' 
			const int a = 0; // plain Ising model has only one interaction term
			const double coupling = ham.interactions[a]->J;
			const double prob = -std::expm1(+2.0*beta*coupling);
			
			// loop over stack as long as non-empty
			while (q>=0)
			{
				// extract last sv in stack
				const int currentidx = cstack[q];
				q--;
			
				// get its neighbours
				const auto nbrs = grid.nbrs(a, currentidx);
	
				// loop over neighbours
				for (std::size_t i = 0; i < nbrs.size(); ++i)
				{
					// extract corresponding sv
					const auto currentnbr = nbrs[i];
					StateVector& candidate = statespace[currentnbr];
	
					// test whether site is added to the cluster
					if (candidate[0] == val)
					{
						if (rng.real() < prob)
						{
							q++;
							cstack[q] = currentnbr;
							clustersize++;
							ham.wolff_flip(candidate);
						}
					}
				}
			}
	
		return clustersize;
	    }
	};

}

*/


#endif
