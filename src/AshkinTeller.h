#ifndef ASHKINTELLER_H
#define ASHKINTELLER_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"
#include "metropolis.h"

// the 3-color Ashkin-Teller model

// numerically treated as embedded Ising models (compare Zhu et. al 2015)
// make sure you set up a the EMCS properly


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class AshkinTellerMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int N = grid.size();

			double mag = 0.0;

			for (int i=0; i<N; i++)
			{
				mag += statespace[i][0];
				mag += statespace[i][1];
				mag += statespace[i][2];
			}

			return std::abs(mag)/double(3*N);
		}
		AshkinTellerMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class AshkinTeller_interaction : public Interaction<StateVector> 
{
	public:
		AshkinTeller_interaction(double J)
		{
			this->J = J;
		}
		StateVector operator() (const StateVector& phi) {return phi;};
};

template <class StateVector, class RNG>
class AshkinTeller_Initializer
{
	public:
		AshkinTeller_Initializer(RNG&) {}
		StateVector newsv(const StateVector& svold) 	{cout << "This should not have happened!" << endl;};
};

// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = int>
class AshkinTeller
{
	public:
		int J, K;
		constexpr static int SymD = 3;
		typedef std::array<SpinType, SymD> StateVector;

		template <typename RNG>
		using MetroInitializer = AshkinTeller_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		AshkinTeller(double J, double K) : J(J), K(K)
		{	
			interactions[0] = new AshkinTeller_interaction<StateVector>(J); 
		}
		
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		AshkinTellerMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}

};


// ----------------------- UPDATE SPECIALIZATIONS ----------------------

namespace MARQOV 
{
	// some typedefs
	typedef AshkinTeller<int> myHamiltonian;
	typedef typename myHamiltonian::StateVector StateVector;
	typedef int ReducedStateVector;
	


	// Wolff
	template <class Lattice>
	struct Wolff<myHamiltonian, Lattice>
	{

		// Wolff coupling
		static inline double wolff_coupling(StateVector& sv1, StateVector& sv2, int color, myHamiltonian& ham) 
		{
			if (sv1[color] == sv2[color]) return 0.0;
			else
			{
				switch (color)
				{
					case 0: return ham.J + ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
					case 1: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
					case 2: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
				}
			}
		}
	
		// Wolff flip
		static inline void wolff_flip(StateVector& sv, const int color) {sv[color] *= -1;}
	
	
		// The actual Wolff step
		template <class DirType, class RNG, class StateSpace>
		static inline int move(myHamiltonian& ham, Lattice& grid, StateSpace& statespace, 
						   RNG& rng, double beta, int rsite, const DirType&)
		{
			const int ncolors = 3;
			int clustersize_sum = 0;
			for (int color=0; color<ncolors; color++) // better select a random color
			{
	
				// prepare stack
	     		std::vector<int> cstack(grid.size(), 0);
	
	     		// add initial site and flip it
	     		int q = 0;
	     		cstack[q] = rsite;
	
	     		wolff_flip(statespace[rsite], color);
	     		int clustersize = 1;
	
	     		// loop over stack as long as non-empty
	     		while (q>=0)
	     		{
	         			 // extract last sv in stack
	         			 const int currentidx = cstack[q];
		    			 StateVector& currentsv = statespace[currentidx];
	         			 q--;
	
	         			 // get its neighbours
					const int a = 0; // spin family (hard-coded)
	         			 const auto nbrs = grid.getnbrs(a, currentidx);
	
	         			 // loop over neighbours
	         			 for (std::size_t i = 0; i < nbrs.size(); ++i)
	         			 {
	         			      // extract corresponding sv
	         			      const auto currentnbr = nbrs[i];
	         			      StateVector& candidate = statespace[currentnbr];
	
						const double coupling  = - wolff_coupling(currentsv, candidate, color, ham);
	
	         			      // test whether site is added to the cluster
						if (coupling > 0)
	         			      {
	         			           if (rng.d() < -std::expm1(-2.0*beta*coupling))
	         			           {
	         			                q++;
	         			                cstack[q] = currentnbr;
	         			                clustersize++;
	         			                wolff_flip(candidate, color);
	         			           }
	         			      }
	         			 }
	         		}
				clustersize_sum += clustersize;
			}
			return clustersize_sum/ncolors;
		}
	};
	
	
	
	
	
	
	// Metropolis
	template <class Lattice>
	struct Metropolis<myHamiltonian, Lattice>
	{

		// coupling of the embedded Ising model
		static inline double metro_coupling(StateVector& sv1, StateVector& sv2, int color, myHamiltonian& ham)
		{
			switch (color)
			{
				case 0: return ham.J + ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
				case 1: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
				case 2: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
			}
		}
	
	
		// propose new configuration
		static inline ReducedStateVector metro_newconf(ReducedStateVector& rsv)
		{
			ReducedStateVector retval(rsv);
			retval = -rsv;
			return retval;
		}
	
	
		// flip action
		static inline void metro_flip(StateVector& sv, const int color)
		{
			sv[color] *= -1;
		}
	
	
		// the actual Metropolis move attempt
		template <class StateSpace, class M, class RNG>
		static inline int move(AshkinTeller<int>& ham, Lattice& grid, StateSpace& statespace, 
						   M& metro, RNG& rng, double beta, int rsite)
		{
			int retval = 0;
			for (int color=0; color<3; color++) // better select a random color
			{
		    		// old state vector at rsite
				StateVector&        svold = statespace[rsite];
				ReducedStateVector rsvold = svold[color];
	
				// propose new configuration
				ReducedStateVector rsvnew = metro_newconf(rsvold);
				
				// interaction part
				double interactionenergydiff = 0;
	
				// set interaction family
				const int a = 0;
	
				// extract neighbours
				auto nbrs = grid.getnbrs(a, rsite);
	
				// sum over neighbours
				double averagevector = 0; 
				for (std::size_t i=0; i<nbrs.size(); ++i)
				{
					// neighbour index
					auto idx = nbrs[i];
					// full neighbour
					auto nbr = ham.interactions[a]->operator()(statespace[idx]);
					// reduced neighbour
					auto rnbr = nbr[color];
					// coupling
					auto cpl = metro_coupling(svold, nbr, color, ham);
					// sum
					averagevector = averagevector + mult(cpl,rnbr);
				}
	
				// energy difference
				const double dE = dot(rsvnew - rsvold, averagevector);
				
				if ( dE <= 0 )
				{
					metro_flip(svold,color);
					retval++;
				}
				else if (rng.d() < exp(-beta*dE))
				{
					metro_flip(svold,color);
					retval++;
				}
			}
			return retval;
		}
	};	
}



#endif
