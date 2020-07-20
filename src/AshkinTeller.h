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

			double mag1 = 0.0;
			double mag2 = 0.0;
			double mag3 = 0.0;

			for (int i=0; i<N; i++)
			{
				mag1 += statespace[i][0];
				mag2 += statespace[i][1];
				mag3 += statespace[i][2];
			}

			return (std::abs(mag1)+std::abs(mag2)+std::abs(mag3))/double(3*N);
		}
		AshkinTellerMag() : name("m") {}
};


// ----------------------------------------------------------------------

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
		double J, K;
		constexpr static int SymD = 3;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		typedef AshkinTeller<int> myHamiltonian;
		template <typename RNG>
		using MetroInitializer = AshkinTeller_Initializer<StateVector, RNG>;

		AshkinTeller(double J, double K) : J(J), K(K), name("AshkinTeller") {}
		
	
		// instantiate and choose observables
		AshkinTellerMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
		}


		// init
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (int i=0; i<grid.size(); i++)
			{
				for (int j=0; j<SymD; j++)
				{
					if (rng.real() > 0.5) statespace[i][j] = 1;
					else statespace[i][j] = -1;
				}
			}
		}

};


// ----------------------- UPDATE SPECIALIZATIONS ----------------------

namespace MARQOV 
{

	// Wolff
	template <class Lattice>
	struct Wolff<AshkinTeller<int>, Lattice>
	{

		// some typedefs
		typedef typename AshkinTeller<int>::StateVector StateVector;
		typedef int ReducedStateVector;

		// Wolff coupling
		static inline double wolff_coupling(StateVector& sv1, StateVector& sv2, int color, AshkinTeller<int>& ham) 
		{
			if (sv1[color] == sv2[color]) return 0.0;
			else
			{
				switch (color)
				{
					case 0: return ham.J + ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
					case 1: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
					case 2: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
					default: throw std::invalid_argument("invalid color!");
				}
			}
			return 0;
		}
	
		// Wolff flip
		static inline void wolff_flip(StateVector& sv, const int color) {sv[color] *= -1;}
	
	
		// The actual Wolff step
		template <class DirType, class RNG, class StateSpace>
		static inline int move(AshkinTeller<int>& ham, Lattice& grid, StateSpace& statespace, 
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
	         			           if (rng.real() < -std::expm1(-2.0*beta*coupling))
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
	struct Metropolis<AshkinTeller<int>, Lattice>
	{

		// some typedefs
		typedef typename AshkinTeller<int>::StateVector StateVector;
		typedef int ReducedStateVector;


		// coupling of the embedded Ising model
		static inline double metro_coupling(StateVector& sv1, StateVector& sv2, int color, AshkinTeller<int>& ham)
		{
			switch (color)
			{
				case 0: return ham.J + ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
				case 1: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
				case 2: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
				default: throw std::invalid_argument("invalid color!");
			}
			return 0;
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
					auto nbr = statespace[idx];
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
				else if (rng.real() < exp(-beta*dE))
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
