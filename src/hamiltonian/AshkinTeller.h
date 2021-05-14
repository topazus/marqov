
/**
 * @file AshkinTeller.h
 * @brief Three-Color Ashkin-Teller Hamiltonian.
 *
 * numerically, this model is treated as embedded Ising models
 *
 * the Wolff algorithm can only be used for K=<0.5!
 *
 * compare Zhu et. al, PRB 91, 224201 (2015) for more details
 *
 * 
 * Note: This Hamiltonian presents an "extrem" use case of MARQOV, where the generic term structure is not used at all! Therefore the general Metropolis and Wolff update algorithms can not be used! Instead the interactions of the model are explicitely coded in the specialized update algorithms below
 */


#ifndef ASHKINTELLER_H
#define ASHKINTELLER_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
//#include "../hamparts.h" // not needed, see above
#include "../metropolis.h"



// ------------------------------ OBSERVABLES ---------------------------

/**
* @brief Magnetization of the Three-Color Ashkin-Teller model
*/
class AshkinTellerMag
{
	public:
		std::string name, desc;
		/**
		* @brief Perform a measurement of the Magnetization
		*
		* @tparam StateSpace	the type of the state space 
		* @tparam Grid the type of the lattice
		* @param statespace the statespace
		* @param grid the lattice
		*
		* @return A scalar value, the magnetization per site
		*/
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const auto N = grid.size();

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
		/**
		* @brief Constructor of the Ashkin-Teller magnetization
		*/
		AshkinTellerMag() : name("m"), desc("Magnetization of the Three-color Ashkin-Teller model") {}
};


// ----------------------------------------------------------------------

/** Dummy function which prints an error if one attempts to use the general Metropolis algorithm.
 *
 * @todo Find a way that this function is not needed in the first place
 *
 * @tparam StateVector the type of the state vector
 * @tparam RNG the type of the random number generator
 */
template <class StateVector, class RNG>
class AshkinTeller_Initializer
{
	public:
		/** Constructor */
		AshkinTeller_Initializer(RNG&) {}

		/** Usually specifies how a random new state vector is generated. 
		  * In this model it is just a place holder and will not be needed! 
		  */
		StateVector newsv(const StateVector& svold) 	
		{
			cout << "This should not have happened!" << endl;
			cout << "This Hamiltonian is not supposed to work with the generic Metropolis algorithm." << endl;
			cout << "Use a specialized implementation instead!" << endl;
		};
};




// ------------------------------ HAMILTONIAN ---------------------------


/**
 * Three-Color Ashkin-Teller Hamiltonian
 *
 * @tparam SpinType the type in which to store the binary magnetization values.
 */

template <typename SpinType = int>
class AshkinTeller
{
	public:
		//  ----  Parameters  ----

		double J; ///< Ising interaction strength
		double K; ///< Four-spin interaction strength

		static constexpr int SymD = 3; ///< use SymD to encode the three colors of the model
		const std::string name;



		//  ---- Definitions  -----

		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = AshkinTeller_Initializer<StateVector, RNG>;

		/** Constructor
		 *  
		 * @param J the Ising interaction
		 * @param K the Four-spin interaction
		 */
		AshkinTeller(double J, double K) : J(J), K(K), name("ThreeColorAshkinTeller"), observables(obs_m) {}
		

	
		//  ----  Observables ----

		AshkinTellerMag       obs_m;
		decltype(std::make_tuple(obs_m)) observables = {std::make_tuple(obs_m)};
	



		//  ----  Initializer  ----

		/** Specifies how the state space is initialized
		*
		* @tparam StateSpace 	the type of the state space
		* @tparam Lattice 		the type of the latticie
		* @tparam RNG				the type of the random number generator
		* @param statespace		the state space
		* @param grid				the lattice
		* @param rng				the random number generator
		*/
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for (decltype(grid.size()) i = 0; i < grid.size(); i++)
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

	/** Specialized Wolff algorithm for the Ashkin-Teller Hamiltonian
	* 
	*
	* @tparam Lattice the type of the lattice
	*/
	template <class Lattice>
	struct Wolff<AshkinTeller<int>, Lattice>
	{

		// some typedefs
		typedef typename AshkinTeller<int>::StateVector StateVector;
		typedef int ReducedStateVector;

		// Wolff coupling
		static inline double wolff_coupling(StateVector& sv1, StateVector& sv2, int color, AshkinTeller<int>& ham) 
		{
			double retval = 0.0;
			if (sv1[color] != sv2[color])
			{
				switch (color)
				{
					case 0: retval = ham.J - ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]); break;
					case 1: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]); break;
					case 2: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]); break;
					default: cout << "invalid color!" << color << endl;
//					default: throw std::invalid_argument("invalid color!"); // catch me!
				}
			}
			return retval;
		}
	
		// Wolff flip
		static inline void wolff_flip(StateVector& sv, const int color) {sv[color] *= -1;}
	
	
		// The actual Wolff step
		template <class DirType, class RNG, class StateSpace>
		static inline int move(AshkinTeller<int>& ham, 
						   Lattice& grid, 
						   StateSpace& statespace, 
						   RNG& rng, 
						   double beta, 
						   int rsite, 
						   const DirType&)
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
	         			 const auto nbrs = grid.nbrs(a, currentidx);
	
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
	
	
	
	
	
	
	/** Specialized Metropolis algorithm for the Ashkin-Teller Hamiltonian
	 *
	 * @tparam Lattice type of the lattice
	 */
	template <class Lattice>
	struct Metropolis<AshkinTeller<int>, Lattice>
	{

		// some typedefs
		typedef typename AshkinTeller<int>::StateVector StateVector;
		typedef int ReducedStateVector;


		// coupling of the embedded Ising model
		static inline double metro_coupling(StateVector& sv1, StateVector& sv2, int color, AshkinTeller<int>& ham)
		{
			double retval = 0.0;
			switch (color)
			{
				case 0: retval = ham.J - ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]); break;
				case 1: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]); break;
				case 2: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]); break;
				default: cout << "invalid color!" << color << endl;
//				default: throw std::invalid_argument("invalid color!"); // catch me!
			}
			return retval;
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
		static inline int move(AshkinTeller<int>& ham, 
						   Lattice& grid, 
						   StateSpace& statespace, 
						   M& metro, 
						   RNG& rng, 
						   double beta, 
						   int rsite)
		{
			int retval = 0;
			for (int color=0; color<3; color++) // better select a random color
			{
		    		// old state vector at rsite
				StateVector&        svold = statespace[rsite];
				ReducedStateVector rsvold = svold[color];
	
				// propose new configuration
				ReducedStateVector rsvnew = metro_newconf(rsvold);
				
				// set interaction family
				const int a = 0;
	
				// extract neighbours
				auto nbrs = grid.nbrs(a, rsite);
	
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
