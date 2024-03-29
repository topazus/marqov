
/**
 * @file AshkinTellerThreeColor.h
 * @brief Three-Color Ashkin-Teller Hamiltonian.
 *
 * Numerically, this model is treated as embedded Ising models
 *
 * The Wolff algorithm can only be used for K=<0.5!
 *
 * Compare Zhu et. al, PRB 91, 224201 (2015) for more details
 *
 * 
 * Note: This Hamiltonian presents an "extreme" use case of MARQOV, where the generic term structure is not used at all! Therefore the general Metropolis and Wolff update algorithms can not be used! Instead the interactions of the model are explicitely coded in the specialized update algorithms below
 */


#ifndef ASHKINTELLER_H
#define ASHKINTELLER_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
//#include "util/hamparts.h" // not needed, see above
#include "util/termcollection.h"
#include "../libmarqov/metropolis.h"



// ------------------------------ OBSERVABLES ---------------------------

/**
* @brief Magnetization of the Three-Color Ashkin-Teller model
*/
class AshkinTellerThreeColorMag
{
	public:
		std::string name{"m"};///< We call ourselves m.
		std::string desc{"Magnetization of the Three-color Ashkin-Teller model"};///< An extended description of the observable.
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

			for (std::size_t i=0; i<N; i++)
			{
				mag1 += statespace[i][0];
				mag2 += statespace[i][1];
				mag3 += statespace[i][2];
			}

			return (std::abs(mag1)+std::abs(mag2)+std::abs(mag3))/double(3*N);
		}
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
			std::cout << "This should not have happened!" << std::endl;
			std::cout << "This Hamiltonian is not supposed to work with the generic Metropolis algorithm." << std::endl;
			std::cout << "Use a specialized implementation instead!" << std::endl;
		};
};




// ------------------------------ HAMILTONIAN ---------------------------

/** Three-Color Ashkin-Teller Hamiltonian
 */
class AshkinTellerThreeColor
{
	public:
		//  ----  Parameters  ----

		double J; ///< Ising interaction strength
		double K; ///< Four-spin interaction strength

		static constexpr int SymD = 3; ///< use SymD to encode the three colors of the model
		const std::string name;

		typedef std::array<int, SymD> StateVector;


		//  ---- Hamiltonian  -----

		/** Constructor
		 *  
		 * @param J the Ising interaction
		 * @param K the Four-spin interaction
		 */
		AshkinTellerThreeColor(double J, double K) : J(J), K(K), name("ThreeColorAshkinTeller"), observables(obs_m) {}
		

		// we need this only for the update algorithms to access ham.interactions[0]->J 
		std::array<Standard_Interaction<StateVector>*, 1> interactions = {new Standard_Interaction<StateVector>(J)};
		// can be avoided by providing an explicit specialization of the Wolff algorithm for this model


	
		//  ----  Observables ----

		AshkinTellerThreeColorMag       obs_m;
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

template <>
class Initializer<AshkinTellerThreeColor>
{
typedef typename AshkinTellerThreeColor::StateVector StateVector;
template <class RNGCache>
static StateVector newsv(const StateVector&, RNGCache&) {return StateVector();}
};



// ----------------------- UPDATE SPECIALIZATIONS ----------------------


namespace MARQOV 
{

	template <class Lattice>
	class Embedder<AshkinTellerThreeColor, Lattice>
	{
		typedef AshkinTellerThreeColor Hamiltonian;
		typedef typename Hamiltonian::StateVector StateVector;
		typedef Space<StateVector, Lattice> StateSpace;
		static constexpr int SymD = Hamiltonian::SymD;

		private:

			const Hamiltonian& ham;
			const Lattice& lat;
			const StateSpace& statespace;

			int rcolor;

		public:

			Embedder(const Hamiltonian& ham, const Lattice& lat, StateSpace& statespace) : ham(ham), lat(lat), statespace(statespace) {};


			template <class RNG>
			void draw(RNG& rng, StateVector& sv)	{ rcolor = rng.integer(3); }


			double coupling(int pos1, int pos2) const
			{
				double retval = 0.0;
				const auto sv1 = statespace[pos1];
				const auto sv2 = statespace[pos2];

				if (sv1[rcolor] != sv2[rcolor])
				{
					switch (rcolor)
					{
						case 0: retval = ham.J - ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]); break;
						case 1: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]); break;
						case 2: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]); break;
                        default: std::cout << "invalid color!" << rcolor << std::endl;
//						default: throw std::invalid_argument("invalid color!"); // catch me!
					}
				}
				return retval;
			}


			void flip(StateVector& sv) { sv[rcolor] *= -1; }
		};
				

	
	
	
	/** Specialized Metropolis algorithm for the Ashkin-Teller Hamiltonian
	 *
	 * @tparam Lattice type of the lattice
	 */
	template <class Lattice>
	struct Metropolis<AshkinTellerThreeColor, Lattice>
	{

		// some typedefs
		typedef typename AshkinTellerThreeColor::StateVector StateVector;
		typedef int ReducedStateVector;


		// coupling of the embedded Ising model
		static inline double metro_coupling(StateVector& sv1, StateVector& sv2, int color, const AshkinTellerThreeColor& ham)
		{
			double retval = 0.0;
			switch (color)
			{
				case 0: retval = ham.J - ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]); break;
				case 1: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]); break;
				case 2: retval = ham.J - ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]); break;
                default: std::cout << "invalid color!" << color << std::endl;
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
		template <class StateSpace,  class RNG>
		static inline int move(const AshkinTellerThreeColor& ham, 
						   const Lattice& grid, 
						   StateSpace& statespace, 
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
