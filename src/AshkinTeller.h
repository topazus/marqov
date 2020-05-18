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
//				cout << statespace[i][0] << "\t" << statespace[i][1] << "\t" << statespace[i][2] << endl;
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



		// using the general Wolff cluster algorithm requires to implement
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class color = int>
		inline double wolff_coupling(StateVector& sv1, StateVector& sv2, const color a=0) const
		{
			if (sv1[a] == sv2[a]) return 0.0;
			else
			{
				switch (a)
				{
					case 0: return J + K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
					case 1: return J + K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
					case 2: return J + K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
				}
			}
		}

		template <class color = int>
		inline void wolff_flip(StateVector& sv, const color a=0) const
		{
			sv[a] *= -1;
		}

};


namespace MARQOV 
{

template <class Lattice>
struct Metropolis<AshkinTeller<int>, Lattice>
{
	typedef typename AshkinTeller<int>::StateVector StateVector;
	typedef int ReducedStateVector;

	static inline double metro_coupling(StateVector& sv1, StateVector& sv2, const int color, const AshkinTeller<int>& ham)
	{
		switch (color)
		{
			case 0: return ham.J + ham.K * (sv1[1]*sv2[1] + sv1[2]*sv2[2]);
			case 1: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[2]*sv2[2]);
			case 2: return ham.J + ham.K * (sv1[0]*sv2[0] + sv1[1]*sv2[1]);
		}
	}


	static inline ReducedStateVector metro_newconf(ReducedStateVector& rsv)
	{
		ReducedStateVector retval(rsv);
		retval = -rsv;
		return retval;
	}


	static inline void metro_flip(StateVector& sv, const int color)
	{
		sv[color] *= -1;
	}


	template <class StateSpace, class M, class RNG>
	static inline int move(	const AshkinTeller<int>& ham, 
						const Lattice& grid, 
						StateSpace& statespace, 
						M& metro,
						RNG& rng, 
						double beta, 
						int rsite)
	{
		int retval = 0;
		for (int color=0; color<3; color++)
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

			double averagevector = 0; 

			// sum over neighbours
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
