#ifndef HEISENBERG_H
#define HEISENBERG_H
#include <array>
#include <cmath>
#include <string>
#include <functional>
#include "vectorhelpers.h"
#include "hamiltonianparts.h"


// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class HeisenbergMag
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			constexpr static int SymD = 3;	// improve me
			const     int N = grid.size();

			std::vector<double> mag(SymD,0) ;

			for (int i=0; i<N; i++)
			{
				for (int j=0; j<SymD; j++)
				{
					mag[j] += statespace[i][j];
				}
			}
			
			double retval = 0;
			for (int j=0; j<SymD; j++)
			{
				retval += mag[j]*mag[j];
			}
			return sqrt(retval)/double(N);
		}
		HeisenbergMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector, class RNG>
class Heisenberg_Initializer
{
	public:
		// provide the spin dimension as a compile-time constant expression
		constexpr static int SymD = std::tuple_size<StateVector>::value;

		// constructors
		Heisenberg_Initializer()   {}
		Heisenberg_Initializer(RNG& rn) : rng(rn) {}

		// generate new statevector
		StateVector newsv(const StateVector&) 
		{
			return rnddir<RNG, double, SymD>(rng);
		};

	private:
		RNG& rng;
};



template <class StateVector>
class Heisenberg_interaction : public Interaction<StateVector> 
{
	public:
		Heisenberg_interaction(double J)
		{
	 		this->J = -J;
		}
		StateVector operator() (const StateVector& phi) {return phi;};
};


// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType, typename MyFPType>
class Heisenberg
{
	public:

		double J;
		constexpr static int SymD = 3;
		typedef MyFPType FPType;
		typedef std::array<SpinType, SymD> StateVector;
		
		template <typename RNG>
		using MetroInitializer =  Heisenberg_Initializer<StateVector, RNG>; 
		// this construction allows to specify a number of template arguments
		// while leaving others open (C++11 feature)

		
		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;

		// requires pointers
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, FPType>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];

		Heisenberg(double J) : J(J) {interactions[0] = new Heisenberg_interaction<StateVector>(J);}
		

		// instantiate and choose observables
		HeisenbergMag obs_m;
		auto getobs() { return std::make_tuple(obs_m); }
		

		// state space initializer
		template <class StateSpace, class Lattice, class RNG>
		void initstatespace(StateSpace& statespace, Lattice& grid, RNG& rng) const
		{
			for(int i=0; i<grid.size(); ++i)
			{
				statespace[i] = rnddir<RNG, typename StateVector::value_type, SymD>(rng);
			}
		}
		

		// using the Wolff cluster algorithm requires to implement 
		// the functions 'wolff_coupling' and 'wolff_flip'

		template <class A> 
		inline auto wolff_coupling(StateVector& sv1, StateVector& sv2, const A a) const
		{
			return dot(sv1, a) * dot(sv2, a);
		}

		template <class A>
		inline void wolff_flip(StateVector& sv, const A a) const
		{
			const double dotp = dot(sv, a);
			for (int i=0; i<SymD; i++) sv[i] -= 2*dotp*a[i];

			normalize(sv);  // necessary?
		}

};

namespace MARQOV
{
    template <class Lattice, class SpinType, class FPType>
    struct Wolff<Heisenberg<SpinType, FPType>, Lattice>
    {
        template <class DirType, class RNG, class StateSpace>
        static inline int move(const Heisenberg<SpinType, FPType>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite, const DirType& rdir)
        {
            typedef typename Heisenberg<SpinType, FPType>::StateVector StateVector;
            std::vector<int> cstack(grid.size(), 0);
            int q = 0;
            ham.wolff_flip(statespace[rsite], rdir);
            cstack[q] = rsite;
            
            int clustersize = 1;
            int current = 0;
            
            while (q>=0)
            {
                current = cstack[q];
                q--;

		// plain Heisenberg model has only one interaction term
		auto coupling = ham.interactions[0]->J; 
		const auto proj1 = coupling*dot(statespace[current], rdir);

		const auto nbrs = grid.getnbrs(0, current);
		for (std::size_t i = 0; i < nbrs.size(); ++i)
		{
			const auto currentidx = nbrs[i];
			StateVector& candidate = statespace[currentidx];

			const auto proj2 = dot(candidate, rdir);

			if (proj1*proj2 > 0)
			{
				const auto prob = -std::expm1(-2.0*beta*proj1*proj2);

				if (rng.real() < prob)
				{
					q++;
					cstack[q] = currentidx;
					clustersize++;
					ham.wolff_flip(candidate, rdir);
				}
			}
		}
	}
	return clustersize;
        }
    };
}
#endif
