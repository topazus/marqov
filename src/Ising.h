#ifndef ISING_H
#define ISING_H
#include <array>
#include <tuple>
#include <string>
#include <functional>
#include "hamiltonianparts.h"
#include "metropolis.h"

// ------------------------------ OBSERVABLES ---------------------------

// Magnetization
class IsingMag
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
			}

			return std::abs(mag)/double(N);
		}
		IsingMag() : name("m") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class Ising_interaction : public Interaction<StateVector> 
{
public:
	Ising_interaction(double J)
	{
		this->J = J;
	}
	StateVector operator() (const StateVector& phi) {return phi;};
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
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = Ising_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		Ising(double J) : J(J) {	interactions[0] = new Ising_interaction<StateVector>(J); }
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		IsingMag       obs_m;
		auto getobs()
		{
			return std::make_tuple(obs_m);
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


		// think about this ... (see update.h)
		template <typename bond_type>
		inline double wolff_scalarize(const std::vector<bond_type>& bond)
		{
			return bond[0];
		}

	
};

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
	for(int a = 0; a < ham.Nalpha; ++a)
	{
		auto nbrs = grid.getnbrs(a, rsite);
		typedef decltype(ham.interactions[a]->operator()(statespace[0])) InteractionType;
		typedef decltype(MARQOV::callbonds<Lattice>(grid, a, rsite, 0, ham.interactions[a]->operator()(statespace[0]))) BondType;
        
		typename MARQOV::Promote_Array<InteractionType, BondType>::CommonArray averagevector = {0};
		// sum over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			auto idx = nbrs[i];
			auto nbr = ham.interactions[a]->operator()(statespace[idx]);
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
    else if (rng.d() < exp(-beta*dE))
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
		StateVector& currentsv = statespace[currentidx];
		q--;
	
		// get its neighbours
		const auto nbrs = grid.getnbrs(a, currentidx);

		// loop over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			// extract corresponding sv
			const auto currentnbr = nbrs[i];
			StateVector& candidate = statespace[currentnbr];

			// test whether site is added to the cluster
			if (candidate[0] == val)
			{
				if (rng.d() < prob)
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
#endif
