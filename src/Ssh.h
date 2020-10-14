#ifndef SSH_H
#define SSH_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include "hamiltonianparts.h"
#include "metropolis.h"


// ------------------------------ OBSERVABLES ---------------------------


// ----------------------------------------------------------------------

class SSHMag
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid> 
		double measure(const StateSpace& statespace, const Grid& grid) 
		{
			const int N = grid.size();
			double mag = 0.0;

			for (int i=0; i<N; i++) 
			{
				 mag += statespace[i][0];
				}
			return mag/double(N);
//			return std::abs(mag)/double(N);
			}
			SSHMag() : name("m"), desc("The Magnetization of the SSH Modell") {}
};


class SSHMagSq
{
	public:
		std::string name, desc;
		template <class StateSpace, class Grid> 
		double measure(const StateSpace& statespace, const Grid& grid) 
		{
			const int N = grid.size();
			double mag = 0.0;

			for (int i=0; i<N; i++) 
			{
				 mag += pow(statespace[i][0],2);
				}
			return mag/double(N);
//			return std::abs(mag)/double(N);
			}
			SSHMagSq() : name("msq"), desc("...") {}
};

template <class StateVector>
class SSH_interaction : public Interaction<StateVector> 
{
public:
	SSH_interaction(double m, double dtau)
	{
//		this->J = m*dtau/4;
		this->J = -m/dtau;
	}
	StateVector get (const StateVector& phi) {return phi;};
};




template <class StateVector>
class SSH_onsite : public OnSite<StateVector, double> 
{
	public:
		SSH_onsite(double m, double k, double dtau)
		{
			this->h = m/dtau + k/2;
//			this->h = m*dtau/4 + k/2;
		}
		double get (const StateVector& phi) {return dot(phi,phi);}; 
};




template <class StateVector, class RNG>
class SSH_Initializer
{
	private:
		RNG& rng;
	public:
		SSH_Initializer()   {}
		SSH_Initializer(RNG& rn) : rng(rn) {}

		// specifies how a random new state vector is generated
		// in this case a simple spin flip
		StateVector newsv(const StateVector& svold) 
		{
			StateVector retval(svold); 
			double amp = 0.2;
			double diff = rng.real(-amp, amp);
			retval[0] += diff;
			return retval;
		};
};



// ------------------------------ HAMILTONIAN ---------------------------

template <typename SpinType = double>
class SSH
{
	public:
		double m, k, dtau;
		constexpr static int SymD = 1;
		const std::string name;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = SSH_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 1;
		static constexpr uint Ngamma = 0;
		
		SSH(double m, double k, double dtau) : m(m), k(k), dtau(dtau), name("SSH")
		{
			interactions[0] = new SSH_interaction<StateVector>(m, dtau); 
			onsite[0] = new SSH_onsite<StateVector>(m, k, dtau);
		}
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, double>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		SSHMag       obs_m;
		SSHMagSq       obs_msq;
		auto getobs()	{return std::make_tuple(obs_m, obs_msq);}


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




/*

namespace MARQOV {

    //Work around GCC Bug for specializations of things in namespaces:
    //https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56480
    // still occurs on 6.4.0

	template <class Lattice>
	struct Metropolis<SSH<double>, Lattice>
	{
	    template <class StateSpace, class M, class RNG>
	    static int move(const SSH<double>& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta, int rsite)
	    {
			typedef typename SSH<double>::StateVector StateVector;
		    	// old state vector at rsite
			StateVector& svold = statespace[rsite];
			// propose new configuration
			StateVector svnew = metro.newsv(svold);
			
			// interaction part
			double interactionenergydiff = 0;
			for(typename std::remove_cv<decltype(ham.Nalpha)> ::type a = 0; a < ham.Nalpha; ++a)
			{
				auto nbrs = grid.getnbrs(a, rsite);
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

			cout << "test" << endl;
			 // onsite energy part
			 double onsiteenergydiff = 0;
			 for (typename std::remove_cv<decltype(ham.Nbeta)>::type b=0; b<ham.Nbeta; ++b)
			 {
			 	// compute the difference
				auto diff = ham.onsite[b]->get(svnew) - ham.onsite[b]->get(svold);
				// multiply the constant
				onsiteenergydiff += dot(ham.onsite[b]->h, diff);
			}





	    		// sum up energy differences
	    		double dE 	= interactionenergydiff + onsiteenergydiff;
	
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
	struct Wolff<SSH<int>, Lattice>
	{
		template <class DirType, class RNG, class StateSpace>
		static inline int move(const SSH<double>& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite, const DirType&)
		{
			typedef typename SSH<int>::StateVector StateVector;
			// prepare stack
			std::vector<int> cstack(grid.size(), 0);
		
			// add initial site and flip it
			int q = 0;
			cstack[q] = rsite;
			const int val = statespace[rsite][0];
			ham.wolff_flip(statespace[rsite]);
			int clustersize = 1;
		
			// compute 'Wolff probability' 
			const int a = 0; // plain SSH model has only one interaction term
			const double coupling = ham.interactions[a]->J;
			const double prob = -std::expm1(+2.0*beta*coupling);
			
			// loop over stack as long as non-empty
			while (q>=0)
			{
				// extract last sv in stack
				const int currentidx = cstack[q];
				q--;
			
				// get its neighbours
				const auto nbrs = grid.getnbrs(a, currentidx);
	
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

		return 0;
	    }
	};

}

*/
#endif
