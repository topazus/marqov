#ifndef EDWARDSANDERSONISING_H
#define EDWARDSANDERSONISING_H
#include <array>
#include <tuple>
#include <string>
#include <complex>
#include <functional>
#include "hamiltonianparts.h"
#include "metropolis.h"


// ------------------------------ OBSERVABLES ---------------------------

// the traditional EA order parameter
class EdwardsAndersonOrderParameter
{
	public:
		int counter = 0;
		std::string name;
		std::vector<int> local_sum;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();

			if (local_sum.size() == 0) 
			{
				local_sum.resize(size);
				for (int i=0; i<size; i++) { local_sum[i] = 0; }
			}

			double retval = 0;

			counter++;
			for (int i=0; i<size; i++)
			{
				local_sum[i] += statespace[i][0];
				retval += pow(local_sum[i],2);
			}

			return retval / double(size) / double(counter) / double(counter);
		}

		EdwardsAndersonOrderParameter() : name("qEA") {}
};


class Susceptibility
{
	private:
		double kx;
	public:
		int counter = 0;
		std::string name;
		std::vector<int> sum_i;
		std::vector<std::vector<int>> sum_ij;

		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			const int size = grid.size();
			std::complex<double> jj(0,1); 

			const double norml = 1. /  double(size) / double(size) / double(counter) / double(counter);

			if (sum_ij.size() == 0) 
			{
				sum_ij.resize(size);
				for (int i=0; i<size; i++)
				{
					sum_ij[i].resize(size);
					for (int j=0; j<size; j++)
					{
						sum_ij[i][j] = 0;
					}
				}
			}

			if (sum_i.size() == 0) 
			{
				sum_i.resize(size);
				for (int i=0; i<size; i++) sum_i[i] = 0;
			}

			std::complex<double> retval = 0;

			counter++;

			for (int i=0; i<size; i++)
			{
				sum_i[i] += statespace[i][0];

				for (int j=0; j<size; j++)
				{
					sum_ij[i][j] += statespace[i][0]*statespace[j][0];
				}
			}


			for (int i=0; i<size; i++)
			{
				for (int j=0; j<size; j++)
				{
					const int dir = 0;

					const std::vector<double> xi = {grid.getcrds(i)[dir]};
					const std::vector<double> xj = {grid.getcrds(j)[dir]};

					auto diff = xi[0] - xj[0];

					if (fabs(diff)>0.5) diff = 1.0 - fabs(diff); // account for PBC

					std::complex<double> phase = std::exp(kx*diff*jj);
//					std::complex<double> phase = std::exp(2.0*M_PI*fabs(diff)*jj);

					/*
					auto indi = IndexOf(i, grid.dim, grid.len);
					auto indj = IndexOf(j, grid.dim, grid.len);

					cout << i;
					for (int k=0; k<indi.size(); k++) cout << indi[k] << "  ";
					cout << endl;
					cout << j;
					for (int k=0; k<indj.size(); k++) cout << indj[k] << "  ";
					cout << endl;

					cout << std::fixed << std::setprecision(0) <<  phase << endl << endl;
					*/

					retval += pow(sum_ij[i][j],2) * phase;
				}
			}

			return norml * std::abs(retval); // or real part? or what ...?

			// open questions:
			// - should the distance vector account for PBC?
			// - correct normalization
			// - order of averages correct?
			// - what to return? absolute value, real part, ...?
		}

		Susceptibility(double kx, std::string name) : kx(kx), name(name) {}
};

// Scalar overlap
class ScalarOverlap
{
	public:
		std::string name;
		template <class StateSpace, class Grid>
		double measure(const StateSpace& statespace, const Grid& grid)
		{
			//
			//
			//
			//
			return 0;
		}
		ScalarOverlap() : name("q") {}
};


// ----------------------------------------------------------------------

template <class StateVector>
class EdwardsAndersonIsing_interaction : public Interaction<StateVector> 
{
public:
	EdwardsAndersonIsing_interaction(double J)
	{
		this->J = J;
	}
	StateVector get (const StateVector& phi) {return phi;};
};


template <class StateVector, class RNG>
class EdwardsAndersonIsing_Initializer
{
	public:
		EdwardsAndersonIsing_Initializer()   {}
		EdwardsAndersonIsing_Initializer(RNG&) {}

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
class EdwardsAndersonIsing
{
	public:
		double J;
		constexpr static int SymD = 1;
		typedef std::array<SpinType, SymD> StateVector;
		template <typename RNG>
		using MetroInitializer = EdwardsAndersonIsing_Initializer<StateVector, RNG>;

		static constexpr uint Nalpha = 1;
		static constexpr uint Nbeta = 0;
		static constexpr uint Ngamma = 0;
		
		EdwardsAndersonIsing(double J) : J(J), obs_chi(0, "chi") , obs_chiKmin(2.0*M_PI, "chiKmin")
		{
			interactions[0] = new EdwardsAndersonIsing_interaction<StateVector>(J); 
		}
		
		// instantiate interaction terms (requires pointers)
		Interaction<StateVector>* interactions[Nalpha];
		OnSite<StateVector, int>* onsite[Nbeta];
		MultiSite<StateVector*,  StateVector>* multisite[Ngamma];
	
		// instantiate and choose observables
		EdwardsAndersonOrderParameter      obs_qEA;
		ScalarOverlap					obs_q;
		Susceptibility					obs_chi;
		Susceptibility					obs_chiKmin;
//		auto getobs()	{return std::make_tuple(obs_qEA, obs_q, obs_chi);}
		auto getobs()	{return std::make_tuple(obs_qEA, obs_chi, obs_chiKmin);}


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






namespace MARQOV {

    //Work around GCC Bug for specializations of things in namespaces:
    //https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56480
    // still occurs on 6.4.0

	template <class Lattice>
	struct Metropolis<EdwardsAndersonIsing<int>, Lattice>
	{
	    template <class StateSpace, class M, class RNG>
	    static int move(const EdwardsAndersonIsing<int>& ham, const Lattice& grid, StateSpace& statespace, M& metro, RNG& rng, double beta, int rsite)
	    {
			typedef typename EdwardsAndersonIsing<int>::StateVector StateVector;
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

//					auto temp = MARQOV::callbonds<Lattice>(grid, a, rsite, i, nbr);
//					cout << temp[0] << endl; 

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

}
#endif
