#ifndef MARQOV_H
#define MARQOV_H

#include <array>
#include <vector>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::flush;


template <class Grid, class Hamiltonian>
class Marqov 
{
	public:
		typedef typename Hamiltonian::StateVector StateVector;
		typedef StateVector* StateSpace;

		// Constructor
		Marqov(Grid& lattice) : ham(),  grid(lattice), rng(0, 1), metro(rng) 
		{
		  	rng.seed(42);
		  	rng.seed(time(NULL));
		  	rng.set_integer_range(lattice.size());
		  	statespace = new typename Hamiltonian::StateVector[lattice.size()];
		}

		// Definition of an EMCS
	    void elementaryMCstep()
	    {
	        for(int i = 0; i < grid.size(); ++i)
	        {
                const int rsite = rng.i(); // choose random site -> Moved one level above
	            metropolisstep(rsite);
	        }
	    }
	    
	    void gameloop()
	    {
	        for (int i = 0; i < nstep; ++i)
	            elementaryMCstep();
	        
	        for(int j = 0; j < nobs; ++j)
	            obs[j]->measure(statespace);
					//improve me: consider that there might be reuse across observables!
	    }
	
	
		void visualize_state_2d(int dim=0, double threshold=0.4)
		{
			cout << "_";
			for(int i = 0; i < grid.length; ++i) cout << " _";
			cout << endl;
			for(int i = 0; i < grid.length; ++i)
			{
				cout << "|";
				for(int j = 0; j < grid.length; ++j)
				{
					double current = statespace[grid.length*i+j][dim];

					if (current > threshold) cout << "O ";
					else if (current < -threshold) cout << "  ";
					else if (current > 0) cout << "o ";
					else if (current < 0) cout << ". ";
				}
				cout << "|" << endl;
			}
			cout << "‾";
			for(int i = 0; i < grid.length; ++i) cout << " ‾";
			cout << endl << endl;
		}


	
		 void init_cold()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i][0] = 1;
			}
		 }

		 void init_hot()
		 {
			for(int i = 0; i < grid.size(); ++i)
			{
				statespace[i] = rnddir<RND, double, 3>(rng);
			}
		 }
	
	private:


	// Single Metropolis update step statevectors on a lattice
	// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
	inline int metropolisstep(int rsite)
	{
	    StateVector& svold = statespace[rsite];
	    StateVector svnew = metro.newsv(svold);

	    // interaction part
	    double interactionenergydiff = 0;
	    for(int a = 0; a < ham.Nalpha; ++a)
	    {
	        auto nbrs = grid.getnbrs(a, rsite);
	        StateVector averagevector = {0};
	
	        for (int i = 0; i < nbrs.size(); ++i)
	        {
	            auto mynbr = nbrs[i];
	            auto myvec = ham.interactions[a]->operator()(statespace[mynbr]);
	            averagevector = averagevector + myvec;
	        }
	        interactionenergydiff += ham.interactions[a]->J * (dot(svnew - svold, averagevector));
	    }
	
	    // onsite energy part
	    double onsiteenergydiff = 0;
	    for (int b = 0; b < ham.Nbeta; ++b)
	    {
            // compute the difference
            auto diff = ham.onsite[b]->operator()(svnew) - ham.onsite[b]->operator()(svold);
            // multiply the constant
            onsiteenergydiff += dot(ham.onsite[b]->h, diff);
	    }
        
	    
	    // multi-site energy
	    double multisiteenergyold = 0;
	    double multisiteenergynew = 0;
	    for (int g = 0; g < ham.Ngamma; ++g)
	    {
	        multisiteenergynew += ham.multisite[g]->operator()(svnew, rsite, statespace);//FIXME: think about this...
	        multisiteenergyold += ham.multisite[g]->operator()(svold, rsite, statespace);//FIXME: think about this...
	        //forgot k_gamma
	    }
	    
	    double dE 	= interactionenergydiff + onsiteenergydiff + (multisiteenergynew - multisiteenergyold);
	
	
	    int retval = 0;
	    if ( dE <= 0 )
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    else if (rng.d() < exp(-ham.beta*dE))
	    {
	        svold = svnew;
	        retval = 1;
	    }
	    
	    return retval;
	}


	void wolff()
	{
	}


	StateSpace statespace;
	Hamiltonian ham;
	Grid& grid;
	RND rng;

	//Get the MetroInitializer from the user, It's required to have one template argument left, the RNG.
	typename Hamiltonian::template MetroInitializer<RND> metro;//C++11

	// number of observables
	static constexpr uint nobs = 0;
	Observable<StateSpace>* obs[5];

	// number of EMCS
	static constexpr int nstep = 2500;
};

#endif