#ifndef UPDATE_H
#define UPDATE_H

#include <iomanip>
#include <vector>

// todo: what about the alpha-loop? currently alpha=0 hard-coded
// implement me: does not support locally fluctating (e.g. random) interaction strengths yet


template <class Grid, class Hamiltonian> 
template <class DirType>
inline int Marqov<Grid, Hamiltonian>::wolffstep_general(int rsite, const DirType& rdir)
{
	// prepare stack
	std::vector<int> cstack(grid.size(), 0);

	// add initial site and flip it
	int q = 0;
	cstack[q] = rsite;

	ham.wolff_flip(statespace[rsite], rdir);
	int clustersize = 1;
	
	// loop over stack as long as non-empty
	while (q>=0)
	{
		// extract last sv in stack
		const int currentidx = cstack[q];
		StateVector& currentsv = statespace[currentidx];
		q--;
	
		// get its neighbours
		int a = 0; // to be replaced by loop over Nalpha
		const auto nbrs = grid.getnbrs(a, currentidx);

		// loop over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			// extract corresponding sv
			const auto currentnbr = nbrs[i];
			StateVector& candidate = statespace[currentnbr];

			// compute 'Wolff coupling'
			const double global_coupling = ham.interactions[a]->J;
			const double local_coupling = 1.0; // not yet implemented: grid.getcoupling(...)
			const double coupling = global_coupling
								* local_coupling
								* ham.wolff_coupling(currentsv, candidate, rdir);

			// test whether site is added to the cluster
			if (coupling > 0)
			{
				if (rng.d() < 1.0-std::exp(-2.0*ham.beta*coupling))
				{
					q++;
					cstack[q] = currentnbr;
					clustersize++;
					ham.wolff_flip(candidate, rdir);
				}
			}
		}
	}
	return clustersize;
}




// in the plain Ising model, the Wolff couling is a constant, which can be
// exploited for optimization

template <class Grid, class Hamiltonian> 
inline int Marqov<Grid, Hamiltonian>::wolffstep_Ising(int rsite)
{
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
	const double prob = 1.0-std::exp(+2.0*ham.beta*coupling);
	
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








template <class Grid, class Hamiltonian> 
inline int Marqov<Grid, Hamiltonian>::wolffstep_Heisenberg(int rsite, const StateVector& rdir)
{
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
		
		int a = 0; // plain Heisenberg model has only one interaction term
		const double coupling = ham.interactions[a]->J; 
		const auto proj1 = coupling*dot(statespace[current], rdir);

		const auto nbrs = grid.getnbrs(a, current);
		for (int i = 0; i < nbrs.size(); ++i)
		{
			const auto currentidx = nbrs[i];
			StateVector& candidate = statespace[currentidx];

			const auto proj2 = dot(candidate, rdir);

			if (proj1*proj2 > 0)
			{
				const auto prob = 1.0 - std::exp(-2.0*ham.beta*proj1*proj2);

				if (rng.d() < prob)
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









// Single Metropolis update step statevectors on a lattice
// returns an integer which encodes whether the flip attempt was successful (1) or not (0)

// todo: does not support locally fluctating (e.g. random) interaction strengths

template <class Grid, class Hamiltonian> 
inline int Marqov<Grid, Hamiltonian>::metropolisstep(int rsite)
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

#endif
