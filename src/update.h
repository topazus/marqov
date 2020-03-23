#ifndef UPDATE_H
#define UPDATE_H

template <class Grid, class Hamiltonian> 
inline int Marqov<Grid, Hamiltonian>::wolffstep(int rsite, StateVector rdir)
{
	std::vector<int> cstack(grid.size(), 0);

	int q = 0;

	reflect(statespace[rsite], rdir);

	cstack[q] = rsite;

	int clustersize = 1;
	int current = 0;

	while (q>=0)
	{
		current = cstack[q];
		q--;
		
		const double proj1 = dot(statespace[current], rdir);

		int a = 0; // to be replaced by loop over Nalpha
		const auto nbrs = grid.getnbrs(a, current);

		for (int i = 0; i < nbrs.size(); ++i)
		{
			const auto mynbr = nbrs[i];
			//auto myvec = ham.interactions[a]->operator()(statespace[mynbr]);
			StateVector& myvec = statespace[mynbr];

			const double proj2 = dot(myvec, rdir);

			if (proj1*proj2 < 0)
			{
				const double prob = 1.0 - exp(2.0*ham.beta*proj1*proj2);

				if (rng.d() < prob)
				{
					q++;

					cstack[q] = mynbr;
					clustersize++;

					reflect(myvec, rdir);

					normalize(myvec); // necessary?
				}
			}
		}
	}

	return clustersize;
}


// Single Metropolis update step statevectors on a lattice
// returns an integer which encodes whether the flip attempt was successful (1) or not (0)
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
