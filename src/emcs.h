#ifndef EMCS_H
#define EMCS_H

// Defines the Elementary Monte Carlo Step (EMCS)


template <class Grid, class Hamiltonian>
double Marqov<Grid, Hamiltonian>::elementaryMCstep(const int ncluster, const int nsweeps)
{
	const int SymD = std::tuple_size<StateVector>::value;
	const int ncolors = 3;

	// cluster updates
	double avgclustersize = 0;
	for (int j=0; j<ncluster; j++)
	{
		const int rsite = rng.i();

		// random direction
//		const auto rdir = rnddir<RND, typename StateVector::value_type, SymD>(rng);
//		avgclustersize += wolffstep_general(rsite, rdir);

		// loop colors
		for (int i=0; i<ncolors; i++) avgclustersize += wolffstep_general(rsite, i);
	}



	// Metropolis sweeps
	for (int j=0; j<nsweeps; j++)
	{
		// loop colors
		for (int s = 0; s<ncolors; s++)
		{
			// loop sites
			for(int i = 0; i < grid.size(); ++i)
			{
				const int rsite = rng.i();
//				metropolisstep(rsite);
				metropolisstep(rsite, reduce_ref, reduce_cpy, s);
			}
		}
	}

	return avgclustersize/ncluster/ncolors;
}
		

#endif
