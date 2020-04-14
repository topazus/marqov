#ifndef EMCS_H
#define EMCS_H

// Defines the Elementary Monte Carlo Step (EMCS)


template <class Grid, class Hamiltonian>
double Marqov<Grid, Hamiltonian>::elementaryMCstep(const int ncluster, const int nsweeps)
{
	const int SymD = std::tuple_size<StateVector>::value;

	// cluster updates
	double avgclustersize = 0;
	for (int j=0; j<ncluster; j++)
	{
		const int rsite = rng.i();
		const auto rdir = rnddir<RND, typename StateVector::value_type, SymD>(rng);
		avgclustersize += wolffstep_general(rsite, rdir);
	}


	// Metropolis sweeps
	for (int j=0; j<nsweeps; j++)
	{
		for(int i = 0; i < grid.size(); ++i)
		{
			const int rsite = rng.i();
			metropolisstep(rsite);
		}
	}

	return avgclustersize/ncluster;
}
		

#endif
