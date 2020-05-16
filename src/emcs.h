#ifndef EMCS_H
#define EMCS_H
#include "metropolis.h"
#include "wolff.h"
// Defines the Elementary Monte Carlo Step (EMCS)


template <class Grid, class Hamiltonian, template<class> class RefType>
double Marqov<Grid, Hamiltonian, RefType>::elementaryMCstep()
{
	const int SymD = std::tuple_size<StateVector>::value;
	const int ncolors = 1;

	// cluster updates
	double avgclustersize = 0;
	for (int j=0; j < mcfg.ncluster; j++)
	{
		const int rsite = rng.i();

		// random direction
		const auto rdir = rnddir<RND, typename StateVector::value_type, SymD>(rng);
		avgclustersize += wolffstep(rsite, rdir);
	}


	// Metropolis sweeps
	for (int j=0; j<mcfg.nsweeps; j++)
	{
		// loop colors
		for (int s = 0; s<ncolors; s++)
		{
			// loop sites
			for(int i = 0; i < this->grid.size(); ++i)
			{
				const int rsite = rng.i();
				metropolisstep(rsite);
			}
		}
	}

	return avgclustersize/mcfg.ncluster/ncolors;
}
#endif
