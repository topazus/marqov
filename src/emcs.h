#ifndef EMCS_H
#define EMCS_H
#include "metropolis.h"
#include "wolff.h"
// Defines the Elementary Monte Carlo Step (EMCS)

namespace MARQOV
{
template <class Grid, class Hamiltonian, template<class> class RefType>
double Core<Grid, Hamiltonian, RefType>::elementaryMCstep()
{
	const int SymD = std::tuple_size<StateVector>::value;

	// cluster updates
	double avgclustersize = 0;
	for (int j=0; j < mcfg.ncluster; j++)
	{
		const int rsite = rngcache.integer(this->grid.size());

		// Heisenberg; random direction
		const auto rdir = rnddir<RNGCache<RNGType>, typename StateVector::value_type, SymD>(rngcache);

		avgclustersize += wolffstep(rsite, rdir);
	}

	
	int metrocounter = 0;


	// Metropolis sweeps
	for (int j=0; j<mcfg.nsweeps; j++)
	{
		// loop sites
		for(decltype(this->grid.size()) i = 0; i < this->grid.size(); ++i)
		{
			const int rsite = rngcache.integer(this->grid.size());
			metrocounter += metropolisstep(rsite);
		}
	}

	return metrocounter/double(mcfg.nsweeps)/this->grid.size();
}
};
#endif
