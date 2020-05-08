#ifndef UPDATE_H
#define UPDATE_H
#include <vector>
#include <type_traits>
#include <cmath>
#include "metropolis.h"
// todo: what about the alpha-loop? currently alpha=0 hard-coded
// implement me: does not support locally fluctating (e.g. random) interaction strengths yet


template <class Grid, class Hamiltonian, template<class> class RefType>
template <class DirType>
inline int Marqov<Grid, Hamiltonian, RefType>::wolffstep_general(int rsite, const DirType& rdir)
{
	// prepare stack
	std::vector<int> cstack(this->grid.size(), 0);

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
		const auto nbrs = this->grid.getnbrs(a, currentidx);

		// loop over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			// extract corresponding sv
			const auto currentnbr = nbrs[i];
			StateVector& candidate = statespace[currentnbr];

			// compute 'Wolff coupling'
			const double global_coupling = ham.interactions[a]->J;
			const auto   local_coupling  = 1;

			/* under construction
			const auto   local_coupling  = grid.getbnds(a, currentnbr)[0];

			// !!!  Wolff and Swendsen-Wang cluster algorithms are only valid if 
			// !!!  all interactions are ferromagnetic, i.e., if all J_ij > 0 (Zhu et. al 2015)


			// even more general would be somthing like that:
			// const auto   local_coupling  = ham.wolff_scalarize(grid.getbnds(a, currentnbr));
			// overkill, or even necessary?
			*/

			const double wolff_coupling  = ham.wolff_coupling(currentsv, candidate, rdir);
			const double coupling = global_coupling * local_coupling * wolff_coupling;

			// test whether site is added to the cluster
			if (coupling > 0)
			{
				if (rng.d() < -std::expm1(-2.0*beta*coupling))
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




// in the plain Ising model, the Wolff coupling is a constant, which can be
// exploited for optimization

template <class Grid, class Hamiltonian, template<class> class RefType>
inline int Marqov<Grid, Hamiltonian, RefType>::wolffstep_Ising(int rsite)
{
	// prepare stack
	std::vector<int> cstack(this->grid.size(), 0);

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
		const auto nbrs = this->grid.getnbrs(a, currentidx);

		// loop over neighbours
		for (int i = 0; i < this->nbrs.size(); ++i)
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








template <class Grid, class Hamiltonian, template<class> class RefType>
inline int Marqov<Grid, Hamiltonian, RefType>::wolffstep_Heisenberg(int rsite, const StateVector& rdir)
{
	std::vector<int> cstack(this->grid.size(), 0);
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

		const auto nbrs = this->grid.getnbrs(a, current);
		for (int i = 0; i < this->nbrs.size(); ++i)
		{
			const auto currentidx = nbrs[i];
			StateVector& candidate = statespace[currentidx];

			const auto proj2 = dot(candidate, rdir);

			if (proj1*proj2 > 0)
			{
				const auto prob = -std::expm1(-2.0*beta*proj1*proj2);

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
/*
// filtered Metropolis prototype ....

// takes a function which takes a StateVector and returns a reduced StateVector

template <class Grid, class Hamiltonian>
template <typename callable1, typename callable2>
inline int Marqov<Grid, Hamiltonian>::metropolisstep(int rsite, callable1 filter_ref, callable2 filter_cpy, int comp)
{
	// old state vector at rsite
	StateVector&     svold = statespace[rsite];
	redStateVector& rsvold = filter_ref(svold, comp);

	// propose new configuration
	redStateVector  rsvnew = metro.newsv(rsvold);
	
	// interaction part
	double interactionenergydiff = 0;
	for(int a = 0; a < ham.Nalpha; ++a)
	{
		// extract neighbours
		auto nbrs = grid.getnbrs(a, rsite);
		double averagevector = {0};
		
		// sum over neighbours
		for (int i = 0; i < nbrs.size(); ++i)
		{
			// neighbour index
			auto idx = nbrs[i];
			// full neighbour
			auto nbr  = ham.interactions[a]->operator()(statespace[idx]);
			// reduced neighbour
			auto rnbr = filter_cpy(nbr, comp);
			
			// coupling in the embedded model
			double cpl = ham.metro_coupling(svold, nbr, comp);

			// compute weighted sum of neighbours
			averagevector = averagevector + mult(cpl,rnbr);
		}

		interactionenergydiff += ham.interactions[a]->J * (dot(rsvnew-rsvold, averagevector));
	}

    // (...)

    double dE 	= interactionenergydiff; // + ... + ...

    int retval = 0;
    if ( dE <= 0 )
    {
        rsvold = rsvnew;
        retval = 1;
    }
    else if (rng.d() < exp(-beta*dE))
    {
        rsvold = rsvnew;
        retval = 1;
    }

    return retval;
}

*/

#endif
