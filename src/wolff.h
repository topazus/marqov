#ifndef WOLFF_H
#define WOLFF_H
#include <vector>
#include <type_traits>
#include <cmath>
// todo: what about the alpha-loop? currently alpha=0 hard-coded

template <class Hamiltonian, class Lattice>
struct Wolff
{
    template <class DirType, class RNG, class StateSpace>
    static inline int move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite, const DirType& rdir);
};

template <class Hamiltonian, class Lattice>
template <class DirType, class RNG, class StateSpace>
int Wolff<Hamiltonian, Lattice>::move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNG& rng, double beta, int rsite, const DirType& rdir)
{
    typedef typename Hamiltonian::StateVector StateVector;
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
		for (std::size_t i = 0; i < nbrs.size(); ++i)
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
				if (rng.real() < -std::expm1(-2.0*beta*coupling))
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

template <class Grid, class Hamiltonian, template<class> class RefType>
template <class DirType>
inline int Core<Grid, Hamiltonian, RefType>::wolffstep(int rsite, const DirType& rdir)
{
    return Wolff<Hamiltonian, Grid>::move(this->ham, this->grid, this->statespace, this->rngcache, this->beta, rsite, rdir);
}
#endif
