/* This file is part of MARQOV:
 * A modern framework for classical spin models on general topologies
 * Copyright (C) 2020-2021, The MARQOV Project
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef WOLFF_H
#define WOLFF_H
#include <vector>
#include <cmath>
// todo: what about the alpha-loop? currently alpha=0 hard-coded

namespace MARQOV
{
    /**
     * This class serves as an entry poi8nt for easily defining your own
     * specializations of the Wolff Algorithm of your Hamiltonians.
     * To that end it has the two prototypical template parameters:
     * @tparam Hamiltonian The Hamiltonian that the Wolff algo will use.
     * @tparam Lattice The Lattice, that the Wolff algo should use.
     */
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

        const double global_coupling = ham.interactions[a]->J;
        
        // add initial site and flip it
        int q = 0;
        cstack[q] = rsite;

		ham.wolff_embedd_init(); // the random direction/color whatever is stored in the Hamiltonian! this function resets it
        
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
            const auto nbrs = grid.nbrs(a, currentidx);
            const auto bnds = grid.bnds(a, currentidx); // important: specify what to do if function is not there

            const auto cpl1 = ham.wolff_coupling(currentsv, nbrs); 
			// rdir is implicitely stored in the Hamiltonian! (actually this function should be named "wolff_embedding")
            const auto cpl2 = ham.wolff_scalarize(grid.nbrs(a, currentidx)); 
			// specifiy what to do with the bond strength; 
			// in practice has to communicate with wolff_coupling via rdir which is stored in the lattice
			// todo: if function is not there, return array of 1's

            // loop over neighbours
            for (std::size_t i = 0; i < nbrs.size(); ++i)
            {
                // extract corresponding sv
                const auto currentnbr = nbrs[i];
                StateVector& candidate = statespace[currentnbr];
                
                const double coupling = global_coupling * cpl1[i] * cpl2[i];
                
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
};


template <class Hamiltonian, class = void>
struct has_wolff_embedding : std::false_type{};

template <class Hamiltonian, class = void>
struct has_wolff_embedding<Hamiltonian, MARQOV::detail::type_sink_<
	decltype( ...
	>> : std::true_type{};


template <class Hamilonian, class BondType>
auto wolff_embedding_helper(Hamiltonian& ham, BondType bond, std::true_type)
{
	return ham.wolff_embedding(bond);
}


template <class Hamilonian, class BondType>
auto wolff_embedding_helper(Hamiltonian& ham, BondType bond, std::false)
{
	return bond[0];
}


template <class Hamilonian, class BondType>
auto wolff_embedding(Hamiltonian& ham, BondType bond)
{
	return wolff_embedding_helper(ham, bond, has_wolff_embedding<Hamiltonian>{});
}


#endif
