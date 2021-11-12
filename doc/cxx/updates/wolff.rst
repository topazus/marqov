.. Copyright (c) 2021, Manuel Schrauth, Florian Goth

Wolff Cluster Algorithm
============================
In contrast to the Metropolis algorithm, this 
is our implemementation of a as general as possible
Cluster update according to the Wollf Cluster algorithm.
In case your needs are not covered for your particular model,
it is possible to provide a partial class specialization.

Writing your own specialization for your newly written Hamiltonian "MYHamiltonian" that still works on all Lattices, is easy, you can start with just copying 
the default Wolff template from wolff.h and empty it from everything you don't need
and you end up with something like this:

.. code-block:: cpp

    namespace MARQOV
    {
        template <Lattice>
        struct Wolff<class MYHamiltonian, class Lattice>
        {
            template <class RNGType, class StateSpace>
            inline int move(const Hamiltonian& ham, const Lattice& grid, StateSpace& statespace, RNGCache<RNGType>& rng, double beta, int rsite);
            {
                int clustersize = 0;
            
                /*INSERT YOUR HIGHLY ADAPTED WOLFF UPDATE HERE!*/
                
                return clustersize;
            }
            std::vector<int> cstack = std::vector<int>(4096/sizeof(int), 0);///< the size of the stack is meant to be preserved across different cluster processes.
        };
    }

The Wolff update uses a full class which enables you to preserve state across the runs of a Wolff update.
MARQOV by default uses this to save reallocations of the cstack array as you see in the previous example.
Another use-case for specializing the Wolff update is if your modell allows heavy optimizations that can be exploited, as e.g. in the Ising Modell.

.. doxygenstruct:: MARQOV::Wolff
   :project: Marqov
   :members:
