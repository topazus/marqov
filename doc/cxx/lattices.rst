.. Copyright (c) 2021, Manuel Schrauth, Florian Goth
Lattices
============
Lattices in MARQOV are user defined classes that expose 
a certain interface such that we can derive a Monte Carlo procedure from them.
As all of MARQOV's components, lattice have required parts and optional parts.
there are two required functions of a lattice, those being

.. code-block:: cpp

   std::vector<int> Lattice::nbrs(const int alpha, const int i) const

which enables us to infer neighbour relations in your graph, since the returned vector should contain the indices
of the neighbours of the site i.
The second function is

.. code-block:: cpp

   std::size_t size() const

which gives us the information how many sites there are in your lattice.

Then there are a few optional functions, that, if present, enable additional 
functionality of MARQOV

.. code-block:: cpp

   std::vector<int> Lattice::crds(const int i) const
   std::vector<int> Lattice::flexnbrs(const int alpha, const int i) const
   
   void writelat(H5::Group& h5loc, const Lattice& l)

.. toctree::
   :caption: Lattices
   :maxdepth: 2
   
   lattices/regular_hypercubic
   lattices/constant_coordination
   lattices/erdos_renyj
   lattices/graph_from_csv
   lattices/regular_random_bond
   lattices/simple_bipartite
   lattices/super_chaos
