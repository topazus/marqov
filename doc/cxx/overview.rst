.. Copyright (c) 2022, Manuel Schrauth, Florian Goth

Architectural Overview
========================
MARQOV consists of quite some moving parts and involves a lot of user definable concepts
that are defined via their interfaces:

.. uml:: marqov.puml

The parts closest to the user are :doc:`MARQOV::Core <core>` part which binds everything together
and :doc:`MARQOV::Scheduler <scheduler>` which enables parallelism and, if multiple simulations are scheduled,
parallel tempering between these simulations.
:doc:`MARQOV::Core <core>` can be adapted in variuos ways. A simple one is the change of the Random Number Generator(RNG) and hence
the RNGConcept can be fulfilled by various implementations. We expect them to follow the interface of the C++11 RNGs. 
See :any:`RNGCache <utilfunctions>` for the expected interface. A different RNG can be selected via :doc:`MARQOV::Config <config>` .
Below that we find the concept of an Elementary Monte Carlo Step (EMCS). We have a default implementation but it is fully user
customizable. See :doc:`MARQOV::EMCS <updates/emcs>` for more.
Next come the concept of our moves, the spin-flip metropolis update in :doc:`MARQOV::Metropolis <updates/metropolis>` and the Wolff cluster update
in :doc:`MARQOV::Wolff <updates/wolff>`.
The Wolff update can be further customized using an :doc:`Embedder <updates/wolff>`
These are the classes that most directly use 
the concept of a Hamiltonian and a Lattice. An example implementing the concept of a Hamiltonian is e.g. the
familiar :doc:`Ising Modell <hamiltonians/ising>` modell that implements a very small part of the expected interface and more
examples can be found in :doc:`Hamiltonians <hamiltonians>`.
Examples for the various lattices can be found in :doc:`Lattices <lattices>`.
The final concept that of course has to be there are the physical observables. Examples of some reusable default observables
can be found in :doc:`Observables <observables>`
