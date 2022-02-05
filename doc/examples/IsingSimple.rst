.. Copyright (c) 2022, Manuel Schrauth, Florian Goth

The Ising Simple example
===========================
This is probably the simplest self-defined model.
It shows how to define a Model, use a predefined lattice, and a predefined observable.
It sets up marqov as a single threaded simulation of a single parameter set without any scheduling.

.. literalinclude:: ../../demo/IsingSimple.cpp
   :language: c++
   :linenos:

