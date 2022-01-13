.. Copyright (c) 2022, Manuel Schrauth, Florian Goth

Initializers
================

Initializers can be used to specify how the Metropolis Algorithm generates a local update 
to a State Vector.
We do provide default behaviour, but is also possible to specify the behaviour for your Hamiltonian.
To that end you should make a declaration in the following form visible before MARQOV::Core gets declared:

.. codeblock:: cpp
  template <>
  class Initializer<MyHamiltonian> : public MyInitializer
  {};


.. doxygenfile:: initializers.h

