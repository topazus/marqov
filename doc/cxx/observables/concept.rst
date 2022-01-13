.. Copyright (c) 2022, Manuel Schrauth, Florian Goth

The concept of an observable
==============================
An observable can be any object that has a measure function that returns something and takes a state space and a lattice as argument. For slightly more verbose I/O we require the presence 
of a name, and optionally an extended description.
MARQOV takes care of generating suitable I/O routines from the signature of this measure function.
Hence the bare minimum of an observable would be

.. code-block:: cpp

    class FirstObs
    {
        public:
            std::string name{"FirstObs"}; ///< The name of the observable
            std::string desc{"My first obs"}; ///< A helpful description that will be used in the HDF5 output files.

            template <class StateSpace, class Grid>
            double measure(const StateSpace& statespace, const Grid& grid)
            {
                return 42;
            }
    };
