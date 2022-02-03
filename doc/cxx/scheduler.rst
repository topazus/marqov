.. Copyright (c) 2021, Manuel Schrauth, Florian Goth

MARQOVScheduler
================
Marqov provides a multithreading scheduler with optional support for MPI.

.. doxygenclass:: MARQOV::CXX11Scheduler
   :members:

If MARQOV is compiled with support for MPI the MPI scheduler is utilized.

.. doxygenclass:: MARQOV::MPIScheduler
   :members: