.. Copyright (c) 2021, Manuel Schrauth, Florian Goth

Prerequisites
=============
MARQOV is light on dependencies:

  * A C++14 compiler
  * HDF5
  * optionally MPI

and cmake for executing the build and configure step.
If you intend to build the documentation doxygen, breathe and sphinx are required.

Compiling the Examples
======================
Assuming you have retrieved the most recent version of `MARQOV <https://git.physik.uni-wuerzburg.de/marqov/marqov>`_
the default way to compile is with the help of `cmake <https://cmake.org/>`_ .

The build step follows the basic cmake procedure:

.. code-block:: shell
   :linenos:
   
   cmake -E make_directory build
   cd build
   cmake -DCMAKE_BUILD_TYPE=Release ..
   cmake --build . --target all --config Release

and you should find compiled examples in src.