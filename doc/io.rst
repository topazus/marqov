===========
I/O
===========



Structure
===========
The Layout of our HDF5 containers is like this

  * root

    * step1
  
      * Environment
      * Config
    
        * MarqovConfig
        * Lattice
        * Hamiltonian

      * Observables
      * State (for restart. final statespace + RNG)

    * step2(Every restart generates a new step)
    * step3
    * ...

Hence the root is a series of steps. Each step consists of a simulation, hence the Config, the Observables, and a State-Dump.
The Config has Marqovconfig + LatticeConfig + HamiltonianConfig.
The Observables are a sequence of observable time series. the State-Part is a state-dump that would be consumed if it serves as input for the next step.

It is also possible to describe the HDF5 layout used by MARQOV in Extendend Backus-Naur Form:

MARQOV HDF5 File layout::

  root = steps;
  steps = step | steps;
  step = "/step", number, stepcontent;
  stepcontent = stepconfig, stepenvironment, stepobservables, stepstate;
  stepconfig = "/config", hamiltonianconfig, latticeconfig, marqovconfig;
  stepenvironment = "/environment", usedcode, hostname, startingdate,...;
  stepobservables = "/observables", observables;
  stepstate = "/state", RNG, rngstate, hamiltonianstatespace;
  hamiltonianconfig = "/hamiltonian", hamiltonianname, hamiltonianparams;
  latticeconfig = "/lattice", latticename, latticeparams;
  marqovconfig = "/marqovconfig", marqovparams;


Config
======
Everything that is required for the simulation and is not deemed to be some part of statespace.

Marqovconfig
============
All parameters from the MarqovConfig Object.

Lattice
=========
Information about the lattice.

Hamiltonian
==============
The parameters of the Hamiltonian. Using a paramnames function in the Hamiltonian the keys can optionally be specified by the user.
Else generic names in the form of *param<i>*, are used

Environment
============
Information about the environment of the host system.


  * Used Code + Version
  * Host
  * Time

Observables
============
Here the time series of the observables follow. Every Observable has its own time series.

State
======
The final state space of the system. This consititutes the Monte Carlo state space as well as the state of the used RNG.

Design rationale
=================
The intent was to have a series of immutable steps where an independent party is able to verify the 
performed calculations.
Modelling a calculation as a series of steps offers the possibility that different steps can be done by different programs.
A code with restarts is hence a calculation where the same program is used with the output of the previous step as input.
additional things can be added in the future to the file without affecting the parseability for older codes.
