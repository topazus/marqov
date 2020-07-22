# I/O
## Structure
The Layout of our HDF5 containers is like this
- root
  - step1
    - Environment
    - Config
      - MarqovConfig
      - Lattice
      - Hamiltonian
    - Observables
    - State (for restart. final statespace + RNG)
  - step2(Every restart generates a new step)
  - step3
  - ...

Hence the root is a series of step. Each step consists of an simulation, hence the Config, the Observables, and a State-Dump.
the Config has Marqovconfig + LatticeConfig + HamiltonianConfig. the Observables are a series of observable time series. the State-Part is a state-dump that would be consumed if it serves as input for the next step.

### Environment
Information about the Environment of the host system.

### Config
Everything that is required for the simulation and is not deemed to be some part of statespace

#### Marqovconfig
All parameters from the MarqovConfig Object

#### Lattice
Information about the lattice

#### Hamiltonian
The parameters of the Hamiltonian. using a paramnames function in the Hamiltonian the keys are modifiable.
Else generic names, param<i>, are used

### Observables
Here the time series of the observables follow. Every Observables has its own time series.

### State
The final state space of the system. THis consititutes the Monte Carlo state space as well as the state of the used RNG.

## Design rationale
The intent was to have a series of immutable steps where an independent party is able to verify the 
performed calculations.
