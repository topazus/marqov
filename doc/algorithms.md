## General notes

* Wolff and SW cluster for RBIM and RFIM discussed in `dotsenko1991` and `newman1996`
* Cluster representations and the Wolff algorithm in arbitrary external fields; includes a C++ header library for most common models: `kent2018`
* Detailed analysis of Hybrid Algorithms (i.e. Metropolis+Wolff) for Spin-1/2 and Spin-3/2 Ising model: `plascak2002`
* Mixed cluster algorithm for 3-state Ising model: `bouabci1996`
* Demon updates `creutz1983`, `rummukainen1993`
* Population Annealing `shchur2018`
 



# Invaded cluster algorithm

#### Key facts:

* not a general purpose algorithm (can only be used to simulate the model at the critical point)
* locates the critical point all on its own and very fast ("negative feedback")
* variant of the Swendsen-Wang algorithm, where clusters are built according to a cluster invasion percolation scheme
* does not sample the Boltzmann distribution exactly

#### Prerequisites

* understand [invasion percolation](http://www.physics.purdue.edu/flow/MMproject/Wilkinson1983.pdf)
* understand [Swendsen-Wang cluster flipping mechanics](https://en.wikipedia.org/wiki/Swendsen%E2%80%93Wang_algorithm)

#### Main literature


* [Seminal paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.2792) by J. Machta, Y. S. Choi, A. Lucke, T. Schweizer, and L. V. Chayes (1995)
* Introduction can also be found [here](https://arxiv.org/pdf/cond-mat/9703179.pdf)
* [Invaded cluster algorithm for Potts model](https://www.math.ucla.edu/~lchayes/lchayes_website_folder/old_publications_folder/ic_potts_96.pdf)
* [XY modell](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.026702)
* [Diluted Potts model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.011103)
* [Frustrated Ising model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.57.88)

#### Parallelization and optimization

* [Parallel Invaded Cluster Algorithm for the Ising Model](https://arxiv.org/pdf/cond-mat/9806127.pdf)
* [Fast cluster invasion](https://hal.archives-ouvertes.fr/hal-01653926/document)


#### Prerequisites

* understand [invasion percolation](http://www.physics.purdue.edu/flow/MMproject/Wilkinson1983.pdf)
* understand [Swendsen-Wang cluster flipping mechanics](https://en.wikipedia.org/wiki/Swendsen%E2%80%93Wang_algorithm)

#### Main literature


* [Seminal paper](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.75.2792) by J. Machta, Y. S. Choi, A. Lucke, T. Schweizer, and L. V. Chayes (1995)
* Introduction can also be found [here](https://arxiv.org/pdf/cond-mat/9703179.pdf)
* [Invaded cluster algorithm for Potts model](https://www.math.ucla.edu/~lchayes/lchayes_website_folder/old_publications_folder/ic_potts_96.pdf)
* [XY modell](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.026702)
* [Diluted Potts model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.76.011103)
* [Frustrated Ising model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.57.88)

#### Parallelization and optimization

* [Parallel Invaded Cluster Algorithm for the Ising Model](https://arxiv.org/pdf/cond-mat/9806127.pdf)
* [Fast cluster invasion](https://hal.archives-ouvertes.fr/hal-01653926/document)




# Worm algorithm

The Worm algorithm does not operate on the space of spin configuration but is rather a diagrammatic sampling method, operating on the space of bond (in the spirit of a high-temperature expansion).

#### Prerequisities

Understand high-temperature expansion (HTE) and (associated) diagrammatic loop counting

* [Lecture on Series Expansions](https://www.youtube.com/watch?v=bMnpf0s-mAk)
* [Corresponding notes](https://ocw.mit.edu/courses/physics/8-334-statistical-mechanics-ii-statistical-physics-of-fields-spring-2014/lecture-notes/)
* [Graphical Enumeration Techniques: Series Expansions and Animal Problems](https://www.researchgate.net/publication/318566213_Graphical_Enumeration_Techniques_Series_Expansions_and_Animal_Problems)
* Most standard textbooks


#### Literature

* [Original paper by Prokof’ev and Boris Svistunov (2001)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.87.160601)
* [WA for the CP^n-1 model (also nice introduction for the Ising case)](https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/231782/1/vetter_worm_cpn-1.pdf)
* [New developments (2018: Lifted worm algorithm for the Ising model](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042126)

#### Code

* [WA implementation for the 2D Ising model (Python)](https://github.com/saforem2/worm_algorithm)
* [WA implementation for the J-Current model (C++)](http://mcwa.csi.cuny.edu/umass/jcurrent.html)
* [WA visualization for the J-Current model (Mathematica)](https://demonstrations.wolfram.com/WormAlgorithmForJCurrentModel/)

#### Disorder

* [High-temperature series expansion for spin glasses. I. Derivation of the series](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.36.546)
* [High-temperature series expansion for spin glasses. II. Analysis of the series](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.36.559)

# HMC
Hamiltonian/Hybrid Monte Carlo
 - For the SO(3) Heisenberg model EOMs are known that work directly on the [phase space of spins](https://arxiv.org/pdf/1402.4114.pdf) Hence no need for auxiliary momenta.