.. Copyright (c) 2021-2022, Manuel Schrauth, Florian Goth



The Problem
============

The general Hamiltonian builds the centerpiece of MARQOV, and allows to cover a large number of physical models and general geometries. It governs the physical interactions in the system and is formally split into three parts
$$
\mathcal{H} = \mathcal{H}\text{int} + \mathcal{H}\text{self} + \mathcal{V}_\text{pot},
$$
interaction terms, on-site single particle contributions and contributions due to a potential, which read
$$
\begin{aligned}
\mathcal{H}\text{int} &= \sum\limits{\alpha=1}^{N_\alpha}J^{(\alpha)} \sum\limits_{\alpha\text{-bonds: }i,j}  \phi_i , D^{(\alpha)}{ij}(\phi_i,\phi_j), \phi_j\
\mathcal{H}\text{self} &= \sum\limits_{\beta=1}^{N_\beta}h^{(\beta)} \sum\limits_{\beta\text{-sites: }i}  g^{(\beta)} (\phi_i)\
\mathcal{V}\text{pot} &= \sum\limits{\gamma=1}^{N_\gamma}k^{(\gamma)} \sum\limits_{\gamma\text{-sites: }i}  v^{(\gamma)} \big(\phi_i,\ldots\big)
\end{aligned}
$$
The configuration (or state) of the individual spins, which are placed on discrete lattice positions \(i\), is denoted by the symbol \( \phi\)1. This covers, e.g. Ising spins where \(\phi\in\pm 1\) and \(N\)-vector models, where \(\phi\in\mathbb{R}^N\). The set of all individual spin configurations hence constitutes the state space of the physical system.
Note that the first term in the above Hamiltonian does only support interactions which can be cast into the given matrix-vector form. More general two-site interactions \(f(\phi_i,\phi_j)\) are not natively covered. However, we provide convenient interfaces where users can nevertheless implement more general interaction terms, though, at the cost of not being able to access the optimized set of standard Monte Carlo update algorithms that come along with MARQOV.
The on-site term allows for full flexibility in implementing self-interactions of any functional form, with the only requirement that the dimension of the target set of global pre-factor \(g^{(\beta)}\) needs to match the dimension of the global constant \(h^{(\beta)}\) in order to provide a meaningful scalar product. Compare the XXZ Heisenberg model below for an example.
The potential term is used to account for any physical interactions which go beyond two-site and self-interactions. Typical examples are long-range Coulomb interactions in a spin ice system or plaquette interactions.

One particular strength of MARQOV is that it is designed for supporting irregular geometries.
To that end we can utilize arbitrary lattice classes given by the user.

Finally we have observables. These are user-defined classes model functions 
on the state space of the Monte Carlo simulation.
