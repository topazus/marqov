# The General Hamiltonian
Our General Hamiltonian looks like this

```math
\begin{aligned}
 \mathcal{H} &= \sum\limits_{\alpha = 1}^{N_\alpha} \mathcal{H}^{(\alpha)} 
              + \sum\limits_{\beta = 1}^{N_\beta} \mathcal{H}^{(\beta)} \\
             &= \sum\limits_{\alpha = 1}^{N_\alpha}J^{(\alpha)} \sum\limits_{\alpha\text{-bonds: }\langle i,j\rangle}  f^{(\alpha)} (S_i,S_j,D_{ij})
              + \sum\limits_{\beta = 1}^{N_\beta}h^{(\beta)} \sum\limits_{\beta\text{-sites: }i}  g^{(\beta)} (S_i)
\end{aligned}
```
# Examples

The following Hamiltonians fit into this framework with the respective substitutions given:


Example 1: Ising model in homogeneous external field

```math
N_\alpha=N_\beta = 1, \quad S_i\in\{-1,+1\}, \quad J^{(1)} = J, \quad h^{(1)} = h
```

```math
f^{(1)} = S_i S_j, \quad g^{(1)} = S_i 
```

Example 2: Diluted Heisenberg model 

```math
N_\alpha=1 \quad S_i\in\mathbb{R}^3, \quad |S_i|^2=1, \quad J^{(1)} = 1
```

```math
f^{(1)} = D_{ij} S_i^T S_j \quad\text{where } D_{ij} \text{ are tabulated scalar values}
```

Example 3: XXZ Heisenberg in Y external field

```math
N_\alpha=N_\beta = 1, \quad S_i\in\mathbb{R}^3, \quad |S_i|^2=1
```
```math
J^{(1)} = (1,1,1)^T, \quad h^{(1)} = (0, h_y, 0)^T
```
```math
f^{(1)} = S_i^T D_{ij} S_j, \quad\text{where}\quad D_{ij}=\begin{pmatrix} J_x & 0 & 0 \\ 0 & J_x & 0 \\ 0 & 0 & J_z \end{pmatrix} \\
\quad g^{(1)} = S_i 

```

