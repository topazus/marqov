# The General Hamiltonian
Our General Hamiltonian looks like this

```math
\begin{aligned}
 \mathcal{H} &= \sum\limits_{\alpha = 1}^{N_\alpha} \mathcal{H}^{(\alpha)} 
              + \sum\limits_{\beta = 1}^{N_\beta} \mathcal{H}^{(\beta)} \\
             &= \sum\limits_{\alpha = 1}^{N_\alpha}J^{(\alpha)} \sum\limits_{\alpha\text{-bonds: }\langle i,j\rangle}  f^{(\alpha)} (\phi_i,\phi_j,D_{ij})
              + \sum\limits_{\beta = 1}^{N_\beta}h^{(\beta)} \sum\limits_{\beta\text{-sites: }i}  g^{(\beta)} (\phi_i)
\end{aligned}
```

it is split in interaction terms and on-site single particle contributions.
# Examples

The following Hamiltonians fit into this framework with the respective substitutions given:


* Example 1: **Ferromagnetic Ising model in homogeneous external field**

    ```math
    \mathcal{H} = -J\sum\limits_{\langle i,j\rangle} \phi_i\phi_j + h\sum\limits_{i} \phi_i
    ```

    ```math
    N_\alpha=N_\beta = 1, \quad \phi_i\in\{-1,+1\}, \quad J^{(1)} = -J, \quad h^{(1)} = h
    ```
    
    ```math
    f^{(1)} = \phi_i \phi_j, \quad g^{(1)} = \phi_i 
    ```

* Example 2: **Diluted Heisenberg antiferromagnet**
    
    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle} \epsilon_i \epsilon_j \vec{\phi_i}\vec{\phi_j}
                    \quad\text{where}\quad \epsilon_i,\epsilon_j\in\{0,1\}
    ```
    Note that there is some redundancy in the notation. If the set of bonds {\langle i,j\rangle} were to include only bonds that actually exist, the \epsilon's can be discarded.

    ```math
    N_\alpha=1 \quad \phi_i\in\mathbb{R}^3, \quad |\phi_i|^2=1, \quad J^{(1)} = (J,J,J)^T
    ```
    
    ```math
    f^{(1)} = \phi_i^T \phi_j
    ```
    
    

* Example 3: **Prototypocal spin glass**
    
    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle} \epsilon_i \epsilon_j \phi_i \phi_j
                    \quad\text{where}\quad \epsilon_i,\epsilon_j\in\left[-1,1\right]
    ```

    ```math
    N_\alpha=1 \quad \phi_i\in\mathbb{R}^N, \quad |\phi_i|^2=1, \quad J^{(1)} = (1,1,\ldots,1)^T
    ```
    
    ```math
    f^{(1)} = D_{ij} \phi_i^T \phi_j \quad\text{where } D_{ij} \text{ are tabulated scalar values}
    ```

* Example 4: **XXZ Heisenberg in Y external field**


    ```math
    N_\alpha=N_\beta = 1, \quad \phi_i\in\mathbb{R}^3, \quad |\phi_i|^2=1
    ```
    ```math
    J^{(1)} = (1,1,1)^T, \quad h^{(1)} = (0, h_y, 0)^T
    ```
    ```math
    f^{(1)} = \phi_i^T D_{ij} \phi_j, \quad\text{where}\quad D_{ij}=\begin{pmatrix} J_x & 0 & 0 \\ 0 & J_x & 0 \\ 0 & 0 & J_z \end{pmatrix} \\
    \quad g^{(1)} = \phi_i 
    
    ```

* Example 5: **Blume-Capel model**
 
    ```math
    \mathcal{H} = -J\sum\limits_{\langle i,j\rangle} \phi_i\phi_j + D\sum\limits_{i} \phi_i^2
    ```
    
    ```math
    N_\alpha=N_\beta = 1, \quad \phi_i\in\{-1,0,+1\}, \quad J^{(1)} = J, \quad h^{(1)} = D
    ```
    
    ```math
    f^{(1)} = \phi_i \phi_j, \quad g^{(1)} = \phi_i^2 
    ```

* Example 6: **Improved O(N) lattice field theory**

    ```math
    \mathcal{H} = -\beta\sum\limits_{\langle i,j\rangle} \phi_i \phi_j + m\sum\limits_i \phi_i^2 + \lambda\sum\limits_i (\phi_i^2-1)^2
    ```

    ```math
    N_\alpha=1, \quad N_\beta=2, \quad \phi_i\in\mathbb{R}^N
    ```
    ```math    
    \quad J^{(1)} = -\beta, \quad h^{(1)}=m, \quad h^{(2)}=\lambda
    ```    
    ```math
    f^{(1)} = \phi_i^T \phi_j, \quad g^{(1)} = \phi_i^2, \quad g^{(2)} = (\phi_i^2-1)^2 
    ```    
    all sites in the lattice are included for \beta=1 und \beta=2
    

* Example 7: **Clock model**

    ```math
    N_\alpha= 1, \quad \phi_i\in\{0,1,\ldots,p-1\}, \quad J^{(1)} = J
    ```
    ```math
    f^{(1)} = \cos\Big(2\pi(\phi_i- \phi_j)/p\Big)
    ```
    
* Example 8: **Disordered Baxter-Wu model**

    ```math
    \mathcal{H} = -J\sum\limits_{\langle i,j,k\rangle} \epsilon_i \epsilon_j \epsilon_k \phi_i\phi_j\phi_k, 
                    \quad\text{where}\quad \epsilon_i,\epsilon_j,\epsilon_k\in\{0,1\}
    ```

    ...
    the third interaction partner might be sneaked in via the interaction D?
