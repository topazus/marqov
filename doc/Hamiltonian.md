
# The General Hamiltonian
Our General Hamiltonian is of the form


```math
\begin{aligned}
 \mathcal{H} &= \mathcal{H}_\text{int}
              + \mathcal{H}_\text{self}
              + \mathcal{V}_\text{pot}  \\
             &= \sum\limits_{\alpha=1}^{N_\alpha} \mathcal{H}^{(\alpha)}
              + \sum\limits_{\beta=1}^{N_\beta}   \mathcal{H}^{(\beta)}
              + \sum\limits_{\gamma=1}^{N_\gamma} \mathcal{V}^{(\gamma)} \\
             &= \sum\limits_{\alpha=1}^{N_\alpha}J^{(\alpha)} \sum\limits_{\alpha\text{-bonds: }\langle i,j\rangle}  \phi_i \, D^{(\alpha)}_{ij}(\phi_i,\phi_j)\, \phi_j
              + \sum\limits_{\beta=1}^{N_\beta}h^{(\beta)} \sum\limits_{\beta\text{-sites: }i}  g^{(\beta)} (\phi_i)
              + \sum\limits_{\gamma=1}^{N_\gamma}k^{(\gamma)} \sum\limits_{\gamma\text{-sites: }i}  v^{(\gamma)} \big(\phi_i,\ldots\big)
\end{aligned}
```
It is split into three parts: interaction terms, on-site single particle contributions and contributions due to a potential. The greek indices describe distinct *families* of bonds ($`\alpha`$) and spins ($`\beta`$, $`\gamma`$), respectively. The following figure describes this concept in some more detail. In panel (a) there is only one type (family) of sites (open circles), but two families of bonds, in this case referring to nearest (black lines) and next-nearest neighbour connections (blue lines). In the second panel (b) we see three different species of sites (red, blue, green), translating to three $`\beta`$-families in our framework. Finally, panel (c) shows a typical example of a disordered lattice, with only one family of either bonds and sites.


<p align="center">
<img src="grids.svg" width="40%" alt="Grids" class="center"></p>


Technical notes:
* The first term does only support interactions which can be cast into the given matrix-vector form and therefore *does not* support more general two-site interactions $`f(\phi_i,\phi
  )`$.

* The dimension of the target set of $`g^{(\beta)}`$ needs to match the dimension of the global constant $`h^{(\beta)}`$ in order to be able to perform a meaningful scalar product. Compare the XXZ Heisenberg model below for an example.


# Examples

The following Hamiltonians fit into this framework with the respective substitutions given:


* Example 1: **Ferromagnetic Ising model in homogeneous external field**

    In the celebrated Ising model, individual spins are binary variables on a regular lattice in the presence of an external field. It represents the most basic toy model for ferromagnetic behaviour if the interaction constant $`J`$ is positive.

    ```math
    \mathcal{H} = -J\sum\limits_{\langle i,j\rangle} \phi_i\phi_j + h\sum\limits_{i} \phi_i \qquad\text{where}\quad \phi_i\in\{-1,+1\}
    ```
    In our framework this translates to one interaction term and one on-site term, which exhibit rather simple structure:

    ```math
    N_\alpha=N_\beta = 1, \quad J^{(1)} = -J, \quad h^{(1)} = h
    ```

    ```math
    D_{ij} = 1, \quad g^{(1)} = \phi_i
    ```

* Example 2: **Diluted Heisenberg model**


    In this model, spins are three-dimensional unit vectors, placed on a regular lattice with quenched random vacancies or impurities. The system obeys a $`O(3)`$ symmetry, which means that the Hamiltonian is invariant under three-dimensional rotations in spin space. It is formally given by:

    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle} \epsilon_i \epsilon_j \vec{\phi_i}\vec{\phi_j}
                    \quad\text{where}\quad \epsilon_i,\epsilon_j\in\{0,1\} \quad\text{and}\quad \vec{\phi}_i\in\mathbb{R}^3, \quad |\vec{\phi}_i|^2=1,
    ```
    Here, binary variables $`\epsilon`$ are employed to encode whether the corresponding site is magnetic ($`\epsilon=1`$) or a non-interacting impurity ($`\epsilon=0`$). Note that there is some redundancy in this notation. If the lattice (which is the set of sites $`i`$ and the set of bonds $`\langle i,j\rangle`$) does not include impurity sites in the first place, the $`\epsilon`$'s can be discarded.

    ```math
    N_\alpha=1,  \quad J^{(1)} = J, \quad D_{ij} =1
    ```


* Example 3: **Potts model**

    The $`q`$-state Potts model represents a generalization of the Ising model, where spins can take on $`q`$ values, rather than only two. The Hamiltonian is given by
    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle}\delta_{\phi_i,\phi_j} \qquad\text{where}\quad \phi_i\in\{0,1,\ldots,q-1\}
    ```
    where $`\delta`$ represtents the Kronecker delta. The 2D ferromagnetic Potts model is known to exhibit a second-order phase transition for $`q\leq 4`$ and a first-order phase transition for $`q\geq 5`$.<br><br>

* Example 4: **Clock model**

    The $`q`$-state Clock model reads
    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle}\cos({\phi_i-\phi_j}) \qquad\text{where}\quad \theta_i = 2\pi k/q\quad \text{and} \quad k \in\{0,1,\ldots,q-1\}
    ```
    In 2D it experiences a BKT transition for $`q\geq 5`$, whereas the clock model for $`q=4`$ comprises two sets of the Ising model and the 3-state clock model is equivalent to the 3-state Potts model. The clock model for $`q=2`$ is simply the Ising model. Furthermore, sending $`q\to\infty`$ we arrive at the XY or O(2) model.

    Note that the model matches the form of the interaction in our general Hamiltonian, since by means of trigonometric identities, it can be rewritten as
    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle}  \vec{\phi_i}\vec{\phi_j}
                    \quad\text{where}\quad \phi_i=(\cos\theta_i, \sin\theta_i)
    ```



* Example 5: **Spin glasses**

    Altough seemingly related to the standard Ising or $`O(N)`$ model, the following system is neither ferro- nor antiferromagnetic but rather becomes a prototypical spin glass, as the $`\epsilon_{ij}`$ take on random positive *and* negative values.

    ```math
    \mathcal{H} = J\sum\limits_{\langle i,j\rangle} \epsilon_{ij} \phi_i \phi_j
                    \quad\text{where}\quad \epsilon_i,\epsilon_j\in\left[-1,1\right] \text{(random)}
    ```



    ```math
    N_\alpha=1 \quad \phi_i\in\mathbb{R}^N, \quad |\phi_i|^2=1, \quad J^{(1)} = 1
    ```
    Note that for *quenched disorder* the interaction strengths are read from a table:    
    ```math
    D_{ij} = \mathrm{diag} \, D(i,j) \quad \text{tabulated}
    ```


* Example 6: **Spin ice with dipolar interactions**

    The term *spin ice* refers to substances, in which the disorder of the magnetic moments at low temperatures is analogous to the proton disorder in water ice. The minimal model is the following:


    ```math
    \mathcal{H} =   - J \sum\limits_{\langle i,j\rangle}\phi_i \phi_j
                    + Da^3 \sum\limits_{i>j}\Bigg[\frac{\phi_i\phi_j}{|\vec{r}_{ij}|^3}
                    - 3 \frac{(\vec{\phi}_i\cdot \vec{r}_{ij})(\vec{\phi}_j\cdot \vec{r}_{ij})}{|\vec{r}_{ij}|^5}\Bigg]
    ```
    This Hamiltonian is typically defined on a *pyrochlore* lattice. The first term represents a regular nearest neighbour Ising interaction. The second term, however, introduces long-range dipolar interactions and enters via a suitable potential energy term in our framework. Note that for this kind of interactions the sum can be efficiently computed using the Ewald summation technique.


    ```math
    N_\gamma=1, \quad k^{(1)} = Da^3
    ```
    ```math
    v^{(1)} = v^{(1)}  \Big(  \phi_i,\{\phi_{j\neq i}\},\{\vec{r}_{ij}\}  \Big)
            = \frac{1}{2}\sum\limits_{j\neq i}\Bigg[\frac{\phi_i\phi_j}{|\vec{r}_{ij}|^3}
                    - 3 \frac{(\vec{\phi}_i\cdot \vec{r}_{ij})(\vec{\phi}_j\cdot \vec{r}_{ij})}{|\vec{r}_{ij}|^5}\Bigg]
    ```



* Example 7: **XXZ Heisenberg antiferromagnet in Y external field**


    ```math
    N_\alpha=N_\beta = 1, \quad \phi_i\in\mathbb{R}^3, \quad |\phi_i|^2=1
    ```
    This example demonstrates how direction-dependent couplings can be introduced.


    ```math
    J^{(1)} = -1, \qquad h^{(1)} = (0, h_y, 0)^T, \qquad D_{ij}=\begin{pmatrix} J_x & 0 & 0 \\ 0 & J_x & 0 \\ 0 & 0 & J_z \end{pmatrix}
    ```

    ```math    
    g^{(1)}:\phi_i\mapsto\mathbb{R}^3, \qquad g^{(1)}(\phi_i) = \phi_i

    ```

* Example 8: **Blume-Capel model**

    ```math
    \mathcal{H} = -J\sum\limits_{\langle i,j\rangle} \phi_i\phi_j + D\sum\limits_{i} \phi_i^2\qquad\text{where}\qquad \phi_i\in\{-1,0,+1\}
    ```

    ```math
    N_\alpha=N_\beta = 1, \quad , \quad J^{(1)} = J, \quad h^{(1)} = D, \quad
    D_{ij}=1 \quad g^{(1)} = \phi_i^2
    ```

* Example 9: **Improved O(N) lattice field theory**

    A model with standard interaction and *two* on-site terms

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
    g^{(1)} = \phi_i^2, \quad g^{(2)} = (\phi_i^2-1)^2
    ```    
    Note that the families $`\beta=1`$ und $`\beta=2`$ may of course contain the same sites.<br><br>


* Example 10: **N-color Ashkin-Teller model**
    ```math
    \mathcal{H} =   - J\sum\limits_{\alpha=1}^N\sum\limits_{\langle i,j\rangle} \phi^\alpha_i\phi^\alpha_j
                    - K \sum\limits_{\alpha<\beta}\sum\limits_{\langle i,j\rangle} \phi^\alpha_i\phi^\alpha_j \phi^\beta_i\phi^\beta_j
    ```  
    At a first glance, an implementation of this Hamiltonian is not straightforwardly possible based on the above general form, as it requires an interaction term that couples spins of *different* families. However, one can use standard embedding techniques to numerical handle the Ashkin-Teller as $N$ Ising models, coupled via their interaction strengths.  <br><br>


* Example 11: **Disordered Baxter-Wu model**

        ```math
        \mathcal{H} = -J\sum\limits_{\langle i,j,k\rangle} \epsilon_i \epsilon_j \epsilon_k \phi_i\phi_j\phi_k,
                        \quad\text{where}\quad \epsilon_i,\epsilon_j,\epsilon_k\in\{0,1\}
        ```
        defined on a triangular two-dimensional lattice, where particles are e.g. regular Ising spins...


        *Is this model covered? The third interaction partner might be sneaked in via the interaction D?*
