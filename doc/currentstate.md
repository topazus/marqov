# Already implemented

* Spin-1/2 Ising model: Implemented in `Ising.h`; Metropolis, specialized and general Wolff cluster algorithm available, all three work fine; tested in 2D and 3D via critical scaling collapses

* O(3) model: Implemented in `Heisenberg.h`; Metropolis, specialized and general Wolff cluster algorithm available, all three work fine; tested in 3D via critical scaling collapses; no phase transition in 2D; reference values can be found in `campostrini2002` and `hasenbusch2011`

* O(3) with $`\phi^4`$ interaction: Implemented in `Phi4.h`; correct critical exponents, but critical temperature is off (see corresponding Issue for details); main references are `hasenbusch2001`, `campostrini2002`, `hasenbusch2011`

* Blume-Capel model: Implemented in `BlumeCapel.h`; Extension to integer spin values $`\{-1,0,1\}`$, which calls for a hybrid MC scheme, as the Wolff algorithm does not affect zero-spins; current implementation of the Metropolis and Wolff algorithms are only obvious guesses at the moment, i.e. as it stand it is unclear whether updates obey the detailed balance condition strictly; nonetheless, tested in 3D via critical scaling collapses for a specific parameter choice; results match expected numbers from `hasenbusch2010`, `hasenbusch2018`

# Up next

* something from the XXZ antiferromagnet family