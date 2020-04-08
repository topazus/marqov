# Already implemented

* Spin-1/2 Ising model: Implemented in `Ising.h`; Metropolis, specialized and general Wolff cluster algorithm available, all three work fine; tested in 2D and 3D via critical scaling collapses

* O(3) model: Implemented in `Heisenberg.h`; Metropolis, specialized and general Wolff cluster algorithm available, all three work fine; tested in 3D via critical scaling collapses; no phase transition in 2D; reference values can be found in `campostrini2002` and `hasenbusch2011`

* O(3) with $`\phi^4`$ interaction: Implemented in `Phi4.h`; correct critical exponents, but critical temperature is off (see corresponding Issue for details); main references are `hasenbusch2001`, `campostrini2002`, `hasenbusch2011`

* Blume-Capel model: Implemented in `BlumeCapel.h`; Extension to integer spin values $`\{-1,0,1\}`$, which calls for a hybrid MC scheme, as the Wolff algorithm does not affect zero-spins; current implementation of the Metropolis and Wolff algorithms are only obvious guesses at the moment, i.e. as it stand it is unclear whether updates obey the detailed balance condition strictly; nonetheless, tested in 3D via critical scaling collapses for a specific parameter choice; results match expected numbers from `hasenbusch2010`, `hasenbusch2018`

* XXZ antiferromagnet: Implemented in `XXZAntiferro.h`; Tested at two points in the T/H phase diagram: Found second-order transition from the AF into the paramagnetic phase at $`\beta=0.6384`$, $`H=0`$ and from the SF into the paramagnetic phase at $`\beta=0.98`$, $`H=4.2`$, consistent with recent literature, such as `selke2009`, `selke2011` and `hu2014`; due to lack of precision critical exponents could not be confirmed so far; note: make sure that the Wolff cluster algorithm is only used at zero external field


# Up next

