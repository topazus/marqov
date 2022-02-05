.. Copyright (c) 2022, Manuel Schrauth, Florian Goth

The Heisenberg Example. But threaded!
======================================
The O(3) Heisenberg model is a slightly more advanced model on a continuous state space.
This demo also gives an example of how to inject a self-defined Wolff Cluster update for your self-defined model into MARQOV.
Additionally we show in this example how to use the MARQOV scheduler in order to have automatic scheduling of the jobs
derived from your parameter space sampling.

.. literalinclude:: ../../demo/HeisenbergSimpleThreaded.cpp
   :language: c++
   :linenos:
