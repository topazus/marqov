Getting Started
===============

We load the pyMARQOV core module with this command::

   from pymarqov.core import monte_carlo_run

This allows us to create a analysis object::

   mc = monte_carlo_run("/home/schrauth/marqov/marqov-demo")

and auto-import Monte Carlo data::

   myobs = [Susceptibility(), BinderFourth()]
   mc.auto_import(myobs)


.. warning::
   This is a paragraph which contains some particularly important information.


Let us have a look what the auto import really does

.. currentmodule:: pymarqov.core.monte_carlo_run
.. autofunction:: auto_import
