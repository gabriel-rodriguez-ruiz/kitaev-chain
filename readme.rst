Welcome

This is my readme file for a project on Kitaev chains and their topological properties.
The idea is to solve the toy model with periodic boundary conditions and with a finite length of the chain.

The tight binding model is simulated with the kwant package.
The finite chain follows the following course:
https://topocondmat.org/w1_topointro/1D.html

The typical workflow with Kwant is as follows:

1) Create an “empty” tight binding system.
2) Set its matrix elements and hoppings.
3) Attach leads (tight binding systems with translational symmetry).
4) Pass the finalized system to a solver.
