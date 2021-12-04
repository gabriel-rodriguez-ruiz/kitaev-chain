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

Technical detail:

When building the system, only one direction of hopping is given, i.e. syst[lat(i, j), lat(i, j-1)] = ... and not also syst[lat(i, j-1), lat(i, j)] = .... The reason is that Builder automatically adds the other direction of the hopping such that the resulting system is Hermitian.

However, it does not hurt to define the opposite direction of hopping as well:

syst[lat(1, 0), lat(0, 0)] = -t
syst[lat(0, 0), lat(1, 0)] = -t.conj()
(assuming that t is complex) is perfectly fine. However, be aware that also

syst[lat(1, 0), lat(0, 0)] = -1
syst[lat(0, 0), lat(1, 0)] = -2
is valid code. In the latter case, the hopping syst[lat(1, 0), lat(0, 0)] is overwritten by the last line and also equals to -2.