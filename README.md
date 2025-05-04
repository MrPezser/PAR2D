# PAR2D

A parallel block-structured solver for the 2D Euler equations, with the option to have axysymmetric fluxes.

Finite volume and some rough DGP1 capabilities

Euler time stepping

Solves the conserved equations using density, velocity, and temperature as primatives

LDFSS (~AUSM/VanLeer - type) flux splitting scheme

Termally Perfect Gas using curve fits for finding enthalpy.

Utility programs for a) generating a structured mesh, and b) decomposing a block-structured mesh into smaller structured blocks.

MPI cpu parallelization.


Plans for future development:
- Stability improvements for higher order method
- Improved I/O including an input file and having a distinct working dir
- Show off some old solutions
- Viscous Fluxes
- Thermal nonequilibrium
- GUI
