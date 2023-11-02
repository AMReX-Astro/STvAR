# Space-Time Variable Code Generator and Solver

This code package contains a Sympy-based code generator that creates C++ code for solving PDEs with finite-differencing.

The integration of these PDEs in time is implemented using the Method of Lines with a Runge-Kutta integrator.

We discretize in space using finite differencing with options for various orders. We implement domain decomposition using the AMReX package for adaptive mesh refinement and parallel processing.

A PDE solver generated using the STvAR package will run on CPUs and GPUs using the AMReX performance portability backend.

We include examples based on the wave equation and numerical general relativity.
