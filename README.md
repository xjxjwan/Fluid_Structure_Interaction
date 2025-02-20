# Fluid-structure interaction using the Ghost Fluid Method
This project involves solving the two-dimensional compressible unsteady Euler equations, both as a stand-alone system, and interacting with a rigid body. Following successful completion of the project, the student will have a two-dimensional code capable of simulating flow around either a stationary or moving rigid body.
Structure:
1. 1D_RiemannGFM: 1D Riemann problem-based GFM, with visualize.py.
2. 2D_RiemannGFM: 2D Riemann problem-based GFM, with visualize2D.py for simulation results visualization.
3. 2D_MUSCL_Hancock: 2D second-order Riemann problem-based solver (single material).
4. 2D_SLIC_FORCE: 2D second-order (SLIC+FORCE) solver (single material).
5. Rigid_Body: 2D Riemann problem-based GFM-like method for simulating fluid interaction with rigid body.
6. In each folder, the 'res' folder contains simulation results, including simulation data and figures.
