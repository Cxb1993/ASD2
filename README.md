# ASD2
A Simulator of Droplet (Unfinished)

# Current Status

* 2D : Seems to be working
* 2D Axisymmetric : It's running but it's not right answer
* 3D : Failed to run

# Methods
 
## Convection term / Level Set : WENO method
## Level Set Reinitialization

* Sussman's work (1999)
* Subcell Fix 

## Viscous term 

* 2D case : second-order differencing
* 2D Axisymmetric : Use factorization technique
* 3D : Use poisson equation to solve viscous term

## Poisson's Solver

* Conjugate Gradient Method : Use Intel MKL BLAS
* BiConjugate Gradient Stabilized Method : Use Intel MKL BLAS

## Dealing with jump between fluids
* Viscous term : use Heaviside function to make implicit solver
* Pressure term : use Ghost Fluid Method

# References
 
## WENO method

* Hu, Changqing, and Chi-Wang Shu. "Weighted essentially non-oscillatory schemes on triangular meshes." Journal of Computational Physics 150.1 (1999): 97-127.
 
## Level Set Reinitialization

* Sussman's Method : Sussman, Mark, Peter Smereka, and Stanley Osher. "A level set approach for computing solutions to incompressible two-phase flow." Journal of Computational physics 114.1 (1994): 146-159.
* Subell Fix (1) :  Min, Chohong, and Frédéric Gibou. "A second order accurate level set method on non-graded adaptive cartesian grids." Journal of Computational Physics 225.1 (2007): 300-321.
* Subell Fix (2) :  Min, Chohong. "On reinitializing level set functions." Journal of computational physics 229.8 (2010): 2764-2772.

## Momentum equation
* 2D Axisymmetric : Kang, Myungjoo, Hyeseon Shim, and Stanley Osher. "Level set based simulations of two-phase oil–water flows in pipes." Journal of Scientific Computing 31.1-2 (2007): 153-184.
* 3D : Lee, Byungjoon, and Myungjoo Kang. "Full 3D Simulations of Two-Phase Core–Annular Flow in Horizontal Pipe Using Level Set Method." Journal of Scientific Computing (2015): 1-27.

## Ghost Fluid Method
* Kang, Myungjoo, Ronald P. Fedkiw, and Xu-Dong Liu. "A boundary condition capturing method for multiphase incompressible flow." Journal of Scientific Computing 15.3 (2000): 323-360.
* Liu, Xu-Dong, Ronald P. Fedkiw, and Myungjoo Kang. "A boundary condition capturing method for Poisson's equation on irregular domains." Journal of computational Physics 160.1 (2000): 151-178.
* Hong, Jeong-Mo, and Chang-Hun Kim. "Discontinuous fluids." ACM Transactions on Graphics (TOG) 24.3 (2005): 915-920.
* Ng, Yen Ting, et al. "Guidelines for Poisson solvers on irregular domains with Dirichlet boundary conditions using the ghost fluid method." Journal of Scientific Computing 41.2 (2009): 300-320.
* Yang, Jianming, and Frederick Stern. "Sharp interface immersed-boundary/level-set method for wave–body interactions." Journal of Computational Physics 228.17 (2009): 6590-6616.

## Curvature Discretization (Future Work)
* Ervik, Åsmund, Karl Yngve Lervåg, and Svend Tollak Munkejord. "A robust method for calculating interface curvature and normal vectors using an extracted local level set." Journal of Computational Physics 257 (2014): 259-277.
* Lervåg, Karl Yngve, Bernhard Müller, and Svend Tollak Munkejord. "Calculation of the interface curvature and normal vector with the level-set method." Computers & Fluids 84 (2013): 218-230.
* Macklin, Paul, and John Lowengrub. "Evolving interfaces via gradients of geometry-dependent interior Poisson problems: application to tumor growth." Journal of Computational Physics 203.1 (2005): 191-220.
* Macklin, Paul, and John S. Lowengrub. "A new ghost cell/level set method for moving boundary problems."

 
