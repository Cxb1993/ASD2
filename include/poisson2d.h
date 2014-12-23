#ifndef __POISSON_H_
#define __POISSON_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <unordered_map>

// MKL for Poisson equation solver
#include <mkl.h>
// for data allocation
#include <mkl_rci.h>
// for poisson equation
#include <mkl_poisson.h>
#include <mkl_dfti.h>
// for matrix and vector operation
#include <mkl_pblas.h>

#include "common.h"
#include "data.h"
#include "bc.h"
class PoissonSolver {
public:
	int GS_1FUniform_3D(double ***phi, double ***rhs, int nx, int ny, int nz, int num_bc_grid,
		double dx, double dy, double dz, std::shared_ptr<BoundaryCondition> PBC);
	int MKL_1FUniform_3D(double ***phi, double ***rhs,
		int nx, int ny, int nz, int num_bc_grid,
		double lenX, double lenY, double lenZ,
		double dx, double dy, double dz, std::shared_ptr<BoundaryCondition> PBC);
	int CG_1FUniform_3D(double ***phi, double ***rhs,
		int nx, int ny, int nz, int num_bc_grid,
		double lenX, double lenY, double lenZ,
		double dx, double dy, double dz, std::shared_ptr<BoundaryCondition> PBC);
	int RCICG_1FUniform_3D(double ***phi, double ***rhs,
		int nx, int ny, int nz, int num_bc_grid,
		double lenX, double lenY, double lenZ,
		double dx, double dy, double dz, std::shared_ptr<BoundaryCondition> PBC);
};

#endif __POISSON_H_