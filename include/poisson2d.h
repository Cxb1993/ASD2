#ifndef __POISSON2D_H_
#define __POISSON2D_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

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
#include "bc2d.h"

class PoissonSolver2D {
	int kNx, kNy, kNumBCGrid;
public:
	PoissonSolver2D(int nx, int ny, int num_bc_grid);

	int GS_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
		std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
		const double lenX, const double lenY, const double dx, const double dy,
		const std::shared_ptr<BoundaryCondition2D>& PBC, const int maxiter);
	int MKL_2FUniform_2D(std::vector<double>& phi, const std::vector<double>& rhs,
		double lenX, double lenY, double dx, double dy, std::shared_ptr<BoundaryCondition2D> PBC);
	
	int CG_2FUniformP_2D(std::vector<double>& ps, const std::vector<double>& rhs,
		std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
		std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
		const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxIter);
	int CG_2FUniform_2D(double *x, double *b,
		std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
		std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
		const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxIter);

	int BiCGStab_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
		std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
		std::vector<double>& M, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx, 
		const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxiter);

	bool IsSymmetric(std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx);

	std::vector<double> InvertMatrixDiagonal(const std::vector<double>& M);
	int idx(int i, int j);
};

#endif __POISSON2D_H_