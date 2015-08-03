#ifndef __POISSON3D_H_
#define __POISSON3D_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <cstdint>
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
#include "bc3d.h"

class PoissonSolver3D {
    int64_t kNx, kNy, kNz, kNumBCGrid;
public:
    PoissonSolver3D(int nx, int ny, int nz, int num_bc_grid);

    int CG_2FUniformU_3D(std::vector<double>& uhat, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& UBC, const int maxIter);
    int CG_2FUniformV_3D(std::vector<double>& vhat, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& VBC, const int maxIter);
    int CG_2FUniformW_3D(std::vector<double>& what, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& WBC, const int maxIter);
    int CG_2FUniformP_3D(std::vector<double>& ps, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& PBC, const int maxIter);
    
    int CG_2FUniform_3D(double *x, double *b,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx, int64_t Anrows, int64_t size, const int maxIter);

    int BiCGStab_2FUniformU_3D(std::vector<double>& uhat, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& UBC, const int maxIter);
    int BiCGStab_2FUniformV_3D(std::vector<double>& vhat, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& VBC, const int maxIter);
    int BiCGStab_2FUniformW_3D(std::vector<double>& what, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& WBC, const int maxIter);
    int BiCGStab_2FUniformP_3D(std::vector<double>& ps, const std::vector<double>& rhs,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
        const std::shared_ptr<BoundaryCondition3D>& PBC, const int maxIter);
    
    int BiCGStab_2FUniform_3D(double *x, double *b,
        std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
        std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx, int64_t Anrows, int64_t size, const int maxIter);

    std::vector<double> InvertMatrixDiagonal(const std::vector<double>& M);
    int idx(int i, int j, int k);
};

#endif __POISSON3D_H_