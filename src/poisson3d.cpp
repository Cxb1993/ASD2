#include "poisson3d.h"

PoissonSolver3D::PoissonSolver3D(int nx, int ny, int nz, int num_bc_grid) :
	kNx(nx), kNy(ny), kNz(nz), kNumBCGrid(num_bc_grid) {
}

int PoissonSolver3D::CG_2FUniformU_3D(std::vector<double>& uhat, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& UBC, const int maxIter) {

	MKL_INT Anrows = (kNx - 1) * kNy * kNz, Ancols = (kNx - 1) * kNy * kNz;
	MKL_INT size = (kNx - 1) * kNy * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);
	
	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "Uhat(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx - 1; i++) {
		b[i + j * (kNx - 1) + k * (kNx - 1) * kNy] = rhs[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)];
		x[i + j * (kNx - 1) + k * (kNx - 1) * kNy] = 0.0;
	}

	CG_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx - 1; i++) {
		uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] = x[i + j * (kNx - 1) + k * (kNx - 1) * kNy];
		assert(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] == uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)]);
		if (std::isnan(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)]) || std::isinf(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)])) {
			std::cout << "Uhat(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid + 1 << " " << j + kNumBCGrid << " " << k + kNumBCGrid << " "
				<< " " << uhat [idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	UBC->ApplyBC_U_3D(uhat);

	return 0;
}

int PoissonSolver3D::CG_2FUniformV_3D(std::vector<double>& vhat, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& VBC, const int maxIter) {

	MKL_INT Anrows = kNx * (kNy - 1) * kNz, Ancols = kNx * (kNy - 1) * kNz;
	MKL_INT size = kNx * (kNy - 1) * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "Vhat(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy - 1; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * (kNy - 1)] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)];
		x[i + j * kNx + k * kNx * (kNy - 1)] = 0.0;
	}

	CG_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy - 1; j++)
	for (int i = 0; i < kNx; i++) {
		vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)] = x[i + j * kNx + k * kNx * (kNy - 1)];
		assert(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)] == vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)]);
		if (std::isnan(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)]) || std::isinf(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)])) {
			std::cout << "Vhat(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid + 1 << " " << k + kNumBCGrid << " "
				<< " " << vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	VBC->ApplyBC_V_3D(vhat);

	return 0;
}

int PoissonSolver3D::CG_2FUniformW_3D(std::vector<double>& what, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& WBC, const int maxIter) {

	MKL_INT Anrows = kNx * kNy * (kNz - 1), Ancols = kNx * kNy * (kNz - 1);
	MKL_INT size = kNx * kNy * (kNz - 1);

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "What(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz - 1; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * kNy] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)];
		x[i + j * kNx + k * kNx * kNy] = 0.0;
	}

	CG_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz - 1; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] = x[i + j * kNx + k * kNx * kNy];
		assert(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] == what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)]);
		if (std::isnan(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)]) || std::isinf(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)])) {
			std::cout << "What(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid << " " << k + kNumBCGrid + 1 << " "
				<< " " << what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] << std::endl;
			exit(1);
		}
	}

	WBC->ApplyBC_U_3D(what);

	return 0;
}

int PoissonSolver3D::CG_2FUniformP_3D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& PBC, const int maxIter) {

	MKL_INT Anrows = kNx * kNy * kNz, Ancols = kNx * kNy * kNz;
	MKL_INT size = kNx * kNy * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "Pressure(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * kNy] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)];
		x[i + j * kNx + k * kNx * kNy] = 0.0;
	}

	CG_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] = x[i + j * kNx + k * kNx * kNy];
		assert(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] == ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)]);
		if (std::isnan(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)]) || std::isinf(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)])) {
			std::cout << "Pressure(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid << " " << k + kNumBCGrid << " "
				<< " " << ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	PBC->ApplyBC_U_3D(ps);

	return 0;
}

int PoissonSolver3D::CG_2FUniform_3D(double *x, double *b,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx, int64_t Anrows, int64_t size, const int maxIter) {

	double *Ax = Data::Allocate1Dd(size);
	double *p = Data::Allocate1Dd(size);
	double *q = Data::Allocate1Dd(size);
	double *r = Data::Allocate1Dd(size);
	double *z = Data::Allocate1Dd(size);

	double alpha = 1.0, beta = 1.0;
	double rho = 1.0, rho1 = 1.0;
	double rnorm2 = 0.0, bnorm2 = 0.0, resid = 0.0;

	const double err_tol = 1.0e-6;
	int iter = 0;
	bool isConverged = false;

	// get Ax(=A*x), using upper triangular matrix (Sparse BLAS)
	// https://software.intel.com/en-us/node/468560

	char transa = 'n';
	mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), x, Ax);
	// r = b - Ax, initial residual
	// r = b
	// https://software.intel.com/en-us/node/468396
	cblas_dcopy(size, b, 1, r, 1);
	// r = r - Ax
	// https://software.intel.com/en-us/node/468394
	cblas_daxpy(size, -1.0, Ax, 1, r, 1);

	// p_0 = r_0
	cblas_dcopy(size, r, 1, p, 1);
	
	bnorm2 = cblas_dnrm2(size, b, 1);
	if (bnorm2 == 0.0)
		bnorm2 = 1.0;
	rnorm2 = cblas_dnrm2(size, r, 1);

	if ((resid = rnorm2 / bnorm2) <= err_tol) {
		isConverged = true;
	}

	std::vector<double> MInvVals = InvertMatrixDiagonal(MVals);
	
	while (iter < maxIter && isConverged == false) {
		// z = M^-1 r
		mkl_cspblas_dcsrgemv(&transa, &Anrows, MInvVals.data(), MRowIdx.data(), MCols.data(), r, z);
		rho = cblas_ddot(size, r, 1, z, 1);

		if (iter == 1) {
			cblas_dcopy(size, z, 1, p, 1);
		}
		else {
			// beta = (rho_i / rho_{i - 1}) (\alpha / omega_{i - 1})
			beta = rho / rho1;
			// p_{i} = r_{i-1} + \beta * (p_{i-1} - \omega_{i - 1}v_{i -  1})
			// p = p - \omega_{i - 1}v_{i -  1} = p - \omega_{i - 1}v_{i -  1}
			cblas_dscal(size, beta, p, 1);
			// p = z + beta * p
			cblas_daxpy(size, 1.0, z, 1, p, 1);
			if (std::isnan(beta) || std::isinf(beta)) {
				std::cout << "poisson equation nan/inf error(beta) : " << iter << " " << beta << std::endl;
				exit(1);
			}
			assert(beta == beta);
		}

		// q = A * p_k, will be resued below
		// https://software.intel.com/en-us/node/468560
		mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), p, q);
		// get alpha (r^T_k * r_k) / (d^T_k q) (= (r^T_k * r_k) / (d^T_k A d_k))
		// https://software.intel.com/en-us/node/468398#D4E53C70-D8FA-4095-A800-4203CAFE64FE
		alpha = rho / (cblas_ddot(size, p, 1, q, 1) + err_tol * err_tol);
		
		// x_k+1 = x_k + alpha * p_k
		// https://software.intel.com/en-us/node/468394
		cblas_daxpy(size, alpha, p, 1, x, 1);
		if (std::isnan(alpha) || std::isinf(alpha)) {
			std::cout << "poisson equation nan/inf error(alpha) : " << iter << " " << alpha << std::endl;
			exit(1);
		}
		assert(alpha == alpha);
		
		// Update r
		// r_k+1 = -alpha * A * p_k + r_k = -alpha * q + r_k
		cblas_daxpy(size, -alpha, q, 1, r, 1);
		
		rnorm2 = cblas_dnrm2(size, r, 1);
		if ((resid = rnorm2 / bnorm2) <= err_tol) {
			isConverged = true;
		}
		
		rho1 = rho;
		
		iter++;
	}

	// std::cout << "CG : " << iter << " " << maxIter << " Err : " << rnorm2 / bnorm2 << std::endl;
	Data::Deallocate1Dd(Ax);
	Data::Deallocate1Dd(p);
	Data::Deallocate1Dd(q);
	Data::Deallocate1Dd(r);
	Data::Deallocate1Dd(z);

	return 0;
}

int PoissonSolver3D::BiCGStab_2FUniformU_3D(std::vector<double>& uhat, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& UBC, const int maxIter) {

	MKL_INT Anrows = (kNx - 1) * kNy * kNz, Ancols = (kNx - 1) * kNy * kNz;
	MKL_INT size = (kNx - 1) * kNy * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "uhat(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx - 1; i++) {
		b[i + j * (kNx - 1) + k * (kNx - 1) * kNy] = rhs[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)];
		x[i + j * (kNx - 1) + k * (kNx - 1) * kNy] = 0.0;
	}

	BiCGStab_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx - 1; i++) {
		uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] = x[i + j * (kNx - 1) + k * (kNx - 1) * kNy];
		assert(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] == uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)]);
		if (std::isnan(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)]) || std::isinf(uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)])) {
			std::cout << "uhat(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid + 1 << " " << j + kNumBCGrid << " " << k + kNumBCGrid << " "
				<< " " << uhat[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	UBC->ApplyBC_U_3D(uhat);

	return 0;
}

int PoissonSolver3D::BiCGStab_2FUniformV_3D(std::vector<double>& vhat, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& VBC, const int maxIter) {

	MKL_INT Anrows = kNx * (kNy - 1) * kNz, Ancols = kNx * (kNy - 1) * kNz;
	MKL_INT size = kNx * (kNy - 1) * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "vhat(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy - 1; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * (kNy - 1)] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)];
		x[i + j * kNx + k * kNx * (kNy - 1)] = 0.0;
	}

	BiCGStab_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy - 1; j++)
	for (int i = 0; i < kNx; i++) {
		vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)] = x[i + j * kNx + k * kNx * (kNy - 1)];
		assert(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)] == vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)]);
		if (std::isnan(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)]) || std::isinf(vhat[idx(i + kNumBCGrid, j + kNumBCGrid + 1, k + kNumBCGrid)])) {
			std::cout << "vhat(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid + 1 << " " << k + kNumBCGrid << " "
				<< " " << vhat[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	VBC->ApplyBC_V_3D(vhat);

	return 0;
}

int PoissonSolver3D::BiCGStab_2FUniformW_3D(std::vector<double>& what, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& WBC, const int maxIter) {

	MKL_INT Anrows = kNx * kNy * (kNz - 1), Ancols = kNx * kNy * (kNz - 1);
	MKL_INT size = kNx * kNy * (kNz - 1);

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "what(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz - 1; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * kNy] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)];
		x[i + j * kNx + k * kNx * kNy] = 0.0;
	}

	BiCGStab_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz - 1; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] = x[i + j * kNx + k * kNx * kNy];
		assert(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] == what[idx(i + kNumBCGrid + 1, j + kNumBCGrid, k + kNumBCGrid + 1)]);
		if (std::isnan(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)]) || std::isinf(what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)])) {
			std::cout << "what(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid << " " << k + kNumBCGrid + 1 << " "
				<< " " << what[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid + 1)] << std::endl;
			exit(1);
		}
	}

	WBC->ApplyBC_W_3D(what);

	return 0;
}


int PoissonSolver3D::BiCGStab_2FUniformP_3D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition3D>& PBC, const int maxIter) {

	MKL_INT Anrows = kNx * kNy * kNz, Ancols = kNx * kNy * kNz;
	MKL_INT size = kNx * kNy * kNz;

	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "Pressure(Poisson equation) Error : the # of rows is invalid!" << std::endl;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx + k * kNx * kNy] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)];
		x[i + j * kNx + k * kNx * kNy] = 0.0;
	}

	BiCGStab_2FUniform_3D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, Anrows, size, maxIter);

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] = x[i + j * kNx + k * kNx * kNy];
		assert(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] == ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)]);
		if (std::isnan(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)]) || std::isinf(ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)])) {
			std::cout << "Pressure(Poisson equation) Error : nan/inf occurred : "
				<< i + kNumBCGrid << " " << j + kNumBCGrid << " " << k + kNumBCGrid << " "
				<< " " << ps[idx(i + kNumBCGrid, j + kNumBCGrid, k + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	PBC->ApplyBC_U_3D(ps);

	return 0;
}

int PoissonSolver3D::BiCGStab_2FUniform_3D(double *x, double *b,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx, int64_t Anrows, int64_t size, const int maxIter) { 

	// http://math.nist.gov/iml++/bicgstab.h.txt
	
	double *Ax = Data::Allocate1Dd(size);
	double *r = Data::Allocate1Dd(size);
	double *rtilde = Data::Allocate1Dd(size);
	double *v = Data::Allocate1Dd(size);
	double *p = Data::Allocate1Dd(size);
	double *phat = Data::Allocate1Dd(size);
	double *s = Data::Allocate1Dd(size);
	double *shat = Data::Allocate1Dd(size);
	double *t = Data::Allocate1Dd(size);
	double bnorm2 = 0.0, rnorm2 = 0.0;

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "the # of rows is invalid!" << std::endl;

	// declare coefficients
	double alpha = 1.0, beta = 1.0, omega = 1.0, resid;
	double rho1 = 1.0, rho2 = 1.0;
	const double err_tol = 1.0e-6;

	for (int k = 0; k < kNz; k++)
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		Ax[i + j * kNx + k * kNx * kNy] = 0.0;
		x[i + j * kNx + k * kNx * kNy] = 0.0;
		r[i + j * kNx + k * kNx * kNy] = 0.0;
		v[i + j * kNx + k * kNx * kNy] = 0.0;
		p[i + j * kNx + k * kNx * kNy] = 0.0;
		phat[i + j * kNx + k * kNx * kNy] = 0.0;
		s[i + j * kNx + k * kNx * kNy] = 0.0;
		shat[i + j * kNx + k * kNx * kNy] = 0.0;
		t[i + j * kNx + k * kNx * kNy] = 0.0;
	}

	int iter = 0;
	bool isConverged = false;

	// get Ax(=A*x), using upper triangular matrix (Sparse BLAS)
	// https://software.intel.com/en-us/node/468560
	char transa = 'n';
	mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), x, Ax);
	// r = b - Ax, initial residual
	// r = b
	// https://software.intel.com/en-us/node/468396
	cblas_dcopy(size, b, 1, r, 1);
	// r = r - Ax
	// https://software.intel.com/en-us/node/468394
	cblas_daxpy(size, -1.0, Ax, 1, r, 1);

	// rtilde = r
	cblas_dcopy(size, r, 1, rtilde, 1);
	
	// norm of b matrix
	bnorm2 = cblas_dnrm2(size, b, 1);

	if (bnorm2 == 0.0)
		bnorm2 = 1.0;

	rnorm2 = cblas_dnrm2(size, r, 1);

	if ((resid = rnorm2 / bnorm2) <= err_tol) {
		// no need to do BiCGStab
		isConverged = true;
	}
	
	std::vector<double> MInvVals = InvertMatrixDiagonal(MVals);

	while (iter < maxIter && isConverged == false) {
		// direction vector
		// rho_i = (\hat{r0}, r_(i-1))
		rho1 = cblas_ddot(size, rtilde, 1, r, 1);
		if (rho1 == 0.0) {
			break;
		}
		
		if (iter == 1) {
			cblas_dcopy(size, r, 1, p, 1);
		}
		else {
			// beta = (rho_i / rho_{i - 1}) (\alpha / omega_{i - 1})
			beta = (rho1 / rho2) * (alpha / omega);
			// p_{i} = r_{i-1} + \beta * (p_{i-1} - \omega_{i - 1}v_{i -  1})
			// p = p - \omega_{i - 1}v_{i -  1} = p - \omega_{i - 1}v_{i -  1}
			cblas_daxpy(size, -omega, v, 1, p, 1);
			// p = (beta - 1.0) * p + p = beta * p
			cblas_dscal(size, beta, p, 1);
			// p = 1.0 * r + p = 1.0 * r + \beta * (p_{i-1} - \omega_{i - 1}v_{i -  1})
			cblas_daxpy(size, 1.0, r, 1, p, 1);
		}
		
		// phat = M^-1 p
		mkl_cspblas_dcsrgemv(&transa, &Anrows, MInvVals.data(), MRowIdx.data(), MCols.data(), p, phat);

		// v  = A * p_i
		// https://software.intel.com/en-us/node/468560
		mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), phat, v);
		
		// v = rho_i / (rilde, v)
		alpha = rho1 / cblas_ddot(size, rtilde, 1, v, 1);

		// s = r_{i - 1}
		cblas_dcopy(size, r, 1, s, 1);
		// s = s - \alpha nu_i = r_{i - 1} - \alpha nu_i
		cblas_daxpy(size, -alpha, v, 1, s, 1);

		if ((resid = cblas_dnrm2(size, s, 1) / bnorm2) < err_tol) {
			cblas_daxpy(size, alpha, phat, 1, x, 1);
			isConverged = true;
			break;
		}

		// shat = M^-1 s
		mkl_cspblas_dcsrgemv(&transa, &Anrows, MInvVals.data(), MRowIdx.data(), MCols.data(), s, shat);
		// t = A shat
		mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), shat, t);

		// omega = (t, s) / (t, t)
		omega = cblas_ddot(size, t, 1, s, 1) / cblas_ddot(size, t, 1, t, 1);

		// x_i = x_{i - 1} + \alpha p_i + \omega s
		// x_i = x_{i - 1} + \alpha p_i
		cblas_daxpy(size, alpha, phat, 1, x, 1);
		// x_i = x_{i} + \omega s = x_{i - 1} + \alpha p_i + \omega s
		cblas_daxpy(size, omega, shat, 1, x, 1);

		// r_i =  s_i
		cblas_dcopy(size, s, 1, r, 1);
		// r = r - \omega t
		cblas_daxpy(size, -omega, t, 1, r, 1);
			
		rnorm2 = cblas_dnrm2(size, r, 1);
		
		if (rnorm2 / bnorm2 <= err_tol)
			isConverged = true;
		
		iter++;
	}
	
	// std::cout << "BiCG : " << iter << " " << maxIter << " Err : " << rnorm2 / bnorm2 << std::endl;
	
	Data::Deallocate1Dd(Ax);
	Data::Deallocate1Dd(r);
	Data::Deallocate1Dd(rtilde);
	Data::Deallocate1Dd(v);
	Data::Deallocate1Dd(p);
	Data::Deallocate1Dd(phat);
	Data::Deallocate1Dd(s);
	Data::Deallocate1Dd(shat);
	Data::Deallocate1Dd(t);

	return 0;
}

std::vector<double> PoissonSolver3D::InvertMatrixDiagonal(const std::vector<double>& M) {
	std::vector<double> MInv(M.size(), 0.0);

	for (int i = 0; i < M.size(); i++) {
		MInv[i] = 1.0 / M[i];
	}

	return MInv;
}

inline int PoissonSolver3D::idx(int i, int j, int k) {
	return (k + (kNz + 2 * kNumBCGrid) * (j + (kNy + 2 * kNumBCGrid) * i));
}