#include "poisson2d.h"

PoissonSolver2D::PoissonSolver2D(int nx, int ny, int num_bc_grid) :
	kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid) {
}

int PoissonSolver2D::GS_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	const double lenX, const double lenY, const double dx, const double dy,
	const std::shared_ptr<BoundaryCondition2D>& PBC, const int maxIter) {
	
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;

	std::vector<double> oldPs(size, 0.0);
	std::vector<double> tmpPs(size, 0.0);
	std::vector<double> b(size, 0.0);

	const double eps = 1.0e-6, omega = 1.2;
	double err = 1.0, err_upper = 0.0, err_lower = 0.0;

	for (int i = 0; i < kNx; i++)
	for (int j = 0; j < kNy; j++) {
		b[i + j * kNx] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid)];
		oldPs[i + j * kNx] = ps[idx(i + kNumBCGrid, j + kNumBCGrid)];
		tmpPs[i + j * kNx] = ps[idx(i + kNumBCGrid, j + kNumBCGrid)];
	}

	long int j = 0;
	double diag = 0.0;
	int iter = 0;
	while (err > eps && iter < maxIter) {
		// i varies from 0 to kNx * kNy
		for (long int i = 0; i < size; i++) {
			tmpPs[i] = b[i];

			// p is a index of AVals(size x size) and ACols(size x size) size = kNx * kNy
			// ARowIdx[i] is the where ACols and AVals start, (i, ACols[p]) element location of A matrix
			// tmpPs[j] get values from ith row
			for (long int p = ARowIdx[i]; p < ARowIdx[i + 1]; p++) {
				j = ACols[p];

				if (i == j)
					diag = AVals[p];
				else
					tmpPs[i] -= AVals[p] * tmpPs[j];
			}
			
			if (diag == 0.0) {
				perror("zero diagonal term in Gauss Seidel Method");
			}

			tmpPs[i] /= diag;
		}

		for (long int i = 0; i < size; i++) {
			tmpPs[i] = oldPs[i] + omega * (tmpPs[i] - oldPs[i]);
		}

		err_upper = 0.0;
		err_lower = 0.0;
		for (long int i = 0; i < size; i++) {
			err_upper += (tmpPs[i] - oldPs[i]) * (tmpPs[i] - oldPs[i]);
			err_lower += oldPs[i] * oldPs[i];
			oldPs[i] = tmpPs[i];
		}

		err = err_upper / (err_lower + 1.0e-100);
		err = std::sqrt(err);
		iter++;
	}

	double rnorm2 = cblas_dnrm2(size, tmpPs.data(), 1);
	// std::cout << iter << " GS Norm : " << rnorm2 << std::endl;

	for (int i = 0; i < kNx; i++)
	for (int j = 0; j < kNy; j++) {
		ps[idx(i + kNumBCGrid, j + kNumBCGrid)] = tmpPs[i + j * kNx];
	}
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)

		if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
			std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	
	return 0;
}

int PoissonSolver2D::MKL_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	double lenX, double lenY, double dx, double dy, std::shared_ptr<BoundaryCondition2D> PBC) {

	double q;
	double ax, bx, ay, by, az, bz;
	MKL_INT ix, iy, iz, i, stat;
	// nx, ny, nz = # of cells
	// MKLnx, MKLny, MKLnz = # of mesh intervals incluidng one extra boundary node of each side,
	// e.g. 0(boundary),1,2, ... , MKLnx(boundary) : index of node
	MKL_INT MKLnx = kNx + 1, MKLny = kNy + 1;
	MKL_INT ipar[128];
	DFTI_DESCRIPTOR_HANDLE xhandle = 0;
	char *BCtype = new char[5];
	if (PBC->m_BC_PW == BC2D::DIRICHLET)
		BCtype[0] = 'D';
	else if (PBC->m_BC_PW == BC2D::NEUMANN)
		BCtype[0] = 'N';
	else if (PBC->m_BC_PW == BC2D::PERIODIC)
		BCtype[0] = 'P';

	if (PBC->m_BC_PE == BC2D::DIRICHLET)
		BCtype[1] = 'D';
	else if (PBC->m_BC_PE == BC2D::NEUMANN)
		BCtype[1] = 'N';
	else if (PBC->m_BC_PE == BC2D::PERIODIC)
		BCtype[1] = 'P';

	if (PBC->m_BC_PS == BC2D::DIRICHLET)
		BCtype[2] = 'D';
	else if (PBC->m_BC_PS == BC2D::NEUMANN)
		BCtype[2] = 'N';
	else if (PBC->m_BC_PS == BC2D::PERIODIC)
		BCtype[2] = 'P';

	if (PBC->m_BC_PN == BC2D::DIRICHLET)
		BCtype[3] = 'D';
	else if (PBC->m_BC_PN == BC2D::NEUMANN)
		BCtype[3] = 'N';
	else if (PBC->m_BC_PN == BC2D::PERIODIC)
		BCtype[3] = 'P';

	// insert null character 
	BCtype[4] = '\0';

	std::vector<double> bd_ax(MKLny + 1);
	std::vector<double> bd_bx(MKLny + 1);
	std::vector<double> bd_ay(MKLnx + 1);
	std::vector<double> bd_by(MKLnx + 1);
	std::vector<double> f((MKLnx + 1) * (MKLny + 1));
	std::vector<double> dpar(5 * MKLnx / 2 + 7);

	for (i = 0; i < 128; i++) {
		ipar[i] = 0;
	}

	/* Defining the rectangular domain 0<x<1, 0<y<1 for 3D Poisson Solver */
	ax = -0.5 * dx;
	bx = lenX + 0.5 * dx;
	ay = -0.5 * dy;
	by = lenY + 0.5 * dy;
	q = 0.0;

	/* Setting the values of the boundary function g(x,y) that is equal to
	the normal derivative of the TRUE solution in the mesh points laying on
	Neumann boundaries */
	for (MKL_INT j = 0; j <= MKLny; j++) {
		bd_ax[j] = 0.0;
		bd_bx[j] = 0.0;

		if (PBC->m_BC_PW == BC2D::DIRICHLET)
			bd_ax[j] = 2.0 * PBC->m_BC_DirichletConstantPW
			- ps[idx(kNumBCGrid, j + kNumBCGrid)];
		if (PBC->m_BC_PE == BC2D::DIRICHLET)
			bd_bx[j] = 2.0 * PBC->m_BC_DirichletConstantPE
			- ps[idx(kNx + kNumBCGrid, j + kNumBCGrid)];
	}

	/* Setting the values of the boundary function g(x,y) that is equal to
	the normal derivative of the TRUE solution in the mesh points laying on
	Neumann boundaries */
	for (MKL_INT i = 0; i <= MKLnx; i++) {
		bd_ay[i] = 0.0;
		bd_by[i] = 0.0;

		if (PBC->m_BC_PS == BC2D::DIRICHLET)
			bd_ay[i] = 2.0 * PBC->m_BC_DirichletConstantPS
			- ps[idx(i + kNumBCGrid, kNumBCGrid)];
		if (PBC->m_BC_PN == BC2D::DIRICHLET)
			bd_by[i] = 2.0 * PBC->m_BC_DirichletConstantPN
			- ps[idx(i + kNumBCGrid, kNy + kNumBCGrid)];
	}
	
	for (MKL_INT j = 0; j <= MKLny; j++)
	for (MKL_INT i = 0; i <= MKLnx; i++) {	
		f[i + (MKLnx + 1) * j]
			= rhs[idx(i + kNumBCGrid - 1, j + kNumBCGrid - 1)];
	}

	d_init_Helmholtz_2D(&ax, &bx, &ay, &by, &MKLnx, &MKLny, BCtype, &q, ipar, dpar.data(), &stat);
	if (stat != 0){
		perror("MKL Error:d_init_Helmholtz_3D");
		exit(1);
	}

	/* Contains error messaging options:, print to screen (default) */
	ipar[1] = 1;
	/* Contains warning messaging options, ignore warning message */
	ipar[2] = 0;
	ipar[21] = 1;

	d_commit_Helmholtz_2D(f.data(), bd_ax.data(), bd_bx.data(), bd_ay.data(), bd_by.data(), &xhandle, ipar, dpar.data(), &stat);
	if (stat != 0) {
		perror("MKL Error:d_commit_Helmholtz_3D");
		exit(1);
	}

	d_Helmholtz_2D(f.data(), bd_ax.data(), bd_bx.data(), bd_ay.data(), bd_by.data(), &xhandle, ipar, dpar.data(), &stat);
	if (stat != 0 && stat != 1) {
		perror("MKL Error:d_Helmholtz_3D");
		exit(1);
	}

	free_Helmholtz_2D(&xhandle, ipar, &stat);
	if (stat != 0) {
		perror("MKL Error:free_Helmholtz_3D");
		exit(1);
	}

	for (MKL_INT j = 0; j <= MKLny; j++)
	for (MKL_INT i = 0; i <= MKLnx; i++) {
		ps[idx(i + kNumBCGrid - 1, j + kNumBCGrid - 1)]
			= f[i + (MKLnx + 1) * j];
	}

	delete[] BCtype;
	MKL_Free_Buffers();

	return 0;
}

int PoissonSolver2D::CG_2FUniformP_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxIter) {
	
	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);
	int64_t idxArr = 0;

	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		idxArr = i + j * kNx;
		if ((PBC->m_BC_PW == BC2D::NEUMANN || PBC->m_BC_PW == BC2D::AXISYM || PBC->m_BC_PW == BC2D::WALL) &&
			(PBC->m_BC_PE == BC2D::NEUMANN || PBC->m_BC_PE == BC2D::AXISYM || PBC->m_BC_PE == BC2D::WALL) &&
			(PBC->m_BC_PS == BC2D::NEUMANN || PBC->m_BC_PS == BC2D::AXISYM || PBC->m_BC_PS == BC2D::WALL) &&
			(PBC->m_BC_PN == BC2D::NEUMANN || PBC->m_BC_PN == BC2D::AXISYM || PBC->m_BC_PN == BC2D::WALL)) {

			if (i == 0 && j == 0)
				continue;
			idxArr--;
		}

		b[idxArr] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid)];
	}
	
	CG_2FUniform_2D(x, b, AVals, ACols, ARowIdx, MVals, MCols, MRowIdx, PBC, size, maxIter);

	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		idxArr = i + j * kNx;
		if ((PBC->m_BC_PW == BC2D::NEUMANN || PBC->m_BC_PW == BC2D::AXISYM || PBC->m_BC_PW == BC2D::WALL) &&
			(PBC->m_BC_PE == BC2D::NEUMANN || PBC->m_BC_PE == BC2D::AXISYM || PBC->m_BC_PE == BC2D::WALL) &&
			(PBC->m_BC_PS == BC2D::NEUMANN || PBC->m_BC_PS == BC2D::AXISYM || PBC->m_BC_PS == BC2D::WALL) &&
			(PBC->m_BC_PN == BC2D::NEUMANN || PBC->m_BC_PN == BC2D::AXISYM || PBC->m_BC_PN == BC2D::WALL)) {
			if (i == 0 && j == 0)
				continue;
			idxArr--;
		}
		
		ps[idx(i + kNumBCGrid, j + kNumBCGrid)] = x[idxArr];
		assert(ps[idx(i + kNumBCGrid, j + kNumBCGrid)] == ps[idx(i + kNumBCGrid, j + kNumBCGrid)]);
		if (std::isnan(ps[idx(i + kNumBCGrid, j + kNumBCGrid)]) || std::isinf(ps[idx(i + kNumBCGrid, j + kNumBCGrid)])) {
			std::cout << "poisson equation nan/inf error : " << i + kNumBCGrid << " " << j + kNumBCGrid
				 << " " << ps[idx(i + kNumBCGrid, j + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}
	
	Data::Deallocate1Dd(b);
	Data::Deallocate1Dd(x);

	return 0;
}


int PoissonSolver2D::CG_2FUniform_2D(double *x, double *b,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxIter) {

	MKL_INT Anrows = size, Ancols = size;
	int64_t idxArr = 0;

	double *Ax = Data::Allocate1Dd(size);
	double *p = Data::Allocate1Dd(size);
	double *q = Data::Allocate1Dd(size);
	double *r = Data::Allocate1Dd(size);
	double *z = Data::Allocate1Dd(size);

	for (int i = 0; i < size; i++) {
		x[i] = 0.0;
		p[i] = 0.0;
		r[i] = 0.0;
		q[i] = 0.0;
		z[i] = 0.0;
	}

	double alpha = 1.0, beta = 1.0;
	double rho = 1.0, rho1 = 1.0;
	double rnorm2 = 0.0, bnorm2 = 0.0, resid = 0.0;

	const double err_tol = 1.0e-6;
	int iter = 0;
	bool isConverged = false;

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "the # of rows is invalid!" << std::endl;

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

	std::cout << "CG : " << iter << " " << maxIter << " Err : " << rnorm2 / bnorm2 << std::endl;
	return 0;

}

int PoissonSolver2D::BiCGStab_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	std::vector<double>& MVals, std::vector<MKL_INT>& MCols, std::vector<MKL_INT>& MRowIdx,
	const std::shared_ptr<BoundaryCondition2D>& PBC, int64_t size, const int maxIter) {
	
	// http://math.nist.gov/iml++/bicgstab.h.txt
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	double *b = Data::Allocate1Dd(size);
	double *x = Data::Allocate1Dd(size);
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

	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx] = rhs[idx(i + kNumBCGrid, j + kNumBCGrid)];
		Ax[i + j * kNx] = 0.0;
		x[i + j * kNx] = 0.0;
		r[i + j * kNx] = 0.0;
		v[i + j * kNx] = 0.0;
		p[i + j * kNx] = 0.0;
		phat[i + j * kNx] = 0.0;
		s[i + j * kNx] = 0.0;
		shat[i + j * kNx] = 0.0;
		t[i + j * kNx] = 0.0;
	}
	
	// declare coefficients
	double alpha = 1.0, beta = 1.0, omega = 1.0, resid;
	double rho1 = 1.0, rho2 = 1.0;
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
	
	std::cout << "BiCG : " << iter << " " << maxIter << " Err : " << rnorm2 / bnorm2 << std::endl;
	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		ps[idx(i + kNumBCGrid, j + kNumBCGrid)] = x[i + j * kNx];
		assert(ps[idx(i + kNumBCGrid, j + kNumBCGrid)] == ps[idx(i + kNumBCGrid, j + kNumBCGrid)]);
		if (std::isnan(ps[idx(i + kNumBCGrid, j + kNumBCGrid)]) || std::isinf(ps[idx(i + kNumBCGrid, j + kNumBCGrid)])) {
			std::cout << "poisson equation nan/inf error : " << i + kNumBCGrid << " " << j + kNumBCGrid 
				<< " " << ps[idx(i + kNumBCGrid, j + kNumBCGrid)] << std::endl;
			exit(1);
		}
	}

	Data::Deallocate1Dd(b);
	Data::Deallocate1Dd(x);
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

std::vector<double> PoissonSolver2D::InvertMatrixDiagonal(const std::vector<double>& M) {
	std::vector<double> MInv(M.size(), 0.0);

	for (int i = 0; i < M.size(); i++) {
		MInv[i] = 1.0 / M[i];
	}

	return MInv;
}

inline int PoissonSolver2D::idx(int i, int j) {
	return (j + (kNy + 2 * kNumBCGrid) * (i));
}