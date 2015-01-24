#include "poisson2d.h"

PoissonSolver2D::PoissonSolver2D(int nx, int ny, int num_bc_grid) :
	kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid) {
}

int PoissonSolver2D::GS_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	double dx, double dy, std::shared_ptr<BoundaryCondition2D> PBC) {
	std::vector<double> oldPs((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> tmpPs((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	const double eps = 1.0e-6, omega = 1.5;
	double err = 1.0, err_upper = 0.0, err_lower = 0.0;

	std::copy(ps.begin(), ps.end(), oldPs.begin());
	std::copy(ps.begin(), ps.end(), tmpPs.begin());

	double d = std::min(dx, dy);

	while (err > eps) {
		PBC->ApplyBC_P_2D(tmpPs);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			tmpPs[idx(i, j)] = 0.25
					* (tmpPs[idx(i - 1, j)] + tmpPs[idx(i + 1, j)] + tmpPs[idx(i, j - 1)] + tmpPs[idx(i, j + 1)]
					- d * d * rhs[idx(i, j)]);
		}

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			tmpPs[idx(i, j)] = oldPs[idx(i, j)] + omega * (tmpPs[idx(i, j)] - oldPs[idx(i, j)]);
		}

		PBC->ApplyBC_P_2D(tmpPs);

		err_upper = 0.0;
		err_lower = 0.0;
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			err_upper += (tmpPs[idx(i, j)] - oldPs[idx(i, j)]) * (tmpPs[idx(i, j)] - oldPs[idx(i, j)]);
			err_lower += oldPs[idx(i, j)] * oldPs[idx(i, j)];
		}

		err = err_upper / (err_lower + 1.0e-40);
		err = std::sqrt(err);
		
		std::copy(tmpPs.begin(), tmpPs.end(), oldPs.begin());
	}

	std::copy(tmpPs.begin(), tmpPs.end(), ps.begin());
	
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
	if (PBC->m_BC_PW == BC::DIRICHLET)
		BCtype[0] = 'D';
	else if (PBC->m_BC_PW == BC::NEUMANN)
		BCtype[0] = 'N';
	else if (PBC->m_BC_PW == BC::PERIODIC)
		BCtype[0] = 'P';

	if (PBC->m_BC_PE == BC::DIRICHLET)
		BCtype[1] = 'D';
	else if (PBC->m_BC_PE == BC::NEUMANN)
		BCtype[1] = 'N';
	else if (PBC->m_BC_PE == BC::PERIODIC)
		BCtype[1] = 'P';

	if (PBC->m_BC_PS == BC::DIRICHLET)
		BCtype[2] = 'D';
	else if (PBC->m_BC_PS == BC::NEUMANN)
		BCtype[2] = 'N';
	else if (PBC->m_BC_PS == BC::PERIODIC)
		BCtype[2] = 'P';

	if (PBC->m_BC_PN == BC::DIRICHLET)
		BCtype[3] = 'D';
	else if (PBC->m_BC_PN == BC::NEUMANN)
		BCtype[3] = 'N';
	else if (PBC->m_BC_PN == BC::PERIODIC)
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

		if (PBC->m_BC_PW == BC::DIRICHLET)
			bd_ax[j] = 2.0 * PBC->m_BC_DirichletConstantPW
			- ps[idx(kNumBCGrid, j + kNumBCGrid)];
		if (PBC->m_BC_PE == BC::DIRICHLET)
			bd_bx[j] = 2.0 * PBC->m_BC_DirichletConstantPE
			- ps[idx(kNx + kNumBCGrid, j + kNumBCGrid)];
	}

	/* Setting the values of the boundary function g(x,y) that is equal to
	the normal derivative of the TRUE solution in the mesh points laying on
	Neumann boundaries */
	for (MKL_INT i = 0; i <= MKLnx; i++) {
		bd_ay[i] = 0.0;
		bd_by[i] = 0.0;

		if (PBC->m_BC_PS == BC::DIRICHLET)
			bd_ay[i] = 2.0 * PBC->m_BC_DirichletConstantPS
			- ps[idx(i + kNumBCGrid, kNumBCGrid)];
		if (PBC->m_BC_PN == BC::DIRICHLET)
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

int PoissonSolver2D::ICPCG_2FUniform_2D(std::vector<double>& ps, const std::vector<double>& rhs,
	std::vector<double>& AVals, std::vector<MKL_INT>& ACols, std::vector<MKL_INT>& ARowIdx,
	double lenX, double lenY, double dx, double dy, std::shared_ptr<BoundaryCondition2D> PBC) {
	
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	std::vector<double> b(kNx * kNy), x(kNx * kNy), Ax(kNx * kNy), r(kNx * kNy), d(kNx * kNy), q(kNx * kNy);

	// check values
	if (Anrows != ARowIdx.size() - 1)
		std::cout << "the # of rows is invalid!" << std::endl;

	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		b[i + j * kNx] = rhs[idx(i, j)];
		x[i + j * kNx] = 0.0;
	}
	// get Ax(=A*x), using upper triangular matrix (Sparse BLAS)
	// https://software.intel.com/en-us/node/468560
	char transa = 'n';
	mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), x.data(), Ax.data());
	// r = b - Ax, initial residual
	// r = b
	// https://software.intel.com/en-us/node/468396
	cblas_dcopy(size, b.data(), 1, r.data(), 1);
	// r = r - Ax
	// https://software.intel.com/en-us/node/468394
	cblas_daxpy(size, -1.0, Ax.data(), 1, r.data(), 1);

	// d_0 = r_0
	cblas_dcopy(size, r.data(), 1, d.data(), 1);

	// declare coefficients
	double alpha = 0.0, alpha1 = 0.0, beta = 0.0, beta1 = 0.0;
	double delta_new = 0.0, delta_old = 0.0, delta0;

	// delta = r * r^T
	delta_old = cblas_ddot(size, r.data(), 1, r.data(), 1);
	delta0 = delta_old;
	delta_new = delta_old;
	const int maxiter = 200;
	const double err_tol = 1.0e-6;
	int iter = 0;
	bool isConverged = false;
	while (iter < maxiter && isConverged == false) {
		// q = A * d_k, reuse every iteration
		// https://software.intel.com/en-us/node/468560
		mkl_cspblas_dcsrgemv(&transa, &Anrows, AVals.data(), ARowIdx.data(), ACols.data(), d.data(), q.data());
		// get alpha (r^T_k * r_k) / (d^T_k q) (= (r^T_k * r_k) / (d^T_k A d_k))
		// https://software.intel.com/en-us/node/468398#D4E53C70-D8FA-4095-A800-4203CAFE64FE
		alpha = delta_new / (cblas_ddot(size, q.data(), 1, d.data(), 1) + err_tol * err_tol);

		// x_k+1 = x_k + alpha * d_k
		// https://software.intel.com/en-us/node/468394
		cblas_daxpy(size, alpha, d.data(), 1, x.data(), 1);

		// r_k+1 = r_k -alpha * A * d_k
		cblas_daxpy(size, -alpha, q.data(), 1, r.data(), 1);

		// delta_old = delta_new
		delta_old = delta_new;
		// delta_mew = r^T r
		delta_new = cblas_ddot(size, r.data(), 1, r.data(), 1);

		beta = delta_new / (delta_old + err_tol * err_tol);

		// d_k+1 = r_k+1 + beta * d_k
		// d = d * beta
		cblas_daxpy(size, beta, d.data(), 1, d.data(), 1);
		// d = d + r_k+1
		cblas_daxpy(size, 1.0, r.data(), 1, d.data(), 1);

		if (delta_new <= err_tol * err_tol * delta0)
			isConverged = true;
		iter++;
	}

	for (int j = 0; j < kNy; j++)
	for (int i = 0; i < kNx; i++) {
		ps[idx(i, j)] = d[i + j * kNx];
	}

	return 0;
}

inline int PoissonSolver2D::idx(int i, int j) {
	return i + (kNx + 2 * kNumBCGrid) * j;
}