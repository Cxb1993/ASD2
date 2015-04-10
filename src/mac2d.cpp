#include "mac2d.h"

MACSolver2D::MACSolver2D(double Re, double We, double Fr,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoI, double muI,
	int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	TimeOrderEnum timeOrder, double cfl, int maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoI(rhoI), kRhoO(rhoI * densityRatio), kRhoRatio(densityRatio),
	kMuI(muI), kMuO(muI * viscosityRatio), kMuRatio(viscosityRatio),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK) {

	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::MACSolver2D(double rhoI, double rhoO, double muI, double muO, double gConstant,
	double L, double U, double sigma, int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	TimeOrderEnum timeOrder, double cfl, int maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoI), kMuScale(muI), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoI(rhoI), kRhoO(rhoO), kMuI(muI), kMuO(muO), kRhoRatio(rhoI / rhoO), kMuRatio(muI / muO),
	kRe(rhoI * L * U / muI), kWe(rhoI * L * U * U / sigma), kFr(U * U / (gConstant * L)),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK) {

	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::~MACSolver2D() {
	if (m_BC)
		m_BC.reset();
}

// Deallocated automatically
int MACSolver2D::AllocateVariables() {
	m_u = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_v = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	
	m_rho = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_kappa = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_mu = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	
	m_ps = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_p = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);

	m_J11 = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_J12 = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_J21 = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_J22 = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);

	return 0;
}

int MACSolver2D::UpdateKappa(const std::vector<double>& ls) {
	std::vector<double> dLSdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dLSdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LSSize((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		
		LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);

		// if size is zero, use forward or backward differencing rather than central differencing
		if (LSSize[idx(i, j)] < 1.0e-100) {
			dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			if (LSSize[idx(i, j)] < 1.0e-100) {
				dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
				LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			}
		}
		if (LSSize[idx(i, j)] < 1.0e-100)
			perror("Div/0 Err in computing kappa");
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		m_kappa[idx(i, j)]
			= -(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] * (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDy * kDy) // phi^2_x \phi_yy
			- 2.0 * dLSdX[idx(i, j)] * dLSdY[idx(i, j)] * (dLSdX[idx(i, j + 1)] - dLSdX[idx(i, j - 1)]) / (2.0 * kDy) //2 \phi_x \phi_y \phi_xy
			+ dLSdY[idx(i, j)] * dLSdY[idx(i, j)] * (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / (kDx * kDx)) // phi^2_y \phi_xx
				/ std::pow(LSSize[idx(i, j)], 3.0);

		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j)] = std::min(m_kappa[idx(i, j)], 1.0 / std::min(kDx, kDy));

		assert(m_kappa[idx(i, j)] == m_kappa[idx(i, j)]);
	}

	return 0;
}

int MACSolver2D::UpdateJumpCond(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls) {
	const double eps = 1.0e-100;
	// derivatives
	double dudX = 0.0, dudY = 0.0;
	double dvdX = 0.0, dvdY = 0.0;
	double dldX = 0.0, dldY = 0.0;
	// normal and tangent vector variable (n, t1, t2)
	// normal vector has X, Y component
	// first tangent vector has only Z component
	double nX = 0.0, nY = 0.0, t1X = 0.0, t1Y = 0.0, t2X = 0.0, t2Y = 0.0, t2Z = 0.0;
	std::vector<double> nXVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> nYVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t1XVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t1YVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t2ZVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	// jump condition originally defined as J11 = [\mu u_x], J12 = [\mu u_y], J21 = [\mu v_x], J22 = [\mu v_y], 
	// but skip [\mu], from Kang, Fedkiw, and Liu's work eq. (30).
	// [\mu] changes with level set, which means I can determine [\mu] later
	
	double nxnx = 0.0, nxny = 0.0, nyny = 0.0;
	double t1xt1x = 0.0, t1xt1y = 0.0, t1yt1y = 0.0;
	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++)
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++) {
		// defined at P grid
		dudX = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		dudY = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
			- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);
		dvdX = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
		dvdY = (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		dldX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dldY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		if (std::fabs(dldX) < eps) {
			dldX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			if (std::fabs(dldX) < eps)
				dldX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		}
		if (std::fabs(dldX) < eps) {
			dldY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			if (std::fabs(dldY) < eps)
				dldY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
		}

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nX = dldX / (std::sqrt(std::pow(std::fabs(dldX), 2.0) + std::pow(std::fabs(dldY), 2.0)) + eps);
		nY = dldY / (std::sqrt(std::pow(std::fabs(dldX), 2.0) + std::pow(std::fabs(dldY), 2.0)) + eps);
		// get first tangent vector
		// smallest magnitude determine what to be performed cross product
		// for 2D, smallest magnitude must be z axis (because it is 2D!)
		t1X = nY / (std::sqrt(std::pow(std::fabs(nX), 2.0) + std::pow(std::fabs(nY), 2.0)) + eps);
		t1Y = -nX / (std::sqrt(std::pow(std::fabs(nX), 2.0) + std::pow(std::fabs(nY), 2.0)) + eps);
		t2Z = (nX * t1Y - nY * t1X);
		nxnx = nX * nX;
		nxny = nX * nY;
		nyny = nY * nY;
		t1xt1x = t1X * t1X;
		t1xt1y = t1X * t1Y;
		t1yt1y = t1Y * t1Y;

		nXVec[idx(i, j)] = nX;
		nYVec[idx(i, j)] = nY;
		t1XVec[idx(i, j)] = t1X;
		t1YVec[idx(i, j)] = t1Y;
		t2ZVec[idx(i, j)] = t2Z;

		// eq. (30) from Kang, Fedkiw, and Liu
		m_J11[idx(i, j)] = dudX * t1xt1x + dudY * t1xt1y 
			+ (dudX * nxnx + dvdX * nxny) * nxnx + (dudY * nxnx + dvdY * nxny) * nxny
			- ((dudX * nxnx + dvdX * nxny) * t1xt1x + (dudY * nxnx + dvdY * nxny) * t1xt1y);
		m_J12[idx(i, j)] = dudX * t1xt1y + dudY * t1yt1y
			+ (dudX * nxnx + dvdX * nxny) * nxny + (dudY * nxnx + dvdY * nxny) * nyny
			- ((dudX * nxny + dvdX * nyny) * t1xt1x + (dudY * nxny + dvdY * nyny) * t1xt1y);
		m_J21[idx(i, j)] = dvdX * t1xt1x + dvdY * t1xt1y
			+ (dudX * nxny + dvdX * nyny) * nxnx + (dudY * nxny + dvdY * nyny) * nxny
			- ((dudX * nxnx + dvdX * nxny) * t1xt1y + (dudY * nxnx + dvdY * nxny) * t1yt1y);
		m_J22[idx(i, j)] = dvdX * t1xt1y + dvdY * t1yt1y
			+ (dudX * nxny + dvdX * nyny) * nxny + (dudY * nxny + dvdY * nyny) * nyny
			- ((dudX * nxny + dvdX * nyny) * t1xt1y + (dudY * nxny + dvdY * nyny) * t1yt1y);
		
		if (std::isnan(m_J11[idx(i, j)]) || std::isinf(m_J11[idx(i, j)])) {
			std::cout << "Jump condition(J11) nan/inf error : " << i << " " << j << " " << m_J11[idx(i, j)] << std::endl;
			exit(1);
		}
		if (std::isnan(m_J12[idx(i, j)]) || std::isinf(m_J12[idx(i, j)])) {
			std::cout << "Jump condition(J12) nan/inf error : " << i << " " << j << " " << m_J12[idx(i, j)] << std::endl;
			exit(1);
		}
		if (std::isnan(m_J21[idx(i, j)]) || std::isinf(m_J21[idx(i, j)])) {
			std::cout << "Jump condition(J21) nan/inf error : " << i << " " << j << " " << m_J21[idx(i, j)] << std::endl;
			exit(1);
		}
		if (std::isnan(m_J22[idx(i, j)]) || std::isinf(m_J22[idx(i, j)])) {
			std::cout << "Jump condition(J22) nan/inf error : " << i << " " << j << " " << m_J22[idx(i, j)] << std::endl;
			exit(1);
		}

		assert(m_J11[idx(i, j)] == m_J11[idx(i, j)]);
		assert(m_J12[idx(i, j)] == m_J12[idx(i, j)]);
		assert(m_J21[idx(i, j)] == m_J21[idx(i, j)]);
		assert(m_J22[idx(i, j)] == m_J22[idx(i, j)]);
	}

	return 0;
}

std::vector<double> MACSolver2D::UpdateFU(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v) {
	
	std::vector<double> cU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> vU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> gU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhsU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v);

	// Viscous term
	vU = this->AddViscosityFU(u, v, ls);

	// Gravity term
	gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
		// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		
		if (lsW >= 0 && lsM >= 0) {
			iRhoEff = 1.0 / kRhoI;
		}
		else if (lsW < 0 && lsM < 0) {
			iRhoEff = 1.0 / kRhoO;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between ls[i - 1,j] and ls[i,j]
			// |(lsW)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			iRhoEff = ((1.0 / kRhoO) * std::fabs(lsW) + (1.0 / kRhoI) * std::fabs(lsM))
				/ (std::fabs(lsW) + std::fabs(lsM));
		}
		else if (lsW < 0 && lsM >= 0) {
			// interface lies between ls[i - 1, j] and ls[i,j]
			// |(lsW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			iRhoEff = ((1.0 / kRhoI) * std::fabs(lsW) + (1.0 / kRhoO) * std::fabs(lsM))
				/ (std::fabs(lsW) + std::fabs(lsM));
		}

		rhsU[idx(i, j)] = -cU[idx(i, j)] + vU[idx(i, j)] * iRhoEff + gU[idx(i, j)];
	}
	
	return rhsU;
}

std::vector<double> MACSolver2D::UpdateFV(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v) {

	std::vector<double> cV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> vV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> gV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhsV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Convection term
	cV = this->AddConvectionFV(u, v);

	// Viscous term
	vV = this->AddViscosityFV(u, v, ls);

	// Gravity term
	gV = this->AddGravityFV();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		
		if (lsS >= 0 && lsM >= 0) {
			iRhoEff = 1.0 / kRhoI;
		}
		else if (lsS < 0 && lsM < 0) {
			iRhoEff = 1.0 / kRhoO;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			iRhoEff = ((1.0 / kRhoO) * std::fabs(lsS) + (1.0 / kRhoI) * std::fabs(lsM))
				/ (std::fabs(lsS) + std::fabs(lsM));
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			iRhoEff = ((1.0 / kRhoI) * std::fabs(lsS) + (1.0 / kRhoO) * std::fabs(lsM))
				/ (std::fabs(lsS) + std::fabs(lsM));
		}
		
		rhsV[idx(i, j)] = -cV[idx(i, j)] + vV[idx(i, j)] * iRhoEff + gV[idx(i, j)];
	}
	
	return rhsV;
}

std::vector<double> MACSolver2D::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> cU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> tmpV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		tmpV[idx(i, j)] = (v[idx(i, j)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)] + v[idx(i, j + 1)]) * 0.25;

	std::vector<double> vecF_UX(kNx + 2 * kNumBCGrid, 0.0), vecF_UY(kNy + 2 * kNumBCGrid);

	// U : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_UX[i] = u[idx(i, j)] * u[idx(i, j)];
		}

		UnitHJWENO5(vecF_UX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j)] = FXP[i];
			LXM[idx(i, j)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_UX.begin(), vecF_UX.end(), 0.0);
	}
	
	// U : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		for (int j = 0; j < kNy + 2 *kNumBCGrid; j++) {
			vecF_UY[j] = u[idx(i, j)] * v[idx(i, j)];
		}

		UnitHJWENO5(vecF_UY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j)] = FYP[j];
			LYM[idx(i, j)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_UY.begin(), vecF_UY.end(), 0.0);
	}
	
	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, rhoEff = 0.0;

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		alphaX = u[idx(i, j)];
		alphaY = tmpV[idx(i, j)];
		
		cU[idx(i, j)]
				= (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
				+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
				- alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
				- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])));
		
		if (isnan(cU[idx(i, j)]) || isinf(cU[idx(i, j)])) {
			std::cout << (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
				+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])) << " " << -alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
				- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])) << std::endl;
		}
		assert(cU[idx(i, j)] == cU[idx(i, j)]);
		if (std::isnan(cU[idx(i, j)]) || std::isinf(cU[idx(i, j)])) {
			std::cout << "U-convection term nan/inf error : " << i << " " << j << " " << cU[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return cU;
}

std::vector<double> MACSolver2D::AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> cV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> tmpU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		tmpU[idx(i, j)] = (u[idx(i, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)] + u[idx(i + 1, j)]) * 0.25;

	std::vector<double> vecF_VX(kNx + 2 * kNumBCGrid, 0.0), vecF_VY(kNy + 2 * kNumBCGrid);
	
	// V : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_VX[i] = v[idx(i, j)] * u[idx(i, j)];
		}

		UnitHJWENO5(vecF_VX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j)] = FXP[i];
			LXM[idx(i, j)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_VX.begin(), vecF_VX.end(), 0.0);
	}

	// V : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_VY[j] = v[idx(i, j)] * v[idx(i, j)];
		}

		UnitHJWENO5(vecF_VY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j)] = FYP[j];
			LYM[idx(i, j)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_VY.begin(), vecF_VY.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		alphaX = tmpU[idx(i, j)];
		alphaY = v[idx(i, j)];
		
		cV[idx(i, j)]
			= (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
			+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
			- alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
			- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])));

		assert(cV[idx(i, j)] == cV[idx(i, j)]);
		if (std::isnan(cV[idx(i, j)]) || std::isinf(cV[idx(i, j)])) {
			std::cout << "V-convection term nan/inf error : " << i << " " << j << " " << cV[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return cV;
}

int MACSolver2D::UnitHJWENO5(
	const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n) {
	int stat = 0;
	const double kWeightEps = 1e-6;

	double alphaPlus[3], alphaMinus[3];
	double weightPlus[3], weightMinus[3];

	std::vector<double> IS0(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS1(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS2(n + 2 * kNumBCGrid, 0.0);
	
	// \dfrac{\Delta^+ F}{\Delta x}
	std::vector<double> DFPlus = std::vector<double>(n + 2  * kNumBCGrid, 0.0);

	// Compute Smoothness for phi^{-}
	for (int i = 0; i < n + 2 * kNumBCGrid - 1; i++)
		DFPlus[i] = F[i + 1] - F[i];

	for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
		IS0[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i - 3] - 2.0 * DFPlus[i - 2] + DFPlus[i - 1]) / d, 2.0)
			+ 1.0 / 4.0
			* std::pow((DFPlus[i - 3] - 4.0 * DFPlus[i - 2] + 3.0 * DFPlus[i - 1]) / d, 2.0);
		IS1[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i - 2] - 2.0 * DFPlus[i - 1] + DFPlus[i]) / d, 2.0)
			+ 1.0 / 4.0 * std::pow((DFPlus[i - 2] - DFPlus[i]) / d, 2.0);
		IS2[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i - 1] - 2.0 * DFPlus[i] + DFPlus[i + 1]) / d, 2.0)
			+ 1.0 / 4.0
			* std::pow((3.0 * DFPlus[i - 1] - 4.0 * DFPlus[i] + DFPlus[i + 1]) / d, 2.0);
	}

	// phi^{-}
	for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
		alphaMinus[0] = 0.1 / std::pow((kWeightEps + IS0[i]), 2.0);
		alphaMinus[1] = 0.6 / std::pow((kWeightEps + IS1[i]), 2.0);
		alphaMinus[2] = 0.3 / std::pow((kWeightEps + IS2[i]), 2.0);

		weightMinus[0] = alphaMinus[0]
			/ (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);
		weightMinus[1] = alphaMinus[1]
			/ (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);
		weightMinus[2] = alphaMinus[2]
			/ (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);

		FM[i] =
			weightMinus[0]
			* (1.0 / 3.0 * DFPlus[i - 3] - 7.0 / 6.0 * DFPlus[i - 2] + 11.0 / 6.0 * DFPlus[i - 1]) / d +
			weightMinus[1]
			* (-1.0 / 6.0 * DFPlus[i - 2] + 5.0 / 6.0 * DFPlus[i - 1] + 1.0 / 3.0 * DFPlus[i]) / d +
			weightMinus[2]
			* (1.0 / 3.0 * DFPlus[i - 1] + 5.0 / 6.0 * DFPlus[i] - 1.0 / 6.0 * DFPlus[i + 1]) / d;

		assert(FM[i] == FM[i]);
	}

	// Compute Smoothness for phi^{+}
	for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
		IS0[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i + 2] - 2.0 * DFPlus[i + 1] + DFPlus[i]) / d, 2.0)
			+ 1.0 / 4.0
			* std::pow((DFPlus[i + 2] - 4.0 * DFPlus[i + 1] + 3.0 * DFPlus[i]) / d, 2.0);
		IS1[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i + 1] - 2.0 * DFPlus[i] + DFPlus[i - 1]) / d, 2.0)
			+ 1.0 / 4.0	* std::pow((DFPlus[i + 1] - DFPlus[i - 1]) / d, 2.0);
		IS2[i]
			= 13.0 / 12.0
			* std::pow((DFPlus[i] - 2.0 * DFPlus[i - 1] + DFPlus[i - 2]) / d, 2.0)
			+ 1.0 / 4.0
			* std::pow((3.0 * DFPlus[i] - 4.0 * DFPlus[i - 1] + DFPlus[i - 2]) / d, 2.0);
	}

	// phi^{+}
	for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
		alphaPlus[0] = 0.1 / std::pow((kWeightEps + IS0[i]), 2.0);
		alphaPlus[1] = 0.6 / std::pow((kWeightEps + IS1[i]), 2.0);
		alphaPlus[2] = 0.3 / std::pow((kWeightEps + IS2[i]), 2.0);

		weightPlus[0] = alphaPlus[0]
			/ (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);
		weightPlus[1] = alphaPlus[1]
			/ (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);
		weightPlus[2] = alphaPlus[2]
			/ (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);

		FP[i] =
			weightPlus[0]
			* (1.0 / 3.0 * DFPlus[i + 2] - 7.0 / 6.0 * DFPlus[i + 1] + 11.0 / 6.0 * DFPlus[i]) / d +
			weightPlus[1]
			* (-1.0 / 6.0 * DFPlus[i + 1] + 5.0 / 6.0 * DFPlus[i] + 1.0 / 3.0 * DFPlus[i - 1]) / d +
			weightPlus[2]
			* (1.0 / 3.0 * DFPlus[i] + 5.0 / 6.0 * DFPlus[i - 1] - 1.0 / 6.0 * DFPlus[i - 2]) / d;

		assert(FP[i] == FP[i]);
	}

	return 0;
}

std::vector<double> MACSolver2D::AddViscosityFU(const std::vector<double>& u, const std::vector<double>& v, 
	const std::vector<double>& ls) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dU[idx(i, j)] = 0.0;
		
		return dU;
	}

	// level set
	double lsUW = 0.0, lsUE = 0.0, lsUS = 0.0, lsUN = 0.0, lsUM = 0.0;
	double lsVW_N = 0.0, lsVW = 0.0, lsVM = 0.0, lsVN = 0.0;
	/*
	-------------------------------------
	|			|			|			|
	|		  lsUN			|			|
	|			|			|			|
	----lsVW_N-------lsVN----------------
	|			|			|			|
  lsUW		  lsUM  (i,j)  lsUE			|
	|			|			|			|
	----lsVW---------lsVM----------------
	|			|			|			|
	|		  lsUS		  lsUE_S		|
	|			|			|			|
	-------------------------------------
	*/
	// subcell
	double theta = 0.0;
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
	double muV_X_S = 0.0, muV_X_N = 0.0;
	double uW = 0.0, uE = 0.0, uS = 0.0, uN = 0.0, uM = 0.0;
	double vW = 0.0, vW_N = 0.0, vM = 0.0, vN = 0.0;
	double visX = 0.0, visY = 0.0;

	// jump condition
	double JUW = 0.0, JUE = 0.0, JUS = 0.0, JUN = 0.0, JUM = 0.0;
	double JVW = 0.0, JVW_N = 0.0, JVM = 0.0, JVN = 0.0;
	// effective Jump condition, effective u(uEff), and effective mu (muEff)
	// J is a jump condition and defined at P grid
	double JEff = 0.0, JO = 0.0, uEff = 0.0, vEff = 0.0, muEff = 0.0;

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		visX = 0.0, visY = 0.0;
		muU_X_W = 0.0; muU_X_E = 0.0; muU_Y_S = 0.0; muU_Y_N = 0.0;
		muV_X_S = 0.0; muV_X_N = 0.0;
		lsUW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
		lsUE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsUS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsUN = 0.5 * (ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsVW_N = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j + 1)]);
		lsVM = 0.5 * (ls[idx(i, j)] + ls[idx(i, j - 1)]);
		lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		uM = u[idx(i, j)];
		uW = u[idx(i - 1, j)];
		uE = u[idx(i + 1, j)];
		uS = u[idx(i, j - 1)];
		uN = u[idx(i, j + 1)];
		vW = v[idx(i - 1, j)];
		vW = v[idx(i - 1, j)];
		vW_N = v[idx(i - 1, j + 1)];
		vM = v[idx(i, j)];
		vN = v[idx(i, j + 1)];
		
		// U part 

		if (lsUW >= 0 && lsUM >= 0) {
			muU_X_W = kMuI * (uM - uW) / kDx;
		}
		else if (lsUW < 0 && lsUM < 0) {
			muU_X_W = kMuO * (uM - uW) / kDx;
		}
		else if (lsUW > 0 && lsUM <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = fabs(lsUW) / (fabs(lsUW) + fabs(lsUM));
			// |(lsUW)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsUM)|
			// |(lsUW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUW = (kMuO - kMuI) * (m_J11[idx(i - 2, j)] + m_J11[idx(i - 1, j)]) * 0.5;
			JUM = (kMuO - kMuI) * (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUW;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_X_W = muEff * (uM - uW) / kDx - muEff * JEff * theta / kMuI;
		}
		else if (lsUW <= 0 && lsUM > 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = fabs(lsUW) / (fabs(lsUW) + fabs(lsUM));
			// |(lsUW)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsUM)|
			// |(lsUW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUW = (kMuI - kMuO) * (m_J11[idx(i - 2, j)] + m_J11[idx(i - 1, j)]) * 0.5;
			JUM = (kMuI - kMuO) * (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUW;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_X_W = muEff * (uM - uW) / kDx + muEff * JEff * theta / kMuO;
		}
		
		if (lsUM >= 0 && lsUE >= 0) {
			muU_X_E = kMuI * (uE - uM) / kDx;
		}
		else if (lsUM < 0 && lsUE < 0) {
			muU_X_E = kMuO * (uE - uM) / kDx;
		}
		else if (lsUM > 0 && lsUE <= 0) {
			// interface lies between u[i, j] and u[i + 1, j]
			theta = fabs(lsUE) / (fabs(lsUE) + fabs(lsUM));
			// |(lsUM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsUE)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsUE)|
			JUM = (kMuO - kMuI) * (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JUE = (kMuO - kMuI) * (m_J11[idx(i, j)] + m_J11[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUE;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_X_E = muEff * (uE - uM) / kDx + muEff * JEff * theta / kMuO;
		}
		else if (lsUM <= 0 && lsUE > 0) {
			// interface lies between u[i, j] and u[i + 1, j]
			theta = fabs(lsUE) / (fabs(lsUE) + fabs(lsUM));
			// |(lsUM)| ===   outside(-)    === |(interface)| === inside(+) === |(lsUE)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsUE)|
			JUM = (kMuI - kMuO) * (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JUE = (kMuI - kMuO) * (m_J11[idx(i, j)] + m_J11[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUE;
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_X_E = muEff * (uE - uM) / kDx - muEff * JEff * theta / kMuI;
		}
		
		if (lsUS >= 0 && lsUM >= 0) {
			muU_Y_S = kMuI * (uM - uS) / kDy;
		}
		else if (lsUS < 0 && lsUM < 0) {
			muU_Y_S = kMuO * (uM - uS) / kDy;
		}
		else if (lsUS > 0 && lsUM <= 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = fabs(lsUS) / (fabs(lsUS) + fabs(lsUM));
			// |(lsUS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsUM)|
			// |(lsUS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (kMuO - kMuI) * (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (kMuO - kMuI) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (uM - uS) / kDy - muEff * JEff * theta / kMuI;
		}
		else if (lsUS <= 0 && lsUM > 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = fabs(lsUS) / (fabs(lsUS) + fabs(lsUM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JUS = (kMuI - kMuO) * (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (kMuI - kMuO) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (uM - uS) / kDy + muEff * JEff * theta / kMuO;
		}
		
		if (lsUM >= 0 && lsUN >= 0) {
			muU_Y_N = kMuI * (uN - uM) / kDy;
		}
		else if (lsUM < 0 && lsUN < 0) {
			muU_Y_N = kMuO * (uN - uM) / kDy;
		}
		else if (lsUM > 0 && lsUN <= 0) {
			// interface lies between u[i, j] and u[i, j + 1]
			theta = fabs(lsUN) / (fabs(lsUN) + fabs(lsUM));
			// |(lsUM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsUN)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsUN)|
			JUM = (kMuO - kMuI) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JUN = (kMuO - kMuI) * (m_J12[idx(i - 1, j + 1)] + m_J12[idx(i, j + 1)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUN;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_N = muEff * (uN - uM) / kDy + muEff * JEff * theta / kMuO;
		}
		else if (lsUM <= 0 && lsUN > 0) {
			// interface lies between u[i, j] and u[i, j + 1]
			theta = fabs(lsUN) / (fabs(lsUN) + fabs(lsUM));
			// |(lsUM)| ===    outside(-)   === |(interface)| ===   inside(+) === |(lsUN)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsUN)|
			JUM = (kMuI - kMuO) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JUN = (kMuI - kMuO) * (m_J12[idx(i - 1, j + 1)] + m_J12[idx(i, j + 1)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUN;
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_Y_N = muEff * (uN - uM) / kDy - muEff * JEff * theta / kMuI;
		}

		// V part

		if (lsVW >= 0 && lsVM >= 0) {
			muV_X_S = kMuI * (vM - vW) / kDx;
		}
		else if (lsVW < 0 && lsVM < 0) {
			muV_X_S = kMuO * (vM - vW) / kDx;
		}
		else if (lsVW > 0 && lsVM <= 0) {
			// interface lies between v[i - 1, j] and v[i, j]
			theta = fabs(lsVW) / (fabs(lsVW) + fabs(lsVM));
			// |(lsVW)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (kMuO - kMuI) * (m_J21[idx(i - 1, j - 1)] + m_J21[idx(i - 1, j)]) * 0.5;
			JVM = (kMuO - kMuI) * (m_J21[idx(i, j - 1)] + m_J21[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_X_S = muEff * (vM - vW) / kDx - muEff * JEff * theta / kMuI;
		}
		else if (lsVW <= 0 && lsVM > 0) {
			// interface lies between v[i - 1, j] and v[i, j]
			theta = fabs(lsUS) / (fabs(lsUS) + fabs(lsUM));
			// |(lsVW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (kMuI - kMuO) * (m_J21[idx(i - 1, j - 1)] + m_J21[idx(i - 1, j)]) * 0.5;
			JVM = (kMuI - kMuO) * (m_J21[idx(i, j - 1)] + m_J21[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_S = muEff * (vM - vW) / kDx + muEff * JEff * theta / kMuO;
		}

		if (lsVW_N >= 0 && lsVN >= 0) {
			muV_X_N = kMuI * (vN - vW_N) / kDx;
		}
		else if (lsVW_N < 0 && lsVN < 0) {
			muV_X_N = kMuO * (vN - vW_N) / kDx;
		}
		else if (lsVW_N > 0 && lsVN <= 0) {
			// interface lies between v[i - 1, j + 1] and v[i, j + 1]
			theta = fabs(lsVW_N) / (fabs(lsVW_N) + fabs(lsVN));
			// |(lsVW_N)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsVN)|
			// |(lsVW_N)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVN)|
			JVW_N = (kMuO - kMuI) * (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j + 1)]) * 0.5;
			JVN = (kMuO - kMuI) * (m_J21[idx(i, j)] + m_J21[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVN + (1.0 - theta) * JVW_N;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_X_N = muEff * (vW_N - vN) / kDx + muEff * JEff * theta / kMuI;
		}
		else if (lsVW_N <= 0 && lsVN > 0) {
			// interface lies between v[i - 1, j + 1] and v[i, j + 1]
			theta = fabs(lsVW_N) / (fabs(lsVW_N) + fabs(lsVN));
			// |(lsVW_N)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVN)|
			// |(lsVW_N)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVN)|
			JVW_N = (kMuI - kMuO) * (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j + 1)]) * 0.5;
			JVN = (kMuI - kMuO) * (m_J21[idx(i, j)] + m_J21[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVN + (1.0 - theta) * JVW_N;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_N = muEff * (vW_N - vN) / kDx - muEff * JEff * theta / kMuO;
		}
		
		visX = 2.0 * (muU_X_E - muU_X_W) / kDx;
		visY = (muU_Y_N - muU_Y_S) / kDy + (muV_X_N - muV_X_S) / kDy;
		
		dU[idx(i, j)] = visX + visY;
		
		if (std::isnan(dU[idx(i, j)]) || std::isinf(dU[idx(i, j)])) {
			std::cout << "U-viscosity term nan/inf error : " << i << " " << j << " " << dU[idx(i, j)] << std::endl;
			exit(1);
		}
		assert(dU[idx(i, j)] == dU[idx(i, j)]);
	}

	return dU;
}

std::vector<double> MACSolver2D::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dV[idx(i, j)] = 0.0;

		return dV;
	}

	// level set
	/*
	----------------lsVN-----------------
	|			|			|			|
	|		  lsUM  (i,j)  lsUE			|
	|			|			|			|
	-----lsVW-------lsVM--------lsVE-----
	|			|			|			|
	|		  lsUS		  lsUE_S		|
	|			|			|			|
	----------------lsVS-----------------
	*/

	double lsVW = 0.0, lsVE = 0.0, lsVS = 0.0, lsVN = 0.0, lsVM = 0.0;
	double lsUM = 0.0, lsUS = 0.0, lsUE = 0.0, lsUE_S = 0.0;
	// jump condition
	double JVW = 0.0, JVE = 0.0, JVS = 0.0, JVN = 0.0, JVM = 0.0;
	double JUM = 0.0, JUS = 0.0, JUE = 0.0, JUE_S = 0.0;
	double theta = 0.0;
	double uS = 0.0, uM = 0.0, uE = 0.0, uE_S = 0.0;
	double vW = 0.0, vE = 0.0, vS = 0.0, vN = 0.0, vM = 0.0;
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0;
	double muU_Y_W = 0.0, muU_Y_E = 0.0;
	double visX = 0.0, visY = 0.0;

	// effective Jump condition, effective v (vEff), and effective mu (muEff)
	double JEff = 0.0, JO = 0.0, vEff = 0.0, uEff = 0.0, muEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		muV_X_W = 0.0; muV_X_E = 0.0; muV_Y_S = 0.0; muV_Y_N = 0.0;
		muU_Y_W = 0.0; muU_Y_E = 0.0;
		visX = 0.0, visY = 0.0;
		lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsVE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
		lsVM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsVS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
		lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsUS = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j - 1)]);
		lsUE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsUE_S = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i + 1, j - 1)]);
		vW = v[idx(i - 1, j)];
		vE = v[idx(i + 1, j)];
		vS = v[idx(i, j - 1)];
		vN = v[idx(i, j + 1)];
		vM = v[idx(i, j)];
		uM = u[idx(i, j)];
		uE = u[idx(i + 1, j)];
		uS = u[idx(i, j - 1)];
		uE_S = u[idx(i + 1, j - 1)];

		// V part

		if (lsVW >= 0 && lsVM >= 0) {
			muV_X_W = kMuI * (vM - vW) / kDx;
		}
		else if (lsVW < 0 && lsVM < 0) {
			muV_X_W = kMuO * (vM - vW) / kDx;
		}
		else if (lsVW >= 0 && lsVM < 0) {
			// interface lies between v[i, j] and v[i - 1, j]
			theta = fabs(lsVW) / (fabs(lsVW) + fabs(lsVM));
			// |(lsVW)| === inside(+) === |(interface)| ===   outside(-)    === |(lsVM)|
			// |(lsVW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (kMuO - kMuI) * (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j - 1)]) * 0.5;
			JVM = (kMuO - kMuI) * (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_X_W = muEff * (vM - vW) / kDx - muEff * JEff * theta / kMuI;
		}
		else if (lsVW < 0 && lsVM >= 0) {
			// interface lies between v[i, j] and v[i - 1, j]
			theta = fabs(lsVW) / (fabs(lsVW) + fabs(lsVM));
			// |(lsVW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (kMuI - kMuO) * (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j - 1)]) * 0.5;
			JVM = (kMuI - kMuO) * (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_W = muEff * (vM - vW) / kDx + muEff * JEff * theta / kMuO;
		}
			
		if (lsVM >= 0 && lsVE >= 0) {
			muV_X_E = kMuI * (vE - vM) / kDx;
		}
		else if (lsVM < 0 && lsVE < 0) {
			muV_X_E = kMuO * (vE - vM) / kDx;
		}
		else if (lsVM >= 0 && lsVE < 0) {
			// interface lies between v[i, j] and v[i + 1, j]
			theta = fabs(lsVE) / (fabs(lsVE) + fabs(lsVM));
			// |(lsVM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsVE)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsVE)|
			JVM = (kMuO - kMuI) * (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JVE = (kMuO - kMuI) * (m_J21[idx(i + 1, j)] + m_J21[idx(i + 1, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVE;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_E = muEff * (vE - vM) / kDx + muEff * JEff * theta / kMuO;
		}
		else if (lsVM < 0 && lsVE >= 0) {
			// interface lies between v[i, j] and v[i + 1, j]
			theta = fabs(lsVE) / (fabs(lsVE) + fabs(lsVM));
			// |(lsvM)| ===    outside(-)   === |(interface)| === inside(+)  === |(lsvE)|
			// |(lsvM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsvE)|
			JVM = (kMuI - kMuO) * (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JVE = (kMuI - kMuO) * (m_J21[idx(i + 1, j)] + m_J21[idx(i + 1, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVE;
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muV_X_E = muEff * (vE - vM) / kDx - muEff * JEff * theta / kMuI;
		}
		
		if (lsVS >= 0 && lsVM >= 0) {
			muV_Y_S = kMuI * (vM - vS) / kDy;
			muU_Y_W = kMuI * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
		}
		else if (lsVS < 0 && lsVM < 0) {
			muV_Y_S = kMuO * (vM - vS) / kDy;
			muU_Y_W = kMuO * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
		}
		else if (lsVS >= 0 && lsVM < 0) {
			// interface lies between v[i, j] and v[i, j - 1]
			theta = fabs(lsVS) / (fabs(lsVS) + fabs(lsVM));
			// |(lsVS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsVM)|
			// |(lsVS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVS = (kMuO - kMuI) * (m_J22[idx(i, j - 2)] + m_J22[idx(i, j - 1)]) * 0.5;
			JVM = (kMuO - kMuI) * (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVS;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_Y_S = muEff * (vM - vS) / kDy - muEff * JEff * theta / kMuO;
		}
		else if (lsVS < 0 && lsVM >= 0) {
			// interface lies between v[i, j] and v[i,j  - 1]
			theta = fabs(lsVS) / (fabs(lsVS) + fabs(lsVM));
			// |(lsVS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsVM)|
			// |(lsVS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVS = (kMuI - kMuO) * (m_J22[idx(i, j - 2)] + m_J22[idx(i, j - 1)]) * 0.5;
			JVM = (kMuI - kMuO) * (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVS;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_Y_S = muEff * (vM - vS) / kDy + muEff * JEff * theta / kMuI;
		}
		
		if (lsVM >= 0 && lsVN >= 0) {
			muV_Y_N = kMuI * (vN - vM) / kDy;
		}
		else if (lsVM < 0 && lsVN < 0) {
			muV_Y_N = kMuO * (vN - vM) / kDy;
		}
		else if (lsVM >= 0 && lsVN < 0) {
			// interface lies between v[i, j] and v[i, j + 1]
			theta = fabs(lsVN) / (fabs(lsVN) + fabs(lsVM));
			// |(lsVM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsVN)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsVN)|
			JVM = (kMuO - kMuI) * (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JVN = (kMuO - kMuI) * (m_J22[idx(i, j)] + m_J22[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVN;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_Y_N = muEff * (vN - vM) / kDy + muEff * JEff * theta / kMuO;
		}
		else if (lsVM < 0 && lsVN >= 0) {
			// interface lies between v[i, j] and v[i, j + 1]
			theta = fabs(lsVN) / (fabs(lsVN) + fabs(lsVM));
			// |(lsVM)| ===    outside(-)   === |(interface)| === inside(+) === |(lsVN)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsVN)|
			JVM = (kMuI - kMuO) * (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JVN = (kMuI - kMuO) * (m_J22[idx(i, j)] + m_J22[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVN;
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muV_Y_N = muEff * (vN - vM) / kDy - muEff * JEff * theta / kMuI;
		}
		
		// U part
		
		if (lsUS >= 0 && lsUM >= 0) {
			muU_Y_W = kMuI * (uM - uS) / kDy;
		}
		else if (lsUS < 0 && lsUM < 0) {
			muU_Y_W = kMuO * (uM - uS) / kDy;
		}
		else if (lsUS >= 0 && lsUM < 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = fabs(lsUS) / (fabs(lsUS) + fabs(lsUM));
			// |(lsUS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsUM)|
			// |(lsUS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (kMuO - kMuI) * (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (kMuO - kMuI) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_W = muEff * (uM - uS) / kDy - muEff * JEff * theta / kMuO;
		}
		else if (lsUS < 0 && lsUM >= 0) {
			// interface lies between v[i,j] and v[i,j - 1]
			theta = fabs(lsUS) / (fabs(lsUS) + fabs(lsUM));
			// |(lsUS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsUM)|
			// |(lsUS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (kMuI - kMuO) * (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (kMuI - kMuO) * (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_Y_W = muEff * (uM - uS) / kDy + muEff * JEff * theta / kMuI;
		}

		if (lsUE >= 0 && lsUE_S >= 0) {
			muU_Y_E = kMuI * (uE - uE_S) / kDy;
		}
		else if (lsUE < 0 && lsUE_S < 0) {
			muU_Y_E = kMuO * (uE - uE_S) / kDy;
		}
		else if (lsUE >= 0 && lsUE_S < 0) {
			// interface lies between u[i + 1, j] and u[i + 1, j - 1]
			theta = fabs(lsUE) / (fabs(lsUE) + fabs(lsUE_S));
			// |(lsUE)| === inside(+) === |(interface)| ===    outside(-)   === |(lsUE_S)|
			// |(lsUE)| === theta * d === |(interface)| === (1 - theta) * d === |(lsUE_S)|
			JUE_S = (kMuO - kMuI) * (m_J12[idx(i, j - 1)] + m_J12[idx(i + 1, j - 1)]) * 0.5;
			JUE = (kMuO - kMuI) * (m_J12[idx(i, j)] + m_J12[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUE_S + (1.0 - theta) * JUE;
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_E = muEff * (uE_S - uE) / kDy + muEff * JEff * theta / kMuO;
		}
		else if (lsUE < 0 && lsUE_S >= 0) {
			// interface lies between v[i,j] and v[i,j - 1]
			theta = fabs(lsUE) / (fabs(lsUE) + fabs(lsUE_S));
			// |(lsUE)| === outside(-) === |(interface)| ===     inside(+)   === |(lsUE_S)|
			// |(lsUE)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUE_S)|
			JUE_S = (kMuI - kMuO) * (m_J12[idx(i, j - 1)] + m_J12[idx(i + 1, j - 1)]) * 0.5;
			JUE = (kMuI - kMuO) * (m_J12[idx(i, j)] + m_J12[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUE_S + (1.0 - theta) * JUE;
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_Y_E = muEff * (uE_S - uE) / kDy - muEff * JEff * theta / kMuI;
		}
		
		visX = (muU_Y_E - muU_Y_W) / kDx + (muV_X_E - muV_X_W) / kDx;
		visY = 2.0 * (muV_Y_N - muV_Y_S) / kDy;
		
		dV[idx(i, j)] = visX + visY;
		
		assert(dV[idx(i, j)] == dV[idx(i, j)]);
		if (std::isnan(dV[idx(i, j)]) || std::isinf(dV[idx(i, j)])) {
			std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << dV[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return dV;
}

std::vector<double> MACSolver2D::AddGravityFU() {
	std::vector<double> gU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if (kFr == 0 && !isnan(kFr) && !isinf(kFr)) {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			gU[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			gU[idx(i, j)] = 0.0;
		}
	}

	return gU;
}

std::vector<double> MACSolver2D::AddGravityFV() {
	std::vector<double> gV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if (kFr == 0 && !isnan(kFr) && !isinf(kFr)) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			gV[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			gV[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gV;
}

int MACSolver2D::GetIntermediateVel(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& lsB, const std::vector<double>& u, const std::vector<double>& v,
	std::vector<double>& uhat, std::vector<double>& vhat) {
		
	// Update rhs
	if (kTimeOrder == TimeOrderEnum::EULER) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
		std::vector<double> FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		FU1 = UpdateFU(LSolver, lsB, u, v);
		FV1 = UpdateFV(LSolver, lsB, u, v);
		
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			uhat[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			vhat[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
	}
	else if (kTimeOrder == TimeOrderEnum::RK2) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FU2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, lsB, u, v);
		FV1 = UpdateFV(LSolver, lsB, u, v);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : L(u^(1))
		// FU2 & FV2 : u^(2)
		ApplyBC_U_2D(FU2);
		ApplyBC_V_2D(FV2);
		FU1 = UpdateFU(LSolver, lsB, FU2, FV2);
		FV1 = UpdateFV(LSolver, lsB, FU2, FV2);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			uhat[idx(i, j)] = 0.5 * u[idx(i, j)] + 0.5 * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			vhat[idx(i, j)] = 0.5 * v[idx(i, j)] + 0.5 * FV2[idx(i, j)];
		}
	}
	else if (kTimeOrder == TimeOrderEnum::RK3) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FU2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, lsB, u, v);
		FV1 = UpdateFV(LSolver, lsB, u, v);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : L(u^(1))
		// FU2 & FV2 : u^(1) + \delta t L (u^(1))
		ApplyBC_U_2D(FU2);
		ApplyBC_V_2D(FV2);
		FU1 = UpdateFU(LSolver, lsB, FU2, FV2);
		FV1 = UpdateFV(LSolver, lsB, FU2, FV2);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : u^(2)
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU1[idx(i, j)] = 0.75 * u[idx(i, j)] + 0.25 * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV1[idx(i, j)] = 0.75 * v[idx(i, j)] + 0.25 * FV2[idx(i, j)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);

		// FU2 & FV2 : L(u^(2))
		// FU1 & FV1 : u^(2) + \delta t L (u^(2))
		ApplyBC_U_2D(FU1);
		ApplyBC_V_2D(FV1);
		FU2 = UpdateFU(LSolver, lsB, FU1, FV1);
		FV2 = UpdateFV(LSolver, lsB, FU1, FV1);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			FU1[idx(i, j)] = FU1[idx(i, j)] + m_dt * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			FV1[idx(i, j)] = FV1[idx(i, j)] + m_dt * FV2[idx(i, j)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);

		// FU1 & FV1 : u^(2) + \delta t L (u^(2))
		// FU2 & FV2 : doesn't need. set to zero.
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			uhat[idx(i, j)] = 1.0 / 3.0 * u[idx(i, j)] + 2.0 / 3.0 * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			vhat[idx(i, j)] = 1.0 / 3.0 * v[idx(i, j)] + 2.0 / 3.0 * FV1[idx(i, j)];
		}
	}

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		if (std::isnan(uhat[idx(i, j)]) || std::isinf(uhat[idx(i, j)])) {
			std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << uhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		if (std::isnan(vhat[idx(i, j)]) || std::isinf(vhat[idx(i, j)])) {
			std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << vhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	return 0;
}

int MACSolver2D::SetPoissonSolver(POISSONTYPE type) {
	m_PoissonSolverType = type;
	if (!m_Poisson)
		m_Poisson = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);

	return 0;
}

int MACSolver2D::SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
	const std::vector<double>& lsB, const std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const int maxIter) {
	if (!m_Poisson) {
		perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
	}
	std::vector<double> rhs((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> iRhoW((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoE((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoN((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	/*
	Solver \nabla \cdot ((\nabla p^*) / (\rho)) = rhs

	Liu, Xu-Dong, Ronald P. Fedkiw, and Myungjoo Kang.
	"A boundary condition capturing method for Poisson's equation on irregular domains."
	Journal of Computational Physics 160.1 (2000): 151-178.
	for level set jump condition
	[p^*] - 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N = dt \sigma \kappa
	[p^*] = 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa

	[p^*_x \ rho] = [((2 mu u_x)_x  + (mu(u_y + v_x))_y) / rho]
	[p^*_y \ rho] = [((mu(u_y + v_x))_x  + (2 mu v_y)_y  ) / rho]
	However, for solving poisson equation,
	[p^*_x \ rho] = 0.0
	[p^*_y \ rho] = 0.0
	why?
	*/
	double eps = 1.0e-100;
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double FW = 0.0, FE = 0.0, FS = 0.0, FN = 0.0;
	// derivative
	double duWX = 0.0, duWY = 0.0, duEX = 0.0, duEY = 0.0, duSX = 0.0, duSY = 0.0, duNX = 0.0, duNY = 0.0, duMX = 0.0, duMY = 0.0;
	double dvWX = 0.0, dvWY = 0.0, dvEX = 0.0, dvEY = 0.0, dvSX = 0.0, dvSY = 0.0, dvNX = 0.0, dvNY = 0.0, dvMX = 0.0, dvMY = 0.0;
	double dlWX = 0.0, dlWY = 0.0, dlEX = 0.0, dlEY = 0.0, dlSX = 0.0, dlSY = 0.0, dlNX = 0.0, dlNY = 0.0, dlMX = 0.0, dlMY = 0.0;
	// normal and tangent vector variable (n, t1, t2)
	double nXW = 0.0, nYW = 0.0, nXE = 0.0, nYE = 0.0, nXS = 0.0, nYS = 0.0, nXN = 0.0, nYN = 0.0, nXM = 0.0, nYM = 0.0;
	// jump at grid node (if interface is at grid node, jump occurs and aW, aE, aS, aN, aM describe that condition)
	// aW, aE, aS, aN, aM : at P grid
	double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aM = 0.0;
	// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
	double a = 0.0, b = 0.0;
	double theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, bEff = 0.0;
	
	if (kWe != 0.0 && !isnan(kWe) && !isinf(kWe)) {
		UpdateKappa(ls);
		ApplyBC_P_2D(m_kappa);
	}
	// A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
	std::vector<double> AVals;
	std::vector<MKL_INT> ACols, ARowIdx, MCols, MRowIdx;
	MKL_INT rowIdx = 0, tmpRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	
	std::vector<double> dudX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dudY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	// stored coef for A matrix, Dictionary but it is ordered
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
//		ps[idx(i, j)] = 0.0;
		rhs[idx(i, j)] = 0.0;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = lsB[idx(i - 1, j)];
		lsE = lsB[idx(i + 1, j)];
		lsM = lsB[idx(i, j)];
		lsS = lsB[idx(i, j - 1)];
		lsN = lsB[idx(i, j + 1)];

		// At P grid
		dudX[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
		dudY[idx(i, j)] = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
			- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);

		dvdX[idx(i, j)] = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
		dvdY[idx(i, j)] = (v[idx(i, j + 1)] - v[idx(i, j)]) / (kDy);

		dldX[idx(i, j)] = (lsE - lsW) / (2.0 * kDx);
		dldY[idx(i, j)] = (lsN - lsS) / (2.0 * kDy);
		if (std::fabs(dldX[idx(i, j)]) < eps) {
			dldX[idx(i, j)] = (lsE - lsM) / kDx;
			if (std::fabs(dldX[idx(i, j)]) < eps)
				dldX[idx(i, j)] = (lsM - lsW) / kDx;
		}
		if (std::fabs(dldY[idx(i, j)]) < eps) {
			dldY[idx(i, j)] = (lsN - lsM) / kDy;
			if (std::fabs(dldY[idx(i, j)]) < eps)
				dldY[idx(i, j)] = (lsM - lsS) / kDy;
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsE = ls[idx(i + 1, j)];
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		lsN = ls[idx(i, j + 1)];

		FW = 0.0;
		FE = 0.0;
		FS = 0.0;
		FN = 0.0;

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nXW = dldX[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps);
		nYW = dldY[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps);
		nXE = dldX[idx(i + 1, j)] / (std::sqrt(std::pow(dldX[idx(i + 1, j)], 2.0) + std::pow(dldY[idx(i + 1, j)], 2.0)) + eps);
		nYE = dldY[idx(i + 1, j)] / (std::sqrt(std::pow(dldX[idx(i + 1, j)], 2.0) + std::pow(dldY[idx(i + 1, j)], 2.0)) + eps);
		nXS = dldX[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps);
		nYS = dldY[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps);
		nXN = dldX[idx(i, j + 1)] / (std::sqrt(std::pow(dldX[idx(i, j + 1)], 2.0) + std::pow(dldY[idx(i, j + 1)], 2.0)) + eps);
		nYN = dldY[idx(i, j + 1)] / (std::sqrt(std::pow(dldX[idx(i, j + 1)], 2.0) + std::pow(dldY[idx(i, j + 1)], 2.0)) + eps);
		nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);
		nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);

		/*
		// [a]_\Gamma = a^+ - a^- = a^inside - a^outside
		// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// (P_M - P_W)/kDx appears
		// if P_M == P_+, P_W == P_-, a^+ means all terms related to P_+, P_W changed to P_M related terms
		// P_W = P_M -(2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa)
		*/
		// FW
		if (lsW * lsM > 0.0 || (lsW <= 0.0 && lsM <= 0.0)) {
			// one fluid, x direction
			FW = 0.0;
			if (lsW <= 0.0)
				iRhoW[idx(i, j)] = 1.0 / kRhoO;
			else
				iRhoW[idx(i, j)] = 1.0 / kRhoI;
		}
		else if (lsW > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
			// b always zero when solving level set (dealing with surface tension)
			aW = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
			iRhoW[idx(i, j)] = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsW) + std::fabs(lsM))
				/ (1.0 / kRhoO * std::fabs(lsW) + 1.0 / kRhoI * std::fabs(lsM));
			FW = iRhoW[idx(i, j)] * aEff / (kDx * kDx);
		}
		else if (lsW <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
			aW = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
			iRhoW[idx(i, j)] = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsW) + std::fabs(lsM))
				/ (1.0 / kRhoI * std::fabs(lsW) + 1.0 / kRhoO * std::fabs(lsM));
			FW = -iRhoW[idx(i, j)] * aEff / (kDx * kDx);
		}
		
		// FE
		// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
		if (lsM * lsE > 0.0 || (lsM <= 0.0 && lsE <= 0.0)) {
			// one fluid, x direction
			FE = 0.0;
			if (lsE <= 0.0)
				iRhoE[idx(i, j)] = 1.0 / kRhoO;
			else
				iRhoE[idx(i, j)] = 1.0 / kRhoI;
		}
		else if (lsM > 0.0 && lsE <= 0.0) {
			// interface lies between ls[i, j] and ls[i + 1, j]
			// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
				+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = (aM * std::fabs(lsE) + aE * std::fabs(lsM)) / (std::fabs(lsE) + std::fabs(lsM));
			iRhoE[idx(i, j)] = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsE))
				/ (1.0 / kRhoI * std::fabs(lsE) + 1.0 / kRhoO * std::fabs(lsM));
			FE = -iRhoE[idx(i, j)] * aEff / (kDx * kDx);
		}
		else if (lsM <= 0.0 && lsE > 0.0) {
			// interface lies between ls[i, j] and ls[i + 1, j]
			// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
				+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = (aM * std::fabs(lsE) + aE * std::fabs(lsM)) / (std::fabs(lsE) + std::fabs(lsM));
			iRhoE[idx(i, j)] = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsE))
				/ (1.0 / kRhoO * std::fabs(lsE) + 1.0 / kRhoI * std::fabs(lsM));
			FE = iRhoE[idx(i, j)] * aEff / (kDx * kDx);
		}

		// FS
		if (lsS * lsM > 0.0 || (lsS <= 0.0 && lsM <= 0.0)) {
			// one fluid, y direction
			FS = 0.0;
			if (lsS <= 0.0)
				iRhoS[idx(i, j)] = 1.0 / kRhoO;
			else
				iRhoS[idx(i, j)] = 1.0 / kRhoI;
		}
		else if (lsS > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i, j] and ls[i, j - 1]
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsS) + std::fabs(lsM));
			iRhoS[idx(i, j)] = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsS) + std::fabs(lsM))
				/ (1.0 / kRhoO * std::fabs(lsS) + 1.0 / kRhoI * std::fabs(lsM));
			FS = iRhoS[idx(i, j)] * aEff / (kDy * kDy);
		}
		else if (lsS <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i, j - 1]
			// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsS) + std::fabs(lsM));
			iRhoS[idx(i, j)] = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsS) + std::fabs(lsM))
				/ (1.0 / kRhoI * std::fabs(lsS) + 1.0 / kRhoO * std::fabs(lsM));
			FS = -iRhoS[idx(i, j)] * aEff / (kDy * kDy);
		}

		// FN
		if (lsM * lsN > 0.0 || (lsM <= 0.0 && lsN <= 0.0)) {
			// one fluid, y direction
			FN = 0.0;
			if (lsN <= 0.0)
				iRhoN[idx(i, j)] = 1.0 / kRhoO;
			else
				iRhoN[idx(i, j)] = 1.0 / kRhoI;
		}
		else if (lsM > 0.0 && lsN <= 0.0) {
			// interface lies between ls[i, j] and ls[i, j + 1]
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nYN) * nXN
				+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = (aM * std::fabs(lsN) + aN * std::fabs(lsM)) / (std::fabs(lsN) + std::fabs(lsM));
			iRhoN[idx(i, j)] = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsN))
				/ (1.0 / kRhoI * std::fabs(lsN) + 1.0 / kRhoO * std::fabs(lsM));
			FN = -iRhoN[idx(i, j)] * aEff / (kDy * kDy);
		}
		else if (lsM <= 0.0 && lsN > 0.0) {
			// interface lies between ls[i, j] and ls[i, j + 1]
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nYN) * nXN
				+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = (aM * std::fabs(lsN) + aN * std::fabs(lsM)) / (std::fabs(lsN) + std::fabs(lsM));
			iRhoN[idx(i, j)] = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsN))
				/ (1.0 / kRhoO * std::fabs(lsN) + 1.0 / kRhoI * std::fabs(lsM));
			FN = iRhoN[idx(i, j)] * aEff / (kDy * kDy);
		}
		
		// -= due to poisson equation form should be -\beta \nabla p = f
		rhs[idx(i, j)] += FW + FE + FS + FN;
		
		assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
		if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " 
				<< rhs[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	// Original value of RHS
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		rhs[idx(i, j)] += div[idx(i, j)];

	// An order of A matrix coef. is very important, hence reverse j order
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
		// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
		// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
		AValsDic.clear();
		AColsDic.clear();
		tmpRowIdx = 0;
		// Add starting rowIdx
		ARowIdx.push_back(rowIdx);

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["S"] = iRhoS[idx(i, j)] / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
		AValsDic["W"] = iRhoW[idx(i, j)] / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["C"] = -(iRhoW[idx(i, j)] + iRhoE[idx(i, j)]) / (kDx * kDx) - (iRhoS[idx(i, j)] + iRhoN[idx(i, j)]) / (kDy * kDy);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["E"] = iRhoE[idx(i, j)] / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["N"] = iRhoN[idx(i, j)] / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;
		
		if (i == kNumBCGrid && m_BC->m_BC_PW == BC::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += iRhoW[idx(i, j)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= iRhoW[idx(i, j)] / (kDx * kDx);
			rhs[idx(i, j)] -= iRhoW[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPW);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC::PERIODIC) {
			AValsDic["W"] = iRhoW[idx(kNumBCGrid + kNx - 1, j)];
		}
		
		// East boundary
		if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += iRhoE[idx(i, j)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= iRhoE[idx(i, j)] / (kDx * kDx);
			rhs[idx(i, j)] -= iRhoE[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPE);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::PERIODIC) {
			AValsDic["E"] = iRhoE[idx(kNumBCGrid, j)];
		}

		if (j == kNumBCGrid && m_BC->m_BC_PS == BC::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += iRhoS[idx(i, j)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= iRhoS[idx(i, j)] / (kDy * kDy);
			rhs[idx(i, j)] -= iRhoS[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPS);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC::PERIODIC) {
			AValsDic["S"] = iRhoS[idx(i, kNumBCGrid + kNy - 1)];
		}
		
		if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += iRhoN[idx(i, j)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= iRhoN[idx(i, j)] / (kDy * kDy);
			rhs[idx(i, j)] -= iRhoN[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPN);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::PERIODIC) {
			AValsDic["N"] = iRhoE[idx(i, kNumBCGrid + kNy + 1)];
		}
		// add non zero values to AVals and ACols
		// KEEP ORDER OF PUSH_BACK!!

		if (AColsDic["S"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["S"]);
			ACols.push_back(AColsDic["S"]);
		}

		if (AColsDic["W"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["W"]);
			ACols.push_back(AColsDic["W"]);
		}

		// Center, it doens't neeed to check
		tmpRowIdx++;
		AVals.push_back(AValsDic["C"]);
		ACols.push_back(AColsDic["C"]);

		if (AColsDic["E"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["E"]);
			ACols.push_back(AColsDic["E"]);
		}

		if (AColsDic["N"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["N"]);
			ACols.push_back(AColsDic["N"]);
		}
		
		rowIdx += tmpRowIdx;
		
		assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
		if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " 
				<< i << " " << j << " " << rhs[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	ARowIdx.push_back(rowIdx);

	if (m_PoissonSolverType == POISSONTYPE::MKL) {
		// For Level set solver test only, legacy poisson equation solver
		/*
		Solver \nabla \cdot ((\nabla p^*) / (\rho)) = rhs

		Liu, Xu-Dong, Ronald P. Fedkiw, and Myungjoo Kang.
		"A boundary condition capturing method for Poisson's equation on irregular domains."
		Journal of Computational Physics 160.1 (2000): 151-178.
		for level set jump condition
		[p^*] - 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N = dt \sigma \kappa
		[p^*] = 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa

		[p^*_x \ rho] = [((2 mu u_x)_x  + (mu(u_y + v_x))_y) / rho]
		[p^*_y \ rho] = [((mu(u_y + v_x))_x  + (2 mu v_y)_y  ) / rho]
		However, for solving poisson equation, 
		[p^*_x \ rho] = 0.0
		[p^*_y \ rho] = 0.0
		why? 
		*/

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			rhs[idx(i, j)] = -div[idx(i, j)];
		
		m_Poisson->MKL_2FUniform_2D(ps, rhs, kLenX, kLenY, kDx, kDy, m_BC);
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;

	}
	else if (m_PoissonSolverType == POISSONTYPE::CG) {
		// std::cout << "Poisson : CG" << std::endl;
		m_Poisson->CG_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);
	
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
		// std::cout << "Poisson : BiCG" << std::endl;
		m_Poisson->BiCGStab_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::GS) {
		// std::cout << "Poisson : GS" << std::endl;
		m_Poisson->GS_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	
	ApplyBC_P_2D(ps);
	
	return 0;
}

std::vector<double> MACSolver2D::GetDivergence(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> div((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		div[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
						+ (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
	
	return div;
}

int MACSolver2D::UpdateVel(std::vector<double>& u, std::vector<double>& v,
	const std::vector<double>& us, const std::vector<double>& vs, const std::vector<double>& ps, const std::vector<double>& ls) {
	
	// velocity update after solving poisson equation
	// ps = p * dt
	double lsW = 0.0, lsM = 0.0, lsS = 0.0, iRhoEff = 0.0, theta = 0.0;
	double nXW = 0.0, nYW = 0.0, nXS = 0.0, nYS = 0.0, nXM = 0.0, nYM = 0.0;
	double aW = 0.0, aS = 0.0, aM = 0.0, aEff = 0.0;
	const double eps = 1.0e-100;
	std::vector<double> dudX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dudY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		// At P grid
		dudX[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
		dudY[idx(i, j)] = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
			- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);

		dvdX[idx(i, j)] = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
		dvdY[idx(i, j)] = (v[idx(i, j + 1)] - v[idx(i, j)]) / (kDy);

		dldX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dldY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);

		if (std::fabs(dldX[idx(i, j)]) < eps) {
			dldX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			if (std::fabs(dldX[idx(i, j)]) < eps)
				dldX[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		}
		if (std::fabs(dldY[idx(i, j)]) < eps) {
			dldY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			if (std::fabs(dldY[idx(i, j)]) < eps)
				dldY[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
		}
	}

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];

		nXW = dldX[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps);
		nYW = dldY[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps);
		nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);
		nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);

		aEff = 0.0;
		if (lsW >= 0 && lsM >= 0) {
			iRhoEff = 1.0 / kRhoI;
		}
		else if (lsW < 0 && lsM < 0) {
			iRhoEff = 1.0 / kRhoO;
		}
		else if (lsW > 0 && lsM <= 0) {
			// interface lies between ls[i - 1, j] and ls[i,j]
			// |(lsW)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aW = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
			iRhoEff = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsW) + std::fabs(lsM))
				/ (1.0 / kRhoO * std::fabs(lsW) + 1.0 / kRhoI * std::fabs(lsM));
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between ls[i - 1, j] and ls[i,j]
			// |(lsW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aW = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
			iRhoEff = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsW) + std::fabs(lsM))
				/ (1.0 / kRhoI * std::fabs(lsW) + 1.0 / kRhoO * std::fabs(lsM));
		}
		
		u[idx(i, j)] = us[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx + iRhoEff * aEff;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];

		nXS = dldX[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps);
		nYS = dldY[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps);
		nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);
		nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);

		aEff = 0.0;
		if (lsS >= 0 && lsM >= 0) {
			iRhoEff = 1.0 / kRhoI;
		}
		else if (lsS < 0 && lsM < 0) {
			iRhoEff = 1.0 / kRhoO;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuO - kMuI)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsS) + std::fabs(lsM));
			iRhoEff = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsS) + std::fabs(lsM))
				/ (1.0 / kRhoO * std::fabs(lsS) + 1.0 / kRhoI * std::fabs(lsM));
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuI - kMuO)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsS) + std::fabs(lsM));
			iRhoEff = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsS) + std::fabs(lsM))
				/ (1.0 / kRhoI * std::fabs(lsS) + 1.0 / kRhoO * std::fabs(lsM));
		}

		v[idx(i, j)] = vs[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy + iRhoEff * aEff;
	}

	ApplyBC_U_2D(u);
	ApplyBC_V_2D(v);
	
	return 0;
}

double MACSolver2D::UpdateDt(const std::vector<double>& u, const std::vector<double>& v) {
	double uAMax = 0.0;
	double vAMax = 0.0;
	double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0;
	double dt = std::numeric_limits<double>::max();

	// get maximum of absolute value
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		uAMax = std::max(uAMax, std::fabs((u[idx(i + 1, j)] * u[idx(i, j)]) * 0.5));
		vAMax = std::max(vAMax, std::fabs((v[idx(i, j + 1)] * v[idx(i, j)]) * 0.5));
	}
	
	Cefl = uAMax / kDx + vAMax / kDy;
	Vefl = std::max(kMuI / kRhoI, kMuO / kRhoO) * (2.0 / (kDx * kDx) + 2.0 / (kDy * kDy));
	Gefl = std::sqrt(std::fabs(kG) / kDy);
	
	dt = std::min(dt,
		1.0 / (0.5 * (Cefl + Vefl +
		std::sqrt(std::pow(Cefl + Vefl, 2.0) +
		4.0 * Gefl))));

	if (std::isnan(dt) || std::isinf(dt)) {
		std::cout << "dt nan/inf error : Cefl, Vefl, Gefl : " << Cefl << " " << Vefl << " " << Gefl << std::endl;
	}

	return dt;
}

int MACSolver2D::SetBC_U_2D(std::string BC_W, std::string BC_E,	std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2D::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2D::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_P_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2D::ApplyBC_U_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_U_2D(arr);
	return 0;
}

int MACSolver2D::ApplyBC_V_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_V_2D(arr);
	return 0;
}

int MACSolver2D::ApplyBC_P_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_P_2D(arr);
	return 0;
}

void MACSolver2D::SetBCConstantUW(double BC_ConstantW) {
	return m_BC->SetBCConstantUW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantUE(double BC_ConstantE) {
	return m_BC->SetBCConstantUE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantUS(double BC_ConstantS) {
	return m_BC->SetBCConstantUS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantUN(double BC_ConstantN) {
	return m_BC->SetBCConstantUN(BC_ConstantN);
}

void MACSolver2D::SetBCConstantVW(double BC_ConstantW) {
	return m_BC->SetBCConstantVW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantVE(double BC_ConstantE) {
	return m_BC->SetBCConstantVE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantVS(double BC_ConstantS) {
	return m_BC->SetBCConstantVS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantVN(double BC_ConstantN) {
	return m_BC->SetBCConstantVN(BC_ConstantN);
}

void MACSolver2D::SetBCConstantPW(double BC_ConstantW) {
	return m_BC->SetBCConstantPW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantPE(double BC_ConstantE) {
	return m_BC->SetBCConstantPE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantPS(double BC_ConstantS) {
	return m_BC->SetBCConstantPS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantPN(double BC_ConstantN) {
	return m_BC->SetBCConstantPN(BC_ConstantN);
}

inline int MACSolver2D::idx(int i, int j) {
	return (j + (kNy + 2 * kNumBCGrid) * (i));
}

int MACSolver2D::SetPLTType(PLTTYPE type) {
	m_PLTType = type;

	return 0;
}

int MACSolver2D::OutRes(const int iter, const double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ps, const std::vector<double>& ls) {
	if (m_PLTType == PLTTYPE::ASCII || m_PLTType == PLTTYPE::BOTH) {
		std::ofstream outF;
		std::string fname_vel(fname_vel_base + "_ASCII.plt");
		std::string fname_div(fname_div_base + "_ASCII.plt");
		if (m_iter == 0) {
			outF.open(fname_vel.c_str(), std::ios::out);

			outF << "TITLE = VEL" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"PS\" " << std::endl;
			outF.close();

			outF.open(fname_div.c_str(), std::ios::out);
			outF << "TITLE = DIV" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"LS\", \"DIV\", \"PS\" " << std::endl;
			outF.close();
		}

		std::vector<double>
			resU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			resV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		std::vector<double> resDiv = GetDivergence(u, v);
		
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			resU[idx(i, j)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			resV[idx(i, j)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;

		outF.open(fname_vel.c_str(), std::ios::app);
		
		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
				<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
				<< static_cast<double>(resU[idx(i, j)]) << std::string(",")
				<< static_cast<double>(resV[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j)]) << std::endl;

		outF.close();

		outF.open(fname_div.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
				<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
				<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
				<< static_cast<double>(resDiv[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j)]) << std::endl;

		outF.close();
	}

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH) {
		INTEGER4 whichFile = 0, stat = 0;
		std::string fname_vel(fname_vel_base + "_BINARY.szplt");
		std::string fname_div(fname_div_base + "_BINARY.szplt");

		std::vector<double>
			resX(kNx * kNy, 0.0),
			resY(kNx * kNy, 0.0),
			resU(kNx * kNy, 0.0),
			resV(kNx * kNy, 0.0),
			resLS(kNx * kNy, 0.0),
			resPs(kNx * kNy, 0.0);

		std::vector<double> resDiv = GetDivergence(u, v);
		
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)	 {
			resX[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseX + (i + 0.5 - kNumBCGrid) * kDx;
			resY[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseY + (j + 0.5 - kNumBCGrid) * kDy;
			resU[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
			resV[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ls[idx(i, j)];
			resPs[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ps[idx(i, j)];
		}

		INTEGER4 Debug = 1;
		INTEGER4 VIsDouble = 1; /* 0 = Single precision, 1 = Double precision*/
		INTEGER4 FileType = 0; /* 0 = Full, 1 = Grid only, 2 = Solution only*/
		INTEGER4 FileFormat = 1; /* 0 = plt, 1 = szplt*/

		if (m_iter == 0) {
			/*
			* Open the file and write the tecplot datafile
			* header information
			*/
			stat = TECINI142(const_cast<char *>(std::string("VELOCITY").c_str()),  /* Name of the entire dataset.  */
							const_cast<char *>(std::string("X, Y, U, V, LS, PS").c_str()),  
							/* Defines the variables for the data file. Each zone must contain each of the variables listed here. 
							* The order of the variables in the list is used to define the variable number (e.g. X is Var 1).*/
							const_cast<char *>(fname_vel.c_str()),
							const_cast<char *>(std::string(".").c_str()),      /* Scratch Directory */
								&FileFormat, &FileType, &Debug, &VIsDouble);
		}
		else {
			whichFile = 1;
			stat = TECFIL142(&whichFile);
		}

		/* Set the parameters for TecZne */
		INTEGER4 ZoneType = 0; /* sets the zone type to
							   * ordered
							   */
		/* Create an IJ-ordered zone, by using IMax and JMax
		* values that are greater than one, and setting KMax to one.
		*/
		INTEGER4 IMax = kNx, JMax = kNy, KMax = 1;

		double   SolTime = curTime;
		INTEGER4 StrandID = iter + 1; /* StaticZone */
		INTEGER4 ParentZn = 0; /* used for surface streams */

		INTEGER4 ICellMax = 0, JCellMax = 0, KCellMax = 0; /* not used */

		INTEGER4 IsBlock = 1; /* Block */

		INTEGER4 NFConns = 0; /* this example does not use face neighbors */
		INTEGER4 FNMode = 0;
		INTEGER4 TotalNumFaceNodes = 1, TotalNumBndryFaces = 1, TotalNumBndryConn = 1;
		INTEGER4 ShrConn = 0;

		/* Create an Ordered Zone */
		stat = TECZNE142((char*) std::to_string(StrandID).c_str(), &ZoneType,
			&IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax,
			&SolTime, &StrandID, &ParentZn, &IsBlock,
			&NFConns, &FNMode, &TotalNumFaceNodes,
			&TotalNumBndryFaces, &TotalNumBndryConn,
			NULL, NULL, NULL, &ShrConn);
			
		INTEGER4 DIsDouble = 1;  /* set DIsDouble to 0, for float
								 * values.
								 */

		INTEGER4 ARRSIZEVAL = IMax * JMax * KMax;

		stat = TECDAT142(&ARRSIZEVAL, resX.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resY.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resU.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resV.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);

		if (m_iter == 0) {
			stat = TECINI142(const_cast<char *>(std::string("DIVERGENCE").c_str()),  /* Name of the entire dataset.  */
				const_cast<char *>(std::string("X, Y, LS, DIV, PS").c_str()),
				/* Defines the variables for the data file. Each zone must contain each of the variables listed here.
				* The order of the variables in the list is used to define the variable number (e.g. X is Var 1).*/
				const_cast<char *>(fname_div.c_str()),
				const_cast<char *>(std::string(".").c_str()),      /* Scratch Directory */
				&FileFormat, &FileType, &Debug, &VIsDouble);
		}

		whichFile = 2;
		stat = TECFIL142(&whichFile);

		/* Create an Ordered Zone */
		stat = TECZNE142((char*) std::to_string(StrandID).c_str(), &ZoneType,
			&IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax,
			&SolTime, &StrandID, &ParentZn, &IsBlock,
			&NFConns, &FNMode, &TotalNumFaceNodes,
			&TotalNumBndryFaces, &TotalNumBndryConn,
			NULL, NULL, NULL, &ShrConn);

		stat = TECDAT142(&ARRSIZEVAL, resX.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resY.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resDiv.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);
	}

	return 0;
}

int MACSolver2D::OutResClose() {
	INTEGER4 stat, whichFile;
	whichFile = 1;
	stat = TECFIL142(&whichFile);

	// close first file (velocity)
	stat = TECEND142();
	stat = TECEND142();

	return 0;
}
