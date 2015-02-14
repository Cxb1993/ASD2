#include "mac2d.h"

MACSolver2D::MACSolver2D(double Re, double We, double Fr,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoI, double muI,
	int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	double cfl,	int maxtime, int maxiter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoI(rhoI), kRhoO(rhoI / densityRatio), kRhoRatio(densityRatio),
	kMuI(muI), kMuO(muI / viscosityRatio), kMuRatio(viscosityRatio),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxiter), kNIterSkip(niterskip), kNumBCGrid(num_bc_grid),
	kWriteVTK(writeVTK) {

	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::MACSolver2D(double rhoI, double rhoO, double muI, double muO, double gConstant,
	double L, double U, double sigma, int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	double cfl, int maxtime, int maxiter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoI), kMuScale(muI), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoI(rhoI), kRhoO(rhoO), kMuI(muI), kMuO(muO), kRhoRatio(rhoI / rhoO), kMuRatio(muI / muO),
	kRe(rhoI * L * U / muI), kWe(rhoI * L * U * U / sigma), kFr(U * U / (gConstant * L)),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxiter), kNIterSkip(niterskip), kNumBCGrid(num_bc_grid),
	kWriteVTK(writeVTK) {

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
	
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		
		LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
	}

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		if (LSSize[idx(i, j)] == 0.0) {
			dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			if (LSSize[idx(i, j)] == 0.0) {
				dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
				LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			}
		}
		if (LSSize[idx(i, j)] == 0.0)
			perror("Div/0 Err in computing kappa");

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

	// derivatives
	double dudX = 0.0, dudY = 0.0;
	double dvdX = 0.0, dvdY = 0.0;
	double dldX = 0.0, dldY = 0.0;
	// normal and tangent vector variable (n, t1, t2)
	// normal vector has X, Y component
	// first tangent vector has only Z component
	double nX = 0.0, nY = 0.0, t1Z = 0.0, t2X = 0.0, t2Y = 0.0;

	// jump condition originally defined as J11 = [\mu u_x], J12 = [\mu u_y], J21 = [\mu v_x], J22 = [\mu v_y], 
	// but skip [\mu], from Kang, Fedkiw, and Liu's work eq. (30).
	
	double nxnx = 0.0, nxny = 0.0, nyny = 0.0;
	double t2xt2x = 0.0, t2xt2y = 0.0, t2yt2y = 0.0;
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++)
	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++) {
		dudX = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
		dudY = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
			- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);
		dvdX = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
		dvdY = (v[idx(i, j + 1)] - v[idx(i, j)]) / (2.0 * kDy);
		dldX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dldY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		if (std::fabs(dldX) < 1.0e-40) {
			dldX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			if (dldX == 0.0)
				dldX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		}
		if (std::fabs(dldX) < 1.0e-40) {
			dldY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			if (dldY == 0.0)
				dldY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
		}

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nX = dldX / (std::fabs(dldX) + 1.0e-100);
		nY = dldY / (std::fabs(dldY) + 1.0e-100);
		// get first tangent vector
		if (std::fabs(nX) == std::min(std::fabs(nX), std::fabs(nY))) {
			// smallest magnitude determine what to be performed cross product
			// T1 = N X (1, 0) / |N X (1, 0)|
			// T1 = N X (1, 0) = (nX, nY) X (1, 0) = -nY \hat{k} / |nY| = (0, 0, t1Z)
			// T2 = N X T1 = (nX, nY, 0) X (0, 0, t1Z) = nY * t1Z \hat{i} - nX * t1Z \hat{j} 
			t1Z = -nY / std::fabs(nY);
		}
		else {
			// smallest magnitude determine what to be performed cross product
			// T1 = N X (0, 1) / |N X (0, 1)|
			// T1 = N X (0, 1) = (nX, nY) X (0, 1) = nX \hat{k} / |nX| = (0, 0, t1Z)
			// T2 = N X T1= (nX, nY, 0) X (0, 0, t1Z) = nY * t1Z \hat{i} - nX * t1Z \hat{j} 
			t1Z = nX / std::fabs(nX);
		}
		t2X = nY * t1Z;
		t2Y = -nX * t1Z;
		nxnx = nX * nX;
		nxny = nX * nY;
		nyny = nY * nY;
		t2xt2x = t2X * t2X;
		t2xt2y = t2Y * t2Y;
		t2yt2y = t2Y * t2Y;

		// eq. (30) from Kang, Fedkiw, and Liu
		m_J11[idx(i, j)] = dudX * t2xt2x + dudY * t2xt2y 
			+ (dudX * nxnx + dvdX * nxny) * nxnx + (dudY * nxnx + dvdY * nxny) * nxny
			+ (dudX * nxnx + dudY * nxny) * t2xt2x + (dvdX * nxnx + dvdY * nxny) * t2xt2y;
		m_J12[idx(i, j)] = dudX * t2xt2y + dudY * t2yt2y
			+ (dudX * nxnx + dvdX * nxny) * nxny + (dudY * nxnx + dvdY * nxny) * nyny
			+ (dudX * nxny + dudY * nyny) * t2xt2y + (dvdX * nxny + dvdY * nyny) * t2yt2y;
		m_J21[idx(i, j)] = dvdX * t2xt2x + dvdY * t2xt2y
			+ (dudX * nxny + dvdX * nyny) * nxnx + (dudY * nxny + dvdY * nyny) * nxny
			+ (dudX * nxnx + dudX * nxny) * t2xt2x + (dvdX * nxnx + dvdY * nxny) * t2xt2y;
		m_J22[idx(i, j)] = dvdX * t2xt2y + dvdY * t2yt2y
			+ (dudX * nxny + dvdX * nyny) * nxny + (dudY * nxny + dvdY * nyny) * nyny
			+ (dudX * nxny + dudY * nyny) * t2xt2y + (dvdX * nxny + dvdY * nyny) * t2yt2y;

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
	std::vector<double> ls, std::vector<double>u, std::vector<double> v) {
	
	std::vector<double> cU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> vU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> gU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhsU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v);

	// Viscous term
	vU = this->AddViscosityFU(u, v, ls);

	// Gravity term
	// gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
		// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, rhoEff = 0.0;
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		
		if (lsW >= 0 && lsM >= 0) {
			rhoEff = kRhoI;
		}
		else if (lsW < 0 && lsM < 0) {
			rhoEff = kRhoO;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between ls[i - 1,j] and ls[i,j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
		}
		else if (lsW < 0 && lsM >= 0) {
			// interface lies between ls[i - 1, j] and ls[i,j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
		}
		// std::cout << i << " " << j << " " << cU[idx(i, j)] << " " << vU[idx(i, j)] << " " << gU[idx(i, j)] << " " << rhoEff << std::endl;
		rhsU[idx(i, j)] = (-cU[idx(i, j)] + vU[idx(i, j)] / rhoEff + gU[idx(i, j)]) * m_dt;
	}

	return rhsU;
}

std::vector<double> MACSolver2D::UpdateFV(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	std::vector<double> ls, std::vector<double>u, std::vector<double> v) {

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
	double theta = 0.0, rhoEff = 0.0;
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		
		if (lsS >= 0 && lsM >= 0) {
			rhoEff = kRhoI;
		}
		else if (lsS < 0 && lsM < 0) {
			rhoEff = kRhoO;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
		}
		
		rhsV[idx(i, j)] = (-cV[idx(i, j)] + vV[idx(i, j)] / rhoEff + gV[idx(i, j)]) * m_dt;
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

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
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

	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
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
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	// subcell
	double theta = 0.0;
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0, muV_X_W = 0.0, muV_X_E = 0.0;
	double vW = 0.0, vE = 0.0, vM = 0.0;
	double visX = 0.0, visY = 0.0;

	// jump condition
	double JW = 0.0, JE = 0.0, JS = 0.0, JN = 0.0, JM = 0.0;
	// effective Jump condition, effective u(uEff), and effective mu (muEff)
	// J is a jump condition and defined at P grid
	double JEff = 0.0, JO = 0.0, uEff = 0.0, vEff = 0.0, muEff = 0.0;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
		visX = 0.0, visY = 0.0;
		muU_X_W = 0.0; muU_X_E = 0.0; muU_Y_S = 0.0; muU_Y_N = 0.0;
		lsW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsN = 0.5 * (ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		vW = 0.25 * (v[idx(i - 2, j)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)], v[idx(i - 2, j + 1)]);
		vM = 0.25 * (v[idx(i - 1, j)] + v[idx(i, j)] + v[idx(i, j + 1)], v[idx(i - 1, j + 1)]);
		vE = 0.25 * (v[idx(i, j)] + v[idx(i + 1, j)] + v[idx(i + 1, j + 1)], v[idx(i, j + 1)]);
		
		if (lsW >= 0 && lsM >= 0) {
			muU_X_W = kMuI * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muV_X_W = kMuI * (vM - vW) / kDx;
		}
		else if (lsW < 0 && lsM < 0) {
			muU_X_W = kMuO * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muV_X_W = kMuO * (vM - vW) / kDx;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between v[i,j] and v[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = (kMuO - kMuI) * m_J11[idx(i - 1, j)];
			JM = (kMuO - kMuI) * m_J11[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i - 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_X_W = muEff * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx
				- muEff * JEff * theta / kMuI;
			
			JW = (kMuO - kMuI) * m_J21[idx(i - 1, j)];
			JM = (kMuO - kMuI) * m_J21[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			vEff = (kMuO * vM * theta + kMuI * vW * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muV_X_W = muEff * (vM - vW) / kDx
				- muEff * JEff * theta / kMuI;
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between v[i,j] and v[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = (kMuI - kMuO) * m_J11[idx(i - 1, j)];
			JM = (kMuI - kMuO) * m_J11[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			uEff = (kMuI * u[idx(i, j)] * theta + kMuO * u[idx(i - 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1- theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_X_W = muEff * (u[idx(i, j)] - u[idx(i - 1,j)]) / kDx
				+ muEff * JEff * theta / kMuO;

			JW = (kMuI - kMuO) * m_J21[idx(i - 1, j)];
			JM = (kMuI - kMuO) * m_J21[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			vEff = (kMuI * vM * theta + kMuO * vW * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_W = muEff * (vM - vW) / kDx
				+ muEff * JEff * theta / kMuO;
		}
		

		if (lsM >= 0 && lsE >= 0) {
			muU_X_E = kMuI * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
			muV_X_E = kMuI * (vE - vM) / kDx;
		}
		else if (lsM < 0 && lsE < 0) {
			muU_X_E = kMuO * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
			muV_X_E = kMuO * (vE - vM) / kDx;
		}
		else if (lsM >= 0 && lsE < 0) {
			// interface lies between u[i,j] and u[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
			JM = (kMuO - kMuI) * m_J11[idx(i, j)];
			JE = (kMuO - kMuI) * m_J11[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			uEff = (kMuI * u[idx(i, j)] * theta + kMuO * u[idx(i + 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_X_E = muEff * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- muEff * JEff * theta / kMuO;

			JM = (kMuO - kMuI) * m_J21[idx(i, j)];
			JE = (kMuO - kMuI) * m_J21[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			uEff = (kMuI * vM * theta + kMuO * vE * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_X_E = muEff * (vE - vM) / kDx
				- muEff * JEff * theta / kMuO;
		}
		else if (lsM < 0 && lsE >= 0) {
			// interface lies between v[i,j] and v[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===   outside(-)    === |(interface)| === inside(+) === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsE)|
			JM = (kMuI - kMuO) * m_J11[idx(i, j)];
			JE = (kMuI - kMuO) * m_J11[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i + 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_X_E = muEff * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- muEff * JEff * theta / kMuI;

			JM = (kMuI - kMuO) * m_J21[idx(i, j)];
			JE = (kMuI - kMuO) * m_J21[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			uEff = (kMuO * vM * theta + kMuI * vE * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_X_E = muEff * (vE - vM) / kDx
				- muEff * JEff * theta / kMuI;
		}
		
		if (lsS >= 0 && lsM >= 0) {
			muU_Y_S = kMuI * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
		}
		else if (lsS < 0 && lsM < 0) {
			muU_Y_S = kMuO * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = (kMuO - kMuI) * m_J12[idx(i, j - 1)];
			JM = (kMuO - kMuI) * m_J12[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i, j - 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				- muEff * JEff * theta / kMuI;
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = (kMuI - kMuO) * m_J12[idx(i, j - 1)];
			JM = (kMuI - kMuO) * m_J12[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			uEff = (kMuI * u[idx(i, j)] * theta + kMuO * u[idx(i, j - 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				+ muEff * JEff * theta / kMuO;
		}
		
		if (lsM >= 0 && lsN >= 0) {
			muU_Y_N = kMuI * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lsM < 0 && lsN < 0) {
			muU_Y_N = kMuO * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lsM >= 0 && lsM < 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsN)|
			JM = (kMuO - kMuI) * m_J12[idx(i, j)];
			JN = (kMuO - kMuI) * m_J12[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			uEff = (kMuI * u[idx(i, j)] * theta + kMuO * u[idx(i, j + 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_N = muEff * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- muEff * JEff * theta / kMuO;
		}
		else if (lsM < 0 && lsN >= 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    outside(-)   === |(interface)| ===   inside(+) === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsN)|
			JM = (kMuI - kMuO) * m_J12[idx(i, j)];
			JN = (kMuI - kMuO) * m_J12[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i, j + 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_Y_N = muEff * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- muEff * JEff * theta / kMuI;
		}
		
		visX = (muU_X_E - muU_X_W) / kDx;
		visY = (muU_Y_N - muU_Y_S) / kDy + (muV_X_E - muV_X_W) / kDy;
		
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
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	// jump condition
	double JW = 0.0, JE = 0.0, JS = 0.0, JN = 0.0, JM = 0.0;
	double theta = 0.0;
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
	double uS = 0.0, uM = 0.0, uN = 0.0;
	double visX = 0.0, visY = 0.0;
	// effective Jump condition, effective v (vEff), and effective mu (muEff)
	double JEff = 0.0, JO = 0.0, vEff = 0.0, uEff = 0.0, muEff = 0.0;
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		muV_X_W = 0.0; muV_X_E = 0.0; muV_Y_S = 0.0; muV_Y_N = 0.0;
		visX = 0.0, visY = 0.0;
		lsW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
		lsN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		uS = 0.25 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)] + u[idx(i + 1, j - 2)] + u[idx(i, j - 1)]);
		uM = 0.25 * (u[idx(i, j)] + u[idx(i + 1, j)] + u[idx(i + 1, j - 1)] + u[idx(i, j - 1)]);
		uN = 0.25 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)] + u[idx(i + 1, j)] + u[idx(i, j)]);

		if (lsW >= 0 && lsM >= 0) {
			muV_X_W = kMuI * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
		}
		else if (lsW < 0 && lsM < 0) {
			muV_X_W = kMuO * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between v[i,j] and v[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| === inside(+) === |(interface)| ===   outside(-)    === |(lsM)|
			// |(lsW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = (kMuO - kMuI) * m_J21[idx(i - 1, j)];
			JM = (kMuO - kMuI) * m_J21[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			vEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i - 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_X_W = muEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				- muEff * JEff * theta / kMuI;
		}
		else if (lsW < 0 && lsM >= 0) {
			// interface lies between v[i,j] and v[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = (kMuI - kMuO) * m_J21[idx(i - 1, j)];
			JM = (kMuI - kMuO) * m_J21[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JW;
			vEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i - 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_W = muEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				+ muEff * JEff * theta / kMuO;
		}
		
		if (lsM >= 0 && lsE >= 0) {
			muV_X_E = kMuI * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lsM < 0 && lsE < 0) {
			muV_X_E = kMuO * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lsM >= 0 && lsE < 0) {
			// interface lies between v[i,j] and v[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
			JM = (kMuO - kMuI) * m_J21[idx(i, j)];
			JE = (kMuO - kMuI) * m_J21[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			vEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i + 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_X_E = muEff * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- muEff * JEff * theta / kMuO;
		}
		else if (lsM < 0 && lsE >= 0) {
			// interface lies between v[i,j] and v[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    outside(-)   === |(interface)| === inside(+)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
			JM = (kMuI - kMuO) * m_J21[idx(i, j)];
			JE = (kMuI - kMuO) * m_J21[idx(i + 1, j)];
			JEff = theta * JM + (1.0 - theta) * JE;
			vEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i + 1, j)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muV_X_E = muEff * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- muEff * JEff * theta / kMuI;
		}
		
		if (lsS >= 0 && lsM >= 0) {
			muV_Y_S = kMuI * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_Y_S = kMuI * (uM - uS) / kDy;
		}
		else if (lsS < 0 && lsM < 0) {
			muV_Y_S = kMuO * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_Y_S = kMuO * (uM - uS) / kDy;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between v[i,j] and v[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = (kMuO - kMuI) * m_J22[idx(i, j - 1)];
			JM = (kMuO - kMuI) * m_J22[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			vEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i, j - 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_Y_S = muEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				+ muEff * JEff * theta / kMuO;

			JS = (kMuO - kMuI) * m_J12[idx(i, j - 1)];
			JM = (kMuO - kMuI) * m_J12[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			uEff = (kMuI * uM * theta + kMuO * uS * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (uM - uS) / kDy
				+ muEff * JEff * theta / kMuO;
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between v[i,j] and v[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = (kMuI - kMuO) * m_J22[idx(i, j - 1)];
			JM = (kMuI - kMuO) * m_J22[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			vEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i, j - 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muV_Y_S = muEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				- muEff * JEff * theta / kMuI;

			JS = (kMuI - kMuO) * m_J12[idx(i, j - 1)];
			JM = (kMuI - kMuO) * m_J12[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			vEff = (kMuO * uM * theta + kMuI * uS * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuO * (1.0 - theta));
			muU_Y_S = muEff * (uM - uS) / kDy
				- muEff * JEff * theta / kMuI;
		}
		
		if (lsM >= 0 && lsN >= 0) {
			muV_Y_N = kMuI * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
			muU_Y_N = kMuI * (uN - uM) / kDy;
		}
		else if (lsM < 0 && lsN < 0) {
			muV_Y_N = kMuO * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
			muU_Y_N = kMuO * (uN - uM) / kDy;
		}
		else if (lsM >= 0 && lsN < 0) {
			// interface lies between v[i,j] and v[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsN)|
			JM = (kMuO - kMuI) * m_J22[idx(i, j)];
			JN = (kMuO - kMuI) * m_J22[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			vEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i, j + 1)] * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muV_Y_N = muEff * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- muEff * JEff * theta / kMuO;

			JM = (kMuO - kMuI) * m_J12[idx(i, j)];
			JN = (kMuO - kMuI) * m_J12[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			vEff = (kMuI * uM * theta + kMuO * uN * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			muU_Y_N = muEff * (uN - uM) / kDy
				- muEff * JEff * theta / kMuO;
		}
		else if (lsM < 0 && lsN >= 0) {
			// interface lies between v[i,j] and v[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    outside(-)   === |(interface)| === inside(+) === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			JM = (kMuI - kMuO) * m_J22[idx(i, j)];
			JN = (kMuI - kMuO) * m_J22[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			vEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i, j + 1)] * (1.0 - theta)
					- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muV_Y_N = muEff * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- muEff * JEff * theta / kMuI;

			JM = (kMuI - kMuO) * m_J21[idx(i, j)];
			JN = (kMuI - kMuO) * m_J21[idx(i, j + 1)];
			JEff = theta * JM + (1.0 - theta) * JN;
			vEff = (kMuO * uM * theta + kMuI * uN * (1.0 - theta)
				- JEff * theta * (1.0 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1.0 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			muU_Y_N = muEff * (uN - uM) / kDy
				- muEff * JEff * theta / kMuI;
		}
		
		visX = (muV_X_E - muV_X_W) / kDx;
		visY = (muV_Y_N - muV_Y_S) / kDy + (muU_Y_N - muU_Y_S) / kDy;

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

	if (kFr == 0) {
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
				gU[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
			gU[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gU;
}

std::vector<double> MACSolver2D::AddGravityFV() {
	std::vector<double> gV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if (kFr == 0) {
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			gV[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			gV[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gV;
}

std::vector<double> MACSolver2D::GetUhat(const std::vector<double>& u, const std::vector<double>& rhsu) {
	std::vector<double> uhat((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Update rhs
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
		uhat[idx(i, j)] = u[idx(i, j)] + m_dt * rhsu[idx(i, j)];
		if (std::isnan(uhat[idx(i, j)]) || std::isinf(uhat[idx(i, j)])) {
			std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << uhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return uhat;
}

std::vector<double> MACSolver2D::GetVhat(const std::vector<double>& v, const std::vector<double>& rhsv) {
	std::vector<double> vhat((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Update rhs
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		vhat[idx(i, j)] = v[idx(i, j)] + m_dt * rhsv[idx(i, j)];
		if (std::isnan(vhat[idx(i, j)]) || std::isinf(vhat[idx(i, j)])) {
			std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << vhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return vhat;
}

int MACSolver2D::SetPoissonSolver(POISSONTYPE type) {
	m_PoissonSolverType = type;
	if (!m_Poisson)
		m_Poisson = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);

	return 0;
}

int MACSolver2D::SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
	const std::vector<double>& lsB, const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v) {
	if (!m_Poisson) {
		perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
	}
	std::vector<double> rhs((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if (m_PoissonSolverType == POISSONTYPE::MKL) {
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
		
		m_Poisson->MKL_2FUniform_2D(ps, rhs,
			kLenX, kLenY, kDx, kDy, m_BC);
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;

	}
	else if (m_PoissonSolverType == POISSONTYPE::ICPCG) {
		// based on preconditioned Conjugate Gradient Method and Fedkiw's paper

	}
	else if (m_PoissonSolverType == POISSONTYPE::CG) {
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
		// level set
		
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
		// 1 / rho is a coef of poisson equation
		// iRhoW, iRhoE : at U grid
		// iRhoS, iRhoN : at V grid
		double iRhoW = 0.0, iRhoE = 0.0, iRhoS = 0.0, iRhoN = 0.0;
		// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
		double a = 0.0, b = 0.0;
		double theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, bEff = 0.0;
		
		UpdateKappa(ls);
		// A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
		std::vector<double> AVals, DiagVals, MVals;
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

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = lsB[idx(i - 1, j)];
			lsE = lsB[idx(i + 1, j)];
			lsM = lsB[idx(i, j)];
			lsS = lsB[idx(i, j - 1)];
			lsN = lsB[idx(i, j + 1)];

			dudX[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
			dudY[idx(i, j)] = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
				- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);

			dvdX[idx(i, j)] = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
				- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
			dvdY[idx(i, j)] = (v[idx(i, j + 1)] - v[idx(i, j)]) / (2.0 * kDy);

			dldX[idx(i, j)] = (lsE - lsW) / (2.0 * kDx);
			dldY[idx(i, j)] = (lsN - lsS) / (2.0 * kDy);
			if (std::fabs(dldX[idx(i, j)]) < 1.0e-40) {
				dldX[idx(i, j)] = (lsE - lsM) / kDx;
				if (dldX[idx(i, j)] == 0.0)
					dldX[idx(i, j)] = (lsM - lsW) / kDx;
			}
			if (std::fabs(dldY[idx(i, j)]) < 1.0e-40) {
				dldY[idx(i, j)] = (lsN - lsM) / kDy;
				if (dldY[idx(i, j)] == 0.0)
					dldY[idx(i, j)] = (lsM - lsS) / kDy;
			}
		}

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = lsB[idx(i - 1, j)];
			lsE = lsB[idx(i + 1, j)];
			lsM = lsB[idx(i, j)];
			lsS = lsB[idx(i, j - 1)];
			lsN = lsB[idx(i, j + 1)];

			FW = 0.0;
			FE = 0.0;
			FS = 0.0;
			FN = 0.0;

			// normal vector = (\nabla \phi) / |\nabla \phi|
			nXW = dldX[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + 1.0e-100);
			nYW = dldY[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + 1.0e-100); 
			nXE = dldX[idx(i + 1, j)] / (std::sqrt(std::pow(dldX[idx(i + 1, j)], 2.0) + std::pow(dldY[idx(i + 1, j)], 2.0)) + 1.0e-100);
			nYE = dldY[idx(i + 1, j)] / (std::sqrt(std::pow(dldX[idx(i + 1, j)], 2.0) + std::pow(dldY[idx(i + 1, j)], 2.0)) + 1.0e-100);
			nXS = dldX[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + 1.0e-100);
			nYS = dldY[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + 1.0e-100);
			nXN = dldX[idx(i, j + 1)] / (std::sqrt(std::pow(dldX[idx(i, j + 1)], 2.0) + std::pow(dldY[idx(i, j + 1)], 2.0)) + 1.0e-100);
			nYN = dldY[idx(i, j + 1)] / (std::sqrt(std::pow(dldX[idx(i, j + 1)], 2.0) + std::pow(dldY[idx(i, j + 1)], 2.0)) + 1.0e-100);
			nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + 1.0e-100);
			nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + 1.0e-100);
			
			/*
			// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
			// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
			// (P_M - P_W)/kDx appears
			// if P_M == P_+, P_W == P_-, a^+ means all terms related to P_+, P_W changed to P_M related terms
			// P_W = P_M -(2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa)
			*/
			// FW
			if (lsW * lsM >= 0) {
				// one fluid, x direction
				FW = 0.0;
			}
			else if (lsW >= 0 && lsM < 0) {
				// interface lies between ls[i,j] and ls[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				// b always zero when solving level set (dealing with surface tension)
				b = 0.0;
				aW = 2 * m_dt * (kMuO - kMuI) 
					* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
					+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
				aM = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
				iRhoW = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsW) + std::fabs(lsM))
					/ (1.0 / kRhoO * std::fabs(lsW) + 1.0 / kRhoI * std::fabs(lsM));
				FW = iRhoW * aEff / (kDx * kDx) - iRhoW * b * theta / ((1.0 / kRhoI) * kDx);
			}
			else if (lsW < 0 && lsM >= 0) {
				// interface lies between ls[i,j] and ls[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				b = 0.0;
				aW = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
					+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
				aM = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * std::fabs(lsW) + aW * std::fabs(lsM)) / (std::fabs(lsW) + std::fabs(lsM));
				iRhoW = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsW) + std::fabs(lsM))
					/ (1.0 / kRhoI * std::fabs(lsW) + 1.0 / kRhoO * std::fabs(lsM));
				FW = -iRhoW * aEff / (kDx * kDx) + iRhoW * b * theta / ((1.0 / kRhoO) * kDx);
			}

			// FE
			// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
			if (lsM * lsE >= 0) {
				// one fluid, x direction
				FE = 0.0;
			}
			else if (lsM >= 0 && lsE < 0) {
				// interface lies between ls[i,j] and ls[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
				b = 0.0;
				aM = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aE = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
					+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
				aEff = (aM * std::fabs(lsE) + aE * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsE));
				1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsE))
					/ (1.0 / kRhoI * std::fabs(lsE) + 1.0 / kRhoO * std::fabs(lsM));
				FE = -iRhoE * aEff / (kDx * kDx) - iRhoE * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM < 0 && lsE >= 0) {
				// interface lies between ls[i,j] and ls[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
				b = 0.0;
				aM = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aE = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
					+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
				aEff = (aM * std::fabs(lsE) + aE * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsE));
				iRhoE = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsE))
					/ (1.0 / kRhoO * std::fabs(lsE) + 1.0 / kRhoI * std::fabs(lsM));
				FE = iRhoE * aEff / (kDx * kDx) + iRhoE * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FS
			if (lsS * lsM >= 0) {
				// one fluid, y direction
				FS = 0.0;
			}
			else if (lsS >= 0 && lsM < 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
					+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
				aM = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsS));
				iRhoS = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsS) + std::fabs(lsM))
					/ (1.0 / kRhoO * std::fabs(lsS) + 1.0 / kRhoI * std::fabs(lsM));
				FS = iRhoS * aEff / (kDy * kDy) - iRhoS * b * theta / ((1.0 / kRhoI) * kDy);
			}
			else if (lsS < 0 && lsM >= 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nXS) * nXS
					+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
				aM = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * std::fabs(lsS) + aS * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsS));
				iRhoS = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsS) + std::fabs(lsM))
					/ (1.0 / kRhoI * std::fabs(lsS) + 1.0 / kRhoO * std::fabs(lsM));
				FS = -iRhoS * aEff / (kDy * kDy) + iRhoS * b * theta / ((1.0 / kRhoO) * kDy);
			}
			

			// FN
			if (lsM * lsN >= 0) {
				// one fluid, y direction
				FN = 0.0;
			}
			else if (lsM >= 0 && lsN < 0) {
				// interface lies between ls[i,j] and ls[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aN = 2 * m_dt * (kMuO - kMuI)
					* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nXN) * nXN
					+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
				aEff = (aM * std::fabs(lsN) + aN * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsN));
				iRhoN = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsN))
					/ (1.0 / kRhoI * std::fabs(lsN) + 1.0 / kRhoO * std::fabs(lsM));
				FN = -iRhoN * aEff / (kDy * kDy) - iRhoN * b * theta / ((1.0 / kRhoO) * kDy);
			}
			else if (lsM < 0 && lsN >= 0) {
				// interface lies between ls[i,j] and ls[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
					+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aN = 2 * m_dt * (kMuI - kMuO)
					* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nXN) * nXN
					+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
				aEff = (aM * std::fabs(lsN) + aN * std::fabs(lsM)) / (std::fabs(lsM) + std::fabs(lsN));
				iRhoN = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsN))
					/ (1.0 / kRhoO * std::fabs(lsN) + 1.0 / kRhoI * std::fabs(lsM));
				FN = iRhoN * aEff / (kDy * kDy) + iRhoN * b * theta / ((1.0 / kRhoI) * kDy);
			}

			// initially set variable coef. considered RHS
			rhs[idx(i, j)] = FW + FE + FS + FN;
			if (rhs[idx(i, j)] != 0.0)
				std::cout << i << " " << j << " " << FW << " " << FE << " " << FS << " " << FN << std::endl;

			assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
			if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
				std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " << rhs[idx(i, j)] << std::endl;
				exit(1);
			}
		}

		// Original value of RHS
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			rhs[idx(i, j)] += -div[idx(i, j)];

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = lsB[idx(i - 1, j)];
			lsE = lsB[idx(i + 1, j)];
			lsM = lsB[idx(i, j)];
			lsS = lsB[idx(i, j - 1)];
			lsN = lsB[idx(i, j + 1)];

			if (lsW * lsM >= 0) {
				// one fluid, x direction
				if (lsW <= 0)
					iRhoW = 1.0 / kRhoO;
				else
					iRhoW = 1.0 / kRhoI;
			}
			else if (lsW >= 0 && lsM < 0) {
				// interface lies between ls[i,j] and ls[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				// b always zero when solving level set (dealing with surface tension)
				iRhoW = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsW) + std::fabs(lsM))
					/ (1.0 / kRhoO * std::fabs(lsW) + 1.0 / kRhoI * std::fabs(lsM));
			}
			else if (lsW < 0 && lsM >= 0) {
				// interface lies between ls[i,j] and ls[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				iRhoW = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsW) + std::fabs(lsM))
					/ (1.0 / kRhoI * std::fabs(lsW) + 1.0 / kRhoO * std::fabs(lsM));
			}

			if (lsM * lsE >= 0) {
				// one fluid, x direction
				if (lsE <= 0)
					iRhoE = 1.0 / kRhoO;
				else
					iRhoE = 1.0 / kRhoI;
			}
			else if (lsM >= 0 && lsE < 0) {
				// interface lies between ls[i,j] and ls[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
				iRhoE = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsE))
					/ (1.0 / kRhoI * std::fabs(lsE) + 1.0 / kRhoO * std::fabs(lsM));
			}
			else if (lsM < 0 && lsE >= 0) {
				// interface lies between ls[i,j] and ls[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
				iRhoE = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsE))
					/ (1.0 / kRhoO * std::fabs(lsE) + 1.0 / kRhoI * std::fabs(lsM));
			}

			if (lsS * lsM >= 0) {
				// one fluid, y direction
				if (lsS <= 0)
					iRhoS = 1.0 / kRhoO;
				else
					iRhoS = 1.0 / kRhoI;
			}
			else if (lsS >= 0 && lsM < 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				iRhoS = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsS) + std::fabs(lsM))
					/ (1.0 / kRhoO * std::fabs(lsS) + 1.0 / kRhoI * std::fabs(lsM));
			}
			else if (lsS < 0 && lsM >= 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				iRhoS = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsS) + std::fabs(lsM))
					/ (1.0 / kRhoI * std::fabs(lsS) + 1.0 / kRhoO * std::fabs(lsM));
			}

			if (lsM * lsN >= 0) {
				// one fluid, y direction
				if (lsN <= 0)
					iRhoN = 1.0 / kRhoO;
				else
					iRhoN = 1.0 / kRhoI;
			}
			else if (lsM >= 0 && lsN < 0) {
				// interface lies between ls[i,j] and ls[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				iRhoN = 1.0 / (kRhoI * kRhoO) * (std::fabs(lsM) + std::fabs(lsN))
					/ (1.0 / kRhoI * std::fabs(lsN) + 1.0 / kRhoO * std::fabs(lsM));
			}
			else if (lsM < 0 && lsN >= 0) {
				// interface lies between ls[i,j] and ls[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				iRhoN = 1.0 / (kRhoO * kRhoI) * (std::fabs(lsM) + std::fabs(lsN))
					/ (1.0 / kRhoO * std::fabs(lsN) + 1.0 / kRhoI * std::fabs(lsM));
			}

			// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
			// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
			// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
			AValsDic.clear();
			AColsDic.clear();
			tmpRowIdx = 0;
			// Add starting rowIdx
			ARowIdx.push_back(rowIdx);

			// Set default values, if a current pointer is in interior, it will not be changed.
			AValsDic["S"] = -iRhoS / (kDy * kDy);
			AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
			AValsDic["W"] = -iRhoW / (kDx * kDx);
			AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
			AValsDic["C"] = (iRhoW + iRhoE) / (kDx * kDx) + (iRhoS + iRhoN) / (kDy * kDy);
			AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
			AValsDic["E"] = -iRhoE / (kDx * kDx);
			AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
			AValsDic["N"] = -iRhoN / (kDy * kDy);
			AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;

			// West boundary
			if (i == kNumBCGrid && m_BC->m_BC_PW == BC::NEUMANN) {
				AValsDic["W"] = 0.0;
				AColsDic["W"] = -1;
				AValsDic["C"] -= iRhoW / (kDx * kDx);

				if (j == kNumBCGrid && m_BC->m_BC_PS == BC::NEUMANN) {
					AValsDic["S"] = 0.0;
					AColsDic["S"] = -1;
					AValsDic["C"] -= iRhoS / (kDy * kDy);
				}
				else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::NEUMANN) {
					AValsDic["N"] = 0.0;
					AColsDic["N"] = -1;
					AValsDic["C"] -= iRhoN / (kDy * kDy);
				}
			}
			// East boundary
			else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::NEUMANN) {
				AValsDic["E"] = 0.0;
				AColsDic["E"] = -1;
				AValsDic["C"] -= iRhoE / (kDx * kDx);

				if (j == kNumBCGrid && m_BC->m_BC_PS == BC::NEUMANN) {
					AValsDic["S"] = 0.0;
					AColsDic["S"] = -1;
					AValsDic["C"] -= iRhoS / (kDy * kDy);
				}
				else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::NEUMANN) {
					AValsDic["N"] = 0.0;
					AColsDic["N"] = -1;
					AValsDic["C"] -= iRhoN / (kDy * kDy);
				}
			}
			// South boundary
			else if (j == kNumBCGrid && m_BC->m_BC_PS == BC::NEUMANN) {
				AValsDic["S"] = 0.0;
				AColsDic["S"] = -1;
				AValsDic["C"] -= iRhoS / (kDy * kDy);

				if (i == kNumBCGrid && m_BC->m_BC_PW == BC::NEUMANN) {
					AValsDic["W"] = 0.0;
					AColsDic["W"] = -1;
					AValsDic["C"] -= iRhoW / (kDx * kDx);
				}
				else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::NEUMANN) {
					AValsDic["E"] = 0.0;
					AColsDic["E"] = -1;
					AValsDic["C"] -= iRhoE / (kDx * kDx);
				}
			}
			// North boundary
			else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC::NEUMANN) {
				AValsDic["N"] = 0.0;
				AColsDic["N"] = -1;
				AValsDic["C"] -= iRhoN / (kDy * kDy);

				if (i == kNumBCGrid && m_BC->m_BC_PW == BC::NEUMANN) {
					AValsDic["W"] = 0.0;
					AColsDic["W"] = -1;
					AValsDic["C"] -= iRhoW / (kDx * kDx);
				}
				else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC::NEUMANN) {
					AValsDic["E"] = 0.0;
					AColsDic["E"] = -1;
					AValsDic["C"] -= iRhoE / (kDx * kDx);
				}
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

			DiagVals.push_back(AValsDic["C"]);
			rowIdx += tmpRowIdx;
			
			assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
			if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
				std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " << rhs[idx(i, j)] << std::endl;
				exit(1);
			}
		}

		ARowIdx.push_back(rowIdx);
		m_Poisson->CG_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC);
	}
	else if (m_PoissonSolverType == POISSONTYPE::GS) {
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			rhs[idx(i, j)] = div[idx(i, j)];

		m_Poisson->GS_2FUniform_2D(ps, rhs, kDx, kDy, m_BC);

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}

	return 0;
}

std::vector<double> MACSolver2D::GetDivergence(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> div((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		div[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
						+ (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
	
	return div;
}

int MACSolver2D::UpdateVel(std::vector<double>& u, std::vector<double>& v,
	const std::vector<double>& us, const std::vector<double>& vs, const std::vector<double>& ps, const std::vector<double>& ls) {

	// velocity update after solving poisson equation
	// ps = p * dt
	double lsW = 0.0, lsM = 0.0, lsS = 0.0, rhoEff = 0.0, theta = 0.0;
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];

		if (lsW >= 0 && lsM >= 0) {
			rhoEff = kRhoI;
		}
		else if (lsW < 0 && lsM < 0) {
			rhoEff = kRhoO;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between ls[i - 1,j] and ls[i,j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
		}
		else if (lsW < 0 && lsM >= 0) {
			// interface lies between ls[i - 1, j] and ls[i,j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
		}
		u[idx(i, j)] = us[idx(i, j)] - (ps[idx(i, j)] - ps[idx(i - 1, j)]) / (kDx * rhoEff);
	}

	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];

		if (lsS >= 0 && lsM >= 0) {
			rhoEff = kRhoI;
		}
		else if (lsS < 0 && lsM < 0) {
			rhoEff = kRhoO;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
		}

		v[idx(i, j)] = vs[idx(i, j)] - (ps[idx(i, j)] - ps[idx(i, j - 1)]) / (kDy * rhoEff);
	}

	return 0;
}

double MACSolver2D::UpdateDt(const std::vector<double>& u, const std::vector<double>& v) {
	double uAMax = 0.0;
	double vAMax = 0.0;
	double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0;
	double dt = std::numeric_limits<double>::max();

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
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
	return i + (kNx + 2 * kNumBCGrid) * j;
}

int MACSolver2D::SetPLTType(PLTTYPE type) {
	m_PLTType = type;

	return 0;
}

int MACSolver2D::OutRes(int iter, double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& phi, const std::vector<double>& ls) {
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
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			resU[idx(i, j)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			resV[idx(i, j)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;

		outF.open(fname_vel.c_str(), std::ios::app);
		
		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy
			<< std::string(", DATAPACKING=POINT")
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
				<< static_cast<double>(phi[idx(i, j)]) << std::endl;

		outF.close();

		outF.open(fname_div.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy
			<< std::string(", DATAPACKING=POINT")
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
				<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
				<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
				<< static_cast<double>(m_kappa[idx(i, j)]) << std::string(",")
				<< static_cast<double>(phi[idx(i, j)]) << std::endl;

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
			resPhi(kNx * kNy, 0.0);

		std::vector<double> resDiv = GetDivergence(u, v);
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			resX[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseX + (i + 0.5 - kNumBCGrid) * kDx;
			resY[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseY + (j + 0.5 - kNumBCGrid) * kDy;
			resU[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
			resV[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ls[idx(i, j)];
			resPhi[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = phi[idx(i, j)];
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
		stat = TECDAT142(&ARRSIZEVAL, resPhi.data(), &DIsDouble);

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
		stat = TECDAT142(&ARRSIZEVAL, resPhi.data(), &DIsDouble);
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
