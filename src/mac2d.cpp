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
	double duMX = 0.0, duMY = 0.0;
	double dvMX = 0.0, dvMY = 0.0;
	double dlMX = 0.0, dlMY = 0.0;
	// normal and tangent vector variable (n, t1, t2)
	// normal vector has X, Y component
	// first tangent vector has only Z component
	double nX = 0.0, nY = 0.0, t1Z = 0.0, t2X = 0.0, t2Y = 0.0;

	// jump condition originally defined as J11 = [\mu u_x], J12 = [\mu u_y], J21 = [\mu v_x], J22 = [\mu v_y], 
	// but skip [\mu], from Kang, Fedkiw, and Liu's work eq. (30).
	
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++)
	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++) {
		duMX = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
		duMY = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
			- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);
		dvMX = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
		dvMY = (v[idx(i, j + 1)] - v[idx(i, j)]) / (2.0 * kDy);
		dlMX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dlMY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		if (dlMX == 0.0) {
			dlMX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			if (dlMX == 0.0)
				dlMX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		}
		if (dlMY == 0.0) {
			dlMY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			if (dlMY == 0.0)
				dlMY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
		}

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nX = dlMX / (std::fabs(dlMX) + 1.0e-100);
		nY = dlMY / (std::fabs(dlMY) + 1.0e-100);
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

		// eq. (30) from Kang, Fedkiw, and Liu
		m_J11[idx(i, j)] = duMX * (t1Z * t1Z) + (nX * nX + nY * nY) * duMX * (nX * nX + nY * nY)
			- (t1Z * t1Z) * duMX * (nX * nX + nY * nY);
		m_J12[idx(i, j)] = duMY * (t1Z * t1Z) + (nX * nX + nY * nY) * duMY * (nX * nX + nY * nY)
			- (t1Z * t1Z) * dvMX * (nX * nX + nY * nY);
		m_J21[idx(i, j)] = dvMX * (t1Z * t1Z) + (nX * nX + nY * nY) * dvMX * (nX * nX + nY * nY)
			- (t1Z * t1Z) * duMY * (nX * nX + nY * nY);
		m_J22[idx(i, j)] = dvMY * (t1Z * t1Z) + (nX * nX + nY * nY) * dvMY * (nX * nX + nY * nY)
			- (t1Z * t1Z) * dvMY * (nX * nX + nY * nY);

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
	cU = this->AddConvectionFU(u, v, ls);

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
		rhsU[idx(i, j)] = (-cU[idx(i, j)] + vU[idx(i, j)] + gU[idx(i, j)]) * m_dt / rhoEff;
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
	cV = this->AddConvectionFV(u, v, ls);

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
		
		rhsV[idx(i, j)] = (-cV[idx(i, j)] + vV[idx(i, j)] + gV[idx(i, j)]) * m_dt / rhoEff;
	}
	
	return rhsV;
}

std::vector<double> MACSolver2D::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& ls) {
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

std::vector<double> MACSolver2D::AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& ls) {
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
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
	double rhoU_X_W = 0.0, rhoU_X_E = 0.0, rhoU_Y_S = 0.0, rhoU_Y_N = 0.0;
	double visX = 0.0, visY = 0.0;

	// jump condition
	double JW = 0.0, JE = 0.0, JS = 0.0, JN = 0.0, JM = 0.0;
	// effective Jump condition, effective u(uEff), and effective mu (muEff)
	// J is a jump condition and defined at P grid
	double JEff = 0.0, JO = 0.0, uEff = 0.0, muEff = 0.0, rhoEff = 0.0;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++) {
		visX = 0.0, visY = 0.0;
		muU_X_W = 0.0; muU_X_E = 0.0; muU_Y_S = 0.0; muU_Y_N = 0.0;
		lsW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsN = 0.5 * (ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		
		if (lsW >= 0 && lsM >= 0) {
			muU_X_W = kMuI * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
		}
		else if (lsW < 0 && lsM < 0) {
			muU_X_W = kMuO * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_W = muEff * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx
				- muEff * JEff * theta / (kMuI / kRhoI);
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_W = muEff * (u[idx(i, j)] - u[idx(i - 1,j)]) / kDx
				+ muEff * JEff * theta / (kMuO / kRhoO);
		}
		

		if (lsM >= 0 && lsE >= 0) {
			muU_X_E = kMuI * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		}
		else if (lsM < 0 && lsE < 0) {
			muU_X_E = kMuO * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			// muU_X_W = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muU_X_E = (muEff / rhoEff) * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- (muEff / rhoEff) * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			// muU_X_W = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muU_X_E = muEff * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- muEff * JEff * theta / (kMuI / kRhoI);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_Y_S = muEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				- muEff * JEff * theta / (kMuI / kRhoI);
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_Y_S = muEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				+ muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_Y_N = muEff * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_Y_N = muEff * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuI / kRhoI);
		}
		
		visX = (muU_X_E - muU_X_W) / kDx;
		visY = (muU_Y_N - muU_Y_S) / kDy;
		
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
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0;
	double rhoV_X_W = 0.0, rhoV_X_E = 0.0, rhoV_Y_S = 0.0, rhoV_Y_N = 0.0;
	double visX = 0.0, visY = 0.0;
	// effective Jump condition, effective v (vEff), and effective mu (muEff)
	double JEff = 0.0, JO = 0.0, vEff = 0.0, muEff = 0.0, rhoEff = 0.0;
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		muV_X_W = 0.0; muV_X_E = 0.0; muV_Y_S = 0.0; muV_Y_N = 0.0;
		visX = 0.0, visY = 0.0;
		lsW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
		lsN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);

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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muV_X_W = muEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				- muEff * JEff * theta / (kMuI / kRhoI);
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muV_X_W = muEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				+ muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muV_X_E = muEff * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muV_X_E = muEff * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- muEff * JEff * theta / (kMuI / kRhoI);
		}
		
		if (lsS >= 0 && lsM >= 0) {
			muV_Y_S = kMuI * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
		}
		else if (lsS < 0 && lsM < 0) {
			muV_Y_S = kMuO * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between v[i,j] and v[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = (kMuO - kMuI) * m_J22[idx(i, j - 1)];
			JM = (kMuO - kMuI) * m_J22[idx(i, j)];
			JEff = theta * JM + (1.0 - theta) * JS;
			vEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i, j - 1)] * (1.0 - theta) - JEff * theta * (1.0 - theta) * kDy)
				/ (kMuI * theta + kMuO * (1.0 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1.0 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muV_Y_S = muEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				+ muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muV_Y_S = muEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				- muEff * JEff * theta / (kMuI / kRhoI);
		}
		
		if (lsM >= 0 && lsN >= 0) {
			muV_Y_N = kMuI * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		}
		else if (lsM < 0 && lsN < 0) {
			muV_Y_N = kMuO * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
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
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muV_Y_N = muEff * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuO / kRhoO);
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
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muV_Y_N = muEff * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuI / kRhoI);
		}
		
		visX = (muV_X_E - muV_X_W) / kDx;
		visY = (muV_Y_N - muV_Y_S) / kDy;

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
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v) {
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
		// level set
		double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
		double FW = 0.0, FE = 0.0, FS = 0.0, FN = 0.0;
		double duWX = 0.0, duWY = 0.0, duEX = 0.0, duEY = 0.0, duSX = 0.0, duSY = 0.0, duNX = 0.0, duNY = 0.0, duMX = 0.0, duMY = 0.0;
		double dvWX = 0.0, dvWY = 0.0, dvEX = 0.0, dvEY = 0.0, dvSX = 0.0, dvSY = 0.0, dvNX = 0.0, dvNY = 0.0, dvMX = 0.0, dvMY = 0.0;
		double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aM = 0.0;
		double rhoW = 0.0, rhoE = 0.0, rhoS = 0.0, rhoN = 0.0; 
		// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
		double a = 0.0, b = 0.0;
		double theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, bEff = 0.0, rhoEff = 0.0;
		UpdateKappa(ls);

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = ls[idx(i - 1, j)];
			lsE = ls[idx(i + 1, j)];
			lsM = ls[idx(i, j)];
			lsS = ls[idx(i, j - 1)];
			lsN = ls[idx(i, j + 1)];

			duWX = (u[idx(i, j)] - u[idx(i - 2, j)]) / (2.0 * kDx);
			duWY = (u[idx(i - 1, j + 1)] - u[idx(i - 1, j - 1)]) / (2.0 * kDy);
			duEX = (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			duEY = (u[idx(i + 1, j + 1)] - u[idx(i + 1, j - 1)]) / (2.0 * kDy);
			duSX = (u[idx(i + 1, j - 1)] - u[idx(i - 1, j - 1)]) / (2.0 * kDx);
			duSY = (u[idx(i, j)] - u[idx(i, j - 2)]) / (2.0 * kDy);
			duNX = (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			duNY = (u[idx(i + 1, j + 1)] - u[idx(i + 1, j - 1)]) / (2.0 * kDy);
			duMX = (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			duMY = (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);
			
			dvWX = (v[idx(i, j)] - v[idx(i - 2, j)]) / (2.0 * kDx);
			dvWY = (v[idx(i - 1, j + 1)] - v[idx(i - 1, j - 1)]) / (2.0 * kDy);
			dvEX = (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			dvEY = (v[idx(i + 1, j + 1)] - v[idx(i + 1, j - 1)]) / (2.0 * kDy);
			dvSX = (v[idx(i + 1, j - 1)] - v[idx(i - 1, j - 1)]) / (2.0 * kDx);
			dvSY = (v[idx(i, j)] - v[idx(i, j - 2)]) / (2.0 * kDy);
			dvNX = (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			dvNY = (v[idx(i + 1, j + 1)] - v[idx(i + 1, j - 1)]) / (2.0 * kDy);
			dvMX = (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			dvMY = (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);

			FW = 0.0;
			FE = 0.0;
			FS = 0.0;
			FN = 0.0;  
			
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
				// interface lies between u[i,j] and u[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				// b always zero when solving level set (dealing with surface tension)
				b = 0.0;
				aW = (2 * m_dt * (kMuO - kMuI)
					* (duWX + duWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)]);
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoW = kRhoO * kRhoI * (fabs(lsW) + fabs(lsM)) / (kRhoI * fabs(lsW) + kRhoO * fabs(lsM));
				FW = (1.0 / rhoW) * aEff / (kDx * kDx) - (1.0 / rhoW) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsW < 0 && lsM >= 0) {
				// interface lies between u[i,j] and u[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				b = 0.0;
				aW = (2 * m_dt * (kMuI - kMuO)
					* (duWX + duWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)]);
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoW = kRhoI * kRhoO * (fabs(lsW) + fabs(lsM)) / (kRhoO * fabs(lsW) + kRhoI * fabs(lsM));
				FW = -(1.0 / rhoW) * aEff / (kDx * kDx) + (1.0 / rhoW) * b * theta / ((1.0 / kRhoI) * kDx);
			}
			
			// FE
			// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
			if (lsM * lsE >= 0) {
				// one fluid, x direction
				FE = 0.0;
			}
			else if (lsM >= 0 && lsE < 0) {
				// interface lies between u[i,j] and u[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
				b = 0.0;
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuO - kMuI)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)]);
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoE = kRhoI * kRhoO * (fabs(lsE) + fabs(lsM)) / (kRhoI * fabs(lsE) + kRhoO * fabs(lsM));
				FE = (1.0 / rhoE) * aEff / (kDx * kDx) - (1.0 / rhoE) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM < 0 && lsE >= 0) {
				// interface lies between u[i,j] and u[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
				b = 0.0;
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuI - kMuO)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)]);
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoE = kRhoO * kRhoI * (fabs(lsE) + fabs(lsM)) / (kRhoO * fabs(lsE) + kRhoI * fabs(lsM));
				FE = (1.0 / rhoE) * aEff / (kDx * kDx) - (1.0 / rhoE) * b * theta / ((1.0 / kRhoI) * kDx);
			}
			
			// FS
			if (lsS * lsM >= 0) {
				// one fluid, y direction
				FS = 0.0;
			}
			else if (lsM >= 0 && lsS < 0) {
				// interface lies between u[i,j] and u[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = (2 * m_dt * (kMuI - kMuO)
					* (duSX + duSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)]);
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoS = kRhoI * kRhoO * (fabs(lsS) + fabs(lsM)) / (kRhoI * fabs(lsS) + kRhoO * fabs(lsM));
				FS = (1.0 / rhoS) * aEff / (kDx * kDx) - (1.0 / rhoS) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM <= 0 && lsS > 0) {
				// interface lies between u[i,j] and u[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = (2 * m_dt * (kMuO - kMuI)
					* (duSX + duSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)]);
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoS = kRhoO * kRhoI * (fabs(lsS) + fabs(lsM)) / (kRhoO * fabs(lsS) + kRhoI * fabs(lsM));
				FS = (1.0 / rhoS) * aEff / (kDx * kDx) - (1.0 / rhoS) * b * theta / ((1.0 / kRhoI) * kDx);
			}
			
			// FN
			if (lsM * lsN >= 0) {
				// one fluid, y direction
				FN = 0.0;
			}
			else if (lsM >= 0 && lsN < 0) {
				// interface lies between u[i,j] and u[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuO - kMuI)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)]);
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoN = kRhoI * kRhoO * (fabs(lsN) + fabs(lsM)) / (kRhoI * fabs(lsN) + kRhoO * fabs(lsM));
				FN = (1.0 / rhoN) * aEff / (kDx * kDx) - (1.0 / rhoN) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM <= 0 && lsN > 0) {
				// interface lies between u[i,j] and u[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuI - kMuO)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)]);
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoN = kRhoO * kRhoI * (fabs(lsN) + fabs(lsM)) / (kRhoO * fabs(lsN) + kRhoI * fabs(lsM));
				FN = (1.0 / rhoN) * aEff / (kDx * kDx) - (1.0 / rhoN) * b * theta / ((1.0 / kRhoI) * kDx);
			}

			rhs[idx(i, j)] = FW + FE + FS + FN;
		}

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			rhs[idx(i, j)] += -div[idx(i, j)];
		
		m_Poisson->MKL_2FUniform_2D(ps, rhs,
			kLenX, kLenY, kDx, kDy, m_BC);
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;

	}
	else if (m_PoissonSolverType == POISSONTYPE::CG) {
		// based on preconditioned Conjugate Gradient Method and Fedkiw's paper

	}
	else if (m_PoissonSolverType == POISSONTYPE::ICPCG) {
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
		double nX = 0.0, nY = 0.0, t1X = 0.0, t1Y = 0.0, t2X = 0.0, t2Y = 0.0;
		double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aM = 0.0;
		double rhoW = 0.0, rhoE = 0.0, rhoS = 0.0, rhoN = 0.0;
		// 1 / rho is a coef of poisson equation
		double iRhoW = 0.0, iRhoE = 0.0, iRhoS = 0.0, iRhoN = 0.0;
		// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
		double a = 0.0, b = 0.0;
		double theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, bEff = 0.0, rhoEff = 0.0;
		
		UpdateKappa(ls);
		// A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
		std::vector<double> AVals, DiagVals, MVals;
		std::vector<MKL_INT> ACols, ARowIdx, MCols, MRowIdx;
		idx(2, 3);
		MKL_INT rowIdx = 0, tmpRowIdx = 0, colIdx = 0;
		MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
		MKL_INT size = kNx * kNy;
		
		// stored coef for A matrix, Dictionary but it is ordered
		std::map<std::string, double> AValsDic;
		std::map<std::string, MKL_INT> AColsDic;
		
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = ls[idx(i - 1, j)];
			lsE = ls[idx(i + 1, j)];
			lsM = ls[idx(i, j)];
			lsS = ls[idx(i, j - 1)];
			lsN = ls[idx(i, j + 1)];

			duWX = (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			duWY = (0.5 * (u[idx(i - 1, j + 1)] + u[idx(i, j + 1)]) 
				- 0.5 * (u[idx(i - 1, j - 1)] + u[idx(i, j - 1)])) / (2.0 * kDy);
			duEX = (u[idx(i + 2, j)] - u[idx(i + 1, j)]) / kDx;
			duEY = (0.5 * (u[idx(i + 1, j + 1)] + u[idx(i + 2, j + 1)]) 
				- 0.5 * (u[idx(i + 1, j - 1)] + u[idx(i + 2, j - 1)])) / (2.0 * kDy);
			duSX = (u[idx(i + 1, j - 1)] - u[idx(i, j - 1)]) / kDx;
			duSY = (0.5 * (u[idx(i, j)] + u[idx(i + 1, j)])
				- 0.5 * (u[idx(i, j - 2)] + u[idx(i + 1, j - 2)])) / (2.0 * kDy);
			duNX = (u[idx(i + 1, j + 1)] - u[idx(i, j + 1)]) / kDx;
			duNY = (0.5 * (u[idx(i, j)] + u[idx(i + 1, j)])
				- 0.5 * (u[idx(i, j + 2)] + u[idx(i + 1, j + 2)])) / (2.0 * kDy);
			duMX = (u[idx(i + 1, j)] - u[idx(i, j)]) / (kDx);
			duMY = (0.5 * (u[idx(i, j + 1)] + u[idx(i + 1, j + 1)]) 
				- 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / (2.0 * kDy);

			dvWX = (0.5 * (v[idx(i, j)] + v[idx(i, j + 1)])
				- 0.5 * (v[idx(i - 2, j)] + v[idx(i - 2, j + 1)])) / (2.0 * kDx);
			dvWY = (v[idx(i - 1, j + 1)] - v[idx(i - 1, j)]) / kDy;
			dvEX = (0.5 * (v[idx(i + 2, j)] + v[idx(i + 2, j + 1)])
				- 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)])) / (2.0 * kDx);
			dvEY = (v[idx(i + 1, j + 1)] - v[idx(i + 1, j)]) / kDy;
			dvSX = (0.5 * (v[idx(i + 1, j - 1)] + v[idx(i + 1, j)])
				- 0.5 * (v[idx(i - 1, j - 1)] + v[idx(i - 1, j)])) / (2.0 * kDx);
			dvSY = (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			dvNX = (0.5 * (v[idx(i + 1, j + 1)] + v[idx(i + 1, j + 2)])
				- 0.5 * (v[idx(i - 1, j + 1)] + v[idx(i - 1, j + 2)])) / (2.0 * kDx);
			dvNY = (v[idx(i, j + 2)] - v[idx(i, j)]) / kDy;
			dvMX = (0.5 * (v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
				- 0.5 * (v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / (2.0 * kDx);
			dvMY = (v[idx(i, j + 1)] - v[idx(i, j)]) / (2.0 * kDy);

			dlWX = (ls[idx(i, j)] - ls[idx(i - 2, j)]) / (2.0 * kDx);
			dlWY = (ls[idx(i - 1, j + 1)] - ls[idx(i - 1, j - 1)]) / (2.0 * kDy);
			dlEX = (ls[idx(i + 2, j)] - ls[idx(i, j)]) / (2.0 * kDx);
			dlEY = (ls[idx(i + 1, j + 1)] - ls[idx(i + 1, j - 1)]) / (2.0 * kDy);
			dlSX = (ls[idx(i + 1, j - 1)] - ls[idx(i - 1, j - 1)]) / (2.0 * kDx);
			dlSY = (ls[idx(i, j)] - ls[idx(i, j - 2)]) / (2.0 * kDy);
			dlNX = (ls[idx(i + 1, j + 1)] - ls[idx(i - 1, j + 1)]) / (2.0 * kDx);
			dlNY = (ls[idx(i, j + 2)] - ls[idx(i, j)]) / (2.0 * kDy);
			dlMX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
			dlMY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
			
			FW = 0.0;
			FE = 0.0;
			FS = 0.0;
			FN = 0.0;

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
				b = 0.0;
				aW = 2 * m_dt * (kMuO - kMuI) 
					* (duWX * dlWX * dlWX + duWY * dlWY * dlWX + dvWX * dlWX * dlWY * dvWY * dlWY * dlWY)
					/ (dlWX * dlWX + dlWY + dlWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
				aM = 2 * m_dt * (kMuO - kMuI)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoW = kRhoO * kRhoI * (fabs(lsW) + fabs(lsM)) / (kRhoO * fabs(lsW) + kRhoI * fabs(lsM));
				iRhoW = 1.0 / (kRhoO * kRhoI) 
					* (fabs(lsW) + fabs(lsM)) / (1.0 / kRhoO * fabs(lsW) + 1.0 / kRhoI * fabs(lsM));
				FW = iRhoW * aEff / (kDx * kDx) - iRhoW * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsW < 0 && lsM >= 0) {
				// interface lies between ls[i,j] and ls[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				b = 0.0;
				aW = 2 * m_dt * (kMuI - kMuO)
					* (duWX * dlWX * dlWX + duWY * dlWY * dlWX + dvWX * dlWX * dlWY * dvWY * dlWY * dlWY)
					/ (dlWX * dlWX + dlWY + dlWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
				aM = 2 * m_dt * (kMuI - kMuO)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoW = kRhoI * kRhoO * (fabs(lsW) + fabs(lsM)) / (kRhoI * fabs(lsW) + kRhoO * fabs(lsM));
				iRhoW = 1.0 / (kRhoI * kRhoO) 
					* (fabs(lsW) + fabs(lsM)) / (1.0 / kRhoI * fabs(lsW) + 1.0 / kRhoO * fabs(lsM));
				FW = -iRhoW * aEff / (kDx * kDx) + iRhoW * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FE
			// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
			if (lsM * lsE >= 0) {
				// one fluid, x direction
				FE = 0.0;
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
				b = 0.0;
				aM = 2 * m_dt * (kMuO - kMuI)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aE = 2 * m_dt * (kMuO - kMuI)
					* (duEX * dlEX * dlEX + duEY * dlEY * dlEX + dvEX * dlEX * dlEY * dvEY * dlEY * dlEY)
					/ (dlEX * dlEX + dlEY + dlEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoE = kRhoI * kRhoO * (fabs(lsE) + fabs(lsM)) / (kRhoI * fabs(lsE) + kRhoO * fabs(lsM));
				iRhoE = 1.0 / (kRhoI * kRhoO)
					* (fabs(lsE) + fabs(lsM)) / (1.0 / kRhoI * fabs(lsE) + 1.0 / kRhoO * fabs(lsM));
				FE = -iRhoE * aEff / (kDx * kDx) + iRhoE * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM < 0 && lsE >= 0) {
				// interface lies between ls[i,j] and ls[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
				b = 0.0;
				aM = 2 * m_dt * (kMuI - kMuO)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aE = 2 * m_dt * (kMuI - kMuO)
					* (duEX * dlEX * dlEX + duEY * dlEY * dlEX + dvEX * dlEX * dlEY * dvEY * dlEY * dlEY)
					/ (dlEX * dlEX + dlEY + dlEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoE = kRhoO * kRhoI * (fabs(lsE) + fabs(lsM)) / (kRhoO * fabs(lsE) + kRhoI * fabs(lsM));
				iRhoE = 1.0 / (kRhoO * kRhoI)
					* (fabs(lsE) + fabs(lsM)) / (1.0 / kRhoO * fabs(lsE) + 1.0 / kRhoI * fabs(lsM));
				FE = iRhoE * aEff / (kDx * kDx) - iRhoE * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FS
			if (lsS * lsM >= 0) {
				// one fluid, y direction
				FS = 0.0;
				if (lsS <= 0)
					iRhoS = 1.0 / kRhoO;
				else
					iRhoS = 1.0 / kRhoI;
			}
			else if (lsM >= 0 && lsS < 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = 2 * m_dt * (kMuI - kMuO)
					* (duSX * dlSX * dlSX + duSY * dlSY * dlSX + dvSX * dlSX * dlSY * dvSY * dlSY * dlSY)
					/ (dlSX * dlSX + dlSY + dlSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
				aM = 2 * m_dt * (kMuI - kMuO)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoS = kRhoO * kRhoI * (fabs(lsS) + fabs(lsM)) / (kRhoO * fabs(lsS) + kRhoI * fabs(lsM));
				iRhoS = 1.0 / (kRhoO * kRhoI)
					* (fabs(lsS) + fabs(lsM)) / (1.0 / kRhoO * fabs(lsS) + 1.0 / kRhoI * fabs(lsM));
				FS = iRhoS * aEff / (kDy * kDy) - iRhoS * b * theta / ((1.0 / kRhoO) * kDy);
			}
			else if (lsM <= 0 && lsS > 0) {
				// interface lies between ls[i,j] and ls[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = 2 * m_dt * (kMuO - kMuI)
					* (duSX * dlSX * dlSX + duSY * dlSY * dlSX + dvSX * dlSX * dlSY * dvSY * dlSY * dlSY)
					/ (dlSX * dlSX + dlSY + dlSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
				aM = 2 * m_dt * (kMuO - kMuI)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoS = kRhoI * kRhoO * (fabs(lsS) + fabs(lsM)) / (kRhoI * fabs(lsS) + kRhoO * fabs(lsM));
				iRhoS = 1.0 / (kRhoI * kRhoO)
					* (fabs(lsS) + fabs(lsM)) / (1.0 / kRhoI * fabs(lsS) + 1.0 / kRhoO * fabs(lsM));
				FS = -iRhoS * aEff / (kDy * kDy) + iRhoS * b * theta / ((1.0 / kRhoI) * kDy);
			}

			// FN
			if (lsM * lsN >= 0) {
				// one fluid, y direction
				FN = 0.0;
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
				b = 0.0;
				aM = 2 * m_dt * (kMuO - kMuI)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aN = 2 * m_dt * (kMuO - kMuI)
					* (duNX * dlNX * dlNX + duNY * dlNY * dlNX + dvNX * dlNX * dlNY * dvNY * dlNY * dlNY)
					/ (dlNX * dlNX + dlNY + dlNY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoN = kRhoI * kRhoO * (fabs(lsN) + fabs(lsM)) / (kRhoI * fabs(lsN) + kRhoO * fabs(lsM));
				iRhoN = 1.0 / (kRhoI * kRhoO)
					* (fabs(lsN) + fabs(lsM)) / (1.0 / kRhoI * fabs(lsN) + 1.0 / kRhoO * fabs(lsM));
				FN = -iRhoN * aEff / (kDx * kDx) + iRhoN * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM <= 0 && lsN > 0) {
				// interface lies between ls[i,j] and ls[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = 2 * m_dt * (kMuI - kMuO)
					* (duMX * dlMX * dlMX + duMY * dlMY * dlMX + dvMX * dlMX * dlMY * dvMY * dlMY * dlMY)
					/ (dlMX * dlMX + dlMY + dlMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)];
				aN = 2 * m_dt * (kMuI - kMuO)
					* (duNX * dlNX * dlNX + duNY * dlNY * dlNX + dvNX * dlNX * dlNY * dvNY * dlNY * dlNY)
					/ (dlNX * dlNX + dlNY + dlNY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoN = kRhoO * kRhoI * (fabs(lsN) + fabs(lsM)) / (kRhoO * fabs(lsN) + kRhoI * fabs(lsM));
				iRhoN = 1.0 / (kRhoO * kRhoI)
					* (fabs(lsN) + fabs(lsM)) / (1.0 / kRhoO * fabs(lsN) + 1.0 / kRhoI * fabs(lsM));
				FN = iRhoN * aEff / (kDx * kDx) - iRhoN * b * theta / ((1.0 / kRhoI) * kDx);
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
			AValsDic["E"] = -iRhoE/ (kDx * kDx);
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
				AValsDic["C"] -= iRhoS / (kDy * kDy);;

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

			// initially set variable coef. considered RHS
			rhs[idx(i, j)] = FW + FE + FS + FN;
			// if (FW + FE + FS + FN > 0.1)
			// 	std::cout << i << " " << j << " " << FW << " " << FE << " " << FS << " " << FN << std::endl;

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

		ARowIdx.push_back(rowIdx);
		m_Poisson->ICPCG_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC);
	}
	else if (m_PoissonSolverType == POISSONTYPE::GS) {

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
		double duWX = 0.0, duWY = 0.0, duEX = 0.0, duEY = 0.0, duSX = 0.0, duSY = 0.0, duNX = 0.0, duNY = 0.0, duMX = 0.0, duMY = 0.0;
		double dvWX = 0.0, dvWY = 0.0, dvEX = 0.0, dvEY = 0.0, dvSX = 0.0, dvSY = 0.0, dvNX = 0.0, dvNY = 0.0, dvMX = 0.0, dvMY = 0.0;
		double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aM = 0.0;
		double rhoW = 0.0, rhoE = 0.0, rhoS = 0.0, rhoN = 0.0;
		// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
		double a = 0.0, b = 0.0;
		double theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, bEff = 0.0, rhoEff = 0.0;
		UpdateKappa(ls);

		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
			lsW = ls[idx(i - 1, j)];
			lsE = ls[idx(i + 1, j)];
			lsM = ls[idx(i, j)];
			lsS = ls[idx(i, j - 1)];
			lsN = ls[idx(i, j + 1)];

			duWX = (u[idx(i, j)] - u[idx(i - 2, j)]) / (2.0 * kDx);
			duWY = (u[idx(i - 1, j + 1)] - u[idx(i - 1, j - 1)]) / (2.0 * kDy);
			duEX = (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			duEY = (u[idx(i + 1, j + 1)] - u[idx(i + 1, j - 1)]) / (2.0 * kDy);
			duSX = (u[idx(i + 1, j - 1)] - u[idx(i - 1, j - 1)]) / (2.0 * kDx);
			duSY = (u[idx(i, j)] - u[idx(i, j - 2)]) / (2.0 * kDy);
			duNX = (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			duNY = (u[idx(i + 1, j + 1)] - u[idx(i + 1, j - 1)]) / (2.0 * kDy);
			duMX = (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			duMY = (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);

			dvWX = (v[idx(i, j)] - v[idx(i - 2, j)]) / (2.0 * kDx);
			dvWY = (v[idx(i - 1, j + 1)] - v[idx(i - 1, j - 1)]) / (2.0 * kDy);
			dvEX = (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			dvEY = (v[idx(i + 1, j + 1)] - v[idx(i + 1, j - 1)]) / (2.0 * kDy);
			dvSX = (v[idx(i + 1, j - 1)] - v[idx(i - 1, j - 1)]) / (2.0 * kDx);
			dvSY = (v[idx(i, j)] - v[idx(i, j - 2)]) / (2.0 * kDy);
			dvNX = (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			dvNY = (v[idx(i + 1, j + 1)] - v[idx(i + 1, j - 1)]) / (2.0 * kDy);
			dvMX = (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			dvMY = (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);

			FW = 0.0;
			FE = 0.0;
			FS = 0.0;
			FN = 0.0;

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
				// interface lies between u[i,j] and u[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				// b always zero when solving level set (dealing with surface tension)
				b = 0.0;
				aW = (2 * m_dt * (kMuO - kMuI)
					* (duWX + duWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)]);
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoW = kRhoO * kRhoI * (fabs(lsW) + fabs(lsM)) / (kRhoI * fabs(lsW) + kRhoO * fabs(lsM));
				FW = (1.0 / rhoW) * aEff / (kDx * kDx) - (1.0 / rhoW) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsW < 0 && lsM >= 0) {
				// interface lies between u[i,j] and u[i - 1,j]
				theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
				// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
				b = 0.0;
				aW = (2 * m_dt * (kMuI - kMuO)
					* (duWX + duWY)
					+ m_dt * kSigma * m_kappa[idx(i - 1, j)]);
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsW) + aW * fabs(lsM)) / (fabs(lsM) + fabs(lsW));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoW = kRhoI * kRhoO * (fabs(lsW) + fabs(lsM)) / (kRhoO * fabs(lsW) + kRhoI * fabs(lsM));
				FW = -(1.0 / rhoW) * aEff / (kDx * kDx) + (1.0 / rhoW) * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FE
			// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
			if (lsM * lsE >= 0) {
				// one fluid, x direction
				FE = 0.0;
			}
			else if (lsM >= 0 && lsE < 0) {
				// interface lies between u[i,j] and u[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
				b = 0.0;
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuO - kMuI)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)]);
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoE = kRhoI * kRhoO * (fabs(lsE) + fabs(lsM)) / (kRhoI * fabs(lsE) + kRhoO * fabs(lsM));
				FE = (1.0 / rhoE) * aEff / (kDx * kDx) - (1.0 / rhoE) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM < 0 && lsE >= 0) {
				// interface lies between u[i,j] and u[i + 1,j]
				theta = fabs(lsE) / (fabs(lsM) + fabs(lsE));
				// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
				b = 0.0;
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuI - kMuO)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i + 1, j)]);
				aEff = (aM * fabs(lsE) + aE * fabs(lsM)) / (fabs(lsM) + fabs(lsE));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoE = kRhoO * kRhoI * (fabs(lsE) + fabs(lsM)) / (kRhoO * fabs(lsE) + kRhoI * fabs(lsM));
				FE = (1.0 / rhoE) * aEff / (kDx * kDx) - (1.0 / rhoE) * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FS
			if (lsS * lsM >= 0) {
				// one fluid, y direction
				FS = 0.0;
			}
			else if (lsM >= 0 && lsS < 0) {
				// interface lies between u[i,j] and u[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = (2 * m_dt * (kMuI - kMuO)
					* (duSX + duSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)]);
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoS = kRhoI * kRhoO * (fabs(lsS) + fabs(lsM)) / (kRhoI * fabs(lsS) + kRhoO * fabs(lsM));
				FS = (1.0 / rhoS) * aEff / (kDx * kDx) - (1.0 / rhoS) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM <= 0 && lsS > 0) {
				// interface lies between u[i,j] and u[i,j - 1]
				theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
				// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				b = 0.0;
				aS = (2 * m_dt * (kMuO - kMuI)
					* (duSX + duSY)
					+ m_dt * kSigma * m_kappa[idx(i, j - 1)]);
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aEff = (aM * fabs(lsS) + aS * fabs(lsM)) / (fabs(lsM) + fabs(lsS));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoS = kRhoO * kRhoI * (fabs(lsS) + fabs(lsM)) / (kRhoO * fabs(lsS) + kRhoI * fabs(lsM));
				FS = (1.0 / rhoS) * aEff / (kDx * kDx) - (1.0 / rhoS) * b * theta / ((1.0 / kRhoI) * kDx);
			}

			// FN
			if (lsM * lsN >= 0) {
				// one fluid, y direction
				FN = 0.0;
			}
			else if (lsM >= 0 && lsN < 0) {
				// interface lies between u[i,j] and u[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = (2 * m_dt * (kMuO - kMuI)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuO - kMuI)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)]);
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1 - theta));
				rhoN = kRhoI * kRhoO * (fabs(lsN) + fabs(lsM)) / (kRhoI * fabs(lsN) + kRhoO * fabs(lsM));
				FN = (1.0 / rhoN) * aEff / (kDx * kDx) - (1.0 / rhoN) * b * theta / ((1.0 / kRhoO) * kDx);
			}
			else if (lsM <= 0 && lsN > 0) {
				// interface lies between u[i,j] and u[i,j + 1]
				theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
				// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
				// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
				b = 0.0;
				aM = (2 * m_dt * (kMuI - kMuO)
					* (duMX + duMY)
					+ m_dt * kSigma * m_kappa[idx(i, j)]);
				aE = (2 * m_dt * (kMuI - kMuO)
					* (duEX + duEY)
					+ m_dt * kSigma * m_kappa[idx(i, j + 1)]);
				aEff = (aM * fabs(lsN) + aN * fabs(lsM)) / (fabs(lsM) + fabs(lsN));
				rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1 - theta));
				rhoN = kRhoO * kRhoI * (fabs(lsN) + fabs(lsM)) / (kRhoO * fabs(lsN) + kRhoI * fabs(lsM));
				FN = (1.0 / rhoN) * aEff / (kDx * kDx) - (1.0 / rhoN) * b * theta / ((1.0 / kRhoI) * kDx);
			}

			rhs[idx(i, j)] = FW + FE + FS + FN;
		}

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
