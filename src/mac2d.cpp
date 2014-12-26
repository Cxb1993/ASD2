#include "mac2d.h"

MACSolver2D::MACSolver2D(double Re, double We, double Fr,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoI, double muI,
	int nx, int ny, int nz, double lenX, double lenY, double lenZ, double cfl,
	int maxtime, int maxiter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoI(rhoI), kRhoO(rhoI / densityRatio), kRhoRatio(densityRatio),
	kMuI(muI), kMuO(muI / viscosityRatio), kMuRatio(viscosityRatio),
	kNx(nx), kNy(ny), kNz(nz),
	kLenX(lenX), kLenY(lenY), kLenZ(lenZ),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)), kDz(lenZ / static_cast<double>(nz)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxiter), kNIterSkip(niterskip), kNumBCGrid(num_bc_grid),
	kWriteVTK(writeVTK) {

	m_iter = 0;
	m_totTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(kDx, kDy) / kUScale);
}

MACSolver2D::MACSolver2D(double rhoI, double rhoO, double muI, double muO, double gConstant,
	double L, double U, double sigma, int nx, int ny, int nz, double lenX, double lenY, double lenZ, double cfl,
	int maxtime, int maxiter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rho1), kMuScale(mu1), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoI(rhoI), kRhoO(rhoO), kMuI(muI), kMuO(muO), kRhoRatio(rhoI / rhoO), kMuRatio(muI / muO),
	kRe(rhoI * L * U / mu1), kWe(rhoI * L * U * U / sigma), kFr(U * U / (gConstant * L)),
	kNx(nx), kNy(ny), kNz(nz),
	kLenX(lenX), kLenY(lenY), kLenZ(lenZ),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)), kDz(lenZ / static_cast<double>(nz)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxiter), kNIterSkip(niterskip), kNumBCGrid(num_bc_grid),
	kWriteVTK(writeVTK) {

	m_iter = 0;
	m_totTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(kDx, kDy) / kUScale);
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

	return 0;
}

int MACSolver2D::UpdateKappa(const std::vector<double>& ls) {
	std::unique_ptr<double[]> dPhidX(new double((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid)));
	std::unique_ptr<double[]> dPhidY(new double((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid)));
	std::unique_ptr<double[]> Squared(new double((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid)));
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		dPhidX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dPhidX[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDx);
		
		Squared[idx(i, j)] = (dPhidX[(idx(i, j)] * dPhidX[(idx(i, j)] + dPhidY[(idx(i, j)] * dPhidY[(idx(i, j)]);
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		m_kappa[idx(i, j)]
			= (dPhidX[idx(i, j)] * dPhidX[idx(i, j)] * (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDy * kDy) // phi^2_x \phi_yy
			- 2.0 * dPhidX[idx(i, j)] * dPhidY[idx(i, j)] * (dPhidX[i][j + 1] - dPhidX[i][j - 1]) / (2.0 * kDy) //2 \phi_x \phi_y \phi_xy
			+ (dPhidY[idx(i, j)] * dPhidY[idx(i, j)]) * (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / (kDx * kDx) // phi^2_y \phi_xx
				/ std::pow(Squared[idx(i, j)], 1.5);

		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j)] = std::min(m_kappa[idx(i, j)], 1.0 / std::min(kDx, kDy);
	}

	return 0;
}

std::vector<double> MACSolver2D::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> cU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> tmpV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		tmpV[idx(i, j)] = (v[idx(i, j)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)] + v[idx(i, j + 1)]) * 0.25;

	std::vector<double> vecF_UX(kNx + 2 * kNumBCGrid, 0.0), vecF_UY(kNy + 2 * kNumBCGrid);

	// U : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_UX[i] = u[idx(i, j)];
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
			vecF_UY[i] = u[idx(i, j)];
		}

		UnitHJWENO5(vecF_UY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; i++) {
			LYP[idx(i, j)] = FYP[i];
			LYM[idx(i, j)] = FYM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_UY.begin(), vecF_UY.end(), 0.0);
	}
	
	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0;

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		alphaX = u[idx(i, j)];
		alphaY = tmpV[idx(i, j)];

		cU[idx(i, j)]
				= 0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
				+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
				- alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
				- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)]));
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

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		tmpU[idx(i, j)] = (u[idx(i, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)] + u[idx(i + 1, j)]) * 0.25;


	std::vector<double> vecF_VX(kNx + 2 * kNumBCGrid, 0.0), vecF_VY(kNy + 2 * kNumBCGrid);
	
	// V : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_VX[i] = v[idx(i, j)];
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
			vecF_VY[i] = v[idx(i, j)];
		}

		UnitHJWENO5(vecF_VY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; i++) {
			LYP[idx(i, j)] = FYP[i];
			LYM[idx(i, j)] = FYM[i];
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
			= 0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
			+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
			- alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
			- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)]));
	}

	return cV;
}

int MACSolver2D::UnitHJWENO5(
	const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n) {
	int stat = 0;
	const double kWeightEps = 1e-6;

	double alphaPlus[3], alphaMinus[3];
	double weightPlus[3], weightMinus[3];

	std::vector<double> IS0 = std::vector<double>(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS1 = std::vector<double>(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS2 = std::vector<double>(n + 2 * kNumBCGrid, 0.0);
	
	// \dfrac{\Delta^+ F}{\Delta x}
	std::vector<double> DFPlus = std::vector<double>(n + 2 * kNumBCGrid, 0.0);

	// Compute Smoothness for phi^{-}
	for (int i = 0; i < n + 2 * m_num_bc_grid - 1; i++)
		DFPlus[i] = F[i + 1] - F[i];

	for (int i = m_num_bc_grid; i < n + m_num_bc_grid; i++) {
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
	for (int i = m_num_bc_grid; i < n + m_num_bc_grid; i++) {
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
	for (int i = m_num_bc_grid; i < n + m_num_bc_grid; i++) {
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
	for (int i = m_num_bc_grid; i < n + m_num_bc_grid; i++) {
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
	std::vector<double> dU = std::vector((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0));
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dU[idx(i, j)] = 0.0;
		
		return dU;
	}

	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	// jump condition
	double JW = 0.0, JE = 0.0, JS = 0.0, JN = 0.0;
	double theta = 0.0;
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
	double rhoU_X_W = 0.0, rhoU_X_E = 0.0, rhoU_Y_S = 0.0, rhoU_Y_N = 0.0;
	double visX = 0.0, visY = 0.0;
	// effective Jump condition, effective u(uEff), and effective mu (muEff)
	double Jeff = 0.0, JO = 0.0, uEff = 0.0, muEff = 0.0;
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		visX = 0.0, visY = 0.0;
		lsW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsS = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i - 1, j - 1)]);
		lsN = 0.5 * (ls[idx(i, j + 1)] + ls[idx(i - 1, j + 1)]);
		
		if (lsW <= 0 && lsM <= 0 && lsE <= 0) {
			// one fluid(inside, negative levelset), x direction
			muU_X_W = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i - 1,j)]) / kDx;
			muU_X_E = kMuI / kRhoI * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		}
		else if (lsW >= 0 && lsM >= 0 && lsE >= 0) {
			// one fluid(outside, positive levelset), x direction
			muU_X_W = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muU_X_E = kMuO / kRhoO * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| ===  inside   === |(interface)| ===    outside      === |(lsM)|
			// |(lsW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i - 2, j)]) / (2.0 * kDx);
			JM = kMuO / kRhoO * (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JW;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i - 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1- theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_W = muEff / rhoEff * (u[idx(i, j)] - u[idx(i - 1,j)]) / kDx
				+ (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
			muU_X_E = kMuI / kRhoI * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		}
		else if (lwW > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| ===  outside  === |(interface)| ===     inside      === |(lsM)|
			// |(lsW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i - 2, j)]) / (2.0 * kDx);
			JM = kMuI / kRhoI * (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JW;
			uEff = (kMuI * u[idx(i, j)] * theta + kMuO * u[idx(i - 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_W = muEff / rhoEff * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx
				- (muEff / rhoEff) * JEff * theta / (kMuO / kRhoO);
			muU_X_E = kMuI / kRhoI * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		}
		else if (lwE > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsE)|
			JM = kMuI / kRhoI * (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			JE = kMuO / kRhoO * (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuI / kRhoI * u[idx(i, j)] * theta + kMuO / kRhoO * u[idx(i + 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_W = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muU_X_E = (muEff / rhoEff) * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- muEff * JEff * theta / (kMuO / kRhoO);
		}
		else if (lwE <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsE)|
			JM = kMuO / kRhoO * (u[idx(i + 1, j)] - u[idx(i - 1, j)]) / (2.0 * kDx);
			JE = kMuI / kRhoI * (u[idx(i + 2, j)] - u[idx(i, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuO / kRhoO * u[idx(i, j)] * theta + kMuI / kRhoI * u[idx(i + 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuO / kRhoO * theta + kMuI / kRhoI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_W = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i - 1, j)]) / kDx;
			muU_X_E = (muEff / rhoEff) * (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
				- (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
		}
		
		if (lsS <= 0 && lsM <= 0 && lsN <= 0) {
			// one fluid(inside, negative levelset), y direction
			muU_Y_S = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
			muU_Y_N = kMuI / kRhoI * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lsS >= 0 && lsM >= 0 && lsN >= 0) {
			// one fluid(outside, postiive levelset), y direction
			muU_Y_S = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
			muU_Y_N = kMuO / kRhoO * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lsS <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside   === |(interface)| ===    outside      === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i, j - 2)]) / (2.0 * kDy);
			JM = kMuO / kRhoO * (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JS;
			uEff = (kMuO * u[idx(i, j)] * theta + kMuI * u[idx(i, j - 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_S = muEff / rhoEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				+ (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
			muU_X_N = kMuI / kRhoI * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lwS > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  outside  === |(interface)| ===     inside      === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i, j - 2)]) / (2.0 * kDy);
			JM = kMuI / kRhoI * (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JS;
			uEff = (kMuI / kRhoI * u[idx(i, j)] * theta + kMuO / kRhoO * u[idx(i, j - 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_S = muEff / rhoEff * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy
				- (muEff / rhoEff) * JEff * theta / (kMuO / kRhoO);
			muU_X_N = kMuI / kRhoI * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy;
		}
		else if (lwN > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			JM = kMuI / kRhoI * (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);
			JN = kMuO / kRhoO * (u[idx(i, j + 2)] - u[idx(i, j)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuI / kRhoI * u[idx(i, j)] * theta + kMuO / kRhoO * u[idx(i, j + 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_S = kMuI / kRhoI * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
			muU_X_N = (muEff / rhoEff) * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuO / kRhoO);
		}
		else if (lwN <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			JM = kMuO / kRhoO * (u[idx(i, j + 1)] - u[idx(i, j - 1)]) / (2.0 * kDy);
			JN = kMuI / kRhoI * (u[idx(i, j + 2)] - u[idx(i, j)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuO / kRhoO * u[idx(i, j)] * theta + kMuI / kRhoI * u[idx(i, j + 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuO / kRhoO * theta + kMuI / kRhoI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_S = kMuO / kRhoO * (u[idx(i, j)] - u[idx(i, j - 1)]) / kDy;
			muU_X_N = (muEff / rhoEff) * (u[idx(i, j + 1)] - u[idx(i, j)]) / kDy
				- (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
		}
		
		visX = (muU_X_E - muU_X_W) / kDx;
		visY = (muU_X_N - muU_X_S) / kDy;
		
		dU[idx(i, j)] = visX + visY;
	}

	return dU;
}


std::vector<double> MACSolver2D::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV = std::vector((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0));
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dV[idx(i, j)] = 0.0;

		return dV;
	}

	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	// jump condition
	double JW = 0.0, JE = 0.0, JS = 0.0, JN = 0.0;
	double theta = 0.0;
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0;
	double rhoV_X_W = 0.0, rhoV_X_E = 0.0, rhoV_Y_S = 0.0, rhoV_Y_N = 0.0;
	double visX = 0.0, visY = 0.0;
	// effective Jump condition, effective v (vEff), and effective mu (muEff)
	double Jeff = 0.0, JO = 0.0, vEff = 0.0, muEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		visX = 0.0, visY = 0.0;
		lsW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
		lsM = 0.5 * (ls[idx(i, j)] + ls[idx(i, j - 1)]);
		lsS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
		lsN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 2)]);

		if (lsW <= 0 && lsM <= 0 && lsE <= 0) {
			// one fluid(inside, negative levelset), x direction
			muV_X_W = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
			muV_X_E = kMuI / kRhoI * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lsW >= 0 && lsM >= 0 && lsE >= 0) {
			// one fluid(outside, positive levelset), x direction
			muV_X_W = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
			muV_X_E = kMuO / kRhoO * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| ===  inside   === |(interface)| ===    outside      === |(lsM)|
			// |(lsW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i - 2, j)]) / (2.0 * kDx);
			JM = kMuO / kRhoO * (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JW;
			uEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i - 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuO * theta + kMuI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_W = muEff / rhoEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				+ (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
			muU_X_E = kMuI / kRhoI * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lwW > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i - 1,j]
			theta = fabs(lsW) / (fabs(lsW) + fabs(lsM));
			// |(lsW)| ===  outside  === |(interface)| ===     inside      === |(lsM)|
			// |(lsW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JW = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i - 2, j)]) / (2.0 * kDx);
			JM = kMuI / kRhoI * (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JW;
			uEff = (kMuI * v[idx(i, j)] * theta + kMuO * v[idx(i - 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuI * theta + kMuO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_W = muEff / rhoEff * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx
				- (muEff / rhoEff) * JEff * theta / (kMuO / kRhoO);
			muU_X_E = kMuI / kRhoI * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx;
		}
		else if (lwE > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsE)|
			JM = kMuI / kRhoI * (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			JE = kMuO / kRhoO * (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuI / kRhoI * v[idx(i, j)] * theta + kMuO / kRhoO * v[idx(i + 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_W = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
			muU_X_E = (muEff / rhoEff) * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- muEff * JEff * theta / (kMuO / kRhoO);
		}
		else if (lwE <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i + 1,j]
			theta = fabs(lsE) / (fabs(lsE) + fabs(lsM));
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsE)|
			JM = kMuO / kRhoO * (v[idx(i + 1, j)] - v[idx(i - 1, j)]) / (2.0 * kDx);
			JE = kMuI / kRhoI * (v[idx(i + 2, j)] - v[idx(i, j)]) / (2.0 * kDx);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuO / kRhoO * v[idx(i, j)] * theta + kMuI / kRhoI * v[idx(i + 1, j)] - JEff * theta * (1 - theta) * kDx)
				/ (kMuO / kRhoO * theta + kMuI / kRhoI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_W = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i - 1, j)]) / kDx;
			muU_X_E = (muEff / rhoEff) * (v[idx(i + 1, j)] - v[idx(i, j)]) / kDx
				- (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
		}

		if (lsS <= 0 && lsM <= 0 && lsN <= 0) {
			// one fluid(inside, negative levelset), y direction
			muU_Y_S = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_Y_N = kMuI / kRhoI * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		}
		else if (lsS >= 0 && lsM >= 0 && lsN >= 0) {
			// one fluid(outside, postiive levelset), y direction
			muU_Y_S = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_Y_N = kMuO / kRhoO * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		}
		else if (lsS <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  inside   === |(interface)| ===    outside      === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i, j - 2)]) / (2.0 * kDy);
			JM = kMuO / kRhoO * (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JS;
			uEff = (kMuO * v[idx(i, j)] * theta + kMuI * v[idx(i, j - 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuO * theta + kMuI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_S = muEff / rhoEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				+ (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
			muU_X_N = kMuI / kRhoI * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		}
		else if (lwS > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i,j - 1]
			theta = fabs(lsS) / (fabs(lsS) + fabs(lsM));
			// |(lsS)| ===  outside  === |(interface)| ===     inside      === |(lsM)|
			// |(lsS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsM)|
			JS = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i, j - 2)]) / (2.0 * kDy);
			JM = kMuI / kRhoI * (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JS;
			uEff = (kMuI / kRhoI * v[idx(i, j)] * theta + kMuO / kRhoO * v[idx(i, j - 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuI * (1.0 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_S = muEff / rhoEff * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy
				- (muEff / rhoEff) * JEff * theta / (kMuO / kRhoO);
			muU_X_N = kMuI / kRhoI * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		}
		else if (lwN > 0 && lsM <= 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			JM = kMuI / kRhoI * (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);
			JN = kMuO / kRhoO * (v[idx(i, j + 2)] - v[idx(i, j)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuI / kRhoI * v[idx(i, j)] * theta + kMuO / kRhoO * v[idx(i, j + 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuI / kRhoI * theta + kMuO / kRhoO * (1 - theta));
			muEff = kMuI * kMuO / (kMuI * theta + kMuO * (1 - theta));
			rhoEff = kRhoI * kRhoO / (kRhoI * theta + kRhoO * (1.0 - theta));
			muU_X_S = kMuI / kRhoI * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_X_N = (muEff / rhoEff) * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- muEff * JEff * theta / (kMuO / kRhoO);
		}
		else if (lwN <= 0 && lsM > 0) {
			// interface lies between u[i,j] and u[i,j + 1]
			theta = fabs(lsN) / (fabs(lsN) + fabs(lsM));
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			JM = kMuO / kRhoO * (v[idx(i, j + 1)] - v[idx(i, j - 1)]) / (2.0 * kDy);
			JN = kMuI / kRhoI * (v[idx(i, j + 2)] - v[idx(i, j)]) / (2.0 * kDy);
			JEff = theta * JM + (1 - theta) * JE;
			uEff = (kMuO / kRhoO * v[idx(i, j)] * theta + kMuI / kRhoI * v[idx(i, j + 1)] - JEff * theta * (1 - theta) * kDy)
				/ (kMuO / kRhoO * theta + kMuI / kRhoI * (1 - theta));
			muEff = kMuO * kMuI / (kMuO * theta + kMuI * (1 - theta));
			rhoEff = kRhoO * kRhoI / (kRhoO * theta + kRhoI * (1.0 - theta));
			muU_X_S = kMuO / kRhoO * (v[idx(i, j)] - v[idx(i, j - 1)]) / kDy;
			muU_X_N = (muEff / rhoEff) * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy
				- (muEff / rhoEff) * JEff * theta / (kMuI / kRhoI);
		}

		visX = (muV_X_E - muV_X_W) / kDx;
		visY = (muV_X_N - muV_X_S) / kDy;

		dV[idx(i, j)] = visX + visY;
	}

	return dV;
}

std::vector<double> MACSolver2D::AddGravityUF() {
	std::vector<double> gU = std::vector((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0));

	if (kFr == 0) {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
				gU[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			gU[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gU;
}

std::vector<double> MACSolver2D::AddGravityVF() {
	std::vector<double> gV = std::vector((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0));

	if (kFr == 0) {
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

int MACSolver2D::SetPoissonSolver(POISSONTYPE type) {
	m_PoissonSolverType = type;
	if (!m_Poisson)
		m_Poisson = std::make_shared<PoissonSolver>();

	return 0;
}

int MACSolver2D::UpdateVel(std::vector<double>& u, std::vector<double>& v,
	const std::vector<double>& us, const std::vector<double>& vs, const std::vector<double>& ps) {

	// velocity update after solving poisson equation
	// ps = p * dt
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			u[idx(i, j)] = us[idx(i, j)] - (ps[idx(i, j)] - ps[idx(i - 1, j])) / (kDx * m_rho[idx(i, j)]);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			v[idx(i, j)] = vs[idx(i, j)] - (ps[idx(i, j)] - ps[idx(i, j - 1])) / (kDy * m_rho[idx(i, j)]);

	return 0;
}

int MACSolver2D::SetBC_U_2D(std::string BC_W, std::string BC_E,	std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2D::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2D::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition>(kNx, kNy, kNz, kNumBCGrid);
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

int MACSolver2D::ApplyBC_W_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_W_2D(arr);
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

inline int idx(int i, int j) {
	return i + (kNx + 2 * kNumBCGrid) * j;
}