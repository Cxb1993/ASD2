#include "mac2dAxisym.h"

MACSolver2DAxisym::MACSolver2DAxisym(double Re, double We, double Fr, GAXISENUM2DAXISYM GAxis,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoH, double muH,
	int nr, int nz, double baseR, double baseZ, double lenR, double lenZ,
	double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kGAxis(GAxis), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoH(rhoH), kRhoL(rhoH * densityRatio), kRhoRatio(densityRatio),
	kMuH(muH), kMuL(muH * viscosityRatio), kMuRatio(viscosityRatio),
	kNr(nr), kNz(nz), kBaseR(baseR), kBaseZ(baseZ), kLenR(lenR), kLenZ(lenZ),
	kDr(lenR / static_cast<double>(nr)), kDz(lenZ / static_cast<double>(nz)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
	kArrSize(
	static_cast<int64_t>(kNr + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNz + 2 * kNumBCGrid)) {

	// positive level set : inside
	// negative level set : outside
	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenR / static_cast<double>(nr), lenZ / static_cast<double>(nz)) / U);
}

MACSolver2DAxisym::MACSolver2DAxisym(double rhoH, double rhoL, double muH, double muL, double gConstant, GAXISENUM2DAXISYM GAxis,
	double L, double U, double sigma, int nr, int nz, double baseR, double baseZ, double lenR, double lenZ,
	double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoH), kMuScale(muH), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoH(rhoH), kRhoL(rhoL), kRhoRatio(rhoL / rhoH), 
	kMuH(muH), kMuL(muL), kMuRatio(muL / muH),
	kRe(rhoH * L * U / muH), kWe(rhoH * L * U * U / sigma), kFr(U * U / (gConstant * L)), kGAxis(GAxis),
	kNr(nr), kNz(nz), kBaseR(baseR), kBaseZ(baseZ), kLenR(lenR), kLenZ(lenZ),
	kDr(lenR / static_cast<double>(nr)), kDz(lenZ / static_cast<double>(nz)),
	kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
	kArrSize(
	static_cast<int64_t>(kNr + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNz + 2 * kNumBCGrid))  {

	// positive level set : inside
	// negative level set : outside
	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenR / static_cast<double>(nr), lenZ / static_cast<double>(nz)) / U);
}

MACSolver2DAxisym::~MACSolver2DAxisym() {
	if (m_BC)
		m_BC.reset();
}

// Deallocated automatically
int MACSolver2DAxisym::AllocateVariables() {
	m_u = std::vector<double>(kArrSize, 0);
	m_v = std::vector<double>(kArrSize, 0);
	
	m_rho = std::vector<double>(kArrSize, 0);
	m_kappa = std::vector<double>(kArrSize, 0);
	m_mu = std::vector<double>(kArrSize, 0);
	
	m_ps = std::vector<double>(kArrSize, 0);
	m_p = std::vector<double>(kArrSize, 0);

	m_J11 = std::vector<double>(kArrSize, 0);
	m_J12 = std::vector<double>(kArrSize, 0);
	m_J21 = std::vector<double>(kArrSize, 0);
	m_J22 = std::vector<double>(kArrSize, 0);

	m_nr = std::vector<double>(kArrSize, 0);
	m_nz = std::vector<double>(kArrSize, 0);

	m_t1r = std::vector<double>(kArrSize, 0);
	m_t1z = std::vector<double>(kArrSize, 0);
	m_t1t = std::vector<double>(kArrSize, 0);

	m_t2r = std::vector<double>(kArrSize, 0);
	m_t2z = std::vector<double>(kArrSize, 0);
	m_t2t = std::vector<double>(kArrSize, 0);

	return 0;
}

std::vector<double> MACSolver2DAxisym::UpdateHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside(kArrSize, 0.0);
	const double eps = std::min(kDr, kDz) * 1.5;

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		// inside
		if (ls[idx(i, j)] >= 0.0)
			Heaviside[idx(i, j)] = 1.0;
		else
			Heaviside[idx(i, j)] = 0.0;
	}

	return Heaviside;
}

std::vector<double> MACSolver2DAxisym::UpdateSmoothHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside(kArrSize, 0.0);
	const double eps = std::min(kDr, kDz) * 1.5;

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		// High
		if (ls[idx(i, j)] > eps)
			Heaviside[idx(i, j)] = 1.0;
		// Lower
		else if (ls[idx(i, j)] < -eps)
			Heaviside[idx(i, j)] = 0.0;
		else
			Heaviside[idx(i, j)] 
				= 0.5 * (1.0 + ls[idx(i, j)] / eps
					+ 1.0 / M_PI * sin(M_PI * ls[idx(i, j)] / eps));
	}
		
	return Heaviside;
}

int MACSolver2DAxisym::UpdateNTK(const std::shared_ptr<LevelSetSolver2D>& LSolver, const std::vector<double>& ls,
	std::tuple<std::vector<double>&, std::vector<double>&> normalVec,
	std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t1Vec,
	std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t2Vec,
	std::vector<double>& kappa) {
	// Lervåg, Karl Yngve, Bernhard Müller, and Svend Tollak Munkejord.
	// "Calculation of the interface curvature and normal vector with the level-set method."
	//  Computers & Fluids 84 (2013) : 218 - 230.
	//
	// Macklin, Paul, and John S. Lowengrub.
	//  "A new ghost cell/level set method for moving boundary problems: application to tumor growth."
	// Journal of scientific computing 35.2-3 (2008): 266-299.
	const double eta = 0.0001;
	const double eps = 1.0e-100;

	std::vector<double> Q(kArrSize, 0.0),
		absDerivLS(kArrSize, 0.0);

	absDerivLS = LSolver->ENO_DerivAbsLS_2D(ls, ls);

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		Q[idx(i, j)] = std::fabs(1.0 - absDerivLS[idx(i, j)]);
	}
	
	double dLSdR = 0.0, dLSdZ = 0.0;
	// determination function
	int Dr = 0.0, Dz = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		// complex level set, such as droplet merging
		if (Q[idx(i, j)] > eta) {
			// Using Dierctional Direction First
			if (Q[idx(i - 1, j)] < eta && Q[idx(i + 1, j)] >= eta)
				Dr = -1;
			else if (Q[idx(i - 1, j)] >= eta && Q[idx(i + 1, j)] < eta)
				Dr = 1;
			else if (Q[idx(i - 1, j)] < eta && Q[idx(i, j)] < eta && Q[idx(i + 1, j)] < eta)
				Dr = 0;
			else if (Q[idx(i - 1, j)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i + 1, j)] >= eta)
				Dr = 0;
			else
				Dr = 2;

			// determination funciton
			if (Q[idx(i, j - 1)] < eta && Q[idx(i, j + 1)] >= eta)
				Dz = -1;
			else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j + 1)] < eta)
				Dz = 1;
			else if (Q[idx(i, j - 1)] < eta && Q[idx(i, j)] < eta && Q[idx(i, j + 1)] < eta)
				Dz = 0;
			else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i, j + 1)] >= eta)
				Dz = 0;
			else
				// undetermined
				Dz = 2;

			if (Dr == -1)
				dLSdR = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDr;
			else if (Dr == 1)
				dLSdR = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDr;
			else if (Dr == 0)
				dLSdR = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDr);

			if (Dz == -1)
				dLSdZ = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDz;
			else if (Dz == 1)
				dLSdZ = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDz;
			else if (Dz == 0)
				dLSdZ = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDz);

			// Curve Fitting Method
			if (Dr == 2) {
				// Breadfirst Search

				// Curve Fitting

				// Normal Vector

				// Tangent Vector

				// Curvature

			}

		}
		// simple level set
		else {
			dLSdR = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDr);
			dLSdZ = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDz);
			std::get<0>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdR * dLSdR + dLSdZ * dLSdZ + eps) * dLSdR;
			std::get<1>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdR * dLSdR + dLSdZ * dLSdZ + eps) * dLSdZ;
		}
	}

	return 0;
}

int MACSolver2DAxisym::UpdateKappa(const std::vector<double>& ls) {
	std::vector<double> dLSdR(kArrSize, 0.0);
	std::vector<double> dLSdZ(kArrSize, 0.0);
	std::vector<double> LSSize(kArrSize, 0.0);
	const double eps = 1.0e-100;
	double rPM = 0.0, rPW = 0.0, rPE = 0.0;

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		dLSdR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDr);
		dLSdZ[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDz);
		
		// exclude symmetric axis
		if (i == kNumBCGrid)
			dLSdR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDr;

		LSSize[idx(i, j)] = std::sqrt(dLSdR[idx(i, j)] * dLSdR[idx(i, j)] + dLSdZ[idx(i, j)] * dLSdZ[idx(i, j)] + eps);

		if (LSSize[idx(i, j)] < eps)
			perror("Div/0 Err in computing kappa");
	}
	
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;

		m_kappa[idx(i, j)]
			= ((ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j )]) / (kDr * kDr) * dLSdZ[idx(i, j)] * dLSdZ[idx(i, j)]
			- 2.0 * (dLSdZ[idx(i + 1, j)] - dLSdZ[idx(i - 1, j)]) / (2.0 * kDr) * dLSdR[idx(i, j)] * dLSdZ[idx(i, j)]
			+ (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDz * kDz) * dLSdR[idx(i, j)] * dLSdR[idx(i, j)]
			+ dLSdR[idx(i, j)] * std::pow(LSSize[idx(i, j)], 2.0) / rPM)
			/ std::pow(LSSize[idx(i, j)], 3.0);
		
		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j)] = sign(m_kappa[idx(i, j)]) * std::fabs(std::min(std::fabs(m_kappa[idx(i, j)]), 1.0 / std::min(kDr, kDz)));
		
		assert(m_kappa[idx(i, j)] == m_kappa[idx(i, j)]);
	}

	return 0;
}

std::vector<double> MACSolver2DAxisym::GetRHSU(const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& H) {
	
	std::vector<double> cU(kArrSize, 0.0);
	std::vector<double> vU(kArrSize, 0.0);
	std::vector<double> gU(kArrSize, 0.0);
	std::vector<double> rhsU(kArrSize, 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v);

	// Viscous term
	vU = this->AddExternalViscosityFU(u, v, ls, H);

	gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsM = 0.0;
	double thetaH = 0.0, rhoEff = 0.0;

	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rhsU[idx(i, j)] = u[idx(i, j)] + m_dt * (-cU[idx(i, j)] + gU[idx(i, j)] + vU[idx(i, j)]);
	}
	
	ApplyBC_U_2D(rhsU);
	return rhsU;
}

std::vector<double> MACSolver2DAxisym::GetRHSV(const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& uhat,
	const std::vector<double>& v, const std::vector<double>& H) {

	std::vector<double> cV(kArrSize, 0.0);
	std::vector<double> vV(kArrSize, 0.0);
	std::vector<double> gV(kArrSize, 0.0);
	std::vector<double> rhsV(kArrSize, 0.0);

	// Convection term
	cV = this->AddConvectionFV(u, v);

	// Viscous term
	vV = this->AddExternalViscosityFV(u, uhat, v, ls, H);

	gV = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsS = 0.0, lsM = 0.0;
	double thetaH = 0.0, rhoEff = 0.0;

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		rhsV[idx(i, j)] = v[idx(i, j)] + m_dt * (-cV[idx(i, j)] + gV[idx(i, j)] + vV[idx(i, j)]);
	}

	ApplyBC_V_2D(rhsV);
	return rhsV;
}

std::vector<double> MACSolver2DAxisym::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> cU(kArrSize, 0.0);
	std::vector<double> tmpV(kArrSize, 0.0);
	std::vector<double> LRP(kArrSize, 0.0);
	std::vector<double> LRM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0);
	std::vector<double> LZM(kArrSize, 0.0);
	
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		tmpV[idx(i, j)] = (v[idx(i, j)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)] + v[idx(i, j + 1)]) * 0.25;

	std::vector<double> vecF_UR(kNr + 2 * kNumBCGrid, 0.0), vecF_UZ(kNz + 2 * kNumBCGrid);

	// U : X direction
	std::vector<double> FRP(kNr + 2 * kNumBCGrid, 0.0), FRM(kNr + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++) {
			vecF_UR[i] = u[idx(i, j)] * u[idx(i, j)];
		}

		UnitHJWENO5(vecF_UR, FRP, FRM, kDr, kNr);

		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++) {
			LRP[idx(i, j)] = FRP[i];
			LRM[idx(i, j)] = FRM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FRP.begin(), FRP.end(), 0.0);
		std::fill(FRM.begin(), FRM.end(), 0.0);
		std::fill(vecF_UR.begin(), vecF_UR.end(), 0.0);
	}
	
	// U : Y direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++) {
		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
			vecF_UZ[j] = u[idx(i, j)] * tmpV[idx(i, j)];
		}

		UnitHJWENO5(vecF_UZ, FZP, FZM, kDz, kNz);

		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
			LZP[idx(i, j)] = FZP[j];
			LZM[idx(i, j)] = FZM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_UZ.begin(), vecF_UZ.end(), 0.0);
	}
	
	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, rhoEff = 0.0;

	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		alphaX = u[idx(i, j)];
		alphaY = tmpV[idx(i, j)];
		
		cU[idx(i, j)]
				= (0.5 * (LRP[idx(i, j)] + LRM[idx(i, j)])
				+ 0.5 * (LZP[idx(i, j)] + LZM[idx(i, j)])
				- alphaX * (0.5 * (LRP[idx(i, j)] - LRM[idx(i, j)]))
				- alphaY * (0.5 * (LZP[idx(i, j)] - LZM[idx(i, j)])));
		
		if (isnan(cU[idx(i, j)]) || isinf(cU[idx(i, j)])) {
			std::cout << (0.5 * (LRP[idx(i, j)] + LRM[idx(i, j)])
				+ 0.5 * (LZP[idx(i, j)] + LZM[idx(i, j)])) << " " << -alphaX * (0.5 * (LRP[idx(i, j)] - LRM[idx(i, j)]))
				- alphaY * (0.5 * (LZP[idx(i, j)] - LZM[idx(i, j)])) << std::endl;
		}
		assert(cU[idx(i, j)] == cU[idx(i, j)]);
		if (std::isnan(cU[idx(i, j)]) || std::isinf(cU[idx(i, j)])) {
			std::cout << "U-convection term nan/inf error : " << i << " " << j << " " << cU[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return cU;
}

std::vector<double> MACSolver2DAxisym::AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> cV(kArrSize, 0.0);
	std::vector<double> tmpU(kArrSize, 0.0);
	std::vector<double> LRP(kArrSize, 0.0);
	std::vector<double> LRM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0);
	std::vector<double> LZM(kArrSize, 0.0);
	
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		tmpU[idx(i, j)] = (u[idx(i, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)] + u[idx(i + 1, j)]) * 0.25;

	std::vector<double> vecF_VR(kNr + 2 * kNumBCGrid, 0.0), vecF_VZ(kNz + 2 * kNumBCGrid);
	
	// V : X direction
	std::vector<double> FRP(kNr + 2 * kNumBCGrid, 0.0), FRM(kNr + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++) {
			vecF_VR[i] = v[idx(i, j)] * tmpU[idx(i, j)];
		}

		UnitHJWENO5(vecF_VR, FRP, FRM, kDr, kNr);

		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++) {
			LRP[idx(i, j)] = FRP[i];
			LRM[idx(i, j)] = FRM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FRP.begin(), FRP.end(), 0.0);
		std::fill(FRM.begin(), FRM.end(), 0.0);
		std::fill(vecF_VR.begin(), vecF_VR.end(), 0.0);
	}

	// V : Y direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++) {
		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
			vecF_VZ[j] = v[idx(i, j)] * v[idx(i, j)];
		}

		UnitHJWENO5(vecF_VZ, FZP, FZM, kDz, kNz);

		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
			LZP[idx(i, j)] = FZP[j];
			LZM[idx(i, j)] = FZM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_VZ.begin(), vecF_VZ.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0;

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		alphaX = tmpU[idx(i, j)];
		alphaY = v[idx(i, j)];
		
		cV[idx(i, j)]
			= (0.5 * (LRP[idx(i, j)] + LRM[idx(i, j)])
			+ 0.5 * (LZP[idx(i, j)] + LZM[idx(i, j)])
			- alphaX * (0.5 * (LRP[idx(i, j)] - LRM[idx(i, j)]))
			- alphaY * (0.5 * (LZP[idx(i, j)] - LZM[idx(i, j)])));

		assert(cV[idx(i, j)] == cV[idx(i, j)]);
		if (std::isnan(cV[idx(i, j)]) || std::isinf(cV[idx(i, j)])) {
			std::cout << "V-convection term nan/inf error : " << i << " " << j << " " << cV[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	return cV;
}

int MACSolver2DAxisym::UnitHJWENO5(
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

std::vector<double> MACSolver2DAxisym::AddExternalViscosityFU(const std::vector<double>& u, const std::vector<double>& v, 
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	// Kang, Myungjoo, Hyeseon Shim, and Stanley Osher.
	// "Level set based simulations of two-phase oil–water flows in pipes."
	// Journal of Scientific Computing 31.1-2 (2007): 153-184.
	std::vector<double> dU(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);

	if (kRe <= 0.0) {
		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
			dU[idx(i, j)] = 0.0;
		
		return dU;
	}

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

	// subcell
	double theta = 0.0, thetaH = 0.0;

	// need for updating rho
	/*
	-------------------------
	|			|			|
	|			|			|
	|			|			|
	---------lsUNHalf---------
	|			|			|
	|	lsW   lsUM   lsM    |
	|			|			|
	---------lsUSHalf--------
	|			|			|
	|			|			|
	|			|			|
	-------------------------
	*/
	/*
	-------------------------
	|			|			|
	|		  lsUN			|
	|			|			|
	----lsVW_N-------lsVN----
	|			|			|
  lsUW		  lsUM  (i,j)  lsUE
	|			|			|
	----lsVW---------lsVM----
	|			|			|
	|		  lsUS			|
	|			|			|
	-------------------------
	*/
	double lsW = 0.0, lsM = 0.0, lsUNHalf = 0.0, lsUSHalf = 0.0;
	double vW = 0.0, vW_N = 0.0, vM = 0.0, vN = 0.0;
	double muUNHalf = 0.0, muUSHalf = 0.0;

	double visR = 0.0, visZ = 0.0;
	
	double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;
	
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		visR = 0.0, visZ = 0.0;
		
		vW = v[idx(i - 1, j)];
		vW_N = v[idx(i - 1, j + 1)];
		vM = v[idx(i, j)];
		vN = v[idx(i, j + 1)];
		
		lsUNHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		lsUSHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];

		muUNHalf = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)]);
		muUSHalf = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)]);
		
		/*
		-------------------------
		|			|			|
		|			|			|
		|			|			|
		---------lsUNHalf---------
		|			|			|
		|	lsW   lsUM   lsM    |
		|			|			|
		---------lsUSHalf--------
		|			|			|
		|			|			|
		|			|			|
		-------------------------
		*/

		
		if (lsW >= 0 && lsM >= 0) {
			rhoEffWE = kRhoH;
			iRhoEffWE = 1.0 / kRhoH;
		}
		else if (lsW <= 0 && lsM <= 0) {
			rhoEffWE = kRhoL;
			iRhoEffWE = 1.0 / kRhoL;
		}
		else if (lsW >= 0 && lsM < 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsW < 0 && lsM >= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsM)|
			// |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffWE = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}

		if (lsUSHalf >= 0 && lsUNHalf >= 0) {
			rhoEffSN = kRhoH;
			iRhoEffSN = 1.0 / kRhoH;
		}
		else if (lsUSHalf <= 0 && lsUNHalf <= 0) {
			rhoEffSN = kRhoL;
			iRhoEffSN = 1.0 / kRhoL;
		}
		else if (lsUSHalf >= 0 && lsUNHalf < 0) {
			// interface lies between lsUSHalf and lsUNHalf
			thetaH = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			// |(lsUSHalf)| ===   High(+)  === |(interface)| ===     Low(-)      === |(lsUNHalf)|
			// |(lsUSHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsUSHalf < 0 && lsUNHalf >= 0) {
			// interface lies between lsUSHalf and lsUNHalf
			// |(lsUSHalf)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsUNHalf)|
			// |(lsUSHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
			thetaH = std::fabs(lsUNHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}

		visZ = (muUNHalf * (vN - vW_N) / kDr - muUSHalf * (vM - vW) / kDr) / kDz;
		
		dU[idx(i, j)] = visZ / rhoEffSN;
		
		if (std::isnan(dU[idx(i, j)]) || std::isinf(dU[idx(i, j)])) {
			std::cout << "U-viscosity term nan/inf error : " << i << " " << j << " " << dU[idx(i, j)] << std::endl;
			exit(1);
		}
		assert(dU[idx(i, j)] == dU[idx(i, j)]);
	}
	
	return dU;
}

std::vector<double> MACSolver2DAxisym::AddExternalViscosityFV(const std::vector<double>& u, const std::vector<double>& uhat,
	const std::vector<double>& v, const std::vector<double>& ls, const std::vector<double>& H) {
	
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);

	if (kRe <= 0.0) {
		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
			dV[idx(i, j)] = 0.0;

		return dV;
	}

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

	// level set
	double theta = 0.0, thetaH = 0.0;
	/*
	-------------------------------------
	|			|			|			|
	|			|    lsM	|			|
	|			|			|			|
	---------lsVWHalf----lsVEHalf--------
	|			|			|			|
	|			|	 lsS	|			|
	|			|			|			|
	-------------------------------------
	*/
	
	/*
	----------------lsVN-----------------
	|			|			|			|
	|		  lsUM  (i,j)  lsUE			|
	|			|			|			|
	-----lsVW-------lsVM--------lsVE-----
	|			|			|			|
	|		  lsUS		  lsUS_E		|
	|			|			|			|
	----------------lsVS-----------------
	*/
	double lsVWHalf = 0.0, lsVEHalf = 0.0, lsM = 0.0, lsS = 0.0;
	double uS = 0.0, uM = 0.0, uE = 0.0, uS_E = 0.0;
	double muVWHalf = 0.0, muVEHalf = 0.0;

	double visR = 0.0, visZ = 0.0;
	double rPM = 0.0;
	
	// effective Jump condition, effective v (vEff), effective mu (muEff), effective rho (rhoEff),
	double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		visR = 0.0, visZ = 0.0;
		rPM = (i + 0.5 - kNumBCGrid) * kDr;

		uM = u[idx(i, j)];
		uE = u[idx(i + 1, j)];
		uS = u[idx(i, j - 1)];
		uS_E = u[idx(i + 1, j - 1)];

		lsVWHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsVEHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i + 1, j)] + ls[idx(i + 1, j - 1)] + ls[idx(i, j - 1)]);
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];

		muVWHalf = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)]);
		muVEHalf = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)]);
		
		/*
		-------------------------------------
		|			|			|			|
		|			|    lsM	|			|
		|			|			|			|
		---------lsVWHalf----lsVEHalf--------
		|			|			|			|
		|			|	 lsS	|			|
		|			|			|			|
		-------------------------------------
		*/
		
		if (lsVWHalf >= 0 && lsVEHalf >= 0) {
			rhoEffWE = kRhoH;
			iRhoEffWE = 1.0 / kRhoH;
		}
		else if (lsVWHalf <= 0 && lsVEHalf <= 0) {
			rhoEffWE = kRhoL;
			iRhoEffWE = 1.0 / kRhoL;
		}
		else if (lsVWHalf >= 0 && lsVEHalf < 0) {
			// interface lies between lsVWHalf and lsVEHalf
			thetaH = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			// |(lsVWHalf)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsVEHalf)|
			// |(lsVWHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsVWHalf < 0 && lsVEHalf >= 0) {
			// interface lies between lsVWHalf and lsVEHalf
			thetaH = std::fabs(lsVEHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			// |(lsVWHalf)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsVEHalf)|
			// |(lsVWHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffWE = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}

		if (lsS >= 0 && lsM >= 0) {
			rhoEffSN = kRhoH;
			iRhoEffSN = 1.0 / kRhoH;
		}
		else if (lsS <= 0 && lsM <= 0) {
			rhoEffSN = kRhoL;
			iRhoEffSN = 1.0 / kRhoL;
		}
		else if (lsS >= 0 && lsM < 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsS < 0 && lsM >= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsM)|
			// |(lsS)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}

		visR = (muVEHalf * (uE - uS_E) / kDz - muVWHalf * (uM - uS) / kDz) / kDr;
		visZ = 0.5 * (mu[idx(i, j)] + mu[idx(i, j - 1)]) / rPM 
			* (0.5 * (uhat[idx(i, j)] + uhat[idx(i + 1, j)]) - 0.5 * (uhat[idx(i, j - 1)] + uhat[idx(i + 1, j - 1)])) / kDz;

		dV[idx(i, j)] = visR / rhoEffWE + visZ / rhoEffSN;
		
		assert(dV[idx(i, j)] == dV[idx(i, j)]);
		if (std::isnan(dV[idx(i, j)]) || std::isinf(dV[idx(i, j)])) {
			std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << dV[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	return dV;	
}

std::vector<double> MACSolver2DAxisym::AddGravityFU() {
	std::vector<double> gU(kArrSize, 0.0);

	if ((kFr != 0 || kG != 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis == GAXISENUM2DAXISYM::R) {
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			gU[idx(i, j)] = -kG;
		}
	}
	else {
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			gU[idx(i, j)] = 0.0;
		}
	}

	return gU;
}

std::vector<double> MACSolver2DAxisym::AddGravityFV() {
	std::vector<double> gV(kArrSize, 0.0);

	if ((kFr != 0 || kG != 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis == GAXISENUM2DAXISYM::Z) {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			gV[idx(i, j)] = -kG;
		}
	}
	else {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			gV[idx(i, j)] = 0.0;
		}
	}

	return gV;
}

std::vector<double> MACSolver2DAxisym::GetUHat(const std::vector<double>& ls, const std::vector<double>& u, 
	const std::vector<double>& rhs, const std::vector<double>& H) {

	const int64_t size = (kNr - 1) * kNz;
	std::vector<double> uhat(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);
	
	// u*, u**, u***
	std::vector<double> us1(size, 0.0), us2(size, 0.0), us3(size, 0.0);
	
	/*
	-------------------------
	|			|			|
	|			|			|
	|			|			|
	---------lsUNHalf---------
	|			|			|
	|	lsW   lsUM   lsM    |
	|			|			|
	---------lsUSHalf--------
	|			|			|
	|			|			|
	|			|			|
	-------------------------
	*/
	double lsW = 0.0, lsM = 0.0, lsUSHalf = 0.0, lsUNHalf = 0.0;
	double muW = 0.0, muM = 0.0, muUSHalf = 0.0, muUNHalf = 0.0;
	double rUM = 0.0, rPW = 0.0, rPM = 0.0;
	double rhoEffWE = 0.0, rhoEffSN = 0.0;
	// theta : lsW, lsE, lsS, lsN portion of grids, is used for linear value
	// thetaH : portion of high density, is used for averaging fluid properties
	double theta = 0.0, thetaH = 0.0;

	int64_t idxArr = 0;
	// Original RHS
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		idxArr = (j - kNumBCGrid) + kNz * (i - 1 - kNumBCGrid);
		
		// first rhs
		us3[idxArr] = rhs[idx(i, j)];
	}

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		// add boundary value due to interpolation
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];
	}

	// first, solve dr * dr factorization
	for (int j = 0; j < kNz; j++) {
		// used for tdma
		std::vector<double> a(kNr - 1, 0.0), b(kNr - 1, 0.0), c(kNr - 1, 0.0), d(kNr - 1, 0.0);

		for (int i = 0; i < kNr - 1; i++) {
			idxArr = j + kNz * i;

			d[i] = us3[idxArr];
		}

		for (int i = 0; i < kNr - 1; i++) {
			rUM = (i + 1.0) * kDr;
			rPW = (i + 0.5) * kDr;
			rPM = (i + 1.5) * kDr;
			
			// i (without boundary array) = i' + kNumBCGrid + 1 (with boundary array)
			// j (without boundary array) = j' + kNumBCGrid (with boundary array)
			lsW = ls[idx(i + kNumBCGrid, j + kNumBCGrid)];
			lsM = ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid)];
			muW = mu[idx(i + kNumBCGrid, j + kNumBCGrid)];
			muM = mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid)];

			if (lsW >= 0.0 && lsM >= 0.0)  {
				thetaH = 1.0; theta = 0.0;
			}
			else if (lsW <= 0.0 && lsM <= 0.0) {
				thetaH = 0.0; theta = 0.0;
			}
			else if (lsW >= 0.0 && lsM < 0.0) {
				// interface lies between u[i - 1, j] and u[i, j]
				// |(lsW)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsM)|
				// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
				theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			}
			else if (lsW < 0.0 && lsM >= 0.0) {
				// interface lies between u[i - 1, j] and u[i, j]
				// |(lsW)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsM)|
				// |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
				theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			}

			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);

			a[i] = -m_dt / (rhoEffWE * rUM) * (2.0 * rPW * muW) / (kDr * kDr);
			b[i] = 1.0 + m_dt / (rhoEffWE * rUM) * 2.0 * (rPW * muW + rPM * muM) / (kDr * kDr);
			c[i] = -m_dt / (rhoEffWE * rUM) * (2.0 * rPM * muM) / (kDr * kDr);

			// Boundary Condition
			if (i == 0 && (m_BC->m_BC_UW == BC2D::NEUMANN || m_BC->m_BC_UW == BC2D::OUTLET)) {
				b[i] += a[i];
				a[i] = 0.0;
			} 
			else if (i == 0 && (m_BC->m_BC_UW == BC2D::DIRICHLET || m_BC->m_BC_UW == BC2D::AXISYM ||
				m_BC->m_BC_UW == BC2D::INLET || m_BC->m_BC_UW == BC2D::WALL)) {
				d[i] -= a[i] * m_BC->m_BC_DirichletConstantUW;
				a[i] = 0.0;
			}
			else if (i == 0 && m_BC->m_BC_UW == BC2D::PERIODIC) {
				d[i] -= a[i] * us3[j + kNz * (kNr - 2)];
				a[i] = 0.0;
			}

			if (i == kNr - 2 && (m_BC->m_BC_UE == BC2D::NEUMANN || m_BC->m_BC_UE == BC2D::OUTLET)) {
				b[i] += c[i];
				c[i] = 0.0;
			}
			else if (i == kNr - 2 && (m_BC->m_BC_UE == BC2D::DIRICHLET || m_BC->m_BC_UE == BC2D::AXISYM ||
				m_BC->m_BC_UE == BC2D::INLET || m_BC->m_BC_UE == BC2D::WALL)) {
				d[i] -= c[i] * m_BC->m_BC_DirichletConstantUE;
				c[i] = 0.0;
			}
			else if (i == kNr - 2 && m_BC->m_BC_UE == BC2D::PERIODIC) {
				d[i] -= c[i] * us3[j + kNz * 0];
				c[i] = 0.0;
			}

		}

		SolveTDMA(a, b, c, d, kNr - 1);

		for (int i = 0; i < kNr - 1; i++) {
			idxArr = j + kNz * i;

			us2[idxArr] = d[i];
		}
	}

	// second, solve dz * dz factorization
	for (int i = 0; i < kNr - 1; i++) {
		std::vector<double> a(kNz, 0.0), b(kNz, 0.0), c(kNz, 0.0), d(kNz, 0.0);

		for (int j = 0; j < kNz; j++) {
			idxArr = j + kNz * i;

			d[j] = us2[idxArr];
		}

		for (int j = 0; j < kNz; j++) {
			// i (without boundary array) = i' + kNumBCGrid + 1 (with boundary array)
			// j (without boundary array) = j' + kNumBCGrid (with boundary array)
			lsUSHalf = 0.25 * (ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid)] + ls[idx(i + kNumBCGrid, j + kNumBCGrid)]
				+ ls[idx(i + kNumBCGrid, j + kNumBCGrid - 1)] + ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid - 1)]);
			lsUNHalf = 0.25 * (ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid)] + ls[idx(i + kNumBCGrid, j + kNumBCGrid)]
				+ ls[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid + 1)]);
			muUSHalf = 0.25 * (mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid)] + mu[idx(i + kNumBCGrid, j + kNumBCGrid)]
				+ mu[idx(i + kNumBCGrid, j + kNumBCGrid - 1)] + mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid - 1)]);
			muUNHalf = 0.25 * (mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid)] + mu[idx(i + kNumBCGrid, j + kNumBCGrid)]
				+ mu[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid + 1)]);

			if (lsUSHalf >= 0 && lsUNHalf >= 0) {
				thetaH = 1.0; theta = 0.0;
			}
			else if (lsUSHalf <= 0 && lsUNHalf <= 0) {
				thetaH = 0.0; theta = 0.0;
			}
			else if (lsUSHalf >= 0 && lsUNHalf < 0) {
				// interface lies between lsUSHalf and lsUNHalf
				// |(lsUSHalf)| ===   High(+)  === |(interface)| ===     Low(-)      === |(lsUNHalf)|
				// |(lsUSHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
				thetaH = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
				theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			}
			else if (lsUSHalf < 0 && lsUNHalf >= 0) {
				// interface lies between lsUSHalf and lsUNHalf
				// |(lsUSHalf)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsUNHalf)|
				// |(lsUSHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
				thetaH = std::fabs(lsUNHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
				theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			}

			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);

			a[j] = -m_dt / rhoEffSN * muUSHalf / (kDz * kDz);
			b[j] = 1.0 + m_dt / rhoEffSN * (muUSHalf + muUNHalf) / (kDz * kDz);
			c[j] = -m_dt / rhoEffSN * muUNHalf / (kDz * kDz);
			
			if (j == 0 && (m_BC->m_BC_US == BC2D::NEUMANN ||
				m_BC->m_BC_US == BC2D::AXISYM || m_BC->m_BC_US == BC2D::OUTLET)) {
				b[j] += a[j];
				a[j] = 0.0;
			}
			else if (j == 0 && (m_BC->m_BC_US == BC2D::DIRICHLET ||
				m_BC->m_BC_US == BC2D::INLET || m_BC->m_BC_US == BC2D::WALL)) {
				d[j] -= a[j] * 2.0 * m_BC->m_BC_DirichletConstantUS;
				a[j] = 0.0;
			}
			else if (j == 0 && m_BC->m_BC_US == BC2D::PERIODIC) {
				d[j] -= a[j] * us2[(kNz - 1) + kNz * i];
				a[j] = 0.0;
			}

			if (j == kNz - 1 && (m_BC->m_BC_UN == BC2D::NEUMANN ||
				m_BC->m_BC_UN == BC2D::AXISYM || m_BC->m_BC_UN == BC2D::OUTLET)) {
				b[j] += c[j];
				c[j] = 0.0;
			}
			else if (j == kNz - 1 && (m_BC->m_BC_UN == BC2D::DIRICHLET ||
				m_BC->m_BC_UN == BC2D::INLET || m_BC->m_BC_UN == BC2D::WALL)) {
				d[j] -= c[j] * 2.0 * m_BC->m_BC_DirichletConstantUN;
				c[j] = 0.0;
			}
			else if (j == kNz - 1 && m_BC->m_BC_UN == BC2D::PERIODIC) {
				d[j] -= c[j] * us2[0 + kNz * i];
				c[j] = 0.0;
			}
		}

		SolveTDMA(a, b, c, d, kNz);

		for (int j = 0; j < kNz; j++) {
			idxArr = j + kNz * i;

			us1[idxArr] = d[j];
		}
	}

	// third, just inverse r * r part
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rUM = (i - kNumBCGrid) * kDr;
		idxArr = (j - kNumBCGrid) + kNz * (i - kNumBCGrid - 1);

		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];

		if (lsW >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsW <= 0.0 && lsM <= 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsW >= 0.0 && lsM < 0.0) {
			// interface lies between u[i - 1, j] and u[i, j]
			// |(lsW)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}
		else if (lsW < 0.0 && lsM >= 0.0) {
			// interface lies between u[i - 1, j] and u[i, j]
			// |(lsW)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsM)|
			// |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}

		rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);

		uhat[idx(i, j)] = us1[idxArr] /
			(1.0 + m_dt / rhoEffWE * (mu[idx(i - 1, j)] + mu[idx(i, j)]) / (rUM * rUM));
	}
	
	// check nan & inf
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		if (std::isnan(uhat[idx(i, j)]) || std::isinf(uhat[idx(i, j)])) {
			std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << uhat[idx(i, j)] 
			<< " isNan? " << isnan(uhat[idx(i, j)]) << " isInf? " << isinf(uhat[idx(i, j)]) << std::endl;
			exit(1);
		}
	}

	return uhat;
}

std::vector<double> MACSolver2DAxisym::GetVHat(const std::vector<double>& ls, const std::vector<double>& v,
	const std::vector<double>& rhs, const std::vector<double>& H) {
	
	// use factorization
	const int64_t size = kNr * (kNz - 1);
	std::vector<double> vhat(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);

	// v*, v**
	std::vector<double> vs0(size, 0.0), vs1(size, 0.0), vs2(size, 0.0);

	/*
	-------------------------------------
	|			|			|			|
	|			|    lsM	|			|
	|			|			|			|
	---------lsVWHalf----lsVEHalf--------
	|			|			|			|
	|			|	 lsS	|			|
	|			|			|			|
	-------------------------------------
	*/

	double lsS = 0.0, lsM = 0.0, lsVWHalf = 0.0, lsVEHalf = 0.0;
	double muS = 0.0, muM = 0.0, muVWHalf = 0.0, muVEHalf = 0.0;
	double rVM = 0.0, rVWHalf = 0.0, rVEHalf = 0.0;
	double rhoEffWE = 0.0, rhoEffSN = 0.0;
	// theta : lsW, lsE, lsS, lsN portion of grids, is used for linear value
	// thetaH : portion of high density, is used for averaging fluid properties
	double theta = 0.0, thetaH = 0.0;

	int64_t idxArr = 0;
		
	// Original RHS
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		idxArr = (j - kNumBCGrid - 1) + (kNz - 1) * (i - kNumBCGrid);

		// first rhs
		vs2[idxArr] = rhs[idx(i, j)];
	}

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		// add boundary value due to interpolation
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];
	}

	for (int j = 0; j < kNz - 1; j++) {
		// used for tdma
		std::vector<double> a(kNr, 0.0), b(kNr, 0.0), c(kNr, 0.0), d(kNr, 0.0);

		for (int i = 0; i < kNr; i++) {
			idxArr = j + (kNz - 1) * i;

			d[i] = vs2[idxArr];
		}

		for (int i = 0; i < kNr; i++) {
			rVM = (i + 0.5) * kDr;
			rVWHalf = i * kDr;
			rVEHalf = (i + 1.0) * kDr;
			
			// i (without boundary array) = i' + kNumBCGrid (with boundary array)
			// j (without boundary array) = j' + kNumBCGrid + 1 (with boundary array)

			lsVWHalf = 0.25 * (ls[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + ls[idx(i + kNumBCGrid - 1, j + kNumBCGrid + 1)]
				+ ls[idx(i + kNumBCGrid, j + kNumBCGrid)] + ls[idx(i + kNumBCGrid - 1, j + kNumBCGrid)]);
			lsVEHalf = 0.25 * (ls[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid + 1)]
				+ ls[idx(i + kNumBCGrid, j + kNumBCGrid)] + ls[idx(i + kNumBCGrid + 1, j + kNumBCGrid)]);
			muVWHalf = 0.25 * (mu[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + mu[idx(i + kNumBCGrid - 1, j + kNumBCGrid + 1)]
				+ mu[idx(i + kNumBCGrid, j + kNumBCGrid)] + mu[idx(i + kNumBCGrid - 1, j + kNumBCGrid)]);
			muVEHalf = 0.25 * (mu[idx(i + kNumBCGrid, j + kNumBCGrid + 1)] + mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid + 1)]
				+ mu[idx(i + kNumBCGrid, j + kNumBCGrid)] + mu[idx(i + kNumBCGrid + 1, j + kNumBCGrid)]);

			if (lsVWHalf >= 0.0 && lsVEHalf >= 0.0)  {
				thetaH = 1.0; theta = 0.0;
			}
			else if (lsVWHalf <= 0.0 && lsVEHalf <= 0.0) {
				thetaH = 0.0; theta = 0.0;
			}
			else if (lsVWHalf >= 0.0 && lsVEHalf < 0.0) {
				// interface lies between lsVWHalf and lsVEHalf
				// |(lsVWHalf)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsVEHalf)|
				// |(lsVWHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
				thetaH = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
				theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			}
			else if (lsVWHalf < 0.0 && lsVEHalf >= 0.0) {
				// interface lies between lsVWHalf and lsVEHalf
				// |(lsVWHalf)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsVEHalf)|
				// |(lsVWHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
				thetaH = std::fabs(lsVEHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
				theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			}

			rhoEffWE = kRhoH * thetaH + kRhoL * (1.0 - thetaH);

			a[i] = -m_dt / (rhoEffWE * rVM) * rVWHalf * muVWHalf / (kDr * kDr);
			b[i] = 1.0 + m_dt / (rhoEffWE * rVM) * (rVWHalf * muVWHalf + rVEHalf * muVEHalf) / (kDr * kDr);
			c[i] = -m_dt / (rhoEffWE * rVM) * rVEHalf * muVEHalf / (kDr * kDr);

			// Boundary Condition
			if (i == 0 && (m_BC->m_BC_VW == BC2D::NEUMANN ||
				m_BC->m_BC_VW == BC2D::AXISYM || m_BC->m_BC_VW == BC2D::OUTLET)) {
				b[i] += a[i];
				a[i] = 0.0;
			}
			else if (i == 0 && (m_BC->m_BC_VW == BC2D::DIRICHLET ||
				m_BC->m_BC_VW == BC2D::INLET || m_BC->m_BC_VW == BC2D::WALL)) {
				d[i] -= a[i] * 2.0 * m_BC->m_BC_DirichletConstantVW;
				a[i] = 0.0;
			}
			else if (i == 0 && m_BC->m_BC_VW == BC2D::PERIODIC) {
				d[i] -= a[i] * vs2[j + (kNz - 1) * (kNr - 1)];
				a[i] = 0.0;
			}

			if (i == kNr - 1 && (m_BC->m_BC_VE == BC2D::NEUMANN ||
				m_BC->m_BC_VE == BC2D::AXISYM || m_BC->m_BC_VE == BC2D::OUTLET)) {
				b[i] += c[i];
				c[i] = 0.0;
			}
			else if (i == kNr - 1 && (m_BC->m_BC_VE == BC2D::DIRICHLET ||
				m_BC->m_BC_VE == BC2D::INLET || m_BC->m_BC_VE == BC2D::WALL)) {
				d[i] -= c[i] * 2.0 * m_BC->m_BC_DirichletConstantVE;
				c[i] = 0.0;
			}
			else if (i == kNr - 1 && m_BC->m_BC_VE == BC2D::PERIODIC) {
				d[i] -= c[i] * vs2[j + (kNz - 1) * 0];
				c[i] = 0.0;
			}

		}

		SolveTDMA(a, b, c, d, kNr);

		for (int i = 0; i < kNr; i++) {
			idxArr = j + (kNz - 1) * i;

			vs1[idxArr] = d[i];
		}
	}

	// second, solve dz * dz factorization
	for (int i = 0; i < kNr; i++) {
		std::vector<double> a(kNz, 0.0), b(kNz, 0.0), c(kNz, 0.0), d(kNz, 0.0);

		for (int j = 0; j < kNz - 1; j++) {
			idxArr = j + (kNz - 1) * i;

			d[j] = vs1[idxArr];
		}

		for (int j = 0; j < kNz - 1; j++) {
			// i (without boundary array) = i' + kNumBCGrid (with boundary array)
			// j (without boundary array) = j' + kNumBCGrid + 1 (with boundary array)

			lsS = ls[idx(i + kNumBCGrid, j + kNumBCGrid)];
			lsM = ls[idx(i + kNumBCGrid, j + kNumBCGrid + 1)];
			muS = mu[idx(i + kNumBCGrid, j + kNumBCGrid)];
			muM = mu[idx(i + kNumBCGrid, j + kNumBCGrid + 1)];

			if (lsS >= 0.0 && lsM >= 0.0)  {
				thetaH = 1.0; theta = 0.0;
			}
			else if (lsS <= 0.0 && lsM <= 0.0) {
				thetaH = 0.0; theta = 0.0;
			}
			else if (lsS >= 0.0 && lsM < 0.0) {
				// interface lies between ls[i, j - 1] and ls[i, j]
				// |(lsS)| ===   High(+)  === |(interface)| ===      Low(-)     === |(lsM)|
				// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
				theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			}
			else if (lsS < 0.0 && lsM >= 0.0) {
				// interface lies between ls[i, j - 1] and ls[i, j]
				// |(lsS)| ===    Low(-)   === |(interface)| ===     High(+)     === |(lsM)|
				// |(lsS)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
				thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
				theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			}

			rhoEffSN = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
			
			a[j] = -m_dt / rhoEffSN * (2.0 * muS) / (kDz * kDz);
			b[j] = 1.0 + m_dt / rhoEffSN * 2.0 * (muS + muM) / (kDz * kDz);
			c[j] = -m_dt / rhoEffSN * (2.0 * muM) / (kDz * kDz);

			if (j == 0 && (m_BC->m_BC_VS == BC2D::NEUMANN || m_BC->m_BC_VS == BC2D::OUTLET)) {
				b[j] += a[j];
				a[j] = 0.0;
			}
			else if (j == 0 && (m_BC->m_BC_VS == BC2D::DIRICHLET || m_BC->m_BC_VS == BC2D::AXISYM || 
				m_BC->m_BC_VS == BC2D::INLET || m_BC->m_BC_VS == BC2D::WALL)) {
				d[j] -= a[j] * m_BC->m_BC_DirichletConstantVS;
				a[j] = 0.0;
			}
			else if (j == 0 && m_BC->m_BC_VS == BC2D::PERIODIC) {
				d[j] -= a[j] * vs1[(kNz - 2) + (kNz - 1) * i];
				a[j] = 0.0;
			}

			if (j == kNz - 2 && (m_BC->m_BC_VN == BC2D::NEUMANN || m_BC->m_BC_VN == BC2D::OUTLET)) {
				b[j] += c[j];
				c[j] = 0.0;
			}
			else if (j == kNz - 2 && (m_BC->m_BC_VN == BC2D::DIRICHLET || m_BC->m_BC_VN == BC2D::AXISYM ||
				m_BC->m_BC_VN == BC2D::INLET || m_BC->m_BC_VN == BC2D::WALL)) {
				d[j] -= c[j] * m_BC->m_BC_DirichletConstantVN;
				c[j] = 0.0;
			}
			else if (j == kNz - 2 && m_BC->m_BC_VN == BC2D::PERIODIC) {
				d[j] -= c[j] * vs1[0 + (kNz - 1) * i];
				c[j] = 0.0;
			}
		}

		SolveTDMA(a, b, c, d, kNz - 1);

		for (int j = 0; j < kNz - 1; j++) {
			idxArr = j + (kNz - 1) * i;

			vs0[idxArr] = d[j];
		}
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		idxArr = (j - kNumBCGrid - 1) + (kNz - 1) * (i - kNumBCGrid);

		vhat[idx(i, j)] = vs0[idxArr];
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		if (std::isnan(vhat[idx(i, j)]) || std::isinf(vhat[idx(i, j)])) {
			std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << vhat[idx(i, j)]
				<< " isNan? " << isnan(vhat[idx(i, j)]) << " isInf? " << isinf(vhat[idx(i, j)]) << std::endl;
			exit(1);
		}
	}

	return vhat;
}

int MACSolver2DAxisym::SetPoissonSolver(POISSONTYPE type) {
	m_PoissonSolverType = type;
	if (!m_Poisson)
		m_Poisson = std::make_shared<PoissonSolver2D>(kNr, kNz, kNumBCGrid);

	return 0;
}

int MACSolver2DAxisym::SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
	const std::vector<double>& ls, const std::vector<double>& lsB,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& H, const int maxIter) {
	if (!m_Poisson) {
		perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
	}
	std::vector<double> rhs(kArrSize, 0.0);
	std::vector<double> rhoW(kArrSize, 0.0),
		rhoE(kArrSize, 0.0),
		rhoS(kArrSize, 0.0),
		rhoN(kArrSize, 0.0);
	std::vector<double> pCoefW(kArrSize, 0.0),
		pCoefE(kArrSize, 0.0),
		pCoefS(kArrSize, 0.0),
		pCoefN(kArrSize, 0.0);
	std::vector<double> FWVec(kArrSize, 0.0),
		FEVec(kArrSize, 0.0),
		FSVec(kArrSize, 0.0),
		FNVec(kArrSize, 0.0);
	std::vector<double> aWVec(kArrSize, 0.0),
		aEVec(kArrSize, 0.0),
		aCVec(kArrSize, 0.0),
		aSVec(kArrSize, 0.0),
		aNVec(kArrSize, 0.0);

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
	
	// theta : lsW, lsE, lsS, lsN portion of grids, is used for linear value
	// thetaH : portion of high density, is used for averaging fluid properties
	double thetaH = 0.0, theta = 0.0;
	// effective values
	double kappaEff = 0.0, rhoEff = 0.0;
	// r axis value
	double rEff = 0.0, rPW = 0.0, rPM = 0.0, rPE = 0.0, rPWHalf = 0.0, rPEHalf = 0.0;
	
	if (kWe != 0.0 && !isnan(kWe) && !isinf(kWe)) {
		UpdateKappa(ls);
		ApplyBC_P_2D(m_kappa);
	}
	// A Matrix is (nr * ny * nz) X (nr * ny * nz) matrix, which is very very huge. hence use sparse blas
	std::vector<double> AVals, DiagVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNr * kNz, Ancols = kNr * kNz;
	MKL_INT size = kNr * kNz;

	// stored coef for A matrix, Dictionary but it is ordered
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		ps[idx(i, j)] = 0.0;
		rhs[idx(i, j)] = 0.0;
	}

	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++) {
		lsW = ls[idx(i - 1, j)];
		lsE = ls[idx(i + 1, j)];
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		lsN = ls[idx(i, j + 1)];
			
		FW = 0.0;
		FE = 0.0;
		FS = 0.0;
		FN = 0.0;
		rEff = 0.0;
		rPW = kBaseR + (i - 0.5 - kNumBCGrid) * kDr;
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rPE = kBaseR + (i + 1.5 - kNumBCGrid) * kDr;
		rPWHalf = kBaseR + (i - kNumBCGrid) * kDr;
		rPEHalf = kBaseR + (i + 1.0 - kNumBCGrid) * kDr;
		
		// [a]_\Gamma = a^+ - a^- = a^inside - a^outside
		// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// (P_M - P_W)/kDr appears
		// if P_M == P_+, P_W == P_-, a^+ means all terms related to P_+, P_W changed to P_M related terms
		// P_W = P_M -(2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa)

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsW >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
			rEff = rPWHalf;
		}
		else if (lsW < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
			rEff = rPWHalf;
		}
		else if (lsW >= 0.0 && lsM < 0.0) {
			thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));	
			rEff = rPM * theta + rPW * (1.0 - theta);
			rEff = std::max(rEff, 0.0);
		}
		else if (lsW < 0.0 && lsM >= 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			rEff = rPM * theta + rPW * (1.0 - theta);
			rEff = std::max(rEff, 0.0);
		}
		
		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);

		// coefficient
		pCoefW[idx(i, j)] = rEff / rhoEff;
		rhoW[idx(i, j)] = rhoEff;

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// + level set(high) : H = 1, - level set(low) : H = 0;
		// Bubble : (+ls)-kRhoH-Liquid, (-ls)-kRhoL-Gas
		FW = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		FW *= pCoefW[idx(i, j)] / (kDr * kDr);
		aWVec[idx(i, j)] = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid cell adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsE >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
			rEff = rPEHalf;
		}
		else if (lsM < 0.0 && lsE < 0.0) {
			thetaH = 0.0; theta = 0.0;
			rEff = rPEHalf;
		}
		else if (lsM >= 0.0 && lsE < 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsM) + std::fabs(lsE));
			theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
			rEff = rPM * theta + rPE * (1.0 - theta);
		}
		else if (lsM < 0.0 && lsE >= 0.0) {
			thetaH = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
			theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
			rEff = rPM * theta + rPE * (1.0 - theta);
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i + 1, j)] * (1.0 - theta);

		// coefficient
		pCoefE[idx(i, j)] = rEff / rhoEff;
		rhoE[idx(i, j)] = rhoEff;
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FE = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);
		FE *= pCoefE[idx(i, j)] / (kDr * kDr);
		aEVec[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsS >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsS < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsS >= 0.0 && lsM < 0.0) {
			thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}
		else if (lsS < 0.0 && lsM >= 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);

		// coefficient
		pCoefS[idx(i, j)] = rPM / rhoEff;
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// + level set(high) : H = 1, - level set(low) : H = 0;
		// Bubble : (+ls)-kRhoH-Liquid, (-ls)-kRhoL-Gas
		FS = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
		FS *= pCoefS[idx(i, j)] / (kDz * kDz);
		aSVec[idx(i, j)] = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsN >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsM < 0.0 && lsN < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsM >= 0.0 && lsN < 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsM) + std::fabs(lsN));
			theta = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
		}
		else if (lsM < 0.0 && lsN >= 0.0) {
			thetaH = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
			theta = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j + 1)] * (1.0 - theta);

		// coefficient
		pCoefN[idx(i, j)] = rPM / rhoEff;
		
		// pressure jump = High - low = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// + level set(high) : H = 1, - level set(low) : H = 0;
		// Bubble : (+ls)-kRhoH-Liquid, (-ls)-kRhoL-Gas
		FN = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1)] - H[idx(i, j)]);
		FN *= pCoefN[idx(i, j)] / (kDz * kDz);
		aNVec[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1)] - H[idx(i, j)]);

		// poisson equation form should be -\beta \nabla p = f
		// pCoefEff has already negative value of rhoEff, then it is just a coefficient.
		// For discretization of pressure gradient, additional term is negative and it goes to RHS
		// Then it is a positive value and don't worrry about the sign
		rhs[idx(i, j)] -= FW + FE + FS + FN;
		FWVec[idx(i, j)] = FW;
		FEVec[idx(i, j)] = FE;
		FSVec[idx(i, j)] = FS;
		FNVec[idx(i, j)] = FN;
		
		// std::cout << i << " " << j << " " << pCoefS[idx(i, j)] << " "	 << pCoefW[idx(i, j)] << " " << pCoefE[idx(i, j)] << " " << pCoefN[idx(i, j)] << std::endl;;
		assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
		if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " 
				<< rhs[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	// Original value of RHS
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		rhs[idx(i, j)] -= div[idx(i, j)];

	// An order of A matrix coef. is very important, hence reverse j order
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++) {
		// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
		// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
		// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
		AValsDic.clear();
		AColsDic.clear();
		tmpRowIdx = 0;
		tmpMRowIdx = 0;
		if ((m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM || m_BC->m_BC_PW == BC2D::WALL) &&
			(m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM || m_BC->m_BC_PE == BC2D::WALL) &&
			(m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM || m_BC->m_BC_PS == BC2D::WALL) &&
			(m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM || m_BC->m_BC_PN == BC2D::WALL) &&
			(i == kNumBCGrid && j == kNumBCGrid))
			continue;

		// Add starting rowIdx
		ARowIdx.push_back(rowIdx);
		DiagRowIdx.push_back(MRowIdx);

		// to make positive definite matrix 
		pCoefS[idx(i, j)] *= -1.0;
		pCoefW[idx(i, j)] *= -1.0;
		pCoefE[idx(i, j)] *= -1.0;
		pCoefN[idx(i, j)] *= -1.0;

		// Set default values, if a current pointer is in interior, it will not be changed.
		// matrix is symmetric
		AValsDic["S"] = pCoefS[idx(i, j)] / (kDz * kDz);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNr;
		AValsDic["W"] = pCoefW[idx(i, j)] / (kDr * kDr);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["C"] = -(pCoefW[idx(i, j)] + pCoefE[idx(i, j)]) / (kDr * kDr)
			- (pCoefS[idx(i, j)] + pCoefN[idx(i, j)]) / (kDz * kDz);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["E"] = pCoefE[idx(i, j)] / (kDr * kDr);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["N"] = pCoefN[idx(i, j)] / (kDz * kDz);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNr;
		
		if (i == kNumBCGrid && (m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM)) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += pCoefW[idx(i, j)] / (kDr * kDr);
		}
		else if (i == kNumBCGrid 
			&& (m_BC->m_BC_PW == BC2D::DIRICHLET || m_BC->m_BC_PW == BC2D::INLET
			|| m_BC->m_BC_PW == BC2D::OUTLET || m_BC->m_BC_PW == BC2D::PRESSURE
			|| m_BC->m_BC_PW == BC2D::WALL)) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			rhs[idx(i, j)] -= pCoefW[idx(i, j)] / (kDr * kDr) 
				* (-ps[idx(i, j)] + 2.0 * m_BC->m_BC_DirichletConstantPW);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::PERIODIC) {
			AValsDic["W"] = pCoefW[idx(kNumBCGrid + kNr - 1, j)];
		}

		// East boundary
		if (i == kNumBCGrid + kNr - 1 && (m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM)) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += pCoefE[idx(i, j)] / (kDr * kDr);
		}
		else if (i == kNumBCGrid + kNr - 1
			&& (m_BC->m_BC_PE == BC2D::DIRICHLET || m_BC->m_BC_PE == BC2D::INLET
			|| m_BC->m_BC_PE == BC2D::OUTLET || m_BC->m_BC_PE == BC2D::PRESSURE
			|| m_BC->m_BC_PE == BC2D::WALL)) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			rhs[idx(i, j)] -= pCoefE[idx(i, j)] / (kDr * kDr) 
				* (-ps[idx(i, j)] + 2.0 * m_BC->m_BC_DirichletConstantPE);
		}
		else if (i == kNumBCGrid + kNr - 1 && m_BC->m_BC_PE == BC2D::PERIODIC) {
			AValsDic["E"] = pCoefE[idx(kNumBCGrid, j)];
		}

		if (j == kNumBCGrid && (m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM)) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += pCoefS[idx(i, j)] / (kDz * kDz);
		}
		else if (j == kNumBCGrid
			&& (m_BC->m_BC_PS == BC2D::DIRICHLET || m_BC->m_BC_PS == BC2D::INLET
			|| m_BC->m_BC_PS == BC2D::OUTLET || m_BC->m_BC_PS == BC2D::PRESSURE
			|| m_BC->m_BC_PS == BC2D::WALL)) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			rhs[idx(i, j)] -= pCoefS[idx(i, j)] / (kDz * kDz) 
				* (-ps[idx(i, j)] + 2.0 * m_BC->m_BC_DirichletConstantPS);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::PERIODIC) {
			AValsDic["S"] = pCoefS[idx(i, kNumBCGrid + kNz - 1)];
		}

		if (j == kNumBCGrid + kNz - 1 && (m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM)) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += pCoefN[idx(i, j)] / (kDz * kDz);
		}
		else if (j == kNumBCGrid + kNz - 1
			&& (m_BC->m_BC_PN == BC2D::DIRICHLET || m_BC->m_BC_PN == BC2D::INLET
			|| m_BC->m_BC_PN == BC2D::OUTLET || m_BC->m_BC_PN == BC2D::PRESSURE
			|| m_BC->m_BC_PN == BC2D::WALL)) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			rhs[idx(i, j)] -= pCoefN[idx(i, j)] / (kDz * kDz) 
				* (-ps[idx(i, j)] + 2.0 * m_BC->m_BC_DirichletConstantPN);
		}
		else if (j == kNumBCGrid + kNz - 1 && m_BC->m_BC_PN == BC2D::PERIODIC) {
			AValsDic["N"] = pCoefN[idx(i, kNumBCGrid)];
		}
		// add non zero values to AVals and ACols
		// KEEP ORDER OF PUSH_BACK!!

		// fix ps[idx(kNumBCGrid, kNumBCGrid)] = 0.0 if all boundary conditions are neumann condition like. 
		if ((m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM || m_BC->m_BC_PW == BC2D::WALL) && 
			(m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM || m_BC->m_BC_PE == BC2D::WALL) && 
			(m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM || m_BC->m_BC_PS == BC2D::WALL) &&
			(m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM || m_BC->m_BC_PN == BC2D::WALL)) {
			if (i == kNumBCGrid + 1 && j == kNumBCGrid) {
				AColsDic["W"] = -1;
				AValsDic["W"] = 0.0;
				rhs[idx(i, j)] -= pCoefW[idx(i, j)] / (kDr * kDr) * (-ps[idx(i, j)]);
			}
			else if (i == kNumBCGrid && j == kNumBCGrid + 1) {
				AColsDic["S"] = -1;
				AValsDic["S"] = 0.0;
				rhs[idx(i, j)] -= pCoefS[idx(i, j)] / (kDz * kDz) * (-ps[idx(i, j)]);
			}

			AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNr - 1;
			AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr - 1; 
			AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr - 1;
			AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNr - 1;

			size = kNr * kNz - 1;
		}

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
		
		tmpMRowIdx++;
		DiagVals.push_back(AValsDic["C"]);
		DiagCols.push_back(AColsDic["C"]);
		
		rowIdx += tmpRowIdx;
		MRowIdx += tmpMRowIdx;

		// std::cout << i << " " << j << " " << AValsDic["S"] << " " << AValsDic["W"] << " " << AValsDic["C"] << " " << AValsDic["E"] << " " << AValsDic["N"] << " " << std::endl;
		if (m_PoissonSolverType == POISSONTYPE::CG && 
			(std::fabs(AValsDic["S"]) + std::fabs(AValsDic["W"]) + std::fabs(AValsDic["E"]) + std::fabs(AValsDic["N"]) - std::fabs(AValsDic["C"]) > 1.0e-6)) {
			std::cout << "Poisson equation Error : Conjugate Gradient Method Failed, Not Diagonally Dominant Matrix "
				<< std::fabs(AValsDic["S"]) + std::fabs(AValsDic["W"]) + std::fabs(AValsDic["E"]) + std::fabs(AValsDic["N"])
				<< " " << std::fabs(AValsDic["C"])
				<< " " << std::fabs(AValsDic["S"]) + std::fabs(AValsDic["W"]) + std::fabs(AValsDic["E"]) + std::fabs(AValsDic["N"]) - std::fabs(AValsDic["C"]) << std::endl;
			exit(1);
		}

 		assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
		if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
			std::cout << "Poisson equation Error : right hand side of poisson equation nan/inf error : " 
				<< i << " " << j << " " << rhs[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);
	
	if (m_PoissonSolverType == POISSONTYPE::CG) {
		// std::cout << "Poisson : CG" << std::endl;
		m_Poisson->CG_2FUniformP_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, m_BC, size, maxIter);
	
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
		// std::cout << "Poisson : BiCG" << std::endl;
		m_Poisson->BiCGStab_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, m_BC, size, maxIter);

		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}

	if ((m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM || m_BC->m_BC_PW == BC2D::WALL) &&
		(m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM || m_BC->m_BC_PE == BC2D::WALL) &&
		(m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM || m_BC->m_BC_PS == BC2D::WALL) &&
		(m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM || m_BC->m_BC_PN == BC2D::WALL))
			ps[idx(kNumBCGrid, kNumBCGrid)] = 0.0;

	ApplyBC_P_2D(ps);

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		m_p[idx(i, j)] = ps[idx(i, j)] / m_dt;

	std::ofstream outF;
	std::string fname("P_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"R\", \"Z\", \"FW\",\"AW\",\"CoefW\",\"FE\",\"AE\",\"CoefE\", \"FS\",\"AS\", \"CoefS\",\"FN\", \"AN\",\"CoefN\", \"F\", \"LS\", \"PS\", \"HatDIV\", \"RHS\",\"KAPPA\"" << std::endl;
		outF.close();
	}
	if (m_iter % kNIterSkip == 0) {
	outF.open(fname.c_str(), std::ios::app);
	
	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNr + 6 << std::string(", J=") << kNz + 6
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	// for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	// for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefW[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FEVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aEVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefE[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefS[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefN[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)] + FEVec[idx(i, j)] + FSVec[idx(i, j)] + FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(div[idx(i, j)]) << std::string(",")
			<< static_cast<double>(rhs[idx(i, j)]) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;

	outF.close();
	}
	return 0;
}

std::vector<double> MACSolver2DAxisym::GetDivergence4Poisson(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> div(kArrSize, 0.0);
	double rPM = 0.0, rPEHalf = 0.0, rPWHalf = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rPEHalf = kBaseR + (i + 1.0 - kNumBCGrid) * kDr;
		rPWHalf = kBaseR + (i - kNumBCGrid) * kDr;

		div[idx(i, j)] = (rPEHalf * u[idx(i + 1, j)] - rPWHalf * u[idx(i, j)]) / kDr
			+ rPM * (v[idx(i, j + 1)] - v[idx(i, j)]) / kDz;
	}

	return div;
}

std::vector<double> MACSolver2DAxisym::GetDivergence(const std::vector<double>& u, const std::vector<double>& v) {
	std::vector<double> div(kArrSize, 0.0);
	double rPM = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		
		div[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]) / rPM + (u[idx(i + 1, j)] - u[idx(i, j)]) / kDr
			+ (v[idx(i, j + 1)] - v[idx(i, j)]) / kDz;
	}
	
	return div;
}

int MACSolver2DAxisym::UpdateVel(std::vector<double>& u, std::vector<double>& v,
	const std::vector<double>& us, const std::vector<double>& vs,
	const std::vector<double>& ps, const std::vector<double>& ls, const std::vector<double>& lsB,
	const std::vector<double>& H) {
	
	// velocity update after solving poisson equation
	// ps = p * dt
	double lsW = 0.0, lsM = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0;
	double rhoEff = 0.0, theta = 0.0, thetaH = 0.0, pCoefEff = 0.0, kappaEff = 0.0;
	const double eps = 1.0e-100;
		
	std::vector<double> aWVec(kArrSize, 0.0),
		aSVec(kArrSize, 0.0);
	std::vector<double> AWW(kArrSize, 0.0),
		AWM(kArrSize, 0.0),
		ASS(kArrSize, 0.0),
		ASM(kArrSize, 0.0);
	std::vector<double> pCoefEffWVec(kArrSize, 0.0),
		pCoefEffSVec(kArrSize, 0.0);

	std::vector<double> aWEff(kArrSize, 0.0),
		aSEff(kArrSize, 0.0);
	
	/*
	// http://ctr.stanford.edu/Summer/SP08/3_1_Moureau.pdf
	// Moreau, V., and O. Desjardins. 
	// "A second-order ghost-fluid method for the primary atomization
	//	 of liquid fuel in air-blast type injectors."
	// Proceedings of the Summer Program. Vol. 143. 2008.
	*/
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsW >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsW < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsW >= 0.0 && lsM < 0.0) {
			thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}
		else if (lsW < 0.0 && lsM >= 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		pCoefEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		u[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		u[idx(i, j)] *= pCoefEff / kDr;

		aWEff[idx(i, j)] = u[idx(i, j)];
		pCoefEffWVec[idx(i, j)] = pCoefEff;
		u[idx(i, j)] += us[idx(i, j)] - pCoefEff * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDr;
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		lsS = ls[idx(i, j - 1)];
		lsM = ls[idx(i, j)];
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsS >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsS < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsS >= 0.0 && lsM < 0.0) {
			thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}
		else if (lsS < 0.0 && lsM >= 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		pCoefEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		v[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
		v[idx(i, j)] *= pCoefEff / kDz; 

		aSEff[idx(i, j)] = v[idx(i, j)];
		pCoefEffSVec[idx(i, j)] = pCoefEff;
		v[idx(i, j)] += vs[idx(i, j)] - pCoefEff * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDz;
	}

	ApplyBC_U_2D(u);
	ApplyBC_V_2D(v);
	
	std::ofstream outF;
	std::string fname("VUpdate_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		// outF << "VARIABLES = \"R\", \"Z\", \"U\", \"V\", \"Uhat\", \"Vhat\", \"LS\", \"PS\", \"pCoefEffW\", \"Rho*AWEff\", \"Rho*Px\", \"Px+AW\",\"pCoefEffS\", \"Rho*ASEff\",\"Rho*Py\", \"Py+AS\", \"UmU\", \"VmV\", \"RealDIV\", \"KAPPA\"" << std::endl;
		outF << "VARIABLES = \"R\", \"Z\", \"U\", \"V\", \"Uhat\", \"Vhat\", \"LS\", \"PS\", \"pCoefEffW\", \"JumpW\", \"Px+AW\",\"pCoefEffS\",\"JumpS\",\"Py+AS\", \"KAPPA\"" << std::endl;
		outF.close();
	}

	if (m_iter % kNIterSkip == 0) {
		outF.open(fname.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << m_iter
			<< std::string("\", I=") << kNr + 6 << std::string(", J=") << kNz + 6
			<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
			<< std::string(", STRANDID=") << m_iter + 1
			<< std::endl;

		double rPM = 0.0, rPEHalf = 0.0, rPWHalf = 0.0;
		// for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		// for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++) {
		for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
		for (int i = 0; i < kNr + 2 * kNumBCGrid; i++) {
			rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
			rPEHalf = kBaseR + (i + 1.0 - kNumBCGrid) * kDr;
			rPWHalf = kBaseR + (i - kNumBCGrid) * kDr;
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			// << static_cast<double>((u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5) << std::string(",")
			// << static_cast<double>((v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5) << std::string(",")
			// << static_cast<double>((us[idx(i, j)] + us[idx(i + 1, j)]) * 0.5) << std::string(",")
			// << static_cast<double>((vs[idx(i, j)] + vs[idx(i, j + 1)]) * 0.5) << std::string(",")
			<< static_cast<double>(u[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)]) << std::string(",")
			<< static_cast<double>(us[idx(i, j)]) << std::string(",")
			<< static_cast<double>(vs[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefEffWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(u[idx(i, j)] - us[idx(i, j)]) << std::string(",")
			<< static_cast<double>(pCoefEffSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)] - vs[idx(i, j)]) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;
		}
		outF.close();
	}
	return 0;
}

double MACSolver2DAxisym::UpdateDt(const std::vector<double>& u, const std::vector<double>& v) {
	double uAMax = 0.0, vAMax = 0.0, kAMax = 0.0;
	double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0, Sefl = 0.0;
	double dt = std::numeric_limits<double>::max();

	// get maximum of absolute value
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		uAMax = std::max(uAMax, std::fabs((u[idx(i + 1, j)] * u[idx(i, j)]) * 0.5));
		vAMax = std::max(vAMax, std::fabs((v[idx(i, j + 1)] * v[idx(i, j)]) * 0.5));
		kAMax = std::max(kAMax, std::fabs(m_kappa[idx(i, j)]));
	}
	
	Cefl = uAMax / kDr + vAMax / kDz;
	Vefl = std::max(kMuH / kRhoH, kMuL / kRhoL) * (2.0 / ((kDr * kDr) + (kDz * kDz)));
	Gefl = std::max(std::sqrt(std::fabs(kG) / kDz), std::sqrt(std::fabs(kG) / kDr));
	Sefl = std::sqrt((kSigma * kAMax)
		/ (std::min(kRhoH, kRhoL) * std::pow(std::min(kDr, kDz), 2.0)));
	
	dt = std::min(dt,
		kCFL / (0.5 * (Cefl +
		std::sqrt(std::pow(Cefl, 2.0) +
		4.0 * std::pow(Gefl, 2.0) + 4.0 * std::pow(Sefl, 2.0)))));

	if (std::isnan(dt) || std::isinf(dt)) {
		std::cout << "dt nan/inf error : Cefl, Vefl, Gefl : " << Cefl << " " << Vefl << " " << Gefl << std::endl;
	}

	return dt;
}

int MACSolver2DAxisym::CheckCompatibility(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& uhat, const std::vector<double>& vhat) {

	std::vector<double> div(kArrSize, 0.0), hatDiv(kArrSize, 0.0);
	double divSum = 0.0, hatDivSum = 0.0;

	div = GetDivergence(u, v);
	hatDiv = GetDivergence(uhat, vhat);
	
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		divSum += div[idx(i, j)];
		hatDivSum += hatDiv[idx(i, j)];
	}

	std::cout << " Sum of divergence of intermediate velocity : " << hatDivSum << std::endl;
	std::cout << " Sum of divergence of velocity : " << divSum << std::endl;

	return 0;
}

int MACSolver2DAxisym::SetBC_U_2D(std::string BC_W, std::string BC_E,	std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kDr, kDz, kNumBCGrid);
	}

	m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2DAxisym::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kDr, kDz, kNumBCGrid);
	}

	m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2DAxisym::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kDr, kDz, kNumBCGrid);
	}

	m_BC->SetBC_P_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2DAxisym::ApplyBC_U_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_U_2D(arr);
	return 0;
}

int MACSolver2DAxisym::ApplyBC_V_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_V_2D(arr);
	return 0;
}

int MACSolver2DAxisym::ApplyBC_P_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_P_2D(arr);
	return 0;
}

void MACSolver2DAxisym::SetBCConstantUW(double BC_ConstantW) {
	return m_BC->SetBCConstantUW(BC_ConstantW);
}

void MACSolver2DAxisym::SetBCConstantUE(double BC_ConstantE) {
	return m_BC->SetBCConstantUE(BC_ConstantE);
}

void MACSolver2DAxisym::SetBCConstantUS(double BC_ConstantS) {
	return m_BC->SetBCConstantUS(BC_ConstantS);
}

void MACSolver2DAxisym::SetBCConstantUN(double BC_ConstantN) {
	return m_BC->SetBCConstantUN(BC_ConstantN);
}

void MACSolver2DAxisym::SetBCConstantVW(double BC_ConstantW) {
	return m_BC->SetBCConstantVW(BC_ConstantW);
}

void MACSolver2DAxisym::SetBCConstantVE(double BC_ConstantE) {
	return m_BC->SetBCConstantVE(BC_ConstantE);
}

void MACSolver2DAxisym::SetBCConstantVS(double BC_ConstantS) {
	return m_BC->SetBCConstantVS(BC_ConstantS);
}

void MACSolver2DAxisym::SetBCConstantVN(double BC_ConstantN) {
	return m_BC->SetBCConstantVN(BC_ConstantN);
}

void MACSolver2DAxisym::SetBCConstantPW(double BC_ConstantW) {
	return m_BC->SetBCConstantPW(BC_ConstantW);
}

void MACSolver2DAxisym::SetBCConstantPE(double BC_ConstantE) {
	return m_BC->SetBCConstantPE(BC_ConstantE);
}

void MACSolver2DAxisym::SetBCConstantPS(double BC_ConstantS) {
	return m_BC->SetBCConstantPS(BC_ConstantS);
}

void MACSolver2DAxisym::SetBCConstantPN(double BC_ConstantN) {
	return m_BC->SetBCConstantPN(BC_ConstantN);
}

void MACSolver2DAxisym::UpdateAmbientPressure(double ambientPressure) {
	// cause I use pseudo-pressure (m_dt * pressure)
	return m_BC->SetAmbientPressure(m_dt * ambientPressure);
}

// https://en.wikibooks.org/wiki/Algorithm_Implementation/Linear_Algebra/Tridiagonal_matrix_algorithm
int MACSolver2DAxisym::SolveTDMA(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n) {
	/*
	// n is the number of unknowns

	|b0 c0 0 ||x0| |d0|
	|a1 b1 c1||x1|=|d1|
	|0  a2 b2||x2| |d2|

	1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

	x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

	2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
	from 1st it.: -| a1x0 + a1g0x1        = a1r0
	-----------------------------
	(b1 - a1g0)x1 + c1x2 = d1 - a1r0

	x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

	3rd iteration:      | a2x1 + b2x2   = d2
	from 2st it. : -| a2x1 + a2g1x2 = a2r2
	-----------------------
	(b2 - a2g1)x2 = d2 - a2r2
	x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
	Finally we have a triangular matrix:
	|1  g0 0 ||x0| |r0|
	|0  1  g1||x1|=|r1|
	|0  0  1 ||x2| |r2|

	Condition: ||bi|| > ||ai|| + ||ci||

	in this version the c matrix reused instead of g
	and             the d matrix reused instead of r and x matrices to report results
	Written by Keivan Moradi, 2014
	*/
	n--; // since we start from x0 (not x1)
	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++) {
		c[i] /= b[i] - a[i] * c[i - 1];
		d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
	}

	d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c[n - 1]);

	for (int i = n; i-- > 0;) {
		d[i] -= c[i] * d[i + 1];
	}

	return 0;
}

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
inline int MACSolver2DAxisym::sign(const double& val) {
	return (val > 0) - (val < 0);
}

inline int MACSolver2DAxisym::idx(int i, int j) {
	return (j + (kNz + 2 * kNumBCGrid) * (i));
}

int MACSolver2DAxisym::SetPLTType(PLTTYPE type) {
	m_PLTType = type;

	return 0;
}

int MACSolver2DAxisym::OutRes(const int iter, const double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ps, const std::vector<double>& ls) {
	if (m_PLTType == PLTTYPE::ASCII || m_PLTType == PLTTYPE::BOTH) {
		std::ofstream outF;
		std::string fname_vel(fname_vel_base + "_ASCII.plt");
		std::string fname_div(fname_div_base + "_ASCII.plt");
		if (m_iter == 0) {
			outF.open(fname_vel.c_str(), std::ios::out);

			outF << "TITLE = VEL" << std::endl;
			outF << "VARIABLES = \"R\", \"Z\", \"U\", \"V\", \"LS\", \"PS\" " << std::endl;
			outF.close();

			outF.open(fname_div.c_str(), std::ios::out);
			outF << "TITLE = DIV" << std::endl;
			outF << "VARIABLES = \"R\", \"Z\", \"LS\", \"DIV\", \"PS\" " << std::endl;
			outF.close();
		}

		std::vector<double>
			resU(kArrSize, 0.0),
			resV(kArrSize, 0.0);

		std::vector<double> resDiv = GetDivergence(u, v);
		
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			resU[idx(i, j)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;

		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			resV[idx(i, j)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
		
		outF.open(fname_vel.c_str(), std::ios::app);
		
		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNr << std::string(", J=") << kNz
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;

		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
				<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
				<< static_cast<double>(resU[idx(i, j)]) << std::string(",")
				<< static_cast<double>(resV[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j)]) << std::endl;

		outF.close();
		
		outF.open(fname_div.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNr << std::string(", J=") << kNz
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;
		
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
				<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
				<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
				<< static_cast<double>(resDiv[idx(i, j)]) << std::string(",")
				// << static_cast<double>(m_kappa[idx(i, j)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j)]) << std::endl;

		outF.close();
	}

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH) {
		INTEGER4 whichFile = 0, stat = 0;
		std::string fname_vel(fname_vel_base + "_BINARY.szplt");
		std::string fname_div(fname_div_base + "_BINARY.szplt");

		std::vector<double>
			resX(kNr * kNz, 0.0), resY(kNr * kNz, 0.0),
			resU(kNr * kNz, 0.0), resV(kNr * kNz, 0.0),
			resDivCompact(kNr * kNz, 0.0), resLS(kNr * kNz, 0.0),
			resPs(kNr * kNz, 0.0);
		std::vector<double> resDivFull = GetDivergence(u, v);
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)	 {
			resX[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
			resY[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = kBaseZ + (j + 0.5 - kNumBCGrid) * kDz;
			resU[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
			resV[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = ls[idx(i, j)];
			resPs[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = ps[idx(i, j)];
			resDivCompact[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = resDivFull[idx(i, j)];
		}

		INTEGER4 Debug = 0;
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
		INTEGER4 IMax = kNr, JMax = kNz, KMax = 1;

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
		stat = TECDAT142(&ARRSIZEVAL, resDivCompact.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);
	}

	return 0;
}

int MACSolver2DAxisym::OutResClose() {
	INTEGER4 stat, whichFile;
	whichFile = 1;
	stat = TECFIL142(&whichFile);

	// close first file (velocity)
	stat = TECEND142();
	stat = TECEND142();

	return 0;
}
