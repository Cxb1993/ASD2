#include "mac2dAxisym.h"

MACSolver2DAxisym::MACSolver2DAxisym(double Re, double We, double Fr, GAXISENUM2DAXISYM GAxis,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoH, double muH,
	int nr, int nz, double baseR, double baseZ, double lenR, double lenZ,
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kGAxis(GAxis), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoH(rhoH), kRhoL(rhoH * densityRatio), kRhoRatio(densityRatio),
	kMuH(muH), kMuL(muH * viscosityRatio), kMuRatio(viscosityRatio),
	kNr(nr), kNz(nz), kBaseR(baseR), kBaseZ(baseZ), kLenR(lenR), kLenZ(lenZ),
	kDr(lenR / static_cast<double>(nr)), kDz(lenZ / static_cast<double>(nz)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
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
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoH), kMuScale(muH), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoH(rhoH), kRhoL(rhoL), kRhoRatio(rhoL / rhoH), 
	kMuH(muH), kMuL(muL), kMuRatio(muL / muH),
	kRe(rhoH * L * U / muH), kWe(rhoH * L * U * U / sigma), kFr(U * U / (gConstant * L)), kGAxis(GAxis),
	kNr(nr), kNz(nz), kBaseR(baseR), kBaseZ(baseZ), kLenR(lenR), kLenZ(lenZ),
	kDr(lenR / static_cast<double>(nr)), kDz(lenZ / static_cast<double>(nz)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
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
		// inside
		if (ls[idx(i, j)] > eps)
			Heaviside[idx(i, j)] = 1.0;
		// outside
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
		
		LSSize[idx(i, j)] = std::sqrt(dLSdR[idx(i, j)] * dLSdR[idx(i, j)] + dLSdZ[idx(i, j)] * dLSdZ[idx(i, j)] + eps);

		if (LSSize[idx(i, j)] < eps)
			perror("Div/0 Err in computing kappa");
	}
	
	ApplyBC_P_2D(dLSdR);
	ApplyBC_P_2D(dLSdZ);
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rPW = kBaseR + (i - 0.5 - kNumBCGrid) * kDr;
		rPE = kBaseR + (i + 1.5 - kNumBCGrid) * kDr;
		m_kappa[idx(i, j)]
			= (dLSdR[idx(i, j)] * dLSdR[idx(i, j)] * (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDz * kDz) // phi^2_x \phi_yy
			- 2.0 * dLSdR[idx(i, j)] * dLSdZ[idx(i, j)] * (dLSdR[idx(i, j + 1)] - dLSdR[idx(i, j - 1)]) / (2.0 * kDz) //2 \phi_x \phi_y \phi_xy
			+ dLSdZ[idx(i, j)] * dLSdZ[idx(i, j)] * (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / (kDr * kDr)) // phi^2_y \phi_xx
			/ std::pow(LSSize[idx(i, j)], 3.0)
			+ dLSdR[idx(i, j)] / (LSSize[idx(i, j)] * rPM);

		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j)] = sign(m_kappa[idx(i, j)]) * std::fabs(std::min(std::fabs(m_kappa[idx(i, j)]), 1.0 / std::min(kDr, kDz)));
		
		assert(m_kappa[idx(i, j)] == m_kappa[idx(i, j)]);
	}

	return 0;
}

int MACSolver2DAxisym::UpdateJumpCond(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls) {
	const double eps = 1.0e-100;
	std::vector<double> U_PGrid(kArrSize, 0.0),
		V_PGrid(kArrSize, 0.0);
	// derivatives
	std::vector<double> dudR(kArrSize, 0.0);
	std::vector<double> dudZ(kArrSize, 0.0);
	std::vector<double> dvdR(kArrSize, 0.0);
	std::vector<double> dvdZ(kArrSize, 0.0);
	std::vector<double> dldR(kArrSize, 0.0);
	std::vector<double> dldZ(kArrSize, 0.0);
	// normal and tangent vector variable (n, t1, t2)
	// normal vector has X, Y component
	// first tangent vector has only Z component
	// double nX = 0.0, nY = 0.0, t1X = 0.0, t1Y = 0.0, t2X = 0.0, t2Y = 0.0, t2Z = 0.0;
	std::vector<double> nR(kArrSize, 0.0);
	std::vector<double> nZ(kArrSize, 0.0);
	std::vector<double> t1R(kArrSize, 0.0);
	std::vector<double> t1Z(kArrSize, 0.0);
	std::vector<double> t2t(kArrSize, 0.0);
	
	// jump condition originally defined as J11 = [\mu u_x], J12 = [\mu u_y], J21 = [\mu v_x], J22 = [\mu v_y], 
	// but skip [\mu], from Kang, Fedkiw, and Liu's work eq. (30).
	// [\mu] changes with level set, which means I can determine [\mu] later
	
	double nrnr = 0.0, nrny = 0.0, nyny = 0.0;
	double t1rt1r = 0.0, t1rt1z = 0.0, t1zt1z = 0.0;
	double theta = 0.0;
	double lsW = 0.0, lsE = 0.0, lsM = 0.0, lsS = 0.0, lsN = 0.0;
	double lsUW = 0.0, lsUM = 0.0, lsUE = 0.0, lsUSHalf = 0.0, lsUNHalf = 0.0;
	double lsVS = 0.0, lsVM = 0.0, lsVN = 0.0, lsVWHalf = 0.0, lsVEHalf = 0.0;

	double uEff = 0.0, vEff = 0.0;
	double dudRW = 0.0, dudRE = 0.0, dudZS = 0.0, dudZN = 0.0;
	double dvdRW = 0.0, dvdRE = 0.0, dvdZS = 0.0, dvdZN = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		// defined at P grid
		dudR[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDr;
		dudZ[idx(i, j)] = (0.25 * (u[idx(i, j)] + u[idx(i + 1, j)] + u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
		- 0.25 * (u[idx(i, j)] + u[idx(i + 1, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / kDz;
		dvdR[idx(i, j)] = (0.25 * (v[idx(i, j)] + v[idx(i, j + 1)] + v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.25 * (v[idx(i, j)] + v[idx(i, j + 1)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / kDr;
		dvdZ[idx(i, j)] = (v[idx(i, j + 1)] - v[idx(i, j)]) / kDz;
		
		dudR[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dudZ[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDz);
		dvdR[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dvdZ[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDz);

		// velocity derivative could have a jump. compensate its derivative
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		lsE = ls[idx(i + 1, j)];
		lsS = ls[idx(i, j - 1)];
		lsN = ls[idx(i, j + 1)];
		lsUW = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsUM = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsUE = 0.5 * (ls[idx(i + 1, j)] + ls[idx(i + 2, j)]);
		
		lsVS = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsVM = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		lsVN = 0.5 * (ls[idx(i, j + 1)] + ls[idx(i, j + 2)]);

		lsUSHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i, j - 1)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)]);
		lsUNHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i, j + 1)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)]);
		lsVWHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i, j - 1)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)]);
		lsVEHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i, j - 1)] + ls[idx(i + 1, j)] + ls[idx(i + 1, j - 1)]);
		
		// Jump occurs when computing dudR
		if (lsW * lsE <= 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudRW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dudRE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDr);
			dudR[idx(i, j)] = dudRE * theta + dudRW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdRW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dvdRE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDr);
			dvdR[idx(i, j)] = dvdRE * theta + dvdRW * (1.0 - theta);
		}

		if (lsS * lsN <= 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudZS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dudZN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDz);
			dudZ[idx(i, j)] = dudZN * theta + dudZS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdZS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dvdZN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDz);
			dvdZ[idx(i, j)] = dvdZN * theta + dvdZS * (1.0 - theta);
		}
		
		dldR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDr);
		dldZ[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDz);
		if (std::fabs(dldR[idx(i, j)]) < eps) {
			dldR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDr;
			if (std::fabs(dldR[idx(i, j)]) < eps)
				dldR[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDr;
		}
		if (std::fabs(dldR[idx(i, j)]) < eps) {
			dldZ[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDz;
			if (std::fabs(dldZ[idx(i, j)]) < eps)
				dldZ[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDz;
		}
		
		// normal vector = (\nabla \phi) / |\nabla \phi|
		nR[idx(i, j)] = dldR[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps);
		nZ[idx(i, j)] = dldZ[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps);
		// get first tangent vector
		// smallest magnitude determine what to be performed cross product
		// for 2D, smallest magnitude must be z axis (because it is 2D!)
		t1R[idx(i, j)] = nZ[idx(i, j)] / (std::sqrt(std::pow(nR[idx(i, j)], 2.0) + std::pow(nZ[idx(i, j)], 2.0)) + eps);
		t1Z[idx(i, j)] = -nR[idx(i, j)] / (std::sqrt(std::pow(nR[idx(i, j)], 2.0) + std::pow(nZ[idx(i, j)], 2.0)) + eps);
		t2t[idx(i, j)] = (nR[idx(i, j)] * t1Z[idx(i, j)] - nZ[idx(i, j)] * t1R[idx(i, j)]);

		nrnr = nR[idx(i, j)] * nR[idx(i, j)];
		nrny = nR[idx(i, j)] * nZ[idx(i, j)];
		nyny = nZ[idx(i, j)] * nZ[idx(i, j)];
		t1rt1r = t1R[idx(i, j)] * t1R[idx(i, j)];
		t1rt1z = t1R[idx(i, j)] * t1Z[idx(i, j)];
		t1zt1z = t1Z[idx(i, j)] * t1Z[idx(i, j)];

		// eq. (30) from Kang, Fedkiw, and Liu
		// u_x
		m_J11[idx(i, j)] = dudR[idx(i, j)] * t1rt1r + dudZ[idx(i, j)] * t1rt1z
			+ (dudR[idx(i, j)] * nrnr + dvdR[idx(i, j)] * nrny) * nrnr + (dudZ[idx(i, j)] * nrnr + dvdZ[idx(i, j)] * nrny) * nrny
			- ((dudR[idx(i, j)] * t1rt1r + dudZ[idx(i, j)] * t1rt1z) * nrnr + (dvdR[idx(i, j)] * t1rt1r + dvdZ[idx(i, j)] * t1rt1z) * nrny);
		m_J11[idx(i, j)] *= (kMuH - kMuL);

		// u_y
		m_J12[idx(i, j)] = dudR[idx(i, j)] * t1rt1z + dudZ[idx(i, j)] * t1zt1z
			+ (dudR[idx(i, j)] * nrnr + dvdR[idx(i, j)] * nrny) * nrny + (dudZ[idx(i, j)] * nrnr + dvdZ[idx(i, j)] * nrny) * nyny
			- ((dudR[idx(i, j)] * t1rt1r + dudZ[idx(i, j)] * t1rt1z) * nrny + (dvdR[idx(i, j)] * t1rt1r + dvdZ[idx(i, j)] * t1rt1z) * nyny);
		m_J12[idx(i, j)] *= (kMuH - kMuL);

		// v_x
		m_J21[idx(i, j)] = dvdR[idx(i, j)] * t1rt1r + dvdZ[idx(i, j)] * t1rt1z
			+ (dudR[idx(i, j)] * nrny + dvdR[idx(i, j)] * nyny) * nrnr + (dudZ[idx(i, j)] * nrny + dvdZ[idx(i, j)] * nyny) * nrny
			- ((dudR[idx(i, j)] * t1rt1z + dudZ[idx(i, j)] * t1zt1z) * nrnr + (dvdR[idx(i, j)] * t1rt1z + dvdZ[idx(i, j)] * t1zt1z) * nrny);
		m_J21[idx(i, j)] *= (kMuH - kMuL);

		// v_y
		m_J22[idx(i, j)] = dvdR[idx(i, j)] * t1rt1z + dvdZ[idx(i, j)] * t1zt1z
			+ (dudR[idx(i, j)] * nrny + dvdR[idx(i, j)] * nyny) * nrny + (dudZ[idx(i, j)] * nrny + dvdZ[idx(i, j)] * nyny) * nyny
			- ((dudR[idx(i, j)] * t1rt1z + dudZ[idx(i, j)] * t1zt1z) * nrny + (dvdR[idx(i, j)] * t1rt1z + dvdZ[idx(i, j)] * t1zt1z) * nyny);
		m_J22[idx(i, j)] *= (kMuH - kMuL);

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

std::vector<double> MACSolver2DAxisym::UpdateFU(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& H) {
	
	std::vector<double> cU(kArrSize, 0.0);
	std::vector<double> vU(kArrSize, 0.0);
	std::vector<double> gU(kArrSize, 0.0);
	std::vector<double> rhsU(kArrSize, 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v);

	// Viscous term
	vU = this->AddViscosityFU(u, v, ls, H);

	gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rhsU[idx(i, j)] = -cU[idx(i, j)] + vU[idx(i, j)] + gU[idx(i, j)];
	}
	
	std::ofstream outF;
	std::string fname("FU_ASCII.plt");
	if (m_iter == 1) {
	outF.open(fname.c_str(), std::ios::out);

	outF << "TITLE = VEL" << std::endl;
	outF << "VARIABLES = \"R\", \"Z\", \"-cU\", \"vV\", \"gU\", \"rhsU\"" << std::endl;
	outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"") << m_iter
	<< std::string("\", I=") << kNr << std::string(", J=") << kNz
	<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
	<< std::string(", STRANDID=") << m_iter + 1
	<< std::endl;

	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
	<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
	<< static_cast<double>(-cU[idx(i, j)]) << std::string(",")
	<< static_cast<double>(vU[idx(i, j)]) << std::string(",")
	<< static_cast<double>(gU[idx(i, j)]) << std::string(",")
	<< static_cast<double>(rhsU[idx(i, j)])
	<< std::endl;

	outF.close();
	
	return rhsU;
}

std::vector<double> MACSolver2DAxisym::UpdateFV(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& H) {

	std::vector<double> cV(kArrSize, 0.0);
	std::vector<double> vV(kArrSize, 0.0);
	std::vector<double> gV(kArrSize, 0.0);
	std::vector<double> rhsV(kArrSize, 0.0);

	// Convection term
	cV = this->AddConvectionFV(u, v);

	// Viscous term
	vV = this->AddViscosityFV(u, v, ls, H);

	gV = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		rhsV[idx(i, j)] = -cV[idx(i, j)] + vV[idx(i, j)] + gV[idx(i, j)];
	}
	
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
		for (int j = 0; j < kNz + 2 *kNumBCGrid; j++) {
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

std::vector<double> MACSolver2DAxisym::AddViscosityFU(const std::vector<double>& u, const std::vector<double>& v, 
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dU(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);

	std::vector<double> visRVec(kArrSize, 0.0),
		visZVec(kArrSize, 0.0),
		UxxVec(kArrSize, 0.0),
		UyyVec(kArrSize, 0.0),
		VxyVec(kArrSize, 0.0),
		VxyNVec(kArrSize, 0.0),
		VxySVec(kArrSize, 0.0);

	std::vector<double> muU_W(kArrSize, 0.0),
		muU_E(kArrSize, 0.0),
		muU_S(kArrSize, 0.0),
		muU_N(kArrSize, 0.0),
		muV_W(kArrSize, 0.0),
		muV_E(kArrSize, 0.0),
		resJ11W(kArrSize, 0.0), 
		resJ11E(kArrSize, 0.0),
		resJ12S(kArrSize, 0.0), 
		resJ12N(kArrSize, 0.0),
		resJ21W(kArrSize, 0.0),
		resJ21E(kArrSize, 0.0),
		checkDiv(kArrSize, 0.0);
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			dU[idx(i, j)] = 0.0;
		
		return dU;
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

	// subcell
	double theta = 0.0;

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
	double lsW = 0.0, lsM = 0.0, lsUNHalf = 0.0, lsUSHalf = 0.0;
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
	// u_xx, u_yy
	double lsUW = 0.0, lsUE = 0.0, lsUS = 0.0, lsUN = 0.0, lsUM = 0.0;
	// v_xy
	double lsVW_N = 0.0, lsVW = 0.0, lsVM = 0.0, lsVN = 0.0;
	double uW = 0.0, uE = 0.0, uS = 0.0, uN = 0.0, uM = 0.0;
	double vW = 0.0, vW_N = 0.0, vM = 0.0, vN = 0.0;
	double muU_R_W = 0.0, muU_R_E = 0.0, muU_Z_S = 0.0, muU_Z_N = 0.0;
	double muV_R_S = 0.0, muV_R_N = 0.0;
	double rhoU_R_W = 0.0, rhoU_R_E = 0.0, rhoU_Z_S = 0.0, rhoU_Z_N = 0.0;
	double rhoV_R_S = 0.0, rhoV_R_N = 0.0;
	double rU_R_W = 0.0, rU_R_E = 0.0, rU = 0.0;
	double visR = 0.0, visZ = 0.0;

	// jump condition
	double JUW = 0.0, JUE = 0.0, JUS = 0.0, JUN = 0.0, JUM = 0.0;
	double JVW = 0.0, JVW_N = 0.0, JVM = 0.0, JVN = 0.0;
	// effective Jump condition, effective u(uEff), effective mu (muEff), and effective rho (rhoEff)
	// J is a jump condition and defined at P grid
	double JEff = 0.0, JO = 0.0, uEff = 0.0, vEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
	const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		visR = 0.0, visZ = 0.0;
		muU_R_W = 0.0; muU_R_E = 0.0; muU_Z_S = 0.0; muU_Z_N = 0.0;
		muV_R_S = 0.0; muV_R_N = 0.0;
		
		lsUW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
		lsUE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
		lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsUS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsUN = 0.5 * (ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsVW_N = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)]);
		lsVM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		uM = u[idx(i, j)];
		uW = u[idx(i - 1, j)];
		uE = u[idx(i + 1, j)];
		uS = u[idx(i, j - 1)];
		uN = u[idx(i, j + 1)];
		vW = v[idx(i - 1, j)];
		vW_N = v[idx(i - 1, j + 1)];
		vM = v[idx(i, j)];
		vN = v[idx(i, j + 1)];
		
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		lsUNHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
		lsUSHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		/*
		// U part 
		if (lsUW >= 0 && lsUM >= 0) {
			muU_R_W = kMuH * (uM - uW) / kDr;
			muU_W[idx(i, j)] = muU_R_W;
		}
		else if (lsUW <= 0 && lsUM <= 0) {
			muU_R_W = kMuL * (uM - uW) / kDr;
			muU_W[idx(i, j)] = muU_R_W;
		}
		else if (lsUW > 0 && lsUM <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsUW) / (std::fabs(lsUW) + std::fabs(lsUM));
			// |(lsUW)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsUM)|
			// |(lsUW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUW = (m_J11[idx(i - 2, j)] + m_J11[idx(i - 1, j)]) * 0.5;
			JUM = (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUW;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_R_W = muEff * (uM - uW) / kDr
				- muEff * JEff * theta / kMuH;
			muU_W[idx(i, j)] = muEff * (uM - uW) / kDr;
			resJ11W[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsUW <= 0 && lsUM > 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsUW) / (std::fabs(lsUW) + std::fabs(lsUM));
			// |(lsUW)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsUM)|
			// |(lsUW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUW = (m_J11[idx(i - 2, j)] + m_J11[idx(i - 1, j)]) * 0.5;
			JUM = (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUW;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_R_W = muEff * (uM - uW) / kDr
				+ muEff * JEff * theta / kMuL;
			muU_W[idx(i, j)] = muEff * (uM - uW) / kDr;
			resJ11W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsUM >= 0 && lsUE >= 0) {
			muU_R_E = kMuH * (uE - uM) / kDr;
			muU_E[idx(i, j)] = muU_R_E;
		}
		else if (lsUM <= 0 && lsUE <= 0) {
			muU_R_E = kMuL * (uE - uM) / kDr;
			muU_E[idx(i, j)] = muU_R_E;
		}
		else if (lsUM > 0 && lsUE <= 0) {
			// interface lies between u[i, j] and u[i + 1, j]
			theta = std::fabs(lsUE) / (std::fabs(lsUE) + std::fabs(lsUM));
			// |(lsUM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsUE)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsUE)|
			JUM = (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JUE = (m_J11[idx(i, j)] + m_J11[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUE;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_R_E = muEff * (uE - uM) / kDr
				+ muEff * JEff * theta / kMuL;
			muU_E[idx(i, j)] = muEff * (uE - uM) / kDr;
			resJ11E[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		else if (lsUM <= 0 && lsUE > 0) {
			// interface lies between u[i, j] and u[i + 1, j]
			theta = std::fabs(lsUE) / (std::fabs(lsUE) + std::fabs(lsUM));
			// |(lsUM)| ===   outside(-)    === |(interface)| === inside(+) === |(lsUE)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsUE)|
			JUM = (m_J11[idx(i - 1, j)] + m_J11[idx(i, j)]) * 0.5;
			JUE = (m_J11[idx(i, j)] + m_J11[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUE;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_R_E = muEff * (uE - uM) / kDr
				- muEff * JEff * theta / kMuH;
			muU_E[idx(i, j)] = muEff * (uE - uM) / kDr;
			resJ11E[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		
		if (lsUS >= 0 && lsUM >= 0) {
			muU_Z_S = kMuH * (uM - uS) / kDz;
			muU_S[idx(i, j)] = muU_Z_S;
		}
		else if (lsUS <= 0 && lsUM <= 0) {
			muU_Z_S = kMuL * (uM - uS) / kDz;
			muU_S[idx(i, j)] = muU_Z_S;
		}
		else if (lsUS > 0 && lsUM <= 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = std::fabs(lsUS) / (std::fabs(lsUS) + std::fabs(lsUM));
			// |(lsUS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsUM)|
			// |(lsUS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_Z_S = muEff * (uM - uS) / kDz
				- muEff * JEff * theta / kMuH;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDz;
			resJ12S[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsUS <= 0 && lsUM > 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = std::fabs(lsUS) / (std::fabs(lsUS) + std::fabs(lsUM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			JUS = (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_Z_S = muEff * (uM - uS) / kDz
				+ muEff * JEff * theta / kMuL;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDz;
			resJ12S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsUM >= 0 && lsUN >= 0) {
			muU_Z_N = kMuH * (uN - uM) / kDz;
			muU_N[idx(i, j)] = muU_Z_N;
		}
		else if (lsUM <= 0 && lsUN <= 0) {
			muU_Z_N = kMuL * (uN - uM) / kDz;
			muU_N[idx(i, j)] = muU_Z_N;
		}
		else if (lsUM > 0 && lsUN <= 0) {
			// interface lies between u[i, j] and u[i, j + 1]
			theta = std::fabs(lsUN) / (std::fabs(lsUN) + std::fabs(lsUM));
			// |(lsUM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsUN)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsUN)|
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JUN = (m_J12[idx(i - 1, j + 1)] + m_J12[idx(i, j + 1)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUN;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_Z_N = muEff * (uN - uM) / kDz 
				+ muEff * JEff * theta / kMuL;
			muU_N[idx(i, j)] = muEff * (uN - uM) / kDz;
			resJ12N[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		else if (lsUM <= 0 && lsUN > 0) {
			// interface lies between u[i, j] and u[i, j + 1]
			theta = std::fabs(lsUN) / (std::fabs(lsUN) + std::fabs(lsUM));
			// |(lsUM)| ===    outside(-)   === |(interface)| ===   inside(+) === |(lsUN)|
			// |(lsUM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsUN)|
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JUN = (m_J12[idx(i - 1, j + 1)] + m_J12[idx(i, j + 1)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUN;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_Z_N = muEff * (uN - uM) / kDz 
				- muEff * JEff * theta / kMuH;
			muU_N[idx(i, j)] = muEff * (uN - uM) / kDz;
			resJ12N[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		*/
		// V part
		/*
		-------------------------------------
		|			|			|			|
		|		  lsUN			|			|
		|			|			|			|
		----lsVW_N-------lsVN----------------
		|			|			|			|
		lsUW	  lsUM  (i,j)  lsUE			|
		|			|			|			|
		----lsVW---------lsVM----------------
		|			|			|			|
		|		  lsUS		  lsUE_S		|
		|			|			|			|
		-------------------------------------
		*/
		/*
		if (lsVW >= 0 && lsVM >= 0) {
			muV_R_S = kMuH * (vM - vW) / kDr;
			muV_W[idx(i, j)] = muV_R_S;
		}
		else if (lsVW <= 0 && lsVM <= 0) {
			muV_R_S = kMuL * (vM - vW) / kDr;
			muV_W[idx(i, j)] = muV_R_S;
		}
		else if (lsVW > 0 && lsVM <= 0) {
			// interface lies between v[i - 1, j] and v[i, j]
			theta = std::fabs(lsVW) / (std::fabs(lsVW) + std::fabs(lsVM));
			// |(lsVW)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (m_J21[idx(i - 1, j - 1)] + m_J21[idx(i - 1, j)]) * 0.5;
			JVM = (m_J21[idx(i, j - 1)] + m_J21[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_R_S = muEff * (vM - vW) / kDr
				- muEff * JEff * theta / kMuH;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDr;
			resJ21W[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsVW <= 0 && lsVM > 0) {
			// interface lies between v[i - 1, j] and v[i, j]
			theta = std::fabs(lsVW) / (std::fabs(lsVW) + std::fabs(lsVM));
			// |(lsVW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (m_J21[idx(i - 1, j - 1)] + m_J21[idx(i - 1, j)]) * 0.5;
			JVM = (m_J21[idx(i, j - 1)] + m_J21[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_R_S = muEff * (vM - vW) / kDr
				+ muEff * JEff * theta / kMuL;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDr;
			resJ21W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}

		if (lsVW_N >= 0 && lsVN >= 0) {
			muV_R_N = kMuH * (vN - vW_N) / kDr;
			muV_E[idx(i, j)] = muV_R_N;
		}
		else if (lsVW_N <= 0 && lsVN <= 0) {
			muV_R_N = kMuL * (vN - vW_N) / kDr;
			muV_E[idx(i, j)] = muV_R_N;
		}
		else if (lsVW_N > 0 && lsVN <= 0) {
			// interface lies between v[i - 1, j + 1] and v[i, j + 1]
			theta = std::fabs(lsVW_N) / (std::fabs(lsVW_N) + std::fabs(lsVN));
			// |(lsVW_N)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsVN)|
			// |(lsVW_N)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVN)|
			JVW_N = (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j + 1)]) * 0.5;
			JVN = (m_J21[idx(i, j)] + m_J21[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVN + (1.0 - theta) * JVW_N;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_R_N = muEff * (vN - vW_N) / kDr
				- muEff * JEff * theta / kMuH;
			muV_E[idx(i, j)] = muEff * (vN - vW_N) / kDr;
			resJ21E[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsVW_N <= 0 && lsVN > 0) {
			// interface lies between v[i - 1, j + 1] and v[i, j + 1]
			theta = std::fabs(lsVW_N) / (std::fabs(lsVW_N) + std::fabs(lsVN));
			// |(lsVW_N)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVN)|
			// |(lsVW_N)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVN)|
			JVW_N = (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j + 1)]) * 0.5;
			JVN = (m_J21[idx(i, j)] + m_J21[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVN + (1.0 - theta) * JVW_N;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_R_N = muEff  * (vN - vW_N) / kDr
				+ muEff * JEff * theta / kMuL;
			muV_E[idx(i, j)] = muEff * (vN - vW_N) / kDr;
			resJ21E[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		*/
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

		double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;
		if (lsW >= 0 && lsM >= 0) {
			rhoEffWE = kRhoH;
			iRhoEffWE = 1.0 / kRhoH;
		}
		else if (lsW <= 0 && lsM <= 0) {
			rhoEffWE = kRhoL;
			iRhoEffWE = 1.0 / kRhoL;
		}
		else if (lsW > 0 && lsM <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			// rhoEffWE = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
			rhoEffWE = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			// rhoEffWE = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
			rhoEffWE = kRhoL * theta + kRhoH * (1.0 - theta);
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
		else if (lsUSHalf > 0 && lsUNHalf <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			// |(lsUSHalf)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsUNHalf)|
			// |(lsUSHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
			// rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
			rhoEffSN = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsUSHalf <= 0 && lsUNHalf > 0) {
			// interface lies between lsUSHalf and lsUNHalf
			theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
			// |(lsUSHalf)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsUNHalf)|
			// |(lsUSHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
			// rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
			rhoEffSN = kRhoL * theta + kRhoH * (1.0 - theta);
			iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}

		muU_R_W = mu[idx(i - 1, j)] * (uM - uW) / kDr;
		muU_R_E = mu[idx(i, j)] * (uE - uM) / kDr;
		muU_Z_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (uM - uS) / kDz;
		muU_Z_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
			* (uN - uM) / kDz;
		muV_R_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (vM - vW) / kDr;
		muV_R_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
			* (vN - vW_N) / kDr;
		
		rU_R_W = kBaseR + (i - 0.5 - kNumBCGrid) * kDr;
		rU_R_E = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rU = kBaseR + (i - kNumBCGrid) * kDr;

		visR = 1.0 / rU * (2.0 * rU_R_E * muU_R_E - 2.0 * rU_R_W * muU_R_W) / kDr;
		visZ = (muU_Z_N - muU_Z_S) / kDz + (muV_R_N - muV_R_S) / kDz;
		visR -= 2.0 * 0.5 * (mu[idx(i - 1, j)] + mu[idx(i, j)]) * u[idx(i, j)] / (rU * rU);
		// checkDiv[idx(i, j)] = (muU_R_E - muU_R_W) / kDr + (muV_R_N - muV_R_S) / kDz;
		// dU[idx(i, j)] = visR / rhoEffWE + visZ / rhoEffSN;
		dU[idx(i, j)] = visR * iRhoEffWE + visZ * iRhoEffSN;
		// dU[idx(i, j)] = visR + visZ;

		if (std::isnan(dU[idx(i, j)]) || std::isinf(dU[idx(i, j)])) {
			std::cout << "U-viscosity term nan/inf error : " << i << " " << j << " " << dU[idx(i, j)] << std::endl;
			exit(1);
		}
		assert(dU[idx(i, j)] == dU[idx(i, j)]);
	}
	/*
	std::ofstream outF;
	std::string fname("VisU_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"dU\", \"UW\", \"UE\", \"UE-UW\",\"J11\",\"US\", \"UN\", \"UN-US\", \"J12\",\"VW\", \"VE\", \"VE-VW\", \"J21\", \"CHKDIV\" " << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNr << std::string(", J=") << kNz
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			<< static_cast<double>(u[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(dU[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_E[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muU_E[idx(i, j)] - muU_W[idx(i, j)]) / kDr) << std::string(",")
			<< static_cast<double>((resJ11E[idx(i, j)] - resJ11W[idx(i, j)]) / (kDr)) << std::string(",")
			<< static_cast<double>(muU_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_N[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muU_N[idx(i, j)] - muU_S[idx(i, j)]) / kDr) << std::string(",")
			<< static_cast<double>((resJ12N[idx(i, j)] - resJ12S[idx(i, j)]) / (kDz)) << std::string(",")
			<< static_cast<double>(muV_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_E[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muV_E[idx(i, j)] - muV_W[idx(i, j)]) / kDz) << std::string(",")
			<< static_cast<double>((resJ21E[idx(i, j)] - resJ21W[idx(i, j)]) / (kDr)) << std::string(",")
			<< static_cast<double>(checkDiv[idx(i, j)]) 
			<< std::endl;

	outF.close();
	*/
	return dU;
}

std::vector<double> MACSolver2DAxisym::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV(kArrSize, 0.0);
	std::vector<double> mu(kArrSize, 0.0);

	std::vector<double> visRVec(kArrSize, 0.0),
		visZVec(kArrSize, 0.0),
		visR2Vec(kArrSize, 0.0),
		visZ2Vec(kArrSize, 0.0),
		UxyVec(kArrSize, 0.0),
		UxyWVec(kArrSize, 0.0),
		UxyEVec(kArrSize, 0.0),
		VxxVec(kArrSize, 0.0),
		VyyVec(kArrSize, 0.0),
		dVConst(kArrSize, 0.0);
	std::vector<double> iRhoHVec(kArrSize, 0.0),
		iRhoVVec(kArrSize, 0.0);

	std::vector<double> muV_W(kArrSize, 0.0),
		muV_E(kArrSize, 0.0),
		muV_S(kArrSize, 0.0),
		muV_N(kArrSize, 0.0),
		muU_S(kArrSize, 0.0),
		muU_N(kArrSize, 0.0),
		checkDiv(kArrSize, 0.0);

	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			dV[idx(i, j)] = 0.0;

		return dV;
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

	// level set
	double theta = 0.0;
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
	double lsVWHalf = 0.0, lsVEHalf = 0.0, lsM = 0.0, lsS = 0.0;
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
	// need for v_xx, v_yy
	double lsVW = 0.0, lsVE = 0.0, lsVS = 0.0, lsVN = 0.0, lsVM = 0.0;
	// need for u_xy
	double lsUM = 0.0, lsUS = 0.0, lsUE = 0.0, lsUE_S = 0.0;
	// jump condition
	double JVW = 0.0, JVE = 0.0, JVS = 0.0, JVN = 0.0, JVM = 0.0;
	double JUM = 0.0, JUS = 0.0, JUE = 0.0, JUE_S = 0.0;
	double uS = 0.0, uM = 0.0, uE = 0.0, uE_S = 0.0;
	double vW = 0.0, vE = 0.0, vS = 0.0, vN = 0.0, vM = 0.0;
	double muV_R_W = 0.0, muV_R_E = 0.0, muV_Z_S = 0.0, muV_Z_N = 0.0;
	double muU_Z_W = 0.0, muU_Z_E = 0.0;
	double rhoV_R_W = 0.0, rhoV_R_E = 0.0, rhoV_Z_S = 0.0, rhoV_Z_N = 0.0;
	double rhoU_Z_W = 0.0, rhoU_Z_E = 0.0;
	double rV_R_W = 0.0, rV_R_E = 0.0, rV = 0.0;
	double visR = 0.0, visZ = 0.0;
	const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

	// effective Jump condition, effective v (vEff), effective mu (muEff), effective rho (rhoEff),
	double JEff = 0.0, JO = 0.0, vEff = 0.0, uEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		muV_R_W = 0.0; muV_R_E = 0.0; muV_Z_S = 0.0; muV_Z_N = 0.0;
		muU_Z_W = 0.0; muU_Z_E = 0.0;
		visR = 0.0, visZ = 0.0;
		
		lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
		lsVE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
		lsVM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
		lsVS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
		lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
		lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
		lsUS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
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
		lsVWHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
		lsVEHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i + 1, j)] + ls[idx(i + 1, j - 1)] + ls[idx(i, j - 1)]);
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		/*
		// V part
		if (lsVW >= 0 && lsVM >= 0) {
			muV_R_W = kMuH * (vM - vW) / kDr;
			muV_W[idx(i, j)] = muV_R_W;
		}
		else if (lsVW <= 0 && lsVM <= 0) {
			muV_R_W = kMuL * (vM - vW) / kDr;
			muV_W[idx(i, j)] = muV_R_W;
		}
		else if (lsVW > 0 && lsVM <= 0) {
			// interface lies between v[i, j] and v[i - 1, j]
			theta = std::fabs(lsVW) / (std::fabs(lsVW) + std::fabs(lsVM));
			// |(lsVW)| === inside(+) === |(interface)| ===   outside(-)    === |(lsVM)|
			// |(lsVW)| === theta * d === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j - 1)]) * 0.5;
			JVM = (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_R_W = muEff * (vM - vW) / kDr
				- muEff * JEff * theta / kMuH;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDr;
			resJ21W[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsVW <= 0 && lsVM > 0) {
			// interface lies between v[i, j] and v[i - 1, j]
			theta = std::fabs(lsVW) / (std::fabs(lsVW) + std::fabs(lsVM));
			// |(lsVW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsVM)|
			// |(lsVW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVW = (m_J21[idx(i - 1, j)] + m_J21[idx(i - 1, j - 1)]) * 0.5;
			JVM = (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVW;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_R_W = muEff * (vM - vW) / kDr
				+ muEff * JEff * theta / kMuL;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDr;
			resJ21W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
			
		if (lsVM >= 0 && lsVE >= 0) {
			muV_R_E = kMuH * (vE - vM) / kDr;
			muV_E[idx(i, j)] = muV_R_E;
		}
		else if (lsVM <= 0 && lsVE <= 0) {
			muV_R_E = kMuL * (vE - vM) / kDr;
			muV_E[idx(i, j)] = muV_R_E;
		}
		else if (lsVM > 0 && lsVE <= 0) {
			// interface lies between v[i, j] and v[i + 1, j]
			theta = std::fabs(lsVE) / (std::fabs(lsVE) + std::fabs(lsVM));
			// |(lsVM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsVE)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsVE)|
			JVM = (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JVE = (m_J21[idx(i + 1, j)] + m_J21[idx(i + 1, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVE;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_R_E = muEff * (vE - vM) / kDr
				+ muEff * JEff * theta / kMuL;
			muV_E[idx(i, j)] = muEff * (vE - vM) / kDr;
			resJ21E[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		else if (lsVM <= 0 && lsVE > 0) {
			// interface lies between v[i, j] and v[i + 1, j]
			theta = std::fabs(lsVE) / (std::fabs(lsVE) + std::fabs(lsVM));
			// |(lsvM)| ===    outside(-)   === |(interface)| === inside(+)  === |(lsvE)|
			// |(lsvM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsvE)|
			JVM = (m_J21[idx(i, j)] + m_J21[idx(i, j - 1)]) * 0.5;
			JVE = (m_J21[idx(i + 1, j)] + m_J21[idx(i + 1, j - 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVE;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_R_E = muEff * (vE - vM) / kDr
				- muEff * JEff * theta / kMuH;
			muV_E[idx(i, j)] = muEff * (vE - vM) / kDr;
			resJ21W[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		
		if (lsVS > 0 && lsVM > 0) {
			muV_Z_S = kMuH * (vM - vS) / kDz;
			muV_S[idx(i, j)] = muV_Z_S;
		}
		else if (lsVS <= 0 && lsVM <= 0) {
			muV_Z_S = kMuL * (vM - vS) / kDz;
			muV_S[idx(i, j)] = muV_Z_S;
		}
		else if (lsVS > 0 && lsVM <= 0) {
			// interface lies between v[i, j] and v[i, j - 1]
			theta = std::fabs(lsVS) / (std::fabs(lsVS) + std::fabs(lsVM));
			// |(lsVS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsVM)|
			// |(lsVS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVS = (m_J22[idx(i, j - 2)] + m_J22[idx(i, j - 1)]) * 0.5;
			JVM = (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVS;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_Z_S = muEff * (vM - vS) / kDz
				- muEff * JEff * theta / kMuH;
			muV_S[idx(i, j)] = muEff * (vM - vS) / kDz;
			resJ22S[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsVS <= 0 && lsVM > 0) {
			// interface lies between v[i, j] and v[i,j  - 1]
			theta = std::fabs(lsVS) / (std::fabs(lsVS) + std::fabs(lsVM));
			// |(lsVS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsVM)|
			// |(lsVS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVM)|
			JVS = (m_J22[idx(i, j - 2)] + m_J22[idx(i, j - 1)]) * 0.5;
			JVM = (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVS;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_Z_S = muEff * (vM - vS) / kDz
				+ muEff * JEff * theta / kMuL;
			muV_S[idx(i, j)] = muEff * (vM - vS) / kDz;
			resJ22S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsVM > 0 && lsVN > 0) {
			muV_Z_N = kMuH * (vN - vM) / kDz;
			muV_N[idx(i, j)] = muV_Z_N;
		}
		else if (lsVM <= 0 && lsVN <= 0) {
			muV_Z_N = kMuL * (vN - vM) / kDz;
			muV_N[idx(i, j)] = muV_Z_N;
		}
		else if (lsVM > 0 && lsVN <= 0) {
			// interface lies between v[i, j] and v[i, j + 1]
			theta = std::fabs(lsVN) / (std::fabs(lsVN) + std::fabs(lsVM));
			// |(lsVM)| ===    inside(+)    === |(interface)| === outside(-) === |(lsVN)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsVN)|
			JVM = (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JVN = (m_J22[idx(i, j)] + m_J22[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVN;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muV_Z_N = muEff * (vN - vM) / kDz
				+ muEff * JEff * theta / kMuL;
			muV_N[idx(i, j)] = muEff * (vN - vM) / kDz;
			resJ22N[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		else if (lsVM <= 0 && lsVN > 0) {
			// interface lies between v[i, j] and v[i, j + 1]
			theta = std::fabs(lsVN) / (std::fabs(lsVN) + std::fabs(lsVM));
			// |(lsVM)| ===    outside(-)   === |(interface)| === inside(+) === |(lsVN)|
			// |(lsVM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsVN)|
			JVM = (m_J22[idx(i, j - 1)] + m_J22[idx(i, j)]) * 0.5;
			JVN = (m_J22[idx(i, j)] + m_J22[idx(i, j + 1)]) * 0.5;
			JEff = theta * JVM + (1.0 - theta) * JVN;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muV_Z_N = muEff * (vN - vM) / kDz
				- muEff * JEff * theta / kMuH;
			muV_N[idx(i, j)] = muEff * (vN - vM) / kDz;
			resJ22N[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		*/
		// U part
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
/*
		if (lsUS > 0 && lsUM > 0) {
			muU_Z_W = kMuH * (uM - uS) / kDz;
			muU_S[idx(i, j)] = muU_Z_W;
		}
		else if (lsUS <= 0 && lsUM <= 0) {
			muU_Z_W = kMuL * (uM - uS) / kDz;
			muU_S[idx(i, j)] = muU_Z_W;
		}
		else if (lsUS > 0 && lsUM <= 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = std::fabs(lsUS) / (std::fabs(lsUS) + std::fabs(lsUM));
			// |(lsUS)| === inside(+) === |(interface)| ===    outside(-)   === |(lsUM)|
			// |(lsUS)| === theta * d === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_Z_W = muEff * (uM - uS) / kDz
				- muEff * JEff * theta / kMuH;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDz;
			resJ12S[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsUS <= 0 && lsUM > 0) {
			// interface lies between u[i, j] and u[i, j - 1]
			theta = std::fabs(lsUS) / (std::fabs(lsUS) + std::fabs(lsUM));
			// |(lsUS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsUM)|
			// |(lsUS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUM)|
			JUS = (m_J12[idx(i - 1, j - 1)] + m_J12[idx(i, j - 1)]) * 0.5;
			JUM = (m_J12[idx(i - 1, j)] + m_J12[idx(i, j)]) * 0.5;
			JEff = theta * JUM + (1.0 - theta) * JUS;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_Z_W = muEff * (uM - uS) / kDz
				+ muEff * JEff * theta / kMuL;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDz;
			resJ12S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}

		if (lsUE_S > 0 && lsUE > 0) {
			muU_Z_E = kMuH * (uE - uE_S) / kDz;
			muU_N[idx(i, j)] = muU_Z_E;
		}
		else if (lsUE_S <= 0 && lsUE <= 0) {
			muU_Z_E = kMuL * (uE - uE_S) / kDz;
			muU_N[idx(i, j)] = muU_Z_E;
		}
		else if (lsUE_S > 0 && lsUE <= 0) {
			// interface lies between u[i + 1, j] and u[i + 1, j - 1]
			theta = std::fabs(lsUE_S) / (std::fabs(lsUE) + std::fabs(lsUE_S));
			// |(lsUE_S)| === outside(-) === |(interface)| ===     inside(+)   === |(lsUE)|
			// |(lsUE_S)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUE)|
			JUE_S = (m_J12[idx(i, j - 1)] + m_J12[idx(i + 1, j - 1)]) * 0.5;
			JUE = (m_J12[idx(i, j)] + m_J12[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUE + (1.0 - theta) * JUE_S;
			muEff = kMuL * kMuH / (kMuL * theta + kMuH * (1.0 - theta));
			muU_Z_E = muEff * (uE - uE_S) / kDz
				- muEff * JEff * theta / kMuH;
			muU_N[idx(i, j)] = muEff * (uE - uE_S) / kDz;
			resJ12N[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		else if (lsUE_S <= 0 && lsUE > 0) {
			// interface lies between u[i + 1, j] and u[i + 1, j - 1]
			theta = std::fabs(lsUE_S) / (std::fabs(lsUE) + std::fabs(lsUE_S));
			// |(lsUE_S)| === inside(+) === |(interface)| ===    outside(-)   === |(lsUE)|
			// |(lsUE_S)| === theta * d === |(interface)| === (1 - theta) * d === |(lsUE)|
			JUE_S = (m_J12[idx(i, j - 1)] + m_J12[idx(i + 1, j - 1)]) * 0.5;
			JUE = (m_J12[idx(i, j)] + m_J12[idx(i + 1, j)]) * 0.5;
			JEff = theta * JUE + (1.0 - theta) * JUE_S;
			muEff = kMuH * kMuL / (kMuH * theta + kMuL * (1.0 - theta));
			muU_Z_E = muEff * (uE - uE_S) / kDz
				+ muEff * JEff * theta / kMuL;
			muU_N[idx(i, j)] = muEff * (uE - uE_S) / kDz;
			resJ12N[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		*/
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
		double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;
		if (lsVWHalf >= 0 && lsVEHalf >= 0) {
			rhoEffWE = kRhoH;
			iRhoEffWE = 1.0 / kRhoH;
		}
		else if (lsVWHalf <= 0 && lsVEHalf <= 0) {
			rhoEffWE = kRhoL;
			iRhoEffWE = 1.0 / kRhoL;
		}
		else if (lsVWHalf > 0 && lsVEHalf <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			// |(lsVWHalf)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsVEHalf)|
			// |(lsVWHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
			// rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
			rhoEffWE = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsVWHalf <= 0 && lsVEHalf > 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
			// |(lsVWHalf)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsVEHalf)|
			// |(lsVWHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
			// rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
			rhoEffWE = kRhoL * theta + kRhoH * (1.0 - theta);
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
		else if (lsS > 0 && lsM <= 0) {
			// interface lies between u[i - 1, j] and u[i, j]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			// rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
			rhoEffSN = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
		}
		else if (lsS <= 0 && lsM > 0) {
			// interface lies between lsS and lsM
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			// rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
			rhoEffSN = kRhoL * theta + kRhoH * (1.0 - theta);
			iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
		}
		
		muU_Z_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (uM - uS) / kDz;
		muU_Z_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
			* (uE - uE_S) / kDz;
		muV_R_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (vM - vW) / kDr;
		muV_R_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
			* (vE - vM) / kDr;
		muV_Z_S = mu[idx(i, j - 1)] * (vM - vS) / kDz;
		muV_Z_N = mu[idx(i, j)] * (vN - vM) / kDz;

		muV_W[idx(i, j)] = muV_R_W;
		muV_E[idx(i, j)] = muV_R_E;
		muV_S[idx(i, j)] = muV_Z_S;
		muV_N[idx(i, j)] = muV_Z_N;
		muU_S[idx(i, j)] = muU_Z_W;
		muU_N[idx(i, j)] = muU_Z_E;

		rV_R_W = kBaseR + (i - kNumBCGrid) * kDr;
		rV_R_E = kBaseR + (i + 1 - kNumBCGrid) * kDr;
		rV = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;

		visR = (muU_Z_E - muU_Z_W) / kDr + 1.0 / rV * (rV_R_E * muV_R_E - rV_R_W * muV_R_W) / kDr;
		visZ = 2.0 * (muV_Z_N - muV_Z_S) / kDz;
		visR -=
			0.5 * (mu[idx(i, j)] + mu[idx(i, j - 1)]) / rV *
			(0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]) - 0.5 * (u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / kDz;
		
		dV[idx(i, j)] = visR * iRhoEffWE + visZ * iRhoEffSN;
		
		iRhoHVec[idx(i, j)] = iRhoEffWE;
		iRhoVVec[idx(i, j)] = iRhoEffSN;
		
		visRVec[idx(i, j)] = visR * iRhoEffWE;
		visZVec[idx(i, j)] = visZ * iRhoEffSN;
		visR2Vec[idx(i, j)] = visR;
		visZ2Vec[idx(i, j)] = visZ;
		assert(dV[idx(i, j)] == dV[idx(i, j)]);
		if (std::isnan(dV[idx(i, j)]) || std::isinf(dV[idx(i, j)])) {
			std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << dV[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	std::ofstream outF;
	std::string fname("VisV_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"KAPPA\", \"mu\", \"dV\",  \"muU_N - muU_S\",  \"muV_E -muV_W\",\"muV_N - muV_S\", \"iRhoEffWE\", \"iRhoEffSN\", \"visR\", \"visZ\", \"visRNoRho\", \"visZNoRho\" " << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNr << std::string(", J=") << kNz
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			<< static_cast<double>(u[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::string(",")
			<< static_cast<double>(mu[idx(i, j)]) << std::string(",")
			<< static_cast<double>(dV[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_N[idx(i, j)] - muU_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_E[idx(i, j)] - muV_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_N[idx(i, j)] - muV_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoHVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoVVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visRVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visZVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visR2Vec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visZ2Vec[idx(i, j)])
			<< std::endl;

	outF.close();
	
	return dV;	
}

std::vector<double> MACSolver2DAxisym::AddGravityFU() {
	std::vector<double> gU(kArrSize, 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM2DAXISYM::R) {
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			gU[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			gU[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gU;
}

std::vector<double> MACSolver2DAxisym::AddGravityFV() {
	std::vector<double> gV(kArrSize, 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM2DAXISYM::Z) {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			gV[idx(i, j)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			gV[idx(i, j)] = -1.0 / kFr;
		}
	}

	return gV;
}

int MACSolver2DAxisym::GetIntermediateVel(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	std::vector<double>& uhat, std::vector<double>& vhat, const std::vector<double>& H) {
		
	// Update rhs
	if (kTimeOrder == TIMEORDERENUM::EULER) {
		std::vector<double> FU1(kArrSize, 0.0);
		std::vector<double> FV1(kArrSize, 0.0);

		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
		
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			uhat[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			vhat[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
	}
	else if (kTimeOrder == TIMEORDERENUM::RK2) {
		std::vector<double> FU1(kArrSize, 0.0),
			FV1(kArrSize, 0.0),
			FU2(kArrSize, 0.0),
			FV2(kArrSize, 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : L(u^(1))
		// FU2 & FV2 : u^(2)
		ApplyBC_U_2D(FU2);
		ApplyBC_V_2D(FV2);
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			uhat[idx(i, j)] = 0.5 * u[idx(i, j)] + 0.5 * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			vhat[idx(i, j)] = 0.5 * v[idx(i, j)] + 0.5 * FV2[idx(i, j)];
		}
	}
	else if (kTimeOrder == TIMEORDERENUM::RK3) {
		std::vector<double> FU1(kArrSize, 0.0),
			FV1(kArrSize, 0.0),
			FU2(kArrSize, 0.0),
			FV2(kArrSize, 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : L(u^(1))
		// FU2 & FV2 : u^(1) + \delta t L (u^(1))
		ApplyBC_U_2D(FU2);
		ApplyBC_V_2D(FV2);
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		
		// FU1 & FV1 : u^(2)
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU1[idx(i, j)] = 0.75 * u[idx(i, j)] + 0.25 * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV1[idx(i, j)] = 0.75 * v[idx(i, j)] + 0.25 * FV2[idx(i, j)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);

		// FU2 & FV2 : L(u^(2))
		// FU1 & FV1 : u^(2) + \delta t L (u^(2))
		ApplyBC_U_2D(FU1);
		ApplyBC_V_2D(FV1);
		FU2 = UpdateFU(LSolver, ls, FU1, FV1, H);
		FV2 = UpdateFV(LSolver, ls, FU1, FV1, H);
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			FU1[idx(i, j)] = FU1[idx(i, j)] + m_dt * FU2[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			FV1[idx(i, j)] = FV1[idx(i, j)] + m_dt * FV2[idx(i, j)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);

		// FU1 & FV1 : u^(2) + \delta t L (u^(2))
		// FU2 & FV2 : doesn't need. set to zero.
		for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
			uhat[idx(i, j)] = 1.0 / 3.0 * u[idx(i, j)] + 2.0 / 3.0 * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
			vhat[idx(i, j)] = 1.0 / 3.0 * v[idx(i, j)] + 2.0 / 3.0 * FV1[idx(i, j)];
		}
	}

	for (int i = kNumBCGrid + 1; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		if (std::isnan(uhat[idx(i, j)]) || std::isinf(uhat[idx(i, j)])) {
			std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << uhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		if (std::isnan(vhat[idx(i, j)]) || std::isinf(vhat[idx(i, j)])) {
			std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << vhat[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	
	return 0;
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
	std::vector<double> iRhoW(kArrSize, 0.0),
		iRhoE(kArrSize, 0.0),
		iRhoS(kArrSize, 0.0),
		iRhoN(kArrSize, 0.0);
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
	double thetaH = 0.0, theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, kappaEff = 0.0, rhoEff = 0.0, nXEff = 0.0, nYEff = 0.0;
	double rPM = 0.0, rEff = 0.0, rPEHalf = 0.0, rPWHalf = 0.0;
	
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
	
	std::vector<double> dudR(kArrSize, 0.0),
		dudZ(kArrSize, 0.0),
		dvdR(kArrSize, 0.0),
		dvdZ(kArrSize, 0.0),
		dldR(kArrSize, 0.0),
		dldZ(kArrSize, 0.0);
	std::vector<double> U_PGrid(kArrSize, 0.0),
		V_PGrid(kArrSize, 0.0);

	double dudRW = 0.0, dudRE = 0.0, dudZS = 0.0, dudZN = 0.0;
	double dvdRW = 0.0, dvdRE = 0.0, dvdZS = 0.0, dvdZN = 0.0;
	// stored coef for A matrix, Dictionary but it is ordered
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++) {
		ps[idx(i, j)] = 0.0;
		rhs[idx(i, j)] = 0.0;
	}
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		lsW = lsB[idx(i - 1, j)];
		lsE = lsB[idx(i + 1, j)];
		lsM = lsB[idx(i, j)];
		lsS = lsB[idx(i, j - 1)];
		lsN = lsB[idx(i, j + 1)];

		// At P grid
		dudR[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dudZ[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDz);
		dvdR[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dvdZ[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDz);

		// Jump occurs when computing dudR
		if (lsW * lsE < 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudRW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dudRE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDr);
			dudR[idx(i, j)] = dudRE * theta + dudRW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdRW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dvdRE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDr);
			dvdR[idx(i, j)] = dvdRE * theta + dvdRW * (1.0 - theta);
		}

		if (lsS * lsN < 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudZS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dudZN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDz);
			dudZ[idx(i, j)] = dudZN * theta + dudZS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdZS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dvdZN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDz);
			dvdZ[idx(i, j)] = dvdZN * theta + dvdZS * (1.0 - theta);
		}
		
		dldR[idx(i, j)] = (lsE - lsW) / (2.0 * kDr);
		dldZ[idx(i, j)] = (lsN - lsS) / (2.0 * kDz);
		if (std::fabs(dldR[idx(i, j)]) < eps) {
			dldR[idx(i, j)] = (lsE - lsM) / kDr;
			if (std::fabs(dldR[idx(i, j)]) < eps)
				dldR[idx(i, j)] = (lsM - lsW) / kDr;
		}
		if (std::fabs(dldZ[idx(i, j)]) < eps) {
			dldZ[idx(i, j)] = (lsN - lsM) / kDz;
			if (std::fabs(dldZ[idx(i, j)]) < eps)
				dldZ[idx(i, j)] = (lsM - lsS) / kDz;
		}
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
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
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rPEHalf = kBaseR + (i + 1.0 - kNumBCGrid) * kDr;
		rPWHalf = kBaseR + (i - kNumBCGrid) * kDr;

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nXW = dldR[idx(i - 1, j)] / (std::sqrt(std::pow(dldR[idx(i - 1, j)], 2.0) + std::pow(dldZ[idx(i - 1, j)], 2.0)) + eps);
		nYW = dldZ[idx(i - 1, j)] / (std::sqrt(std::pow(dldR[idx(i - 1, j)], 2.0) + std::pow(dldZ[idx(i - 1, j)], 2.0)) + eps);
		nXE = dldR[idx(i + 1, j)] / (std::sqrt(std::pow(dldR[idx(i + 1, j)], 2.0) + std::pow(dldZ[idx(i + 1, j)], 2.0)) + eps);
		nYE = dldZ[idx(i + 1, j)] / (std::sqrt(std::pow(dldR[idx(i + 1, j)], 2.0) + std::pow(dldZ[idx(i + 1, j)], 2.0)) + eps);
		nXS = dldR[idx(i, j - 1)] / (std::sqrt(std::pow(dldR[idx(i, j - 1)], 2.0) + std::pow(dldZ[idx(i, j - 1)], 2.0)) + eps);
		nYS = dldZ[idx(i, j - 1)] / (std::sqrt(std::pow(dldR[idx(i, j - 1)], 2.0) + std::pow(dldZ[idx(i, j - 1)], 2.0)) + eps);
		nXN = dldR[idx(i, j + 1)] / (std::sqrt(std::pow(dldR[idx(i, j + 1)], 2.0) + std::pow(dldZ[idx(i, j + 1)], 2.0)) + eps);
		nYN = dldZ[idx(i, j + 1)] / (std::sqrt(std::pow(dldR[idx(i, j + 1)], 2.0) + std::pow(dldZ[idx(i, j + 1)], 2.0)) + eps);
		nXM = dldR[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps);
		nYM = dldZ[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps);
		aEff = 0.0;
		
		// [a]_\Gamma = a^+ - a^- = a^inside - a^outside
		// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// (P_M - P_W)/kDr appears
		// if P_M == P_+, P_W == P_-, a^+ means all terms related to P_+, P_W changed to P_M related terms
		// P_W = P_M -(2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa)
		/*
		// FW
		if (lsW <= 0.0 && lsM <= 0.0) {
			// one fluid, x direction
			FW = 0.0;
			iRhoW[idx(i, j)] = -1.0 / kRhoL;
		}
		else if (lsW > 0.0 && lsM > 0.0) {
			// one fluid, x direction
			FW = 0.0;
			iRhoW[idx(i, j)] = -1.0 / kRhoH;
		}
		else if (lsW > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)      === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
			// b always zero when solving level set (dealing with surface tension)
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i - 1, j)] * nXW + dudZ[idx(i - 1, j)] * nYW) * nXW
				+ (dvdR[idx(i - 1, j)] * nXW + dvdZ[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoW[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FW = aEff * iRhoW[idx(i, j)] / (kDr * kDr);
		}
		else if (lsW <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i - 1, j)] * nXW + dudZ[idx(i - 1, j)] * nYW) * nXW
				+ (dvdR[idx(i - 1, j)] * nXW + dvdZ[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoW[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FW = -aEff * iRhoW[idx(i, j)] / (kDr * kDr);
		}
		aWVec[idx(i, j)] = aEff; aEff = 0.0;
		
		// FE
		// p_E - p_M = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \kSigma \kappa
		if (lsE <= 0.0 && lsM <= 0.0) {
			// one fluid, x direction
			FE = 0.0;
			iRhoE[idx(i, j)] = -1.0 / kRhoL;
		}
		else if (lsE > 0.0 && lsM > 0.0) {
			// one fluid, x direction
			FE = 0.0;
			iRhoE[idx(i, j)] = -1.0 / kRhoH;
		}
		else if (lsE > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i, j] and ls[i + 1, j]
			theta = std::fabs(lsE) / (std::fabs(lsE) + std::fabs(lsM));
			// |(lsM)| ===   outside(-)    === |(interface)| === inside(+)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d  === |(lsE)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i + 1, j)] * nXE + dudZ[idx(i + 1, j)] * nYE) * nXE
				+ (dvdR[idx(i + 1, j)] * nXE + dvdZ[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = aM * theta + aE * (1.0 - theta);
			iRhoE[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FE = aEff * iRhoE[idx(i, j)] / (kDr * kDr);
		}
		else if (lsE <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i + 1, j]
			theta = std::fabs(lsE) / (std::fabs(lsE) + std::fabs(lsM));
			// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i + 1, j)] * nXE + dudZ[idx(i + 1, j)] * nYE) * nXE
				+ (dvdR[idx(i + 1, j)] * nXE + dvdZ[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = aM * theta + aE * (1.0 - theta);
			iRhoE[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FE = -aEff * iRhoE[idx(i, j)] / (kDr * kDr);
		}
		aEVec[idx(i, j)] = aEff; aEff = 0.0;
		// FS
		if (lsS <= 0.0 && lsM <= 0.0) {
			// one fluid, x direction
			FS = 0.0;
			iRhoS[idx(i, j)] = -1.0 / kRhoL;
		}
		else if (lsS > 0.0 && lsM > 0.0) {
			// one fluid, x direction
			FS = 0.0;
			iRhoS[idx(i, j)] = -1.0 / kRhoH;
		}
		else if (lsS > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i, j] and ls[i, j - 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j - 1)] * nXS + dudZ[idx(i, j - 1)] * nYS) * nXS
				+ (dvdR[idx(i, j - 1)] * nXS + dvdZ[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			// rhoS[idx(i, j)] = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoS[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			// FS = aEff / (rhoS[idx(i, j)] * kDz * kDz);
			FS = aEff * iRhoS[idx(i, j)] / (kDz * kDz);
		}
		else if (lsS <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i, j - 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j - 1)] * nXS + dudZ[idx(i, j - 1)] * nYS) * nXS
				+ (dvdR[idx(i, j - 1)] * nXS + dvdZ[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			// rhoS[idx(i, j)] = kRhoL * theta + kRhoH * (1.0 - theta);
			iRhoS[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			// FS = -aEff / (rhoS[idx(i, j)] * kDz * kDz);
			FS = -aEff * iRhoS[idx(i, j)] / (kDz * kDz);
		}
		// if (FS != 0.0)
		// 	std::cout << "FS " << i << " " << j << " " << FS << " " << aEff << " " <<aEff / (m_dt * kSigma) << " " << iRhoS[idx(i, j)] << std::endl;

		aSVec[idx(i, j)] = aEff; aEff = 0.0;
		// FN
		if (lsN <= 0.0 && lsM <= 0.0) {
			// one fluid, x direction
			FN = 0.0;
			iRhoN[idx(i, j)] = -1.0 / kRhoL;
		}
		else if (lsN > 0.0 && lsM > 0.0) {
			// one fluid, x direction
			FN = 0.0;
			iRhoN[idx(i, j)] = -1.0 / kRhoH;
		}
		else if (lsN > 0.0 && lsM <= 0.0) {
			// interface lies between ls[i, j] and ls[i, j + 1]
			theta = std::fabs(lsN) / (std::fabs(lsN) + std::fabs(lsM));
			// |(lsM)| ===    outside      === |(interface)| ===   inside  === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j + 1)] * nXN + dudZ[idx(i, j + 1)] * nYN) * nXN
				+ (dvdR[idx(i, j + 1)] * nXN + dvdZ[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = aM * theta + aN * (1.0 - theta);
			iRhoN[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FN = aEff * iRhoN[idx(i, j)] / (kDz * kDz);
		}
		else if (lsN <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i, j + 1]
			theta = std::fabs(lsN) / (std::fabs(lsN) + std::fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j + 1)] * nXN + dudZ[idx(i, j + 1)] * nYN) * nXN
				+ (dvdR[idx(i, j + 1)] * nXN + dvdZ[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = aM * theta + aN * (1.0 - theta);
			iRhoN[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FN = -aEff * iRhoN[idx(i, j)] / (kDz * kDz);
		}
		// if (FN != 0.0)
		// 	std::cout << "FN " << i << " " << j << " " << FN << " " << aEff << " " <<aEff / (m_dt * kSigma) << " " << iRhoN[idx(i, j)] << std::endl;
		aNVec[idx(i, j)] = aEff; aEff = 0.0;
	*/
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
		if (theta == 0.0)
			rEff = rPWHalf;
		else
			rEff = rPM - kDr * (1.0 - theta);
		// coefficient
		iRhoW[idx(i, j)] = rEff / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0;
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FW = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		FW *= iRhoW[idx(i, j)] / (kDr * kDr);
		aWVec[idx(i, j)] = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid cell adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsE >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsM < 0.0 && lsE < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else if (lsM >= 0.0 && lsE < 0.0) {
			thetaH = std::fabs(lsM) / (std::fabs(lsM) + std::fabs(lsE));
			theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
		}
		else if (lsM < 0.0 && lsE >= 0.0) {
			thetaH = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
			theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		if (theta == 0.0)
			rEff = rPEHalf;
		else
			rEff = rPM + kDr * (1.0 - theta);
		iRhoE[idx(i, j)] = rEff / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i + 1, j)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FE = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);
		FE *= iRhoE[idx(i, j)] / (kDr * kDr);
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
		iRhoS[idx(i, j)] = rPM / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FS = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
		FS *= iRhoS[idx(i, j)] / (kDz * kDz);
		
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
		iRhoN[idx(i, j)] = rPM / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j + 1)] * (1.0 - theta);
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FN = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1)] - H[idx(i, j)]);
		FN *= iRhoN[idx(i, j)] / (kDz * kDz);
		
		// poisson equation form should be -\beta \nabla p = f
		// iRhoEff has already negative value of rhoEff, then it is just a coefficient.
		// For discretization of pressure gradient, additional term is negative and it goes to RHS
		// Then it is a positive value and don't worrry about the sign
		rhs[idx(i, j)] -= FW + FE + FS + FN;
		FWVec[idx(i, j)] = FW;
		FEVec[idx(i, j)] = FE;
		FSVec[idx(i, j)] = FS;
		FNVec[idx(i, j)] = FN;
		
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
		// Add starting rowIdx
		ARowIdx.push_back(rowIdx);
		DiagRowIdx.push_back(MRowIdx);

		iRhoS[idx(i, j)] *= -1.0;
		iRhoW[idx(i, j)] *= -1.0;
		iRhoE[idx(i, j)] *= -1.0;
		iRhoN[idx(i, j)] *= -1.0;

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["S"] = iRhoS[idx(i, j)] / (kDz * kDz);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNr;
		AValsDic["W"] = iRhoW[idx(i, j)] / (kDr * kDr);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["C"] = -(iRhoW[idx(i, j)] + iRhoE[idx(i, j)]) / (kDr * kDr)
			- (iRhoS[idx(i, j)] + iRhoN[idx(i, j)]) / (kDz * kDz);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["E"] = iRhoE[idx(i, j)] / (kDr * kDr);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNr;
		AValsDic["N"] = iRhoN[idx(i, j)] / (kDz * kDz);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNr;
		
		if (i == kNumBCGrid && (m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM)) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += iRhoW[idx(i, j)] / (kDr * kDr);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= iRhoW[idx(i, j)] / (kDr * kDr);
			rhs[idx(i, j)] -= iRhoW[idx(i, j)] / (kDr * kDr) * (m_BC->m_BC_DirichletConstantPW);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::PERIODIC) {
			AValsDic["W"] = iRhoW[idx(kNumBCGrid + kNr - 1, j)];
		}

		// East boundary
		if (i == kNumBCGrid + kNr - 1 && (m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM)) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += iRhoE[idx(i, j)] / (kDr * kDr);
		}
		else if (i == kNumBCGrid + kNr - 1 && m_BC->m_BC_PE == BC2D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= iRhoE[idx(i, j)] / (kDr * kDr);
			rhs[idx(i, j)] -= iRhoE[idx(i, j)] / (kDr * kDr) * (m_BC->m_BC_DirichletConstantPE);
		}
		else if (i == kNumBCGrid + kNr - 1 && m_BC->m_BC_PE == BC2D::PERIODIC) {
			AValsDic["E"] = iRhoE[idx(kNumBCGrid, j)];
		}

		if (j == kNumBCGrid && (m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM)) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += iRhoS[idx(i, j)] / (kDz * kDz);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= iRhoS[idx(i, j)] / (kDz * kDz);
			rhs[idx(i, j)] -= iRhoS[idx(i, j)] / (kDz * kDz) * (m_BC->m_BC_DirichletConstantPS);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::PERIODIC) {
			AValsDic["S"] = iRhoS[idx(i, kNumBCGrid + kNz - 1)];
		}

		if (j == kNumBCGrid + kNz - 1 && (m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM)) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += iRhoN[idx(i, j)] / (kDz * kDz);
		}
		else if (j == kNumBCGrid + kNz - 1 && m_BC->m_BC_PN == BC2D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= iRhoN[idx(i, j)] / (kDz * kDz);
			rhs[idx(i, j)] -= iRhoN[idx(i, j)] / (kDz * kDz) * (m_BC->m_BC_DirichletConstantPN);
		}
		else if (j == kNumBCGrid + kNz - 1 && m_BC->m_BC_PN == BC2D::PERIODIC) {
			AValsDic["N"] = iRhoN[idx(i, kNumBCGrid + kNz + 1)];
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
		
		tmpMRowIdx++;
		DiagVals.push_back(AValsDic["C"]);
		DiagCols.push_back(AColsDic["C"]);
		// std::cout << AValsDic["S"] << " " << AValsDic["W"] << " " << AValsDic["C"] << " " << AValsDic["E"] << " " << AValsDic["N"] << std::endl;
		rowIdx += tmpRowIdx;
		MRowIdx += tmpMRowIdx;
		assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
		if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " 
				<< i << " " << j << " " << rhs[idx(i, j)] << std::endl;
			exit(1);
		}
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);
	
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

		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			rhs[idx(i, j)] = -div[idx(i, j)];
		
		m_Poisson->MKL_2FUniform_2D(ps, rhs, kLenR, kLenZ, kDr, kDz, m_BC);
		
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;

	}
	else if (m_PoissonSolverType == POISSONTYPE::CG) {
		// std::cout << "Poisson : CG" << std::endl;
		m_Poisson->CG_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenR, kLenZ, kDr, kDz, m_BC, maxIter);
	
		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
		// std::cout << "Poisson : BiCG" << std::endl;
		m_Poisson->BiCGStab_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenR, kLenZ, kDr, kDz, m_BC, maxIter);

		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::GS) {
		// std::cout << "Poisson : GS" << std::endl;
		m_Poisson->GS_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenR, kLenZ, kDr, kDz, m_BC, maxIter);

		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	
	ApplyBC_P_2D(ps);

	std::ofstream outF;
	std::string fname("P_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"X\", \"Y\", \"FW\",\"RW\", \"AW\",\"FE\",\"RE\", \"AE\", \"FS\", \"FN\", \"F\",  \"RS\", \"RN\", \"LS\", \"PS\", \"HatDIV\", \"RHS\", \"P_W\", \"P_S\",\"KAPPA\"" << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);
	
	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNr + 6 << std::string(", J=") << kNz + 6
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	// for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	// 	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoW[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FEVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoE[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aEVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)] + FEVec[idx(i, j)] + FSVec[idx(i, j)] + FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoS[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoN[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(div[idx(i, j)]) << std::string(",")
			<< static_cast<double>(rhs[idx(i, j)]) << std::string(",")
			<< static_cast<double>((ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDr) << std::string(",")
			<< static_cast<double>((ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDz) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;

	outF.close();
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
	double rPM = 0.0, rPEHalf = 0.0, rPWHalf = 0.0;
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		rPM = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
		rPEHalf = kBaseR + (i + 1.0 - kNumBCGrid) * kDr;
		rPWHalf = kBaseR + (i - kNumBCGrid) * kDr;
		div[idx(i, j)] = (rPEHalf * u[idx(i + 1, j)] - rPWHalf * u[idx(i, j)]) / (rPM * kDr)
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
	double uEff = 0.0, vEff = 0.0, rhoEff = 0.0, theta = 0.0, thetaH = 0.0, iRhoEff = 0.0;
	double nXEff = 0.0, nYEff = 0.0, kappaEff = 0.0;
	double nXW = 0.0, nYW = 0.0, nXS = 0.0, nYS = 0.0, nXM = 0.0, nYM = 0.0;
	double aW = 0.0, aS = 0.0, aM = 0.0, aEff = 0.0;
	const double eps = 1.0e-100;
	std::vector<double> dudR(kArrSize, 0.0),
		dudZ(kArrSize, 0.0),
		dvdR(kArrSize, 0.0),
		dvdZ(kArrSize, 0.0),
		dldR(kArrSize, 0.0),
		dldZ(kArrSize, 0.0);

	std::vector<double> U_PGrid(kArrSize, 0.0),
		V_PGrid(kArrSize, 0.0);
	double dudRW = 0.0, dudRE = 0.0, dudZS = 0.0, dudZN = 0.0;
	double dvdRW = 0.0, dvdRE = 0.0, dvdZS = 0.0, dvdZN = 0.0;

	std::vector<double> nrVec(kArrSize, 0.0),
		nyVec(kArrSize, 0.0);
	std::vector<double> AWW(kArrSize, 0.0),
		AWM(kArrSize, 0.0),
		ASS(kArrSize, 0.0),
		ASM(kArrSize, 0.0);
	std::vector<double> iRhoEffWVec(kArrSize, 0.0),
		iRhoEffSVec(kArrSize, 0.0);

	std::vector<double> aWEff(kArrSize, 0.0),
		aSEff(kArrSize, 0.0);
	
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsE = ls[idx(i + 1, j)];
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		lsN = ls[idx(i, j + 1)];
		
		// At P grid 
		dudR[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dudZ[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDz);
		dvdR[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDr);
		dvdZ[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDz);

		// Jump occurs when computing dudR
		if (lsW * lsE <= 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudRW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dudRE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDr);
			dudR[idx(i, j)] = dudRE * theta + dudRW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdRW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDr);
			dvdRE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDr);
			dvdR[idx(i, j)] = dvdRE * theta + dvdRW * (1.0 - theta);
		}

		if (lsS * lsN <= 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudZS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dudZN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDz);
			dudZ[idx(i, j)] = dudZN * theta + dudZS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdZS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDz);
			dvdZN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDz);
			dvdZ[idx(i, j)] = dvdZN * theta + dvdZS * (1.0 - theta);
		}
		
		dldR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDr);
		dldZ[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDz);

		if (std::fabs(dldR[idx(i, j)]) < eps) {
			dldR[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDr;
			if (std::fabs(dldR[idx(i, j)]) < eps)
				dldR[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDr;
		}
		if (std::fabs(dldZ[idx(i, j)]) < eps) {
			dldZ[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDz;
			if (std::fabs(dldZ[idx(i, j)]) < eps)
				dldZ[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDz;
		}	
	}
	
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
		aW = 0.0;
		aM = 0.0;
		
		nXW = dldR[idx(i - 1, j)] / (std::sqrt(std::pow(dldR[idx(i - 1, j)], 2.0) + std::pow(dldZ[idx(i - 1, j)], 2.0)) + eps * eps);
		nYW = dldZ[idx(i - 1, j)] / (std::sqrt(std::pow(dldR[idx(i - 1, j)], 2.0) + std::pow(dldZ[idx(i - 1, j)], 2.0)) + eps * eps);
		nXM = dldR[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps * eps);
		nYM = dldZ[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps * eps);
		nrVec[idx(i, j)] = nXM;
		nyVec[idx(i, j)] = nYM;
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
		iRhoEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		u[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		u[idx(i, j)] *= iRhoEff / kDr;

		/*
		aEff = 0.0; iRhoEff = 0.0;
		u[idx(i, j)] = 0.0;
		if (lsW >= 0 && lsM >= 0) {
			// rhoEff = kRhoH;
			iRhoEff = 1.0 / kRhoH;
		}
		else if (lsW <= 0 && lsM <= 0) {
			// rhoEff = kRhoL;
			iRhoEff = 1.0 / kRhoL;
		}
		else if (lsW > 0 && lsM <= 0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i - 1, j)] * nXW + dudZ[idx(i - 1, j)] * nYW) * nXW
				+ (dvdR[idx(i - 1, j)] * nXW + dvdZ[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			u[idx(i, j)] = -aEff * iRhoEff / kDr;
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i - 1, j)] * nXW + dudZ[idx(i - 1, j)] * nYW) * nXW
				+ (dvdR[idx(i - 1, j)] * nXW + dvdZ[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			u[idx(i, j)] = aEff * iRhoEff / kDr;
		}

		*/
		aWEff[idx(i, j)] = u[idx(i, j)];
		iRhoEffWVec[idx(i, j)] = iRhoEff;
		u[idx(i, j)] += us[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDr;
	}

	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNz + kNumBCGrid; j++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		aS = 0.0;
		aM = 0.0;
				
		nXS = dldR[idx(i, j - 1)] / (std::sqrt(std::pow(dldR[idx(i, j - 1)], 2.0) + std::pow(dldZ[idx(i, j - 1)], 2.0)) + eps * eps);
		nYS = dldZ[idx(i, j - 1)] / (std::sqrt(std::pow(dldR[idx(i, j - 1)], 2.0) + std::pow(dldZ[idx(i, j - 1)], 2.0)) + eps * eps);
		nXM = dldR[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps * eps);
		nYM = dldZ[idx(i, j)] / (std::sqrt(std::pow(dldR[idx(i, j)], 2.0) + std::pow(dldZ[idx(i, j)], 2.0)) + eps * eps);
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
		iRhoEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid

		v[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
		v[idx(i, j)] *= iRhoEff / kDz;

		/*
		aEff = 0.0; iRhoEff = 0.0;
		v[idx(i, j)] = 0.0;
		if (lsS >= 0 && lsM >= 0) {
			iRhoEff = 1.0 / kRhoH;
		}
		else if (lsS <= 0 && lsM <= 0) {
			iRhoEff = 1.0 / kRhoL;
		}
		else if (lsS > 0 && lsM <= 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| ===  inside(+) === |(interface)| ===     outside(-)  === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j - 1)] * nXS + dudZ[idx(i, j - 1)] * nYS) * nXS
				+ (dvdR[idx(i, j - 1)] * nXS + dvdZ[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			v[idx(i, j)] = -aEff * iRhoEff / kDz;
		}
		else if (lsS <= 0 && lsM > 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j - 1)] * nXS + dudZ[idx(i, j - 1)] * nYS) * nXS
				+ (dvdR[idx(i, j - 1)] * nXS + dvdZ[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudR[idx(i, j)] * nXM + dudZ[idx(i, j)] * nYM) * nXM
				+ (dvdR[idx(i, j)] * nXM + dvdZ[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			v[idx(i, j)] = aEff * iRhoEff / kDz;
		}
		*/

		aSEff[idx(i, j)] = v[idx(i, j)];
		iRhoEffSVec[idx(i, j)] = iRhoEff;
		v[idx(i, j)] += vs[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDz;
	}
	
	std::ofstream outF;
	std::string fname("VUpdate_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"Uhat\", \"Vhat\", \"LS\", \"PS\", \"iRhoEffW\", \"Rho*AWEff\", \"Px\", \"Rho*Px\", \"Px+AW\",\"iRhoEffS\", \"Rho*ASEff\",\"Py\", \"Rho*Py\", \"Py+AS\", \"UmU\", \"VmV\", \"RealDIV\", \"KAPPA\"" << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNr + 2 * kNumBCGrid << std::string(", J=") << kNz + 2 * kNumBCGrid
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	ApplyBC_U_2D(u);
	ApplyBC_V_2D(v);

	//for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
	//	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = 0; j < kNz + 2 * kNumBCGrid; j++)
	for (int i = 0; i < kNr + 2 * kNumBCGrid; i++)
			outF << kBaseR + static_cast<double>(i + 0.5 - kNumBCGrid) * kDr << std::string(",")
			<< kBaseZ + static_cast<double>(j + 0.5 - kNumBCGrid) * kDz << std::string(",")
			<< static_cast<double>((u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5) << std::string(",")
			<< static_cast<double>((v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5) << std::string(",")
			<< static_cast<double>((us[idx(i, j)] + us[idx(i + 1, j)]) * 0.5) << std::string(",")
			<< static_cast<double>((vs[idx(i, j)] + vs[idx(i, j + 1)]) * 0.5) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoEffWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aWEff[idx(i, j)]) << std::string(",")
			<< static_cast<double>(-(ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDr) << std::string(",")
			<< static_cast<double>(-iRhoEffWVec[idx(i, j)] * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDr) << std::string(",")
			<< static_cast<double>(u[idx(i, j)] - us[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoEffSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aSEff[idx(i, j)]) << std::string(",")
			<< static_cast<double>(-(ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDz) << std::string(",")
			<< static_cast<double>(-iRhoEffSVec[idx(i, j)] * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDz) << std::string(",")
			<< static_cast<double>(v[idx(i, j)] - vs[idx(i, j)]) << std::string(",")
			<< static_cast<double>(u[idx(i + 1, j)] - u[idx(i, j)]) / kDr << std::string(",")
			<< static_cast<double>(v[idx(i, j + 1)] - v[idx(i, j)]) / kDz << std::string(",")
			<< static_cast<double>((u[idx(i + 1, j)] - u[idx(i, j)]) / kDr + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDz) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;

	outF.close();
	
	return 0;
}

double MACSolver2DAxisym::UpdateDt(const std::vector<double>& u, const std::vector<double>& v) {
	double uAMax = 0.0;
	double vAMax = 0.0;
	double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0;
	double dt = std::numeric_limits<double>::max();

	// get maximum of absolute value
	for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++) {
		uAMax = std::max(uAMax, std::fabs((u[idx(i + 1, j)] * u[idx(i, j)]) * 0.5));
		vAMax = std::max(vAMax, std::fabs((v[idx(i, j + 1)] * v[idx(i, j)]) * 0.5));
	}
	
	Cefl = uAMax / kDr + vAMax / kDz;
	Vefl = std::max(kMuH, kMuL) * (2.0 / (kDr * kDr) + 2.0 / (kDz * kDz));
	Gefl = std::max(std::sqrt(std::fabs(kG) / kDz), std::sqrt(std::fabs(kG) / kDr));
	
	dt = std::min(dt,
		1.0 / (0.5 * (Cefl + Vefl +
		std::sqrt(std::pow(Cefl + Vefl, 2.0) +
		4.0 * Gefl))));

	if (std::isnan(dt) || std::isinf(dt)) {
		std::cout << "dt nan/inf error : Cefl, Vefl, Gefl : " << Cefl << " " << Vefl << " " << Gefl << std::endl;
	}

	return dt;
}

int MACSolver2DAxisym::SetBC_U_2D(std::string BC_W, std::string BC_E,	std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kNumBCGrid);
	}

	m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2DAxisym::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kNumBCGrid);
	}

	m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int MACSolver2DAxisym::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNr, kNz, kNumBCGrid);
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
			outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"PS\" " << std::endl;
			outF.close();

			outF.open(fname_div.c_str(), std::ios::out);
			outF << "TITLE = DIV" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"LS\", \"DIV\", \"PS\" " << std::endl;
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
			resX(kNr * kNz, 0.0),
			resY(kNr * kNz, 0.0),
			resU(kNr * kNz, 0.0),
			resV(kNr * kNz, 0.0),
			resDiv(kNr * kNz, 0.0),
			resLS(kNr * kNz, 0.0),
			resPs(kNr * kNz, 0.0);

		for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)	 {
			resX[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = kBaseR + (i + 0.5 - kNumBCGrid) * kDr;
			resY[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = kBaseZ + (j + 0.5 - kNumBCGrid) * kDz;
			resU[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
			resV[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = ls[idx(i, j)];
			resPs[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = ps[idx(i, j)];
			resDiv[(i - kNumBCGrid) + kNr * (j - kNumBCGrid)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDr
			 + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDz;
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
		stat = TECDAT142(&ARRSIZEVAL, resDiv.data(), &DIsDouble);
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
