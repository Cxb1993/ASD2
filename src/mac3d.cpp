#include "mac3d.h"

MACSolver3D::MACSolver3D(double Re, double We, double Fr, GAXISENUM GAxis,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoH, double muH,
	int nx, int ny, int nz, double baseX, double baseY, double baseZ, double lenX, double lenY, double lenZ,
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kGAxis(GAxis), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoH(rhoH), kRhoL(rhoH * densityRatio), kRhoRatio(densityRatio),
	kMuH(muH), kMuL(muH * viscosityRatio), kMuRatio(viscosityRatio),
	kNx(nx), kNy(ny), kNz(nz), kBaseX(baseX), kBaseY(baseY), kBaseZ(baseZ), kLenX(lenX), kLenY(lenY), kLenZ(lenZ),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)), kDz(lenZ / static_cast<double>(nz)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
	kArrSize(
	static_cast<int64_t>(kNx + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNy + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNz + 2 * kNumBCGrid)) {

	// positive level set : inside
	// negative level set : outside
	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * 
		std::min(std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)),
			lenZ / static_cast<double>(nz)) / U);
}

MACSolver3D::MACSolver3D(double rhoH, double rhoL, double muH, double muL, double gConstant, GAXISENUM GAxis,
	double L, double U, double sigma, int nx, int ny, int nz, double baseX, double baseY, double baseZ, double lenX, double lenY, double lenZ,
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoH), kMuScale(muH), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoH(rhoH), kRhoL(rhoL), kRhoRatio(rhoL / rhoH), 
	kMuH(muH), kMuL(muL), kMuRatio(muL / muH),
	kRe(rhoH * L * U / muH), kWe(rhoH * L * U * U / sigma), kFr(U * U / (gConstant * L)), kGAxis(GAxis),
	kNx(nx), kNy(ny), kNz(nz), kBaseX(baseX), kBaseY(baseY), kBaseZ(baseZ), kLenX(lenX), kLenY(lenY), kLenZ(lenZ),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)), kDz(lenZ / static_cast<double>(nz)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
	kArrSize(
	static_cast<int64_t>(kNx + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNy + 2 * kNumBCGrid) *
	static_cast<int64_t>(kNz + 2 * kNumBCGrid)) {

	// positive level set : inside
	// negative level set : outside
	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL *
		std::min(std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)),
		lenZ / static_cast<double>(nz)) / U);
}

MACSolver3D::~MACSolver3D() {
	if (m_BC)
		m_BC.reset();
}

// Deallocated automatically
int MACSolver3D::AllocateVariables() {
	m_u = std::vector<double>(kArrSize, 0.0);
	m_v = std::vector<double>(kArrSize, 0.0);
	m_w = std::vector<double>(kArrSize, 0.0);
	
	m_rho = std::vector<double>(kArrSize, 0.0);
	m_kappa = std::vector<double>(kArrSize, 0.0);
	m_mu = std::vector<double>(kArrSize, 0.0);
	
	m_ps = std::vector<double>(kArrSize, 0.0);
	m_p = std::vector<double>(kArrSize, 0.0);

	m_nx = std::vector<double>(kArrSize, 0.0);
	m_ny = std::vector<double>(kArrSize, 0.0);
	m_nz = std::vector<double>(kArrSize, 0.0);

	m_t1x = std::vector<double>(kArrSize, 0.0);
	m_t1y = std::vector<double>(kArrSize, 0.0);
	m_t1z = std::vector<double>(kArrSize, 0.0);

	m_t2x = std::vector<double>(kArrSize, 0.0);
	m_t2y = std::vector<double>(kArrSize, 0.0);
	m_t2z = std::vector<double>(kArrSize, 0.0);

	return 0;
}

std::vector<double> MACSolver3D::UpdateHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside(kArrSize, 0.0);
	const double eps = std::min(std::min(kDx, kDy), kDz) * 1.5;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
	for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
		// inside
		if (ls[idx(i, j, k)] >= 0.0)
			Heaviside[idx(i, j, k)] = 1.0;
		else
			Heaviside[idx(i, j, k)] = 0.0;
	}

	return Heaviside;
}

std::vector<double> MACSolver3D::UpdateSmoothHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside(kArrSize, 0.0);
	const double eps = std::min(std::min(kDx, kDy), kDz) * 1.5;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
	for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
		// inside
		if (ls[idx(i, j, k)] > eps)
			Heaviside[idx(i, j, k)] = 1.0;
		// outside
		else if (ls[idx(i, j, k)] < -eps)
			Heaviside[idx(i, j, k)] = 0.0;
		else
			Heaviside[idx(i, j, k)]
			= 0.5 * (1.0 + ls[idx(i, j, k)] / eps
			+ 1.0 / M_PI * sin(M_PI * ls[idx(i, j, k)] / eps));
	}
		
	return Heaviside;
}

int MACSolver3D::UpdateNTK(const std::shared_ptr<LevelSetSolver3D>& LSolver, const std::vector<double>& ls,
	std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> normalVec,
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

	absDerivLS = LSolver->ENO_DerivAbsLS_3D(ls, ls);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		Q[idx(i, j, k)] = std::fabs(1.0 - absDerivLS[idx(i, j, k)]);
	}
	
	double dLSdX = 0.0, dLSdY = 0.0, dLSdZ = 0.0;
	// determination function
	int Dx = 0.0, Dy = 0.0, Dz = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// complex level set, such as droplet merging
		if (Q[idx(i, j, k)] > eta) {
			// Using Dierctional Direction First
			if (Q[idx(i - 1, j, k)] < eta && Q[idx(i + 1, j, k)] >= eta)
				Dx = -1;
			else if (Q[idx(i - 1, j, k)] >= eta && Q[idx(i + 1, j, k)] < eta)
				Dx = 1;
			else if (Q[idx(i - 1, j, k)] < eta && Q[idx(i, j, k)] < eta && Q[idx(i + 1, j, k)] < eta)
				Dx = 0;
			else if (Q[idx(i - 1, j, k)] >= eta && Q[idx(i, j, k)] >= eta && Q[idx(i + 1, j, k)] >= eta)
				Dx = 0;
			else
				Dx = 2;

			// determination funciton
			if (Q[idx(i, j - 1, k)] < eta && Q[idx(i, j + 1, k)] >= eta)
				Dy = -1;
			else if (Q[idx(i, j - 1, k)] >= eta && Q[idx(i, j + 1, k)] < eta)
				Dy = 1;
			else if (Q[idx(i, j - 1, k)] < eta && Q[idx(i, j, k)] < eta && Q[idx(i, j + 1, k)] < eta)
				Dy = 0;
			else if (Q[idx(i, j - 1, k)] >= eta && Q[idx(i, j, k)] >= eta && Q[idx(i, j + 1, k)] >= eta)
				Dy = 0;
			else
				// undetermined
				Dy = 2;

			// determination funciton
			if (Q[idx(i, j, k - 1)] < eta && Q[idx(i, j, k + 1)] >= eta)
				Dz = -1;
			else if (Q[idx(i, j, k - 1)] >= eta && Q[idx(i, j, k + 1)] < eta)
				Dz = 1;
			else if (Q[idx(i, j, k - 1)] < eta && Q[idx(i, j, k)] < eta && Q[idx(i, j, k + 1)] < eta)
				Dz = 0;
			else if (Q[idx(i, j, k - 1)] >= eta && Q[idx(i, j, k)] >= eta && Q[idx(i, j, k + 1)] >= eta)
				Dz = 0;
			else
				// undetermined
				Dz = 2;

			if (Dx == -1)
				dLSdX = (ls[idx(i, j, k)] - ls[idx(i - 1, j, k)]) / kDx;
			else if (Dx == 1)
				dLSdX = (ls[idx(i + 1, j, k)] - ls[idx(i, j, k)]) / kDx;
			else if (Dx == 0)
				dLSdX = (ls[idx(i + 1, j, k)] - ls[idx(i - 1, j, k)]) / (2.0 * kDx);

			if (Dy == -1)
				dLSdY = (ls[idx(i, j, k)] - ls[idx(i, j - 1, k)]) / kDy;
			else if (Dy == 1)
				dLSdY = (ls[idx(i, j + 1, k)] - ls[idx(i, j, k)]) / kDy;
			else if (Dy == 0)
				dLSdY = (ls[idx(i, j + 1, k)] - ls[idx(i, j - 1, k)]) / (2.0 * kDy);

			if (Dz == -1)
				dLSdZ = (ls[idx(i, j, k)] - ls[idx(i, j, k - 1)]) / kDz;
			else if (Dz == 1)
				dLSdZ = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k)]) / kDz;
			else if (Dz == 0)
				dLSdZ = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k - 1)]) / (2.0 * kDz);

			// Curve Fitting Method
			if (Dx == 2) {
				// Breadfirst Search

				// Curve Fitting

				// Normal Vector

				// Tangent Vector

				// Curvature

			}

		}
		// simple level set
		else {
			dLSdX = (ls[idx(i + 1, j, k)] - ls[idx(i - 1, j, k)]) / (2.0 * kDx);
			dLSdY = (ls[idx(i, j + 1, k)] - ls[idx(i, j - 1, k)]) / (2.0 * kDy);
			dLSdZ = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k - 1)]) / (2.0 * kDz);
			std::get<0>(normalVec)[idx(i, j, k)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + dLSdZ * dLSdZ + eps) * dLSdX;
			std::get<1>(normalVec)[idx(i, j, k)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + dLSdZ * dLSdZ + eps) * dLSdY;
			std::get<2>(normalVec)[idx(i, j, k)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + dLSdZ * dLSdZ + eps) * dLSdZ;
		}
	}

	return 0;
}

int MACSolver3D::UpdateKappa(const std::vector<double>& ls) {
	std::vector<double> dLSdX(kArrSize, 0.0), dLSdY(kArrSize, 0.0), dLSdZ(kArrSize, 0.0);
	std::vector<double> LSSize(kArrSize, 0.0);
	
	const double eps = 1.0e-100;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		dLSdX[idx(i, j, k)] = (ls[idx(i + 1, j, k)] - ls[idx(i - 1, j, k)]) / (2.0 * kDx);
		dLSdY[idx(i, j, k)] = (ls[idx(i, j + 1, k)] - ls[idx(i, j - 1, k)]) / (2.0 * kDy);
		dLSdZ[idx(i, j, k)] = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k - 1)]) / (2.0 * kDz);

		LSSize[idx(i, j, k)] = std::sqrt(dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)]
			+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)]
			+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)]);

		// if size is zero, use forward or backward differencing rather than central differencing
		if (LSSize[idx(i, j, k)] < eps) {
			dLSdX[idx(i, j, k)] = (ls[idx(i + 1, j, k)] - ls[idx(i, j, k)]) / kDx;
			LSSize[idx(i, j, k)] = std::sqrt(dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)]
				+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)]
				+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)]);
			if (LSSize[idx(i, j, k)] < eps) {
				dLSdY[idx(i, j, k)] = (ls[idx(i, j + 1, k)] - ls[idx(i, j, k)]) / kDy;
				LSSize[idx(i, j, k)] = std::sqrt(dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)]
					+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)]
					+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)]);
				if (LSSize[idx(i, j, k)] < eps) {
					dLSdZ[idx(i, j, k)] = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k)]) / kDz;
					LSSize[idx(i, j, k)] = std::sqrt(dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)]
						+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)]
						+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)]);
				}
			}
		}
		if (LSSize[idx(i, j, k)] < eps)
			perror("Div/0 Err in computing kappa");
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		m_kappa[idx(i, j, k)]
			= -(dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)] * (ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy) // phi^2_x \phi_yy
			- 2.0 * dLSdX[idx(i, j, k)] * dLSdY[idx(i, j, k)] * (ls[idx(i + 1, j + 1, k)] - ls[idx(i - 1, j + 1, k)] - ls[idx(i + 1, j - 1, k)] + ls[idx(i - 1, j - 1, k)]) / (4.0 * kDx * kDy) //2 \phi_x \phi_y \phi_yx
			+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)] * (ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx) // phi^2_y \phi_xx
			+ dLSdX[idx(i, j, k)] * dLSdX[idx(i, j, k)] * (ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz) // phi^2_x \phi_zz
			- 2.0 * dLSdX[idx(i, j, k)] * dLSdZ[idx(i, j, k)] * (ls[idx(i + 1, j, k + 1)] - ls[idx(i - 1, j, k + 1)] - ls[idx(i + 1, j, k - 1)] + ls[idx(i - 1, j, k - 1)]) / (4.0 * kDx * kDz) //2 \phi_x \phi_z \phi_xz
			+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)] * (ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx) // phi^2_y \phi_xx
			+ dLSdY[idx(i, j, k)] * dLSdY[idx(i, j, k)] * (ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz) // phi^2_y \phi_zz
			- 2.0 * dLSdY[idx(i, j, k)] * dLSdZ[idx(i, j, k)] * (ls[idx(i, j + 1, k + 1)] - ls[idx(i, j - 1, k + 1)] - ls[idx(i, j + 1, k - 1)] + ls[idx(i, j - 1, k - 1)]) / (4.0 * kDy * kDz) //2 \phi_y \phi_z \phi_yz
			+ dLSdZ[idx(i, j, k)] * dLSdZ[idx(i, j, k)] * (ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDy * kDy)) // phi^2_z \phi_yy
			/ std::pow(LSSize[idx(i, j, k)], 3.0);

		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j, k)] = std::fabs(std::min(std::fabs(m_kappa[idx(i, j, k)]), 1.0 / std::min(std::min(kDx, kDy), kDz)));

		assert(m_kappa[idx(i, j, k)] == m_kappa[idx(i, j, k)]);
	}

	return 0;
}

std::vector<double> MACSolver3D::UpdateFU(const std::shared_ptr<LevelSetSolver3D>& LSolver,
	const std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& H) {
	
	std::vector<double> cU(kArrSize, 0.0);
	std::vector<double> vU(kArrSize, 0.0);
	std::vector<double> gU(kArrSize, 0.0);
	std::vector<double> rhsU(kArrSize, 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v, w);

	// Viscous term
	vU = this->AddViscosityFU(u, v, w, ls, H);

	gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		rhsU[idx(i, j, k)] = -cU[idx(i, j, k)] + vU[idx(i, j, k)] + gU[idx(i, j, k)];
	}
	
	return rhsU;
}

std::vector<double> MACSolver3D::UpdateFV(const std::shared_ptr<LevelSetSolver3D>& LSolver,
	const std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& H) {

	std::vector<double> cV(kArrSize, 0.0);
	std::vector<double> vV(kArrSize, 0.0);
	std::vector<double> gV(kArrSize, 0.0);
	std::vector<double> rhsV(kArrSize, 0.0);

	// Convection term
	cV = this->AddConvectionFV(u, v, w);

	// Viscous term
	vV = this->AddViscosityFV(u, v, w, ls, H);

	gV = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		rhsV[idx(i, j, k)] = -cV[idx(i, j, k)] + vV[idx(i, j, k)] + gV[idx(i, j, k)];
	}
	
	return rhsV;
}

std::vector<double> MACSolver3D::UpdateFW(const std::shared_ptr<LevelSetSolver3D>& LSolver,
	const std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& H) {

	std::vector<double> cW(kArrSize, 0.0);
	std::vector<double> vW(kArrSize, 0.0);
	std::vector<double> gW(kArrSize, 0.0);
	std::vector<double> rhsW(kArrSize, 0.0);

	// Convection term
	cW = this->AddConvectionFW(u, v, w);

	// Viscous term
	vW = this->AddViscosityFW(u, v, w, ls, H);

	gW = this->AddGravityFW();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		rhsW[idx(i, j, k)] = -cW[idx(i, j, k)] + vW[idx(i, j, k)] + gW[idx(i, j, k)];
	}

	return rhsW;
}

std::vector<double> MACSolver3D::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {
	std::vector<double> cU(kArrSize, 0.0);
	std::vector<double> tmpV(kArrSize, 0.0), tmpW(kArrSize, 0.0);
	std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
	std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0), LZM(kArrSize, 0.0);
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		tmpV[idx(i, j, k)] = (v[idx(i, j, k)] + v[idx(i - 1, j, k)] + v[idx(i - 1, j + 1, k)] + v[idx(i, j + 1, k)]) * 0.25;
		tmpW[idx(i, j, k)] = (w[idx(i, j, k)] + w[idx(i - 1, j, k)] + w[idx(i - 1, j, k + 1)] + w[idx(i, j, k + 1)]) * 0.25;
	}
		
	std::vector<double> vecF_UX(kNx + 2 * kNumBCGrid, 0.0), vecF_UY(kNy + 2 * kNumBCGrid), vecF_UZ(kNz + 2 * kNumBCGrid);

	// U : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_UX[i] = u[idx(i, j, k)] * u[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_UX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j, k)] = FXP[i];
			LXM[idx(i, j, k)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_UX.begin(), vecF_UX.end(), 0.0);
	}
	
	// U : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_UY[j] = u[idx(i, j, k)] * tmpV[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_UY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j, k)] = FYP[j];
			LYM[idx(i, j, k)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_UY.begin(), vecF_UY.end(), 0.0);
	}

	// U : Z direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			vecF_UZ[k] = u[idx(i, j, k)] * tmpW[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_UZ, FZP, FZM, kDz, kNz);

		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			LZP[idx(i, j, k)] = FZP[k];
			LZM[idx(i, j, k)] = FZM[k];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_UZ.begin(), vecF_UZ.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, alphaZ = 0.0;

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		alphaX = u[idx(i, j, k)];
		alphaY = tmpV[idx(i, j, k)];
		alphaZ = tmpW[idx(i, j, k)];
		
		cU[idx(i, j, k)]
				= (0.5 * (LXP[idx(i, j, k)] + LXM[idx(i, j, k)])
				+ 0.5 * (LYP[idx(i, j, k)] + LYM[idx(i, j, k)])
				+ 0.5 * (LZP[idx(i, j, k)] + LZM[idx(i, j, k)])
				- alphaX * (0.5 * (LXP[idx(i, j, k)] - LXM[idx(i, j, k)]))
				- alphaY * (0.5 * (LYP[idx(i, j, k)] - LYM[idx(i, j, k)])));
		
		assert(cU[idx(i, j, k)] == cU[idx(i, j, k)]);
		if (std::isnan(cU[idx(i, j, k)]) || std::isinf(cU[idx(i, j, k)])) {
			std::cout << "U-convection term nan/inf error : " << i << " " << j << " " 
			<< cU[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	return cU;
}

std::vector<double> MACSolver3D::AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {
	std::vector<double> cV(kArrSize, 0.0);
	std::vector<double> tmpU(kArrSize, 0.0), tmpW(kArrSize, 0.0);
	std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
	std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0), LZM(kArrSize, 0.0);
	
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		tmpU[idx(i, j, k)] = (u[idx(i, j, k)] + u[idx(i, j - 1, k)] + u[idx(i + 1, j - 1, k)] + u[idx(i + 1, j, k)]) * 0.25;
		tmpW[idx(i, j, k)] = (w[idx(i, j, k)] + u[idx(i, j - 1, k)] + w[idx(i, j - 1, k + 1)] + w[idx(i, j, k + 1)]) * 0.25;
	}
		
	std::vector<double> vecF_VX(kNx + 2 * kNumBCGrid, 0.0), vecF_VY(kNy + 2 * kNumBCGrid), vecF_VZ(kNz + 2 * kNumBCGrid);
	
	// V : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_VX[i] = v[idx(i, j, k)] * tmpU[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j, k)] = FXP[i];
			LXM[idx(i, j, k)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_VX.begin(), vecF_VX.end(), 0.0);
	}

	// V : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_VY[j] = v[idx(i, j, k)] * v[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j, k)] = FYP[j];
			LYM[idx(i, j, k)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_VY.begin(), vecF_VY.end(), 0.0);
	}
	
	// V : Z direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			vecF_VZ[k] = v[idx(i, j, k)] * tmpW[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VZ, FZP, FZM, kDz, kNz);

		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			LZP[idx(i, j, k)] = FYP[k];
			LZM[idx(i, j, k)] = FYM[k];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_VZ.begin(), vecF_VZ.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, alphaZ = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		alphaX = tmpU[idx(i, j, k)];
		alphaY = v[idx(i, j, k)];
		alphaZ = tmpW[idx(i, j, k)];
		
		cV[idx(i, j, k)]
			= 0.5 * (LXP[idx(i, j, k)] + LXM[idx(i, j, k)])
			+ 0.5 * (LYP[idx(i, j, k)] + LYM[idx(i, j, k)])
			+ 0.5 * (LZP[idx(i, j, k)] + LZM[idx(i, j, k)])
			- alphaX * (0.5 * (LXP[idx(i, j, k)] - LXM[idx(i, j, k)]))
			- alphaY * (0.5 * (LYP[idx(i, j, k)] - LYM[idx(i, j, k)]))
			- alphaZ * (0.5 * (LZP[idx(i, j, k)] - LZM[idx(i, j, k)]));

		assert(cV[idx(i, j, k)] == cV[idx(i, j, k)]);
		if (std::isnan(cV[idx(i, j, k)]) || std::isinf(cV[idx(i, j, k)])) {
			std::cout << "V-convection term nan/inf error : " << i << " " << j << " " 
			<< cV[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	return cV;
}

std::vector<double> MACSolver3D::AddConvectionFW(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {
	std::vector<double> cW(kArrSize, 0.0);
	std::vector<double> tmpU(kArrSize, 0.0), tmpV(kArrSize, 0.0);
	std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
	std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0), LZM(kArrSize, 0.0);


	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		tmpU[idx(i, j, k)] = (u[idx(i, j, k)] + u[idx(i, j, k - 1)] + u[idx(i + 1, j, k - 1)] + u[idx(i + 1, j, k)]) * 0.25;
		tmpV[idx(i, j, k)] = (w[idx(i, j, k)] + u[idx(i, j + 1, k)] + w[idx(i, j + 1, k - 1)] + w[idx(i, j + 1, k)]) * 0.25;
	}

	std::vector<double> vecF_VX(kNx + 2 * kNumBCGrid, 0.0), vecF_VY(kNy + 2 * kNumBCGrid), vecF_VZ(kNz + 2 * kNumBCGrid);

	// V : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_VX[i] = w[idx(i, j, k)] * tmpU[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j, k)] = FXP[i];
			LXM[idx(i, j, k)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_VX.begin(), vecF_VX.end(), 0.0);
	}

	// V : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_VY[j] = w[idx(i, j, k)] * tmpV[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j, k)] = FYP[j];
			LYM[idx(i, j, k)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_VY.begin(), vecF_VY.end(), 0.0);
	}


	// V : Z direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			vecF_VZ[k] = w[idx(i, j, k)] * w[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_VZ, FZP, FZM, kDz, kNz);

		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			LZP[idx(i, j, k)] = FYP[k];
			LZM[idx(i, j, k)] = FYM[k];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_VZ.begin(), vecF_VZ.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, alphaZ = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
		alphaX = tmpU[idx(i, j, k)];
		alphaY = tmpV[idx(i, j, k)];
		alphaZ = w[idx(i, j, k)];

		cW[idx(i, j, k)]
			= 0.5 * (LXP[idx(i, j, k)] + LXM[idx(i, j, k)])
			+ 0.5 * (LYP[idx(i, j, k)] + LYM[idx(i, j, k)])
			+ 0.5 * (LZP[idx(i, j, k)] + LZM[idx(i, j, k)])
			- alphaX * (0.5 * (LXP[idx(i, j, k)] - LXM[idx(i, j, k)]))
			- alphaY * (0.5 * (LYP[idx(i, j, k)] - LYM[idx(i, j, k)]))
			- alphaZ * (0.5 * (LZP[idx(i, j, k)] - LZM[idx(i, j, k)]));

		assert(cW[idx(i, j, k)] == cW[idx(i, j, k)]);
		if (std::isnan(cW[idx(i, j, k)]) || std::isinf(cW[idx(i, j, k)])) {
			std::cout << "V-convection term nan/inf error : " << i << " " << j << " "
				<< cW[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	return cW;
}

int MACSolver3D::UnitHJWENO5(
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

std::vector<double> MACSolver3D::AddViscosityFU(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dU(kArrSize, 0.0), mu(kArrSize, 0.0);
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
			dU[idx(i, j, k)] = 0.0;
		
		return dU;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		mu[idx(i, j, k)] = kMuL + (kMuH - kMuL) * H[idx(i, j, k)];

	// subcell
	// thetaH : portion of liquid
	double theta = 0.0, thetaH = 0.0;
	
	// need for updating rho
	/*
	xy plane (k fixed)					xz plane (j fixed)
	-------------------------			-------------------------
	|			|			|			|			|			|
	|			|			|			|			|			|
	|			|			|			|			|			|
	---------lsUNHalf--------			---------lsUTHalf--------
	|			|			|			|			|			|
	|	lsW   lsUM(*)  lsM  |			|	lsW   lsUM(*)  lsM  |
	|			|			|			|			|			|
	---------lsUSHalf--------			---------lsUBHalf--------
	|			|			|			|			|			|
	|			|			|			|			|			|
	|			|			|			|			|			|
	-------------------------			-------------------------
	*/
	double lsW = 0.0, lsM = 0.0, lsUSHalf = 0.0, lsUNHalf = 0.0, lsUBHalf = 0.0, lsUTHalf = 0.0;
	/*
	xy plane (k fixed)					xz plane (j fixed)
	-------------------------			-------------------------
	|			|			|			|			|			|
	|		  lsUN			|			|		  lsUT			|
	|			|			|			|			|			|
	----lsVW_N-------lsVN----			----lsWW_T-------lsWT----
	|			|			|			|			|			|
  lsUW		  lsUM  (i,j) lsUE			|		  lsUM	(i,k)	|
	|			|			|			|			|			|
	----lsVW---------lsVM----			-----lsWW--------lwWM----
	|			|			|			|			|			|
	|		  lsUS			|			|		  lsUB			|
	|			|			|			|			|			|
	-------------------------			-------------------------
	*/
	// u_xx, u_yy, u_zz
	double lsUW = 0.0, lsUE = 0.0, lsUS = 0.0, lsUN = 0.0, lsUM = 0.0, lsUB = 0.0, lsUT = 0.0;
	// v_xy
	double lsVW_N = 0.0, lsVW = 0.0, lsVM = 0.0, lsVN = 0.0;
	// w_xz
	double lsWW_T = 0.0, lsWW = 0.0, lsWM = 0.0, lsWT = 0.0;
	double uW = 0.0, uE = 0.0, uS = 0.0, uN = 0.0, uM = 0.0, uB = 0.0, uT = 0.0;
	double vW = 0.0, vW_N = 0.0, vM = 0.0, vN = 0.0;
	double wW = 0.0, wW_T = 0.0, wM = 0.0, wT = 0.0;
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0, muU_Z_B = 0.0, muU_Z_T = 0.0;
	double muV_X_S = 0.0, muV_X_N = 0.0, muW_X_B = 0.0, muW_X_T = 0.0;
	double visX = 0.0, visY = 0.0, visZ = 0.0;

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		muU_X_W = 0.0; muU_X_E = 0.0; muU_Y_S = 0.0; muU_Y_N = 0.0, muU_Z_B = 0.0, muU_Z_T = 0.0;
		muV_X_S = 0.0; muV_X_N = 0.0, muW_X_B = 0.0, muW_X_T = 0.0;
		visX = 0.0; visY = 0.0; visZ = 0.0;
		
		lsUM = 0.5 * (ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]);
		lsUW = 0.5 * (ls[idx(i - 2, j, k)] + ls[idx(i - 1, j, k)]);
		lsUE = 0.5 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]);
		lsUS = 0.5 * (ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsUN = 0.5 * (ls[idx(i - 1, j + 1, k)] + ls[idx(i, j + 1, k)]);
		lsUB = 0.5 * (ls[idx(i - 1, j, k - 1)] + ls[idx(i, j, k - 1)]);
		lsUT = 0.5 * (ls[idx(i - 1, j, k + 1)] + ls[idx(i, j, k + 1)]);
		lsVW = 0.5 * (ls[idx(i - 1, j - 1, k)] + ls[idx(i - 1, j, k)]);
		lsVW_N = 0.5 * (ls[idx(i - 1, j, k)] + ls[idx(i - 1, j + 1, k)]);
		lsVM = 0.5 * (ls[idx(i, j - 1, k)] + ls[idx(i, j, k)]);
		lsVN = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]);
		lsWW = 0.5 * (ls[idx(i - 1, j, k - 1)] + ls[idx(i - 1, j, k)]);
		lsWW_T = 0.5 * (ls[idx(i - 1, j, k)] + ls[idx(i - 1, j, k + 1)]);
		lsVM = 0.5 * (ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]);
		lsWT = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]);

		uM = u[idx(i, j, k)];
		uW = u[idx(i - 1, j, k)];
		uE = u[idx(i + 1, j, k)];
		uS = u[idx(i, j - 1, k)];
		uN = u[idx(i, j + 1, k)];
		uB = u[idx(i, j, k - 1)];
		uT = u[idx(i, j, k + 1)];

		vW = v[idx(i - 1, j, k)];
		vW_N = v[idx(i - 1, j + 1, k)];
		vM = v[idx(i, j, k)];
		vN = v[idx(i, j + 1, k)];
		wW = w[idx(i - 1, j, k)];
		wW_T = w[idx(i - 1, j, k + 1)];
		wM = w[idx(i, j, k)];
		wT = w[idx(i, j, k + 1)];
		
		lsW = ls[idx(i - 1, j, k)];
		lsM = ls[idx(i, j, k)];
		lsUSHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsUNHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j + 1, k)] + ls[idx(i, j + 1, k)]);
		lsUBHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j, k - 1)] + ls[idx(i, j, k - 1)]);
		lsUTHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j, k + 1)] + ls[idx(i, j, k + 1)]);
		
		double rhoEffWE = 0.0, rhoEffSN = 0.0, rhoEffBT = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0, iRhoEffBT = 0.0;
		if (lsW >= 0 && lsM >= 0) {
			thetaH = 1.0;
		}
		else if (lsW <= 0 && lsM <= 0) {
			thetaH = 0.0;
		}
		else  {
			thetaH = (std::max(lsW, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsW) + std::fabs(lsM));
		}
		iRhoEffWE = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsUSHalf >= 0 && lsUNHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsUSHalf <= 0 && lsUNHalf <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsUSHalf, 0.0) + std::max(lsUNHalf, 0.0))
				/ (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
		}
		iRhoEffSN = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsUBHalf >= 0 && lsUTHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsUBHalf <= 0 && lsUTHalf <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsUBHalf, 0.0) + std::max(lsUTHalf, 0.0))
				/ (std::fabs(lsUBHalf) + std::fabs(lsUTHalf));
		}
		iRhoEffBT = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		muU_X_W = mu[idx(i - 1, j, k)] * (uM - uW) / kDx;
		muU_X_E = mu[idx(i, j, k)] * (uE - uM) / kDx;
		muU_Y_S = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (uM - uS) / kDy;
		muU_Y_N = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j + 1, k)] + mu[idx(i, j + 1, k)])
			* (uN - uM) / kDy;
		muU_Z_B = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (uM - uB) / kDz;
		muU_Z_T = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k + 1)] + mu[idx(i, j, k + 1)])
			* (uT - uM) / kDz;
		muV_X_S = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (vM - vW) / kDx;
		muV_X_N = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j + 1, k)] + mu[idx(i, j + 1, k)])
			* (vN - vW_N) / kDx;
		muW_X_B = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (wM - wW) / kDx;
		muW_X_T = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k + 1)] + mu[idx(i, j, k + 1)])
			* (wT - wW_T) / kDx;

		visX = 2.0 * (muU_X_E - muU_X_W) / kDx;
		visY = (muU_Y_N - muU_Y_S) / kDy + (muV_X_N - muV_X_S) / kDy;
		visZ = (muU_Z_T - muU_Z_B) / kDz + (muW_X_T - muW_X_T) / kDz;
		dU[idx(i, j, k)] = visX * iRhoEffWE + visY * iRhoEffSN + visZ * iRhoEffBT;

		if (std::isnan(dU[idx(i, j, k)]) || std::isinf(dU[idx(i, j, k)])) {
			std::cout << "U-viscosity term nan/inf error : " << i << " " << j << " " << k << " " << dU[idx(i, j, k)] << std::endl;
			exit(1);
		}
		assert(dU[idx(i, j, k)] == dU[idx(i, j, k)]);
	}
	
	return dU;
}

std::vector<double> MACSolver3D::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV(kArrSize, 0.0), mu(kArrSize, 0.0);

	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
			dV[idx(i, j, k)] = 0.0;

		return dV;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		mu[idx(i, j, k)] = kMuL + (kMuH - kMuL) * H[idx(i, j, k)];

	// subcell
	// thetaH : portion of liquid
	double theta = 0.0, thetaH = 0.0;
	// need for updating rho
	
	/*
	* = lsVM(v[idx(i, j, k)])
	xy plane (k fixed)							yz plane (i fixed)
	-------------------------------------		-------------------------
	|			|			|			|		|			|			|
	|			|    lsM	|			|		|			|			|
	|			|			|			|		|			|			|
	---------lsVWHalf-*--lsVEHalf--------		--------lsVTHalf---------
	|			|			|			|		|			|			|
	|			|	 lsS	|			|		|	lsS		*	 lsM	|
	|			|			|			|		|			|			|
	-------------------------------------		--------lsVBHalf---------
												|			|			|
												|			|			|
												|			|			|
												-------------------------
	*/
	double lsVWHalf = 0.0, lsVEHalf = 0.0, lsM = 0.0, lsS = 0.0, lsVBHalf = 0.0, lsVTHalf = 0.0;
	/*
	xy plane (k fixed)							yz plane (i fixed)
	----------------lsVN-----------------		-------------------------
	|			|			|			|		|			|			|
	|		  lsUM  (i,j)  lsUE			|		|		  lsVT			|
	|			|			|			|		|			|			|
	-----lsVW-------lsVM--------lsVE-----		---lsWS_T-------lsWT-----
	|			|			|			|		|			|			|
	|		  lsUS		  lsUE_S		|	  lsVS	  	  lsVM	(j,k) lsVN
	|			|			|			|		|			|			|
	----------------lsVS-----------------		----lsWS--------lsWM-----
												|			|			|
												|		  lsVB			|
												|			|			|
												-------------------------
	*/
	// v_xx, v_yy, v_zz
	double lsVW = 0.0, lsVE = 0.0, lsVS = 0.0, lsVN = 0.0, lsVM = 0.0, lsVB = 0.0, lsVT = 0.0;
	// u_xy
	double lsUM = 0.0, lsUS = 0.0, lsUE = 0.0, lsUE_S = 0.0;
	// w_yz
	double lsWM = 0.0, lsWS = 0.0, lsWT = 0.0, lsWS_T = 0.0;
	
	// jump condition
	double uS = 0.0, uM = 0.0, uE = 0.0, uE_S = 0.0;
	double vM = 0.0, vW = 0.0, vE = 0.0, vS = 0.0, vN = 0.0, vB = 0.0, vT = 0.0;
	double wS = 0.0, wM = 0.0, wT = 0.0, wS_T = 0.0;
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0, muV_Z_B = 0.0, muV_Z_T = 0.0;
	double muU_Y_W = 0.0, muU_Y_E = 0.0, muW_Y_B = 0.0, muW_Y_T = 0.0;
	double visX = 0.0, visY = 0.0, visZ = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		muV_X_W = 0.0; muV_X_E = 0.0; muV_Y_S = 0.0; muV_Y_N = 0.0; muV_Z_B = 0.0; muV_Z_T = 0.0;
		muU_Y_W = 0.0; muU_Y_E = 0.0; muW_Y_B = 0.0; muW_Y_T = 0.0;
		visX = 0.0; visY = 0.0; visZ = 0.0;
		
		lsVM = 0.5 * (ls[idx(i, j - 1, k)] + ls[idx(i, j, k)]);
		lsVW = 0.5 * (ls[idx(i - 1, j - 1, k)] + ls[idx(i - 1, j, k)]);
		lsVE = 0.5 * (ls[idx(i + 1, j - 1, k)] + ls[idx(i + 1, j, k)]);
		lsVS = 0.5 * (ls[idx(i, j - 2, k)] + ls[idx(i, j - 1, k)]);
		lsVN = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]);
		lsVB = 0.5 * (ls[idx(i, j - 1, k - 1)] + ls[idx(i, j, k - 1)]);
		lsVT = 0.5 * (ls[idx(i, j - 1, k + 1)] + ls[idx(i, j, k + 1)]);
		lsUM = 0.5 * (ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]);
		lsUS = 0.5 * (ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsUE = 0.5 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]);
		lsUE_S = 0.5 * (ls[idx(i, j - 1, k)] + ls[idx(i + 1, j - 1, k)]);
		lsWM = 0.5 * (ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]);
		lsWT = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]);
		lsWS = 0.5 * (ls[idx(i, j - 1, k - 1)] + ls[idx(i, j - 1, k)]);
		lsWS_T = 0.5 * (ls[idx(i, j - 1, k)] + ls[idx(i, j - 1, k + 1)]);

		vM = v[idx(i, j, k)];
		vW = v[idx(i - 1, j, k)];
		vE = v[idx(i + 1, j, k)];
		vS = v[idx(i, j - 1, k)];
		vN = v[idx(i, j + 1, k)];
		vB = v[idx(i, j, k - 1)];
		vT = v[idx(i, j, k + 1)];

		uM = u[idx(i, j, k)];
		uE = u[idx(i + 1, j, k)];
		uS = u[idx(i, j - 1, k)];
		uE_S = u[idx(i + 1, j - 1, k)];
		wM = w[idx(i, j, k)];
		wT = w[idx(i, j, k + 1)];
		wS = w[idx(i - 1, j, k)];
		wS_T = w[idx(i - 1, j, k + 1)];
		
		lsS = ls[idx(i, j - 1, k)];
		lsM = ls[idx(i, j, k)];
		lsVWHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsVEHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)] + ls[idx(i + 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsVBHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j, k - 1)] + ls[idx(i, j - 1, k - 1)] + ls[idx(i, j - 1, k)]);
		lsVTHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j, k + 1)] + ls[idx(i, j - 1, k + 1)] + ls[idx(i, j - 1, k)]);
		
		muV_X_W = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (vM - vW) / kDx;
		muV_X_E = 0.25 * (mu[idx(i, j, k)] + mu[idx(i + 1, j, k)] + mu[idx(i + 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (vE - vM) / kDx;
		muV_Y_S = mu[idx(i, j - 1, k)] * (vM - vS) / kDy;
		muV_Y_N = mu[idx(i, j, k)] * (vN - vM) / kDy;
		muV_Z_B = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j, k - 1)] + mu[idx(i, j - 1, k - 1)] + mu[idx(i, j - 1, k)])
			* (vM - vB) / kDz;
		muV_Z_T = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j, k + 1)] + mu[idx(i, j - 1, k + 1)] + mu[idx(i, j - 1, k)])
			* (vT - vM) / kDz;
		muU_Y_W = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (uM - uS) / kDy;
		muU_Y_E = 0.25 * (mu[idx(i, j, k)] + mu[idx(i + 1, j, k)] + mu[idx(i + 1, j - 1, k)] + mu[idx(i, j - 1, k)])
			* (uE - uE_S) / kDy;
		muW_Y_B = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j, k - 1)] + mu[idx(i, j - 1, k - 1)] + mu[idx(i, j - 1, k)])
			* (wM - wS) / kDy;
		muW_Y_T = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j, k + 1)] + mu[idx(i, j - 1, k + 1)] + mu[idx(i, j - 1, k)])
			* (wT - wS_T) / kDy;

		double rhoEffWE = 0.0, rhoEffSN = 0.0, rhoEffBT = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0, iRhoEffBT = 0.0;
		if (lsVWHalf >= 0 && lsVEHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsVWHalf <= 0 && lsVEHalf <= 0) {
			thetaH = 0.0;
		}
		else  {
			thetaH = (std::max(lsVWHalf, 0.0) + std::max(lsVEHalf, 0.0))
				/ (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
		}
		iRhoEffWE = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsS >= 0 && lsM >= 0) {
			thetaH = 1.0;
		}
		else if (lsS <= 0 && lsM <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsS, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsS) + std::fabs(lsM));
		}
		iRhoEffSN = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsVBHalf >= 0 && lsVTHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsVBHalf <= 0 && lsVTHalf <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsVBHalf, 0.0) + std::max(lsVTHalf, 0.0))
				/ (std::fabs(lsVBHalf) + std::fabs(lsVTHalf));
		}
		iRhoEffBT = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		visX = (muU_Y_E - muU_Y_W) / kDx + (muV_X_E - muV_X_W) / kDx;
		visY = 2.0 * (muV_Y_N - muV_Y_S) / kDy;
		visZ = (muW_Y_T - muW_Y_B) / kDz + (muV_Z_T - muV_Z_B) / kDz;
		
		dV[idx(i, j, k)] = visX * iRhoEffWE + visY * iRhoEffSN + visZ * iRhoEffBT;
		
		assert(dV[idx(i, j, k)] == dV[idx(i, j, k)]);
		if (std::isnan(dV[idx(i, j, k)]) || std::isinf(dV[idx(i, j, k)])) {
			std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << k << " " << dV[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}
	
	return dV;	
}

std::vector<double> MACSolver3D::AddViscosityFW(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dW(kArrSize, 0.0), mu(kArrSize, 0.0);

	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
			dW[idx(i, j, k)] = 0.0;

		return dW;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		mu[idx(i, j, k)] = kMuL + (kMuH - kMuL) * H[idx(i, j, k)];

	// subcell
	// thetaH : portion of liquid
	double theta = 0.0, thetaH = 0.0;
	// need for updating rho

	/*
	* = lsWM(w[idx(i, j, k)])
	xz plane (k fixed)							yz plane (i fixed)
	-------------------------------------		-------------------------------------
	|			|			|			|		|			|			|			|
	|			|    lsM	|			|		|			|	 lsM	|			|
	|			|			|			|		|			|			|			|
	---------lsWWHalf-*--lsWEHalf--------		---------lsWSHalf--*-lsWNHalf--------
	|			|			|			|		|			|			|			|
	|			|	 lsB	|			|		|	  		|	 lsB	|			|
	|			|			|			|		|			|			|			|
	-------------------------------------		------------------------------------
												
	*/
	double lsM = 0.0, lsB = 0.0, lsWWHalf = 0.0, lsWEHalf = 0.0, lsWSHalf = 0.0, lsWNHalf = 0.0;
	/*
	xz plane (j fixed)							yz plane (i fixed)
	----------------lsWT-----------------		-----------------lsWT-----------------
	|			|			|			|		|			|			|			|
	|		  lsUM  (i,k)  lsUE			|		|		  lsVM	(j,k) lsVN			|
	|			|			|			|		|			|			|			|
	-----lsWW-------lsWM--------lsWE-----		-----lsWS--------lsWM-------lsWN-----
	|			|			|			|		|			|			|			|
	|		  lsUB		  lsUE_B		|		|	  	  lsVB		  lsVN_B		|
	|			|			|			|		|			|			|			|
	----------------lsWB-----------------		-----------------lsWB----------------
	*/
	// u_xy
	double lsUM = 0.0, lsUB = 0.0, lsUE = 0.0, lsUE_B = 0.0;
	// v_yz
	double lsVM = 0.0, lsVB = 0.0, lsVN = 0.0, lsVN_B = 0.0;
	// w_xx, w_yy, w_zz
	double lsWW = 0.0, lsWE = 0.0, lsWS = 0.0, lsWN = 0.0, lsWM = 0.0, lsWB = 0.0, lsWT = 0.0;

	// jump condition
	double uB = 0.0, uM = 0.0, uE = 0.0, uE_B = 0.0;
	double vB = 0.0, vM = 0.0, vN = 0.0, vN_B = 0.0;
	double wM = 0.0, wW = 0.0, wE = 0.0, wS = 0.0, wN = 0.0, wB = 0.0, wT = 0.0;
	
	double muW_X_W = 0.0, muW_X_E = 0.0, muW_Y_S = 0.0, muW_Y_N = 0.0, muW_Z_B = 0.0, muW_Z_T = 0.0;
	double muU_Z_W = 0.0, muU_Z_E = 0.0, muV_Z_S = 0.0, muV_Z_N = 0.0;
	double visX = 0.0, visY = 0.0, visZ = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
		muW_X_W = 0.0; muW_X_E = 0.0; muW_Y_S = 0.0; muW_Y_N = 0.0; muW_Z_B = 0.0; muW_Z_T = 0.0;
		muU_Z_W = 0.0; muU_Z_E = 0.0; muV_Z_S = 0.0; muV_Z_N = 0.0;
		visX = 0.0; visY = 0.0; visZ = 0.0;

		lsWM = 0.5 * (ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]);
		lsWW = 0.5 * (ls[idx(i - 1, j, k - 1)] + ls[idx(i - 1, j, k)]);
		lsWE = 0.5 * (ls[idx(i + 1, j, k - 1)] + ls[idx(i + 1, j, k)]);
		lsWS = 0.5 * (ls[idx(i, j - 1, k - 1)] + ls[idx(i, j - 1, k)]);
		lsWN = 0.5 * (ls[idx(i, j + 1, k - 1)] + ls[idx(i, j + 1, k)]);
		lsWB = 0.5 * (ls[idx(i, j, k - 2)] + ls[idx(i, j, k - 1)]);
		lsWT = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]);
		lsUM = 0.5 * (ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]);
		lsUB = 0.5 * (ls[idx(i - 1, j, k - 1)] + ls[idx(i, j, k - 1)]);
		lsUE = 0.5 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]);
		lsUE_B = 0.5 * (ls[idx(i, j - 1, k)] + ls[idx(i + 1, j - 1, k)]);
		lsVM = 0.5 * (ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]);
		lsVB = 0.5 * (ls[idx(i, j - 1, k - 1)] + ls[idx(i, j - 1, k)]);
		lsVN = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]);
		lsVN_B = 0.5 * (ls[idx(i, j + 1, k - 1)] + ls[idx(i, j + 1, k)]);

		wM = w[idx(i, j, k)];
		wW = w[idx(i - 1, j, k)];
		wE = w[idx(i + 1, j, k)];
		wS = w[idx(i, j - 1, k)];
		wN = w[idx(i, j + 1, k)];
		wB = w[idx(i, j, k - 1)];
		wT = w[idx(i, j, k + 1)];

		uM = u[idx(i, j, k)];
		uE = u[idx(i + 1, j, k)];
		uB = u[idx(i, j - 1, k)];
		uE_B = u[idx(i + 1, j - 1, k)];
		vM = v[idx(i, j, k)];
		vB = v[idx(i, j, k - 1)];
		vN = v[idx(i, j + 1, k)];
		vN_B = v[idx(i, j + 1, k - 1)];

		lsB = ls[idx(i, j, k - 1)];
		lsM = ls[idx(i, j, k)];
		lsWWHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)] + ls[idx(i - 1, j, k - 1)] + ls[idx(i, j, k - 1)]);
		lsWEHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)] + ls[idx(i + 1, j, k - 1)] + ls[idx(i, j, k - 1)]);
		lsWSHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j - 1, k)] + ls[idx(i, j - 1, k - 1)] + ls[idx(i, j, k - 1)]);
		lsWNHalf = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)] + ls[idx(i, j + 1, k - 1)] + ls[idx(i, j, k - 1)]);

		muW_X_W = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (wM - wW) / kDx;
		muW_X_E = 0.25 * (mu[idx(i, j, k)] + mu[idx(i + 1, j, k)] + mu[idx(i + 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (wE - wM) / kDx;
		muW_Y_S = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j - 1, k)] + mu[idx(i, j - 1, k - 1)] + mu[idx(i, j, k - 1)])
			* (wM - wS) / kDy;
		muW_Y_N = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j + 1, k)] + mu[idx(i, j + 1, k - 1)] + mu[idx(i, j, k - 1)])
			* (wN - wM) / kDy;
		muW_Z_B = mu[idx(i, j, k - 1)] * (wM - wB) / kDz;
		muW_Z_T = mu[idx(i, j, k)] * (wT - wM) / kDz;
		muU_Z_W = 0.25 * (mu[idx(i, j, k)] + mu[idx(i - 1, j, k)] + mu[idx(i - 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (uM - uB) / kDz;
		muU_Z_E = 0.25 * (mu[idx(i, j, k)] + mu[idx(i + 1, j, k)] + mu[idx(i + 1, j, k - 1)] + mu[idx(i, j, k - 1)])
			* (uE - uE_B) / kDz;
		muV_Z_S = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j - 1, k)] + mu[idx(i, j - 1, k - 1)] + mu[idx(i, j - 1, k - 1)])
			* (vM - vB) / kDz;
		muV_Z_N = 0.25 * (mu[idx(i, j, k)] + mu[idx(i, j + 1, k)] + mu[idx(i, j + 1, k - 1)] + mu[idx(i, j + 1, k - 1)])
			* (vN - vN_B) / kDz;

		double rhoEffWE = 0.0, rhoEffSN = 0.0, rhoEffBT = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0, iRhoEffBT = 0.0;
		if (lsWWHalf >= 0 && lsWEHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsWWHalf <= 0 && lsWEHalf <= 0) {
			thetaH = 0.0;
		}
		else  {
			thetaH = (std::max(lsWWHalf, 0.0) + std::max(lsWEHalf, 0.0))
				/ (std::fabs(lsWWHalf) + std::fabs(lsWEHalf));
		}
		iRhoEffWE = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsWSHalf >= 0 && lsWNHalf >= 0) {
			thetaH = 1.0;
		}
		else if (lsWSHalf <= 0 && lsWNHalf <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsWSHalf, 0.0) + std::max(lsWNHalf, 0.0))
				/ (std::fabs(lsWSHalf) + std::fabs(lsWNHalf));
		}
		iRhoEffSN = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		if (lsB >= 0 && lsM >= 0) {
			thetaH = 1.0;
		}
		else if (lsB <= 0 && lsM <= 0) {
			thetaH = 0.0;
		}
		else {
			thetaH = (std::max(lsB, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsB) + std::fabs(lsM));
		}
		iRhoEffBT = 1.0 / kRhoH * thetaH + 1.0 / kRhoL * (1.0 - thetaH);

		visX = (muU_Z_E - muU_Z_W) / kDx + (muW_X_E - muW_X_W) / kDx;
		visY = (muV_Z_N - muV_Z_S) / kDy + (muW_Y_N - muW_Y_S) / kDy;
		visZ = 2.0 * (muW_Z_T - muW_Z_B) / kDz;

		dW[idx(i, j, k)] = visX * iRhoEffWE + visY * iRhoEffSN + visZ * iRhoEffBT;

		assert(dW[idx(i, j, k)] == dW[idx(i, j, k)]);
		if (std::isnan(dW[idx(i, j, k)]) || std::isinf(dW[idx(i, j, k)])) {
			std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << k << " " << dW[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	return dW;
}

std::vector<double> MACSolver3D::AddGravityFU() {
	std::vector<double> gU(kArrSize, 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM::X) {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			gU[idx(i, j, k)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			gU[idx(i, j, k)] = -1.0 / kFr;
		}
	}

	return gU;
}

std::vector<double> MACSolver3D::AddGravityFV() {
	std::vector<double> gV(kArrSize, 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM::Y) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			gV[idx(i, j, k)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			gV[idx(i, j, k)] = -1.0 / kFr;
		}
	}

	return gV;
}

std::vector<double> MACSolver3D::AddGravityFW() {
	std::vector<double> gW(kArrSize, 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM::Z) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			gW[idx(i, j, k)] = 0.0;
		}
	}
	else {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			gW[idx(i, j, k)] = -1.0 / kFr;
		}
	}

	return gW;
}

int MACSolver3D::GetIntermediateVel(const std::shared_ptr<LevelSetSolver3D>& LSolver, const std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	std::vector<double>& uhat, std::vector<double>& vhat, std::vector<double>& what,
	const std::vector<double>& H) {
		
	// Update rhs
	if (kTimeOrder == TIMEORDERENUM::EULER) {
		std::vector<double> FU1(kArrSize, 0.0);
		std::vector<double> FV1(kArrSize, 0.0);
		std::vector<double> FW1(kArrSize, 0.0);

		FU1 = UpdateFU(LSolver, ls, u, v, w, H);
		FV1 = UpdateFV(LSolver, ls, u, v, w, H);
		FW1 = UpdateFW(LSolver, ls, u, v, w, H);
		
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			uhat[idx(i, j, k)] = u[idx(i, j, k)] + m_dt * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			vhat[idx(i, j, k)] = v[idx(i, j, k)] + m_dt * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			what[idx(i, j, k)] = w[idx(i, j, k)] + m_dt * FW1[idx(i, j, k)];
		}
	}
	else if (kTimeOrder == TIMEORDERENUM::RK2) {
		std::vector<double> FU1(kArrSize, 0.0), FV1(kArrSize, 0.0), FW1(kArrSize, 0.0),
			FU2(kArrSize, 0.0), FV2(kArrSize, 0.0), FW2(kArrSize, 0.0);

		// FU1 & FV1 & FW1: L(u^(0))
		// FU2 & FV2 & FW2: u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, w, H);
		FV1 = UpdateFV(LSolver, ls, u, v, w, H);
		FW1 = UpdateFW(LSolver, ls, u, v, w, H);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU2[idx(i, j, k)] = u[idx(i, j, k)] + m_dt * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV2[idx(i, j, k)] = v[idx(i, j, k)] + m_dt * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW2[idx(i, j, k)] = w[idx(i, j, k)] + m_dt * FW1[idx(i, j, k)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		std::fill(FW1.begin(), FW1.end(), 0.0);
		
		// FU1 & FV1 : L(u^(1))
		// FU2 & FV2 : u^(2)
		ApplyBC_U_3D(FU2);
		ApplyBC_V_3D(FV2);
		ApplyBC_V_3D(FW2);
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, FW2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, FW2, H);
		FW1 = UpdateFW(LSolver, ls, FU2, FV2, FW2, H);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU2[idx(i, j, k)] = FU2[idx(i, j, k)] + m_dt * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV2[idx(i, j, k)] = FV2[idx(i, j, k)] + m_dt * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW2[idx(i, j, k)] = FW2[idx(i, j, k)] + m_dt * FW1[idx(i, j, k)];
		}
		
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			uhat[idx(i, j, k)] = 0.5 * u[idx(i, j, k)] + 0.5 * FU2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			vhat[idx(i, j, k)] = 0.5 * v[idx(i, j, k)] + 0.5 * FV2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			what[idx(i, j, k)] = 0.5 * w[idx(i, j, k)] + 0.5 * FW2[idx(i, j, k)];
		}
	}
	else if (kTimeOrder == TIMEORDERENUM::RK3) {
		std::vector<double> FU1(kArrSize, 0.0), FV1(kArrSize, 0.0), FW1(kArrSize, 0.0),
			FU2(kArrSize, 0.0), FV2(kArrSize, 0.0), FW2(kArrSize, 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, w, H);
		FV1 = UpdateFV(LSolver, ls, u, v, w, H);
		FW1 = UpdateFW(LSolver, ls, u, v, w, H);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU2[idx(i, j, k)] = u[idx(i, j, k)] + m_dt * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV2[idx(i, j, k)] = v[idx(i, j, k)] + m_dt * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW2[idx(i, j, k)] = w[idx(i, j, k)] + m_dt * FW1[idx(i, j, k)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		std::fill(FW1.begin(), FW1.end(), 0.0);
		
		// FU1 & FV1 & FW1: L(u^(1))
		// FU2 & FV2 & FW2: u^(1) + \delta t L (u^(1))
		ApplyBC_U_3D(FU2);
		ApplyBC_V_3D(FV2);
		ApplyBC_W_3D(FW2);
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, FW2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, FW2, H);
		FW1 = UpdateFV(LSolver, ls, FU2, FV2, FW2, H);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU2[idx(i, j, k)] = FU2[idx(i, j, k)] + m_dt * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) 
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV2[idx(i, j, k)] = FV2[idx(i, j, k)] + m_dt * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW2[idx(i, j, k)] = FW2[idx(i, j, k)] + m_dt * FW1[idx(i, j, k)];
		}
		std::fill(FU1.begin(), FU1.end(), 0.0);
		std::fill(FV1.begin(), FV1.end(), 0.0);
		std::fill(FW1.begin(), FW1.end(), 0.0);
		
		// FU1 & FV1 & FW1 : u^(2)
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU1[idx(i, j, k)] = 0.75 * u[idx(i, j, k)] + 0.25 * FU2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV1[idx(i, j, k)] = 0.75 * v[idx(i, j, k)] + 0.25 * FV2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW1[idx(i, j, k)] = 0.75 * w[idx(i, j, k)] + 0.25 * FW2[idx(i, j, k)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);
		std::fill(FW2.begin(), FW2.end(), 0.0);

		// FU2 & FV2 & FW2: L(u^(2))
		// FU1 & FV1 & FW1: u^(2) + \delta t L (u^(2))
		ApplyBC_U_3D(FU1);
		ApplyBC_V_3D(FV1);
		ApplyBC_W_3D(FW1);
		FU2 = UpdateFU(LSolver, ls, FU1, FV1, FW1, H);
		FV2 = UpdateFV(LSolver, ls, FU1, FV1, FW1, H);
		FW2 = UpdateFW(LSolver, ls, FU1, FV1, FW1, H);
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FU1[idx(i, j, k)] = FU1[idx(i, j, k)] + m_dt * FU2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			FV1[idx(i, j, k)] = FV1[idx(i, j, k)] + m_dt * FV2[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
			FW1[idx(i, j, k)] = FW1[idx(i, j, k)] + m_dt * FW2[idx(i, j, k)];
		}
		std::fill(FU2.begin(), FU2.end(), 0.0);
		std::fill(FV2.begin(), FV2.end(), 0.0);
		std::fill(FW2.begin(), FW2.end(), 0.0);

		// FU1 & FV1 : u^(2) + \delta t L (u^(2))
		// FU2 & FV2 : doesn't need. set to zero.
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			uhat[idx(i, j, k)] = 1.0 / 3.0 * u[idx(i, j, k)] + 2.0 / 3.0 * FU1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			vhat[idx(i, j, k)] = 1.0 / 3.0 * v[idx(i, j, k)] + 2.0 / 3.0 * FV1[idx(i, j, k)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			what[idx(i, j, k)] = 1.0 / 3.0 * w[idx(i, j, k)] + 2.0 / 3.0 * FW1[idx(i, j, k)];
		}
	}

	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		if (std::isnan(uhat[idx(i, j, k)]) || std::isinf(uhat[idx(i, j, k)])) {
			std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << k << " " << uhat[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		if (std::isnan(vhat[idx(i, j, k)]) || std::isinf(vhat[idx(i, j, k)])) {
			std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << k << " " << vhat[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
		if (std::isnan(what[idx(i, j, k)]) || std::isinf(what[idx(i, j, k)])) {
			std::cout << "What term nan/inf error : " << i << " " << j << " " << k << " " << what[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}
	
	return 0;
}

int MACSolver3D::SetPoissonSolver(POISSONTYPE type) {
	m_PoissonSolverType = type;
	if (!m_Poisson)
		m_Poisson = std::make_shared<PoissonSolver3D>(kNx, kNy, kNz, kNumBCGrid);

	return 0;
}

int MACSolver3D::SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
	const std::vector<double>& ls, const std::vector<double>& lsBackup,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& H, const int maxIter) {
	if (!m_Poisson) {
		perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
	}
	std::vector<double> rhs(kArrSize, 0.0);
	std::vector<double> iRhoW(kArrSize, 0.0), iRhoE(kArrSize, 0.0),
		iRhoS(kArrSize, 0.0), iRhoN(kArrSize, 0.0), iRhoB(kArrSize, 0.0), iRhoT(kArrSize, 0.0);
	
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
	double lsM = 0.0, lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsB = 0.0, lsT = 0.0;
	double FW = 0.0, FE = 0.0, FS = 0.0, FN = 0.0, FB = 0.0, FT = 0.0;
	// normal and tangent vector variable (n, t1, t2)
	double nXW = 0.0, nYW = 0.0, nZW = 0.0, nXE = 0.0, nYE = 0.0, nZE = 0.0,
		nXS = 0.0, nYS = 0.0, nZS = 0.0, nXN = 0.0, nYN = 0.0, nZN = 0.0,
		nXB = 0.0, nYB = 0.0, nZB = 0.0, nXT = 0.0, nYT = 0.0, nZT = 0.0,
		nXM = 0.0, nYM = 0.0, nZM = 0.0;
	// jump at grid node (if interface is at grid node, jump occurs and aW, aE, aS, aN, aM describe that condition)
	// jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
	double thetaH = 0.0, theta = 0.0, uEff = 0.0, vEff = 0.0, wEff = 0.0,
		kappaEff = 0.0, rhoEff = 0.0, nXEff = 0.0, nYEff = 0.0, nZEff = 0.0;
	
	if (kWe != 0.0 && !isnan(kWe) && !isinf(kWe)) {
		UpdateKappa(ls);
		ApplyBC_P_3D(m_kappa);
	}
	// A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
	std::vector<double> AVals, DiagVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy * kNz, Ancols = kNx * kNy * kNz;
	MKL_INT size = kNx * kNy * kNz;
	
	std::vector<double> dldX(kArrSize, 0.0), dldY(kArrSize, 0.0), dldZ(kArrSize, 0.0);
	
	// stored coef for A matrix, Dictionary but it is ordered
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
	for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
		ps[idx(i, j, k)] = 0.0;
		rhs[idx(i, j, k)] = 0.0;
	}
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		lsM = ls[idx(i, j, k)];
		lsW = ls[idx(i - 1, j, k)];
		lsE = ls[idx(i + 1, j, k)];
		lsS = ls[idx(i, j - 1, k)];
		lsN = ls[idx(i, j + 1, k)];
		lsB = ls[idx(i, j, k - 1)];
		lsT = ls[idx(i, j, k + 1)];

		dldX[idx(i, j, k)] = (lsE - lsW) / (2.0 * kDx);
		dldY[idx(i, j, k)] = (lsN - lsS) / (2.0 * kDy);
		dldZ[idx(i, j, k)] = (lsT - lsB) / (2.0 * kDz);
		if (std::fabs(dldX[idx(i, j, k)]) < eps) {
			dldX[idx(i, j, k)] = (lsE - lsM) / kDx;
			if (std::fabs(dldX[idx(i, j, k)]) < eps)
				dldX[idx(i, j, k)] = (lsM - lsW) / kDx;
		}
		if (std::fabs(dldY[idx(i, j, k)]) < eps) {
			dldY[idx(i, j, k)] = (lsN - lsM) / kDy;
			if (std::fabs(dldY[idx(i, j, k)]) < eps)
				dldY[idx(i, j, k)] = (lsM - lsS) / kDy;
		}
		if (std::fabs(dldZ[idx(i, j, k)]) < eps) {
			dldZ[idx(i, j, k)] = (lsT - lsM) / kDz;
			if (std::fabs(dldZ[idx(i, j, k)]) < eps)
				dldZ[idx(i, j, k)] = (lsM - lsB) / kDz;
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		lsM = ls[idx(i, j, k)];
		lsW = ls[idx(i - 1, j, k)];
		lsE = ls[idx(i + 1, j, k)];
		lsS = ls[idx(i, j - 1, k)];
		lsN = ls[idx(i, j + 1, k)];
		lsB = ls[idx(i, j, k - 1)];
		lsT = ls[idx(i, j, k + 1)];

		FW = 0.0;
		FE = 0.0;
		FS = 0.0;
		FN = 0.0;
		FB = 0.0;
		FT = 0.0;

		// normal vector = (\nabla \phi) / |\nabla \phi|
		nXW = dldX[idx(i - 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i - 1, j, k)], 2.0) + std::pow(dldY[idx(i - 1, j, k)], 2.0) + std::pow(dldZ[idx(i - 1, j, k)], 2.0)) + eps);
		nYW = dldY[idx(i - 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i - 1, j, k)], 2.0) + std::pow(dldY[idx(i - 1, j, k)], 2.0) + std::pow(dldZ[idx(i - 1, j, k)], 2.0)) + eps);
		nZW = dldZ[idx(i - 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i - 1, j, k)], 2.0) + std::pow(dldY[idx(i - 1, j, k)], 2.0) + std::pow(dldZ[idx(i - 1, j, k)], 2.0)) + eps);
		nXE = dldX[idx(i + 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i + 1, j, k)], 2.0) + std::pow(dldY[idx(i + 1, j, k)], 2.0) + std::pow(dldZ[idx(i + 1, j, k)], 2.0)) + eps);
		nYE = dldY[idx(i + 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i + 1, j, k)], 2.0) + std::pow(dldY[idx(i + 1, j, k)], 2.0) + std::pow(dldZ[idx(i + 1, j, k)], 2.0)) + eps);
		nZE = dldZ[idx(i + 1, j, k)] / (std::sqrt(std::pow(dldX[idx(i + 1, j, k)], 2.0) + std::pow(dldY[idx(i + 1, j, k)], 2.0) + std::pow(dldZ[idx(i + 1, j, k)], 2.0)) + eps);
		nXS = dldX[idx(i, j - 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j - 1, k)], 2.0) + std::pow(dldY[idx(i, j - 1, k)], 2.0) + std::pow(dldZ[idx(i, j - 1, k)], 2.0)) + eps);
		nYS = dldY[idx(i, j - 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j - 1, k)], 2.0) + std::pow(dldY[idx(i, j - 1, k)], 2.0) + std::pow(dldZ[idx(i, j - 1, k)], 2.0)) + eps);
		nZS = dldZ[idx(i, j - 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j - 1, k)], 2.0) + std::pow(dldY[idx(i, j - 1, k)], 2.0) + std::pow(dldZ[idx(i, j - 1, k)], 2.0)) + eps);
		nXN = dldX[idx(i, j + 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j + 1, k)], 2.0) + std::pow(dldY[idx(i, j + 1, k)], 2.0) + std::pow(dldZ[idx(i, j + 1, k)], 2.0)) + eps);
		nYN = dldY[idx(i, j + 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j + 1, k)], 2.0) + std::pow(dldY[idx(i, j + 1, k)], 2.0) + std::pow(dldZ[idx(i, j + 1, k)], 2.0)) + eps);
		nZN = dldZ[idx(i, j + 1, k)] / (std::sqrt(std::pow(dldX[idx(i, j + 1, k)], 2.0) + std::pow(dldY[idx(i, j + 1, k)], 2.0) + std::pow(dldZ[idx(i, j + 1, k)], 2.0)) + eps);
		nXB = dldX[idx(i, j, k - 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k - 1)], 2.0) + std::pow(dldY[idx(i, j, k - 1)], 2.0) + std::pow(dldZ[idx(i, j, k - 1)], 2.0)) + eps);
		nYB = dldY[idx(i, j, k - 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k - 1)], 2.0) + std::pow(dldY[idx(i, j, k - 1)], 2.0) + std::pow(dldZ[idx(i, j, k - 1)], 2.0)) + eps);
		nZB = dldZ[idx(i, j, k - 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k - 1)], 2.0) + std::pow(dldY[idx(i, j, k - 1)], 2.0) + std::pow(dldZ[idx(i, j, k - 1)], 2.0)) + eps);
		nXT = dldX[idx(i, j, k + 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k + 1)], 2.0) + std::pow(dldY[idx(i, j, k + 1)], 2.0) + std::pow(dldZ[idx(i, j, k + 1)], 2.0)) + eps);
		nYT = dldY[idx(i, j, k + 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k + 1)], 2.0) + std::pow(dldY[idx(i, j, k + 1)], 2.0) + std::pow(dldZ[idx(i, j, k + 1)], 2.0)) + eps);
		nZT = dldZ[idx(i, j, k + 1)] / (std::sqrt(std::pow(dldX[idx(i, j, k + 1)], 2.0) + std::pow(dldY[idx(i, j, k + 1)], 2.0) + std::pow(dldZ[idx(i, j, k + 1)], 2.0)) + eps);
		nXM = dldX[idx(i, j, k)] / (std::sqrt(std::pow(dldX[idx(i, j, k)], 2.0) + std::pow(dldY[idx(i, j, k)], 2.0) + std::pow(dldZ[idx(i, j, k)], 2.0)) + eps);
		nYM = dldY[idx(i, j, k)] / (std::sqrt(std::pow(dldX[idx(i, j, k)], 2.0) + std::pow(dldY[idx(i, j, k)], 2.0) + std::pow(dldZ[idx(i, j, k)], 2.0)) + eps);
		nZM = dldZ[idx(i, j, k)] / (std::sqrt(std::pow(dldX[idx(i, j, k)], 2.0) + std::pow(dldY[idx(i, j, k)], 2.0) + std::pow(dldZ[idx(i, j, k)], 2.0)) + eps);
		
		// [a]_\Gamma = a^+ - a^- = a^inside - a^outside
		// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// (P_M - P_W)/kDx appears
		// if P_M == P_+, P_W == P_-, a^+ means all terms related to P_+, P_W changed to P_M related terms
		// P_W = P_M -(2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa)

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsW >= 0 && lsM >= 0) {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsW <= 0 && lsM <= 0) {
			thetaH = 0.0; theta = 0.0;
		}
		else  {
			thetaH = (std::max(lsW, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		// coefficient
		iRhoW[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i - 1, j, k)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0;
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FW = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i - 1, j, k)]);
		FW *= iRhoW[idx(i, j, k)] / (kDx * kDx);
	
		// thetaH = portion of kRhoH
		// theta = portion of fluid cell adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsE >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsM < 0.0 && lsE < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsE, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsE) + std::fabs(lsM));
			theta = std::fabs(lsE) / (std::fabs(lsE) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoE[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i + 1, j, k)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FE = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j, k)] - H[idx(i, j, k)]);
		FE *= iRhoE[idx(i, j, k)] / (kDx * kDx);
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsS >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsS < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsS, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}
		
		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoS[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j - 1, k)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FS = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i, j - 1, k)]);
		FS *= iRhoS[idx(i, j, k)] / (kDy * kDy);
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsN >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsM < 0.0 && lsN < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsN, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsN) + std::fabs(lsM));
			theta = std::fabs(lsN) / (std::fabs(lsN) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoN[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j + 1, k)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FN = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1, k)] - H[idx(i, j, k)]);
		FN *= iRhoN[idx(i, j, k)] / (kDy * kDy);

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsB >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsB < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsB, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsB) + std::fabs(lsM));
			theta = std::fabs(lsB) / (std::fabs(lsB) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoB[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j, k - 1)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FB = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i, j, k - 1)]);
		FB *= iRhoB[idx(i, j, k)] / (kDz * kDz);

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsM >= 0.0 && lsT >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsM < 0.0 && lsT < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsT, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsT) + std::fabs(lsM));
			theta = std::fabs(lsT) / (std::fabs(lsT) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoT[idx(i, j, k)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j, k + 1)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FT = m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k + 1)] - H[idx(i, j, k)]);
		FT *= iRhoT[idx(i, j, k)] / (kDz * kDz);
		
		// poisson equation form should be -\beta \nabla p = f
		// iRhoEff has already negative value of rhoEff, then it is just a coefficient.
		// For discretization of pressure gradient, additional term is negative and it goes to RHS
		// Then it is a positive value and don't worrry about the sign
		rhs[idx(i, j, k)] -= FW + FE + FS + FN + FB + FT;
		
		assert(rhs[idx(i, j, k)] == rhs[idx(i, j, k)]);
		if (std::isnan(rhs[idx(i, j, k)]) || std::isinf(rhs[idx(i, j, k)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " << k << " " 
				<< rhs[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}
	
	// Original value of RHS
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		rhs[idx(i, j, k)] -= div[idx(i, j, k)];

	// An order of A matrix coef. is very important, hence reverse j order
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
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

		iRhoB[idx(i, j, k)] *= -1.0;
		iRhoS[idx(i, j, k)] *= -1.0;
		iRhoW[idx(i, j, k)] *= -1.0;
		iRhoE[idx(i, j, k)] *= -1.0;
		iRhoN[idx(i, j, k)] *= -1.0;
		iRhoT[idx(i, j, k)] *= -1.0;

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["B"] = iRhoB[idx(i, j, k)] / (kDz * kDz);
		AColsDic["B"] = (k - 1 - kNumBCGrid) + kNz * (j - kNumBCGrid) + kNz * kNy * (i - kNumBCGrid);
		AValsDic["S"] = iRhoS[idx(i, j, k)] / (kDy * kDy);
		AColsDic["S"] = (k - kNumBCGrid) + kNz * (j - 1 - kNumBCGrid) + kNz * kNy * (i - kNumBCGrid);
		AValsDic["W"] = iRhoW[idx(i, j, k)] / (kDx * kDx);
		AColsDic["W"] = (k - kNumBCGrid) + kNz * (j - kNumBCGrid) + kNz * kNy * (i - 1 - kNumBCGrid);
		AValsDic["C"] = -(iRhoW[idx(i, j, k)] + iRhoE[idx(i, j, k)]) / (kDx * kDx)
			- (iRhoS[idx(i, j, k)] + iRhoN[idx(i, j, k)]) / (kDy * kDy)
			- (iRhoB[idx(i, j, k)] + iRhoT[idx(i, j, k)]) / (kDz * kDz);
		AColsDic["C"] = (k - kNumBCGrid) + kNz * (j - kNumBCGrid) + kNz * kNy * (i - kNumBCGrid);
		AValsDic["E"] = iRhoE[idx(i, j, k)] / (kDx * kDx);
		AColsDic["E"] = (k - kNumBCGrid) + kNz * (j - kNumBCGrid) + kNz * kNy * (i  + 1 - kNumBCGrid);
		AValsDic["N"] = iRhoN[idx(i, j, k)] / (kDy * kDy);
		AColsDic["N"] = (k - kNumBCGrid) + kNz * (j + 1- kNumBCGrid) + kNz * kNy * (i - kNumBCGrid);
		AValsDic["T"] = iRhoT[idx(i, j, k)] / (kDz * kDz);
		AColsDic["T"] = (k + 1 - kNumBCGrid) + kNz * (j - kNumBCGrid) + kNz * kNy * (i - kNumBCGrid);
		
		if (i == kNumBCGrid && m_BC->m_BC_PW == BC3D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += iRhoW[idx(i, j, k)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC3D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= iRhoW[idx(i, j, k)] / (kDx * kDx);
			rhs[idx(i, j, k)] -= iRhoW[idx(i, j, k)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPW);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC3D::PERIODIC) {
			AValsDic["W"] = iRhoW[idx(kNumBCGrid + kNx - 1, j, k)];
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC3D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += iRhoE[idx(i, j, k)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC3D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= iRhoE[idx(i, j, k)] / (kDx * kDx);
			rhs[idx(i, j, k)] -= iRhoE[idx(i, j, k)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPE);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC3D::PERIODIC) {
			AValsDic["E"] = iRhoE[idx(kNumBCGrid, j, k)];
		}

		if (j == kNumBCGrid && m_BC->m_BC_PS == BC3D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += iRhoS[idx(i, j, k)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC3D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= iRhoS[idx(i, j, k)] / (kDy * kDy);
			rhs[idx(i, j, k)] -= iRhoS[idx(i, j, k)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPS);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC3D::PERIODIC) {
			AValsDic["S"] = iRhoS[idx(i, kNumBCGrid + kNy - 1, k)];
		}

		if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC3D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += iRhoN[idx(i, j, k)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC3D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= iRhoN[idx(i, j, k)] / (kDy * kDy);
			rhs[idx(i, j, k)] -= iRhoN[idx(i, j, k)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPN);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC3D::PERIODIC) {
			AValsDic["N"] = iRhoN[idx(i, kNumBCGrid + kNy + 1, k)];
		}

		if (k == kNumBCGrid && m_BC->m_BC_PB == BC3D::NEUMANN) {
			AColsDic["B"] = -1;
			AValsDic["B"] = 0.0;
			AValsDic["C"] += iRhoB[idx(i, j, k)] / (kDz * kDz);
		}
		else if (k == kNumBCGrid && m_BC->m_BC_PB == BC3D::DIRICHLET) {
			AColsDic["B"] = -1;
			AValsDic["B"] = 0.0;
			AValsDic["C"] -= iRhoB[idx(i, j, k)] / (kDz * kDz);
			rhs[idx(i, j, k)] -= iRhoB[idx(i, j, k)] / (kDz * kDz) * (m_BC->m_BC_DirichletConstantPB);
		}
		else if (k == kNumBCGrid && m_BC->m_BC_PB == BC3D::PERIODIC) {
			AValsDic["B"] = iRhoB[idx(i, j, kNumBCGrid + kNz - 1)];
		}

		if (k == kNumBCGrid + kNz - 1 && m_BC->m_BC_PT == BC3D::NEUMANN) {
			AColsDic["T"] = -1;
			AValsDic["T"] = 0.0;
			AValsDic["C"] += iRhoT[idx(i, j, k)] / (kDz * kDz);
		}
		else if (k == kNumBCGrid + kNz - 1 && m_BC->m_BC_PT == BC3D::DIRICHLET) {
			AColsDic["T"] = -1;
			AValsDic["T"] = 0.0;
			AValsDic["C"] -= iRhoT[idx(i, j, k)] / (kDz * kDz);
			rhs[idx(i, j, k)] -= iRhoT[idx(i, j, k)] / (kDz * kDz) * (m_BC->m_BC_DirichletConstantPT);
		}
		else if (k == kNumBCGrid + kNz - 1 && m_BC->m_BC_PT == BC3D::PERIODIC) {
			AValsDic["T"] = iRhoT[idx(i, j, kNumBCGrid + kNz + 1)];
		}

		// add non zero values to AVals and ACols
		// KEEP ORDER OF PUSH_BACK!!
		if (AColsDic["B"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["B"]);
			ACols.push_back(AColsDic["B"]);
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

		if (AColsDic["T"] >= 0) {
			tmpRowIdx++;
			AVals.push_back(AValsDic["T"]);
			ACols.push_back(AColsDic["T"]);
		}
		
		tmpMRowIdx++;
		DiagVals.push_back(AValsDic["C"]);
		DiagCols.push_back(AColsDic["C"]);
		
		rowIdx += tmpRowIdx;
		MRowIdx += tmpMRowIdx;
		assert(rhs[idx(i, j, k)] == rhs[idx(i, j, k)]);
		if (std::isnan(rhs[idx(i, j, k)]) || std::isinf(rhs[idx(i, j, k)])) {
			std::cout << "right hand side of poisson equation nan/inf error : " 
				<< i << " " << j << " " << rhs[idx(i, j, k)] << std::endl;
			exit(1);
		}
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);
	
	if (m_PoissonSolverType == POISSONTYPE::CG) {
		// std::cout << "Poisson : CG" << std::endl;
		m_Poisson->CG_2FUniform_3D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenX, kLenY, kLenZ, kDx, kDy, kDz, m_BC, maxIter);
	
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
			if (std::isnan(ps[idx(i, j, k)]) || std::isinf(ps[idx(i, j, k)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << k << " " << ps[idx(i, j, k)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
		// std::cout << "Poisson : BiCG" << std::endl;
		m_Poisson->BiCGStab_2FUniform_3D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenX, kLenY, kLenZ, kDx, kDy, kDz, m_BC, maxIter);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
			if (std::isnan(ps[idx(i, j, k)]) || std::isinf(ps[idx(i, j, k)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << k << " " << ps[idx(i, j, k)] << std::endl;
	}
	else {
		std::cout << "Poisson Solver not set!" << std::endl;
		exit(1);
	}
	
	ApplyBC_P_3D(ps);
	
	return 0;
}

std::vector<double> MACSolver3D::GetDivergence(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {
	std::vector<double> div(kArrSize, 0.0);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		div[idx(i, j, k)] = (u[idx(i + 1, j, k)] - u[idx(i, j, k)]) / kDx
		+ (v[idx(i, j + 1, k)] - v[idx(i, j, k)]) / kDy
		+ (w[idx(i, j, k + 1)] - w[idx(i, j, k)]) / kDz;
	
	return div;
}

int MACSolver3D::UpdateVel(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w,
	const std::vector<double>& us, const std::vector<double>& vs, const std::vector<double>& ws,
	const std::vector<double>& ps, const std::vector<double>& ls, const std::vector<double>& lsBackup,
	const std::vector<double>& H) {
	
	// velocity update after solving poisson equation
	// ps = p * dt
	double lsW = 0.0, lsM = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsB = 0.0, lsT = 0.0;
	double uEff = 0.0, vEff = 0.0, rhoEff = 0.0, theta = 0.0, thetaH = 0.0, iRhoEff = 0.0;
	double nXEff = 0.0, nYEff = 0.0, nZEff = 0.0, kappaEff = 0.0;
	double nXW = 0.0, nYW = 0.0, nZW = 0.0,
		nXS = 0.0, nYS = 0.0, nZS = 0.0,
		nXB = 0.0, nYB = 0.0, nZB = 0.0,
		nXM = 0.0, nYM = 0.0, nZM = 0.0;
	double aW = 0.0, aS = 0.0, aM = 0.0, aEff = 0.0;
	const double eps = 1.0e-100;
	std::vector<double> dudX(kArrSize, 0.0), dudY(kArrSize, 0.0), dudZ(kArrSize, 0.0),
		dvdX(kArrSize, 0.0), dvdY(kArrSize, 0.0), dvdZ(kArrSize, 0.0),
		dwdX(kArrSize, 0.0), dwdY(kArrSize, 0.0), dwdZ(kArrSize, 0.0),
		dldX(kArrSize, 0.0), dldY(kArrSize, 0.0), dldZ(kArrSize, 0.0);

	std::vector<double> U_PGrid(kArrSize, 0.0),
		V_PGrid(kArrSize, 0.0);
	double dudXW = 0.0, dudXE = 0.0, dudYS = 0.0, dudYN = 0.0;
	double dvdXW = 0.0, dvdXE = 0.0, dvdYS = 0.0, dvdYN = 0.0;
	
	/*
	// http://ctr.stanford.edu/Summer/SP08/3_1_Moureau.pdf
	// Moreau, V., and O. Desjardins. 
	// "A second-order ghost-fluid method for the primary atomization
	//	 of liquid fuel in air-blast type injectors."
	// Proceedings of the Summer Program. Vol. 143. 2008.
	*/
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		lsW = ls[idx(i - 1, j, k)];
		lsM = ls[idx(i, j, k)];
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsW >= 0 && lsM >= 0) {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsW <= 0 && lsM <= 0) {
			thetaH = 0.0; theta = 0.0;
		}
		else  {
			thetaH = (std::max(lsW, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsW) + std::fabs(lsM));
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i - 1, j, k)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		u[idx(i, j, k)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i - 1, j, k)]);
		u[idx(i, j, k)] *= iRhoEff / kDx;

		u[idx(i, j, k)] += us[idx(i, j, k)] - iRhoEff * (ps[idx(i, j, k)] - ps[idx(i - 1, j, k)]) / kDx;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		lsM = ls[idx(i, j, k)];
		lsS = ls[idx(i, j - 1, k)];
		
		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsS >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsS < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsS, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsS) + std::fabs(lsM));
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j - 1, k)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		v[idx(i, j, k)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i, j - 1, k)]);
		v[idx(i, j, k)] *= iRhoEff / kDy;

		v[idx(i, j, k)] += vs[idx(i, j, k)] - iRhoEff * (ps[idx(i, j, k)] - ps[idx(i, j - 1, k)]) / kDy;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++) {
		lsM = ls[idx(i, j, k)];
		lsB = ls[idx(i, j, k - 1)];

		// thetaH = portion of kRhoH
		// theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
		if (lsB >= 0.0 && lsM >= 0.0)  {
			thetaH = 1.0; theta = 0.0;
		}
		else if (lsB < 0.0 && lsM < 0.0) {
			thetaH = 0.0; theta = 0.0;
		}
		else {
			thetaH = (std::max(lsB, 0.0) + std::max(lsM, 0.0))
				/ (std::fabs(lsB) + std::fabs(lsM));
			theta = std::fabs(lsB) / (std::fabs(lsB) + std::fabs(lsM));
		}

		rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
		iRhoEff = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j, k)] * theta + m_kappa[idx(i, j, k - 1)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// + level set : H = 1, - level set : H = 0; 
		// Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
		w[idx(i, j, k)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j, k)] - H[idx(i, j, k - 1)]);
		w[idx(i, j, k)] *= iRhoEff / kDz;

		w[idx(i, j, k)] += ws[idx(i, j, k)] - iRhoEff * (ps[idx(i, j, k)] - ps[idx(i, j, k - 1)]) / kDz;
	}
	ApplyBC_U_3D(u);
	ApplyBC_V_3D(v);
	ApplyBC_W_3D(w);

	return 0;
}

double MACSolver3D::UpdateDt(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {
	double uAMax = 0.0, vAMax = 0.0, wAMax = 0.0;
	double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0;
	double dt = std::numeric_limits<double>::max();

	// get maximum of absolute value
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		uAMax = std::max(uAMax, std::fabs((u[idx(i + 1, j, k)] * u[idx(i, j, k)]) * 0.5));
		vAMax = std::max(vAMax, std::fabs((v[idx(i, j + 1, k)] * v[idx(i, j, k)]) * 0.5));
		wAMax = std::max(wAMax, std::fabs((w[idx(i, j, k + 1)] * w[idx(i, j, k)]) * 0.5));
	}
	
	Cefl = uAMax / kDx + vAMax / kDy + wAMax / kDz;
	Vefl = std::max(kMuH / kRhoH, kMuL / kRhoL) * (2.0 / (kDx * kDx) + 2.0 / (kDy * kDy));
	Gefl = std::max(std::max(std::sqrt(std::fabs(kG) / kDy), std::sqrt(std::fabs(kG) / kDx)), std::sqrt(std::fabs(kG) / kDz));
	
	dt = std::min(dt,
		1.0 / (0.5 * (Cefl + Vefl +
		std::sqrt(std::pow(Cefl + Vefl, 2.0) +
		4.0 * Gefl))));

	if (std::isnan(dt) || std::isinf(dt)) {
		std::cout << "dt nan/inf error : Cefl, Vefl, Gefl : " << Cefl << " " << Vefl << " " << Gefl << std::endl;
	}

	return dt;
}

int MACSolver3D::SetBC_U_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_U_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int MACSolver3D::SetBC_V_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_V_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int MACSolver3D::SetBC_W_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_W_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int MACSolver3D::SetBC_P_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_P_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int MACSolver3D::ApplyBC_U_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_U_3D(arr);
	return 0;
}

int MACSolver3D::ApplyBC_V_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_V_3D(arr);
	return 0;
}

int MACSolver3D::ApplyBC_W_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_W_3D(arr);
	return 0;
}

int MACSolver3D::ApplyBC_P_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_P_3D(arr);
	return 0;
}

void MACSolver3D::SetBCConstantUW(double BC_ConstantW) {
	return m_BC->SetBCConstantUW(BC_ConstantW);
}

void MACSolver3D::SetBCConstantUE(double BC_ConstantE) {
	return m_BC->SetBCConstantUE(BC_ConstantE);
}

void MACSolver3D::SetBCConstantUS(double BC_ConstantS) {
	return m_BC->SetBCConstantUS(BC_ConstantS);
}

void MACSolver3D::SetBCConstantUN(double BC_ConstantN) {
	return m_BC->SetBCConstantUN(BC_ConstantN);
}

void MACSolver3D::SetBCConstantUB(double BC_ConstantB) {
	return m_BC->SetBCConstantUB(BC_ConstantB);
}

void MACSolver3D::SetBCConstantUT(double BC_ConstantT) {
	return m_BC->SetBCConstantUT(BC_ConstantT);
}

void MACSolver3D::SetBCConstantVW(double BC_ConstantW) {
	return m_BC->SetBCConstantVW(BC_ConstantW);
}

void MACSolver3D::SetBCConstantVE(double BC_ConstantE) {
	return m_BC->SetBCConstantVE(BC_ConstantE);
}

void MACSolver3D::SetBCConstantVS(double BC_ConstantS) {
	return m_BC->SetBCConstantVS(BC_ConstantS);
}

void MACSolver3D::SetBCConstantVN(double BC_ConstantN) {
	return m_BC->SetBCConstantVN(BC_ConstantN);
}

void MACSolver3D::SetBCConstantVB(double BC_ConstantB) {
	return m_BC->SetBCConstantVB(BC_ConstantB);
}

void MACSolver3D::SetBCConstantVT(double BC_ConstantT) {
	return m_BC->SetBCConstantVT(BC_ConstantT);
}

void MACSolver3D::SetBCConstantWW(double BC_ConstantW) {
	return m_BC->SetBCConstantWW(BC_ConstantW);
}

void MACSolver3D::SetBCConstantWE(double BC_ConstantE) {
	return m_BC->SetBCConstantWE(BC_ConstantE);
}

void MACSolver3D::SetBCConstantWS(double BC_ConstantS) {
	return m_BC->SetBCConstantWS(BC_ConstantS);
}

void MACSolver3D::SetBCConstantWN(double BC_ConstantN) {
	return m_BC->SetBCConstantWN(BC_ConstantN);
}

void MACSolver3D::SetBCConstantWB(double BC_ConstantB) {
	return m_BC->SetBCConstantWB(BC_ConstantB);
}

void MACSolver3D::SetBCConstantWT(double BC_ConstantT) {
	return m_BC->SetBCConstantWT(BC_ConstantT);
}

void MACSolver3D::SetBCConstantPW(double BC_ConstantW) {
	return m_BC->SetBCConstantPW(BC_ConstantW);
}

void MACSolver3D::SetBCConstantPE(double BC_ConstantE) {
	return m_BC->SetBCConstantPE(BC_ConstantE);
}

void MACSolver3D::SetBCConstantPS(double BC_ConstantS) {
	return m_BC->SetBCConstantPS(BC_ConstantS);
}

void MACSolver3D::SetBCConstantPN(double BC_ConstantN) {
	return m_BC->SetBCConstantPN(BC_ConstantN);
}

void MACSolver3D::SetBCConstantPB(double BC_ConstantB) {
	return m_BC->SetBCConstantPB(BC_ConstantB);
}

void MACSolver3D::SetBCConstantPT(double BC_ConstantT) {
	return m_BC->SetBCConstantPT(BC_ConstantT);
}

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
inline int MACSolver3D::sign(const double& val) {
	return (val > 0) - (val < 0);
}

inline int MACSolver3D::idx(int i, int j, int k) {
	return (k + (kNz + 2 * kNumBCGrid) * j + (kNz + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid) * i);
}

int MACSolver3D::SetPLTType(PLTTYPE type) {
	m_PLTType = type;

	return 0;
}

int MACSolver3D::OutRes(const int iter, const double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
	const std::vector<double>& ps, const std::vector<double>& ls) {
	if (m_PLTType == PLTTYPE::ASCII || m_PLTType == PLTTYPE::BOTH) {
		std::ofstream outF;
		std::string fname_vel(fname_vel_base + "_ASCII.plt");
		std::string fname_div(fname_div_base + "_ASCII.plt");
		if (m_iter == 0) {
			outF.open(fname_vel.c_str(), std::ios::out);

			outF << "TITLE = VEL" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"LS\", \"PS\" " << std::endl;
			outF.close();

			outF.open(fname_div.c_str(), std::ios::out);
			outF << "TITLE = DIV" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"Z\", \"LS\", \"DIV\", \"PS\" " << std::endl;
			outF.close();
		}

		std::vector<double>
			resU(kArrSize, 0.0), resV(kArrSize, 0.0), resW(kArrSize, 0.0);

		std::vector<double> resDiv = GetDivergence(u, v, w);
		
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			resU[idx(i, j, k)] = (u[idx(i, j, k)] + u[idx(i + 1, j, k)]) * 0.5;
			resV[idx(i, j, k)] = (v[idx(i, j, k)] + v[idx(i, j + 1, k)]) * 0.5;
			resW[idx(i, j, k)] = (w[idx(i, j, k)] + w[idx(i, j, k + 1)]) * 0.5;
		}
		

		outF.open(fname_vel.c_str(), std::ios::app);
		
		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy << std::string(", K=") << kNz
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;

		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
				<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
				<< kBaseZ + static_cast<double>(k + 0.5 - kNumBCGrid) * kDz << std::string(",")
				<< static_cast<double>(resU[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(resV[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(resW[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(ls[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j, k)]) << std::endl;

		outF.close();
		
		outF.open(fname_div.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << kNx << std::string(", J=") << kNy << std::string(", K=") << kNz
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter + 1
			<< std::endl;
		
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
				<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
				<< kBaseZ + static_cast<double>(k + 0.5 - kNumBCGrid) * kDz << std::string(",")
				<< static_cast<double>(ls[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(resDiv[idx(i, j, k)]) << std::string(",")
				// << static_cast<double>(m_kappa[idx(i, j, k)]) << std::string(",")
				<< static_cast<double>(ps[idx(i, j, k)]) << std::endl;

		outF.close();
	}

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH) {
		INTEGER4 whichFile = 0, stat = 0;
		std::string fname_vel(fname_vel_base + "_BINARY.szplt");
		std::string fname_div(fname_div_base + "_BINARY.szplt");

		std::vector<double>
			resX(kNx * kNy * kNz, 0.0),
			resY(kNx * kNy * kNz, 0.0),
			resZ(kNx * kNy * kNz, 0.0),
			resU(kNx * kNy * kNz, 0.0),
			resV(kNx * kNy * kNz, 0.0),
			resW(kNx * kNy * kNz, 0.0),
			resDiv(kNx * kNy * kNz, 0.0),
			resLS(kNx * kNy * kNz, 0.0),
			resPs(kNx * kNy * kNz, 0.0);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) 
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			resX[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = kBaseX + (i + 0.5 - kNumBCGrid) * kDx;
			resY[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = kBaseY + (j + 0.5 - kNumBCGrid) * kDy;
			resZ[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = kBaseZ + (k + 0.5 - kNumBCGrid) * kDz;
			resU[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = (u[idx(i, j, k)] + u[idx(i + 1, j, k)]) * 0.5;
			resV[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = (v[idx(i, j, k)] + v[idx(i, j + 1, k)]) * 0.5;
			resW[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = (w[idx(i, j, k)] + w[idx(i, j, k + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = ls[idx(i, j, k)];
			resPs[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = ps[idx(i, j, k)];
			resDiv[(i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid)] = (u[idx(i + 1, j, k)] - u[idx(i, j, k)]) / kDx
				+ (v[idx(i, j + 1, k)] - v[idx(i, j, k)]) / kDy + (w[idx(i, j, k + 1)] - w[idx(i, j, k)]) / kDz;
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
							const_cast<char *>(std::string("X, Y, Z, U, V, W, LS, PS").c_str()),  
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
		INTEGER4 IMax = kNx, JMax = kNy, KMax = kNz;

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
		stat = TECDAT142(&ARRSIZEVAL, resZ.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resU.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resV.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resW.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);

		if (m_iter == 0) {
			stat = TECINI142(const_cast<char *>(std::string("DIVERGENCE").c_str()),  /* Name of the entire dataset.  */
				const_cast<char *>(std::string("X, Y, Y, LS, DIV, PS").c_str()),
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
		stat = TECDAT142(&ARRSIZEVAL, resZ.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resDiv.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);
	}

	return 0;
}

int MACSolver3D::OutResClose() {
	INTEGER4 stat, whichFile;
	whichFile = 1;
	stat = TECFIL142(&whichFile);

	// close first file (velocity)
	stat = TECEND142();
	stat = TECEND142();

	return 0;
}
