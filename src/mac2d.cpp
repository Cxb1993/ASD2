#include "mac2d.h"

MACSolver2D::MACSolver2D(double Re, double We, double Fr, GAXISENUM GAxis,
	double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoH, double muH,
	int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRe(Re), kWe(We), kFr(Fr),
	kLScale(L), kUScale(U), kSigma(sigma),
	kG(kFr * L / (U * U)), kGAxis(GAxis), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
	kRhoH(rhoH), kRhoL(rhoH * densityRatio), kRhoRatio(densityRatio),
	kMuH(muH), kMuL(muH * viscosityRatio), kMuRatio(viscosityRatio),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK) {

	// positive level set : inside
	// negative level set : outside
	m_iter = 0;
	m_curTime = 0.0;
	m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::MACSolver2D(double rhoH, double rhoL, double muH, double muL, double gConstant, GAXISENUM GAxis,
	double L, double U, double sigma, int nx, int ny, double baseX, double baseY, double lenX, double lenY,
	TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
	bool writeVTK) :
	kRhoScale(rhoH), kMuScale(muH), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
	kRhoH(rhoH), kRhoL(rhoL), kRhoRatio(rhoL / rhoH), 
	kMuH(muH), kMuL(muL), kMuRatio(muL / muH),
	kRe(rhoH * L * U / muH), kWe(rhoH * L * U * U / sigma), kFr(U * U / (gConstant * L)), kGAxis(GAxis),
	kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
	kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
	kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
	kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK) {

	// positive level set : inside
	// negative level set : outside
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

	m_nx = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_ny = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);

	m_t1x = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_t1y = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_t1z = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);

	m_t2x = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_t2y = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);
	m_t2z = std::vector<double>((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0);

	return 0;
}

std::vector<double> MACSolver2D::UpdateHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	const double eps = std::min(kDx, kDy) * 1.5;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
		// inside
		if (ls[idx(i, j)] >= 0.0)
			Heaviside[idx(i, j)] = 1.0;
		else
			Heaviside[idx(i, j)] = 0.0;
	}

	return Heaviside;
}

std::vector<double> MACSolver2D::UpdateSmoothHeavisideFunc(const std::vector<double>& ls) {
	std::vector<double> Heaviside((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	const double eps = std::min(kDx, kDy) * 1.5;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
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

int MACSolver2D::UpdateNTK(const std::shared_ptr<LevelSetSolver2D>& LSolver, const std::vector<double>& ls,
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

	std::vector<double> Q((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		absDerivLS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	absDerivLS = LSolver->ENO_DerivAbsLS_2D(ls, ls);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		Q[idx(i, j)] = std::fabs(1.0 - absDerivLS[idx(i, j)]);
	}
	
	double dLSdX = 0.0, dLSdY = 0.0;
	// determination function
	int Dx = 0.0, Dy = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		// complex level set, such as droplet merging
		if (Q[idx(i, j)] > eta) {
			// Using Dierctional Direction First
			if (Q[idx(i - 1, j)] < eta && Q[idx(i + 1, j)] >= eta)
				Dx = -1;
			else if (Q[idx(i - 1, j)] >= eta && Q[idx(i + 1, j)] < eta)
				Dx = 1;
			else if (Q[idx(i - 1, j)] < eta && Q[idx(i, j)] < eta && Q[idx(i + 1, j)] < eta)
				Dx = 0;
			else if (Q[idx(i - 1, j)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i + 1, j)] >= eta)
				Dx = 0;
			else
				Dx = 2;

			// determination funciton
			if (Q[idx(i, j - 1)] < eta && Q[idx(i, j + 1)] >= eta)
				Dy = -1;
			else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j + 1)] < eta)
				Dy = 1;
			else if (Q[idx(i, j - 1)] < eta && Q[idx(i, j)] < eta && Q[idx(i, j + 1)] < eta)
				Dy = 0;
			else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i, j + 1)] >= eta)
				Dy = 0;
			else
				// undetermined
				Dy = 2;

			if (Dx == -1)
				dLSdX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
			else if (Dx == 1)
				dLSdX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			else if (Dx == 0)
				dLSdX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);

			if (Dy == -1)
				dLSdY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
			else if (Dy == 1)
				dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			else if (Dy == 0)
				dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);

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
			dLSdX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
			dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
			std::get<0>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + eps) * dLSdX;
			std::get<1>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + eps) * dLSdY;
		}
	}

	return 0;
}

int MACSolver2D::UpdateKappa(const std::vector<double>& ls) {
	std::vector<double> dLSdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dLSdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LSSize((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	const double eps = 1.0e-100;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		
		LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);

		// if size is zero, use forward or backward differencing rather than central differencing
		if (LSSize[idx(i, j)] < eps) {
			dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			if (LSSize[idx(i, j)] < eps) {
				dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
				LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)]);
			}
		}
		if (LSSize[idx(i, j)] < eps)
			perror("Div/0 Err in computing kappa");
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		m_kappa[idx(i, j)]
			= -(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] * (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDy * kDy) // phi^2_x \phi_yy
			- 2.0 * dLSdX[idx(i, j)] * dLSdY[idx(i, j)] * (ls[idx(i + 1, j + 1)] - ls[idx(i - 1, j + 1)] -ls[idx(i + 1, j - 1)] + ls[idx(i - 1, j - 1)]) / (4.0 * kDx * kDy) //2 \phi_x \phi_y \phi_xy
			// - 2.0 * dLSdX[idx(i, j)] * dLSdY[idx(i, j)] * (dLSdX[idx(i, j + 1)] - dLSdX[idx(i, j - 1)]) / (2.0 * kDy) //2 \phi_x \phi_y \phi_xy
			+ dLSdY[idx(i, j)] * dLSdY[idx(i, j)] * (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / (kDx * kDx)) // phi^2_y \phi_xx
				/ std::pow(LSSize[idx(i, j)], 3.0);

		// curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
		m_kappa[idx(i, j)] = std::fabs(std::min(std::fabs(m_kappa[idx(i, j)]), 1.0 / std::min(kDx, kDy)));
		
		assert(m_kappa[idx(i, j)] == m_kappa[idx(i, j)]);
	}

	return 0;
}

int MACSolver2D::UpdateJumpCond(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls) {
	const double eps = 1.0e-100;
	std::vector<double> U_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		V_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	// derivatives
	std::vector<double> dudX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dudY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dvdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dvdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dldX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> dldY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	// normal and tangent vector variable (n, t1, t2)
	// normal vector has X, Y component
	// first tangent vector has only Z component
	// double nX = 0.0, nY = 0.0, t1X = 0.0, t1Y = 0.0, t2X = 0.0, t2Y = 0.0, t2Z = 0.0;
	std::vector<double> nX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> nY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t1X((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t1Y((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> t2Z((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	// jump condition originally defined as J11 = [\mu u_x], J12 = [\mu u_y], J21 = [\mu v_x], J22 = [\mu v_y], 
	// but skip [\mu], from Kang, Fedkiw, and Liu's work eq. (30).
	// [\mu] changes with level set, which means I can determine [\mu] later
	
	double nxnx = 0.0, nxny = 0.0, nyny = 0.0;
	double t1xt1x = 0.0, t1xt1y = 0.0, t1yt1y = 0.0;
	double theta = 0.0;
	double lsW = 0.0, lsE = 0.0, lsM = 0.0, lsS = 0.0, lsN = 0.0;
	double lsUW = 0.0, lsUM = 0.0, lsUE = 0.0, lsUSHalf = 0.0, lsUNHalf = 0.0;
	double lsVS = 0.0, lsVM = 0.0, lsVN = 0.0, lsVWHalf = 0.0, lsVEHalf = 0.0;

	double uEff = 0.0, vEff = 0.0;
	double dudXW = 0.0, dudXE = 0.0, dudYS = 0.0, dudYN = 0.0;
	double dvdXW = 0.0, dvdXE = 0.0, dvdYS = 0.0, dvdYN = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		// defined at P grid
		dudX[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx;
		dudY[idx(i, j)] = (0.25 * (u[idx(i, j)] + u[idx(i + 1, j)] + u[idx(i, j + 1)] + u[idx(i + 1, j + 1)])
		- 0.25 * (u[idx(i, j)] + u[idx(i + 1, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)])) / kDy;
		dvdX[idx(i, j)] = (0.25 * (v[idx(i, j)] + v[idx(i, j + 1)] + v[idx(i + 1, j)] + v[idx(i + 1, j + 1)])
			- 0.25 * (v[idx(i, j)] + v[idx(i, j + 1)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)])) / kDx;
		dvdY[idx(i, j)] = (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
		
		dudX[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dudY[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDy);
		dvdX[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dvdY[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDy);

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
		
		// Jump occurs when computing dudX
		if (lsW * lsE <= 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudXW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dudXE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDx);
			dudX[idx(i, j)] = dudXE * theta + dudXW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdXW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dvdXE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDx);
			dvdX[idx(i, j)] = dvdXE * theta + dvdXW * (1.0 - theta);
		}

		if (lsS * lsN <= 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudYS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dudYN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDy);
			dudY[idx(i, j)] = dudYN * theta + dudYS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdYS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dvdYN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDy);
			dvdY[idx(i, j)] = dvdYN * theta + dvdYS * (1.0 - theta);
		}
		
		dldX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
		dldY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
		if (std::fabs(dldX[idx(i, j)]) < eps) {
			dldX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
			if (std::fabs(dldX[idx(i, j)]) < eps)
				dldX[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		}
		if (std::fabs(dldX[idx(i, j)]) < eps) {
			dldY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
			if (std::fabs(dldY[idx(i, j)]) < eps)
				dldY[idx(i, j)] = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
		}
		
		// normal vector = (\nabla \phi) / |\nabla \phi|
		nX[idx(i, j)] = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);
		nY[idx(i, j)] = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps);
		// get first tangent vector
		// smallest magnitude determine what to be performed cross product
		// for 2D, smallest magnitude must be z axis (because it is 2D!)
		t1X[idx(i, j)] = nY[idx(i, j)] / (std::sqrt(std::pow(nX[idx(i, j)], 2.0) + std::pow(nY[idx(i, j)], 2.0)) + eps);
		t1Y[idx(i, j)] = -nX[idx(i, j)] / (std::sqrt(std::pow(nX[idx(i, j)], 2.0) + std::pow(nY[idx(i, j)], 2.0)) + eps);
		t2Z[idx(i, j)] = (nX[idx(i, j)] * t1Y[idx(i, j)] - nY[idx(i, j)] * t1X[idx(i, j)]);

		nxnx = nX[idx(i, j)] * nX[idx(i, j)];
		nxny = nX[idx(i, j)] * nY[idx(i, j)];
		nyny = nY[idx(i, j)] * nY[idx(i, j)];
		t1xt1x = t1X[idx(i, j)] * t1X[idx(i, j)];
		t1xt1y = t1X[idx(i, j)] * t1Y[idx(i, j)];
		t1yt1y = t1Y[idx(i, j)] * t1Y[idx(i, j)];

		// eq. (30) from Kang, Fedkiw, and Liu
		// u_x
		m_J11[idx(i, j)] = dudX[idx(i, j)] * t1xt1x + dudY[idx(i, j)] * t1xt1y
			+ (dudX[idx(i, j)] * nxnx + dvdX[idx(i, j)] * nxny) * nxnx + (dudY[idx(i, j)] * nxnx + dvdY[idx(i, j)] * nxny) * nxny
			- ((dudX[idx(i, j)] * t1xt1x + dudY[idx(i, j)] * t1xt1y) * nxnx + (dvdX[idx(i, j)] * t1xt1x + dvdY[idx(i, j)] * t1xt1y) * nxny);
		m_J11[idx(i, j)] *= (kMuH - kMuL);

		// u_y
		m_J12[idx(i, j)] = dudX[idx(i, j)] * t1xt1y + dudY[idx(i, j)] * t1yt1y
			+ (dudX[idx(i, j)] * nxnx + dvdX[idx(i, j)] * nxny) * nxny + (dudY[idx(i, j)] * nxnx + dvdY[idx(i, j)] * nxny) * nyny
			- ((dudX[idx(i, j)] * t1xt1x + dudY[idx(i, j)] * t1xt1y) * nxny + (dvdX[idx(i, j)] * t1xt1x + dvdY[idx(i, j)] * t1xt1y) * nyny);
		m_J12[idx(i, j)] *= (kMuH - kMuL);

		// v_x
		m_J21[idx(i, j)] = dvdX[idx(i, j)] * t1xt1x + dvdY[idx(i, j)] * t1xt1y
			+ (dudX[idx(i, j)] * nxny + dvdX[idx(i, j)] * nyny) * nxnx + (dudY[idx(i, j)] * nxny + dvdY[idx(i, j)] * nyny) * nxny
			- ((dudX[idx(i, j)] * t1xt1y + dudY[idx(i, j)] * t1yt1y) * nxnx + (dvdX[idx(i, j)] * t1xt1y + dvdY[idx(i, j)] * t1yt1y) * nxny);
		m_J21[idx(i, j)] *= (kMuH - kMuL);

		// v_y
		m_J22[idx(i, j)] = dvdX[idx(i, j)] * t1xt1y + dvdY[idx(i, j)] * t1yt1y
			+ (dudX[idx(i, j)] * nxny + dvdX[idx(i, j)] * nyny) * nxny + (dudY[idx(i, j)] * nxny + dvdY[idx(i, j)] * nyny) * nyny
			- ((dudX[idx(i, j)] * t1xt1y + dudY[idx(i, j)] * t1yt1y) * nxny + (dvdX[idx(i, j)] * t1xt1y + dvdY[idx(i, j)] * t1yt1y) * nyny);
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

std::vector<double> MACSolver2D::UpdateFU(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& H) {
	
	std::vector<double> cU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> vU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> gU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhsU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Convection term
	cU = this->AddConvectionFU(u, v);

	// Viscous term
	vU = this->AddViscosityFU(u, v, ls, H);

	gU = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		rhsU[idx(i, j)] = -cU[idx(i, j)] + vU[idx(i, j)] + gU[idx(i, j)];
	}
	
	return rhsU;
}

std::vector<double> MACSolver2D::UpdateFV(const std::shared_ptr<LevelSetSolver2D>& LSolver,
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& H) {

	std::vector<double> cV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> vV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> gV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhsV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	// Convection term
	cV = this->AddConvectionFV(u, v);

	// Viscous term
	vV = this->AddViscosityFV(u, v, ls, H);

	gV = this->AddGravityFU();

	// Get RHS(Right Hand Side)
	// level set
	double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
	double theta = 0.0, iRhoEff = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		rhsV[idx(i, j)] = -cV[idx(i, j)] + vV[idx(i, j)] + gV[idx(i, j)];
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
			vecF_UY[j] = u[idx(i, j)] * tmpV[idx(i, j)];
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
			vecF_VX[i] = v[idx(i, j)] * tmpU[idx(i, j)];
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
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> mu((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> visXVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		visYVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		UxxVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		UyyVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		VxyVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		VxyNVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		VxySVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> muU_W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muU_E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muU_S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muU_N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muV_W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muV_E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ11W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0), 
		resJ11E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ12S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0), 
		resJ12N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ21W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ21E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		checkDiv((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dU[idx(i, j)] = 0.0;
		
		return dU;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
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
	double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
	double muV_X_S = 0.0, muV_X_N = 0.0;
	double rhoU_X_W = 0.0, rhoU_X_E = 0.0, rhoU_Y_S = 0.0, rhoU_Y_N = 0.0;
	double rhoV_X_S = 0.0, rhoV_X_N = 0.0;
	double visX = 0.0, visY = 0.0;

	// jump condition
	double JUW = 0.0, JUE = 0.0, JUS = 0.0, JUN = 0.0, JUM = 0.0;
	double JVW = 0.0, JVW_N = 0.0, JVM = 0.0, JVN = 0.0;
	// effective Jump condition, effective u(uEff), effective mu (muEff), and effective rho (rhoEff)
	// J is a jump condition and defined at P grid
	double JEff = 0.0, JO = 0.0, uEff = 0.0, vEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
	const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

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
			muU_X_W = kMuH * (uM - uW) / kDx;
			muU_W[idx(i, j)] = muU_X_W;
		}
		else if (lsUW <= 0 && lsUM <= 0) {
			muU_X_W = kMuL * (uM - uW) / kDx;
			muU_W[idx(i, j)] = muU_X_W;
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
			muU_X_W = muEff * (uM - uW) / kDx
				- muEff * JEff * theta / kMuH;
			muU_W[idx(i, j)] = muEff * (uM - uW) / kDx;
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
			muU_X_W = muEff * (uM - uW) / kDx
				+ muEff * JEff * theta / kMuL;
			muU_W[idx(i, j)] = muEff * (uM - uW) / kDx;
			resJ11W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsUM >= 0 && lsUE >= 0) {
			muU_X_E = kMuH * (uE - uM) / kDx;
			muU_E[idx(i, j)] = muU_X_E;
		}
		else if (lsUM <= 0 && lsUE <= 0) {
			muU_X_E = kMuL * (uE - uM) / kDx;
			muU_E[idx(i, j)] = muU_X_E;
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
			muU_X_E = muEff * (uE - uM) / kDx
				+ muEff * JEff * theta / kMuL;
			muU_E[idx(i, j)] = muEff * (uE - uM) / kDx;
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
			muU_X_E = muEff * (uE - uM) / kDx
				- muEff * JEff * theta / kMuH;
			muU_E[idx(i, j)] = muEff * (uE - uM) / kDx;
			resJ11E[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		
		if (lsUS >= 0 && lsUM >= 0) {
			muU_Y_S = kMuH * (uM - uS) / kDy;
			muU_S[idx(i, j)] = muU_Y_S;
		}
		else if (lsUS <= 0 && lsUM <= 0) {
			muU_Y_S = kMuL * (uM - uS) / kDy;
			muU_S[idx(i, j)] = muU_Y_S;
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
			muU_Y_S = muEff * (uM - uS) / kDy
				- muEff * JEff * theta / kMuH;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDy;
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
			muU_Y_S = muEff * (uM - uS) / kDy
				+ muEff * JEff * theta / kMuL;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDy;
			resJ12S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsUM >= 0 && lsUN >= 0) {
			muU_Y_N = kMuH * (uN - uM) / kDy;
			muU_N[idx(i, j)] = muU_Y_N;
		}
		else if (lsUM <= 0 && lsUN <= 0) {
			muU_Y_N = kMuL * (uN - uM) / kDy;
			muU_N[idx(i, j)] = muU_Y_N;
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
			muU_Y_N = muEff * (uN - uM) / kDy 
				+ muEff * JEff * theta / kMuL;
			muU_N[idx(i, j)] = muEff * (uN - uM) / kDy;
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
			muU_Y_N = muEff * (uN - uM) / kDy 
				- muEff * JEff * theta / kMuH;
			muU_N[idx(i, j)] = muEff * (uN - uM) / kDy;
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
			muV_X_S = kMuH * (vM - vW) / kDx;
			muV_W[idx(i, j)] = muV_X_S;
		}
		else if (lsVW <= 0 && lsVM <= 0) {
			muV_X_S = kMuL * (vM - vW) / kDx;
			muV_W[idx(i, j)] = muV_X_S;
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
			muV_X_S = muEff * (vM - vW) / kDx
				- muEff * JEff * theta / kMuH;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDx;
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
			muV_X_S = muEff * (vM - vW) / kDx
				+ muEff * JEff * theta / kMuL;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDx;
			resJ21W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}

		if (lsVW_N >= 0 && lsVN >= 0) {
			muV_X_N = kMuH * (vN - vW_N) / kDx;
			muV_E[idx(i, j)] = muV_X_N;
		}
		else if (lsVW_N <= 0 && lsVN <= 0) {
			muV_X_N = kMuL * (vN - vW_N) / kDx;
			muV_E[idx(i, j)] = muV_X_N;
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
			muV_X_N = muEff * (vN - vW_N) / kDx
				- muEff * JEff * theta / kMuH;
			muV_E[idx(i, j)] = muEff * (vN - vW_N) / kDx;
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
			muV_X_N = muEff  * (vN - vW_N) / kDx
				+ muEff * JEff * theta / kMuL;
			muV_E[idx(i, j)] = muEff * (vN - vW_N) / kDx;
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

		muU_X_W = mu[idx(i - 1, j)] * (uM - uW) / kDx;
		muU_X_E = mu[idx(i, j)] * (uE - uM) / kDx;
		muU_Y_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (uM - uS) / kDy;
		muU_Y_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
			* (uN - uM) / kDy;
		muV_X_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (vM - vW) / kDx;
		muV_X_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
			* (vN - vW_N) / kDx;

		visX = 2.0 * (muU_X_E - muU_X_W) / kDx;
		visY = (muU_Y_N - muU_Y_S) / kDy + (muV_X_N - muV_X_S) / kDy;
		// checkDiv[idx(i, j)] = (muU_X_E - muU_X_W) / kDx + (muV_X_N - muV_X_S) / kDy;
		// dU[idx(i, j)] = visX / rhoEffWE + visY / rhoEffSN;
		dU[idx(i, j)] = visX * iRhoEffWE + visY * iRhoEffSN;
		// dU[idx(i, j)] = visX + visY;

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
		<< std::string("\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
			<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
			<< static_cast<double>(u[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(dU[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_E[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muU_E[idx(i, j)] - muU_W[idx(i, j)]) / kDx) << std::string(",")
			<< static_cast<double>((resJ11E[idx(i, j)] - resJ11W[idx(i, j)]) / (kDx)) << std::string(",")
			<< static_cast<double>(muU_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_N[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muU_N[idx(i, j)] - muU_S[idx(i, j)]) / kDx) << std::string(",")
			<< static_cast<double>((resJ12N[idx(i, j)] - resJ12S[idx(i, j)]) / (kDy)) << std::string(",")
			<< static_cast<double>(muV_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_E[idx(i, j)]) << std::string(",")
			<< static_cast<double>((muV_E[idx(i, j)] - muV_W[idx(i, j)]) / kDy) << std::string(",")
			<< static_cast<double>((resJ21E[idx(i, j)] - resJ21W[idx(i, j)]) / (kDx)) << std::string(",")
			<< static_cast<double>(checkDiv[idx(i, j)]) 
			<< std::endl;

	outF.close();
	*/
	return dU;
}

std::vector<double> MACSolver2D::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& ls, const std::vector<double>& H) {
	// This is incompressible viscous flow, which means velocity is CONTINUOUS!
	std::vector<double> dV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> mu((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> visXVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		visYVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		UxyVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		UxyWVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		UxyEVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		VxxVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		VyyVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dVConst((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> iRhoHVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoVVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> muV_W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muV_E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muV_S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muV_N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muU_S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		muU_N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ21W((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ21E((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ22S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ22N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ12S((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		resJ12N((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		checkDiv((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if (kRe <= 0.0) {
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			dV[idx(i, j)] = 0.0;

		return dV;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
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
	double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0;
	double muU_Y_W = 0.0, muU_Y_E = 0.0;
	double rhoV_X_W = 0.0, rhoV_X_E = 0.0, rhoV_Y_S = 0.0, rhoV_Y_N = 0.0;
	double rhoU_Y_W = 0.0, rhoU_Y_E = 0.0;
	double visX = 0.0, visY = 0.0;
	const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

	// effective Jump condition, effective v (vEff), effective mu (muEff), effective rho (rhoEff),
	double JEff = 0.0, JO = 0.0, vEff = 0.0, uEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
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
			muV_X_W = kMuH * (vM - vW) / kDx;
			muV_W[idx(i, j)] = muV_X_W;
		}
		else if (lsVW <= 0 && lsVM <= 0) {
			muV_X_W = kMuL * (vM - vW) / kDx;
			muV_W[idx(i, j)] = muV_X_W;
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
			muV_X_W = muEff * (vM - vW) / kDx
				- muEff * JEff * theta / kMuH;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDx;
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
			muV_X_W = muEff * (vM - vW) / kDx
				+ muEff * JEff * theta / kMuL;
			muV_W[idx(i, j)] = muEff * (vM - vW) / kDx;
			resJ21W[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
			
		if (lsVM >= 0 && lsVE >= 0) {
			muV_X_E = kMuH * (vE - vM) / kDx;
			muV_E[idx(i, j)] = muV_X_E;
		}
		else if (lsVM <= 0 && lsVE <= 0) {
			muV_X_E = kMuL * (vE - vM) / kDx;
			muV_E[idx(i, j)] = muV_X_E;
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
			muV_X_E = muEff * (vE - vM) / kDx
				+ muEff * JEff * theta / kMuL;
			muV_E[idx(i, j)] = muEff * (vE - vM) / kDx;
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
			muV_X_E = muEff * (vE - vM) / kDx
				- muEff * JEff * theta / kMuH;
			muV_E[idx(i, j)] = muEff * (vE - vM) / kDx;
			resJ21W[idx(i, j)] = -muEff * JEff * theta / kMuH;
		}
		
		if (lsVS > 0 && lsVM > 0) {
			muV_Y_S = kMuH * (vM - vS) / kDy;
			muV_S[idx(i, j)] = muV_Y_S;
		}
		else if (lsVS <= 0 && lsVM <= 0) {
			muV_Y_S = kMuL * (vM - vS) / kDy;
			muV_S[idx(i, j)] = muV_Y_S;
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
			muV_Y_S = muEff * (vM - vS) / kDy
				- muEff * JEff * theta / kMuH;
			muV_S[idx(i, j)] = muEff * (vM - vS) / kDy;
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
			muV_Y_S = muEff * (vM - vS) / kDy
				+ muEff * JEff * theta / kMuL;
			muV_S[idx(i, j)] = muEff * (vM - vS) / kDy;
			resJ22S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}
		
		if (lsVM > 0 && lsVN > 0) {
			muV_Y_N = kMuH * (vN - vM) / kDy;
			muV_N[idx(i, j)] = muV_Y_N;
		}
		else if (lsVM <= 0 && lsVN <= 0) {
			muV_Y_N = kMuL * (vN - vM) / kDy;
			muV_N[idx(i, j)] = muV_Y_N;
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
			muV_Y_N = muEff * (vN - vM) / kDy
				+ muEff * JEff * theta / kMuL;
			muV_N[idx(i, j)] = muEff * (vN - vM) / kDy;
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
			muV_Y_N = muEff * (vN - vM) / kDy
				- muEff * JEff * theta / kMuH;
			muV_N[idx(i, j)] = muEff * (vN - vM) / kDy;
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
			muU_Y_W = kMuH * (uM - uS) / kDy;
			muU_S[idx(i, j)] = muU_Y_W;
		}
		else if (lsUS <= 0 && lsUM <= 0) {
			muU_Y_W = kMuL * (uM - uS) / kDy;
			muU_S[idx(i, j)] = muU_Y_W;
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
			muU_Y_W = muEff * (uM - uS) / kDy
				- muEff * JEff * theta / kMuH;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDy;
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
			muU_Y_W = muEff * (uM - uS) / kDy
				+ muEff * JEff * theta / kMuL;
			muU_S[idx(i, j)] = muEff * (uM - uS) / kDy;
			resJ12S[idx(i, j)] = muEff * JEff * theta / kMuL;
		}

		if (lsUE_S > 0 && lsUE > 0) {
			muU_Y_E = kMuH * (uE - uE_S) / kDy;
			muU_N[idx(i, j)] = muU_Y_E;
		}
		else if (lsUE_S <= 0 && lsUE <= 0) {
			muU_Y_E = kMuL * (uE - uE_S) / kDy;
			muU_N[idx(i, j)] = muU_Y_E;
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
			muU_Y_E = muEff * (uE - uE_S) / kDy
				- muEff * JEff * theta / kMuH;
			muU_N[idx(i, j)] = muEff * (uE - uE_S) / kDy;
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
			muU_Y_E = muEff * (uE - uE_S) / kDy
				+ muEff * JEff * theta / kMuL;
			muU_N[idx(i, j)] = muEff * (uE - uE_S) / kDy;
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

		muU_Y_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (uM - uS) / kDy;
		muU_Y_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
			* (uE - uE_S) / kDy;
		muV_X_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
			* (vM - vW) / kDx;
		muV_X_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
			* (vE - vM) / kDx;
		muV_Y_S = mu[idx(i, j - 1)] * (vM - vS) / kDy;
		muV_Y_N = mu[idx(i, j)] * (vN - vM) / kDy;

		muV_W[idx(i, j)] = muV_X_W;
		muV_E[idx(i, j)] = muV_X_E;
		muV_S[idx(i, j)] = muV_Y_S;
		muV_N[idx(i, j)] = muV_Y_N;
		muU_S[idx(i, j)] = muU_Y_W;
		muU_N[idx(i, j)] = muU_Y_E;

		visX = (muU_Y_E - muU_Y_W) / kDx + (muV_X_E - muV_X_W) / kDx;
		visY = 2.0 * (muV_Y_N - muV_Y_S) / kDy;
		checkDiv[idx(i, j)] = (muV_Y_N - muV_Y_S) / kDy + (muU_Y_E - muU_Y_W) / kDx;
		
		// dV[idx(i, j)] = visX / rhoEffWE + visY / rhoEffSN;
		dV[idx(i, j)] = visX * iRhoEffWE + visY * iRhoEffSN;
		// dV[idx(i, j)] = visX + visY;

		iRhoHVec[idx(i, j)] = iRhoEffWE;
		iRhoVVec[idx(i, j)] = iRhoEffSN;
		visXVec[idx(i, j)] = visX * iRhoEffWE;
		visYVec[idx(i, j)] = visY * iRhoEffSN;
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
		outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"mu\", \"dV\", \"muU_S\", \"muU_N\", \"muV_W\", \"muV_E\",\"muV_S\", \"muV_N\", \"visX\", \"visY\" " << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
			<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
			<< static_cast<double>(u[idx(i, j)]) << std::string(",")
			<< static_cast<double>(v[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(mu[idx(i, j)]) << std::string(",")
			<< static_cast<double>(dV[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muU_N[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_W[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_E[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_S[idx(i, j)]) << std::string(",")
			<< static_cast<double>(muV_N[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visXVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(visYVec[idx(i, j)]) << std::string(",")
			<< std::endl;

	outF.close();
	
	return dV;	
}

std::vector<double> MACSolver2D::AddGravityFU() {
	std::vector<double> gU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM::X) {
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

std::vector<double> MACSolver2D::AddGravityFV() {
	std::vector<double> gV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	if ((kFr == 0 || kG == 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis != GAXISENUM::Y) {
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
	const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
	std::vector<double>& uhat, std::vector<double>& vhat, const std::vector<double>& H) {
		
	// Update rhs
	if (kTimeOrder == TIMEORDERENUM::EULER) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
		std::vector<double> FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
		
		for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
			uhat[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
		}
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
			vhat[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
		}
	}
	else if (kTimeOrder == TIMEORDERENUM::RK2) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FU2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
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
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
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
	else if (kTimeOrder == TIMEORDERENUM::RK3) {
		std::vector<double> FU1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FU2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
			FV2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

		// FU1 & FV1 : L(u^(0))
		// FU2 & FV2 : u^(1)
		FU1 = UpdateFU(LSolver, ls, u, v, H);
		FV1 = UpdateFV(LSolver, ls, u, v, H);
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
		FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
		FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
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
		FU2 = UpdateFU(LSolver, ls, FU1, FV1, H);
		FV2 = UpdateFV(LSolver, ls, FU1, FV1, H);
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
	const std::vector<double>& ls, const std::vector<double>& lsB,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& H, const int maxIter) {
	if (!m_Poisson) {
		perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
	}
	std::vector<double> rhs((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> rhoW((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		rhoE((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		rhoS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		rhoN((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> iRhoW((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoE((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoN((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> FWVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		FEVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		FSVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		FNVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> aWVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		aEVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		aCVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		aSVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		aNVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

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
	
	if (kWe != 0.0 && !isnan(kWe) && !isinf(kWe)) {
		UpdateKappa(ls);
		ApplyBC_P_2D(m_kappa);
	}
	// A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
	std::vector<double> AVals, DiagVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	
	std::vector<double> dudX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dudY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> U_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		V_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	double dudXW = 0.0, dudXE = 0.0, dudYS = 0.0, dudYN = 0.0;
	double dvdXW = 0.0, dvdXE = 0.0, dvdYS = 0.0, dvdYN = 0.0;
	// stored coef for A matrix, Dictionary but it is ordered
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
		ps[idx(i, j)] = 0.0;
		rhs[idx(i, j)] = 0.0;
	}
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = lsB[idx(i - 1, j)];
		lsE = lsB[idx(i + 1, j)];
		lsM = lsB[idx(i, j)];
		lsS = lsB[idx(i, j - 1)];
		lsN = lsB[idx(i, j + 1)];

		// At P grid
		dudX[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dudY[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDy);
		dvdX[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dvdY[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDy);

		// Jump occurs when computing dudX
		if (lsW * lsE < 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudXW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dudXE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDx);
			dudX[idx(i, j)] = dudXE * theta + dudXW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdXW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dvdXE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDx);
			dvdX[idx(i, j)] = dvdXE * theta + dvdXW * (1.0 - theta);
		}

		if (lsS * lsN < 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudYS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dudYN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDy);
			dudY[idx(i, j)] = dudYN * theta + dudYS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdYS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dvdYN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDy);
			dvdY[idx(i, j)] = dvdYN * theta + dvdYS * (1.0 - theta);
		}
		
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
		aEff = 0.0;
		
		// [a]_\Gamma = a^+ - a^- = a^inside - a^outside
		// [p^*] = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// p_M - p_W = 2 dt[mu](\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa
		// (P_M - P_W)/kDx appears
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
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoW[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FW = aEff * iRhoW[idx(i, j)] / (kDx * kDx);
		}
		else if (lsW <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| === outside(-) === |(interface)| ===     inside(+)      === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d    === |(lsM)|
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoW[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FW = -aEff * iRhoW[idx(i, j)] / (kDx * kDx);
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
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
				+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = aM * theta + aE * (1.0 - theta);
			iRhoE[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FE = aEff * iRhoE[idx(i, j)] / (kDx * kDx);
		}
		else if (lsE <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i + 1, j]
			theta = std::fabs(lsE) / (std::fabs(lsE) + std::fabs(lsM));
			// |(lsM)| ===   inside(+)     === |(interface)| === outside(-)  === |(lsE)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d   === |(lsE)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aE = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i + 1, j)] * nXE + dudY[idx(i + 1, j)] * nYE) * nXE
				+ (dvdX[idx(i + 1, j)] * nXE + dvdY[idx(i + 1, j)] * nYE) * nYE)
				+ m_dt * kSigma * m_kappa[idx(i + 1, j)];
			aEff = aM * theta + aE * (1.0 - theta);
			iRhoE[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FE = -aEff * iRhoE[idx(i, j)] / (kDx * kDx);
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
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			// rhoS[idx(i, j)] = kRhoH * theta + kRhoL * (1.0 - theta);
			iRhoS[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			// FS = aEff / (rhoS[idx(i, j)] * kDy * kDy);
			FS = aEff * iRhoS[idx(i, j)] / (kDy * kDy);
		}
		else if (lsS <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i, j - 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===     inside(+)   === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			// rhoS[idx(i, j)] = kRhoL * theta + kRhoH * (1.0 - theta);
			iRhoS[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			// FS = -aEff / (rhoS[idx(i, j)] * kDy * kDy);
			FS = -aEff * iRhoS[idx(i, j)] / (kDy * kDy);
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
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nYN) * nXN
				+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = aM * theta + aN * (1.0 - theta);
			iRhoN[idx(i, j)] = -1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			FN = aEff * iRhoN[idx(i, j)] / (kDy * kDy);
		}
		else if (lsN <= 0.0 && lsM > 0.0) {
			// interface lies between ls[i, j] and ls[i, j + 1]
			theta = std::fabs(lsN) / (std::fabs(lsN) + std::fabs(lsM));
			// |(lsM)| ===    inside       === |(interface)| ===   outside === |(lsN)|
			// |(lsM)| === (1 - theta) * d === |(interface)| === theta * d === |(lsN)|
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aN = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j + 1)] * nXN + dudY[idx(i, j + 1)] * nYN) * nXN
				+ (dvdX[idx(i, j + 1)] * nXN + dvdY[idx(i, j + 1)] * nYN) * nYN)
				+ m_dt * kSigma * m_kappa[idx(i, j + 1)];
			aEff = aM * theta + aN * (1.0 - theta);
			iRhoN[idx(i, j)] = -1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			FN = -aEff * iRhoN[idx(i, j)] / (kDy * kDy);
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
		// coefficient
		iRhoW[idx(i, j)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0;
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FW = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
		FW *= iRhoW[idx(i, j)] / (kDx * kDx);
	
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
		iRhoE[idx(i, j)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i + 1, j)] * (1.0 - theta);

		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FE = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);
		FE *= iRhoE[idx(i, j)] / (kDx * kDx);
		
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
		iRhoS[idx(i, j)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
		
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FS = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
		FS *= iRhoS[idx(i, j)] / (kDy * kDy);
		
		// if (FS != 0.0)
		// 	std::cout << "FS " << i << " " << j << " " << FS << " " << FN << " " 
		// 	<< m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)])
		// 	<< " " << iRhoS[idx(i, j)]
		// 	<< " " << (H[idx(i, j)] - H[idx(i, j - 1)])
		// 	<< " " << m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]) << std::endl;
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
		iRhoN[idx(i, j)] = 1.0 / rhoEff;
		kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j + 1)] * (1.0 - theta);
		// pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
		// kRhoH - Liquid, kRhoL - Gas
		// - level set : H = 1, + level set : H = 0; 
		// Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
		FN = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1)] - H[idx(i, j)]);
		FN *= iRhoN[idx(i, j)] / (kDy * kDy);
		
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
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		rhs[idx(i, j)] -= div[idx(i, j)];

	// An order of A matrix coef. is very important, hence reverse j order
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

		iRhoS[idx(i, j)] *= -1.0;
		iRhoW[idx(i, j)] *= -1.0;
		iRhoE[idx(i, j)] *= -1.0;
		iRhoN[idx(i, j)] *= -1.0;

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["S"] = iRhoS[idx(i, j)] / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
		AValsDic["W"] = iRhoW[idx(i, j)] / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["C"] = -(iRhoW[idx(i, j)] + iRhoE[idx(i, j)]) / (kDx * kDx)
			- (iRhoS[idx(i, j)] + iRhoN[idx(i, j)]) / (kDy * kDy);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["E"] = iRhoE[idx(i, j)] / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["N"] = iRhoN[idx(i, j)] / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;
		
		if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += iRhoW[idx(i, j)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= iRhoW[idx(i, j)] / (kDx * kDx);
			rhs[idx(i, j)] -= iRhoW[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPW);
		}
		else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::PERIODIC) {
			AValsDic["W"] = iRhoW[idx(kNumBCGrid + kNx - 1, j)];
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC2D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += iRhoE[idx(i, j)] / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC2D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= iRhoE[idx(i, j)] / (kDx * kDx);
			rhs[idx(i, j)] -= iRhoE[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPE);
		}
		else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC2D::PERIODIC) {
			AValsDic["E"] = iRhoE[idx(kNumBCGrid, j)];
		}

		if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += iRhoS[idx(i, j)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= iRhoS[idx(i, j)] / (kDy * kDy);
			rhs[idx(i, j)] -= iRhoS[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPS);
		}
		else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::PERIODIC) {
			AValsDic["S"] = iRhoS[idx(i, kNumBCGrid + kNy - 1)];
		}

		if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC2D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += iRhoN[idx(i, j)] / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC2D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= iRhoN[idx(i, j)] / (kDy * kDy);
			rhs[idx(i, j)] -= iRhoN[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPN);
		}
		else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC2D::PERIODIC) {
			AValsDic["N"] = iRhoN[idx(i, kNumBCGrid + kNy + 1)];
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
		m_Poisson->CG_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);
	
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
		// std::cout << "Poisson : BiCG" << std::endl;
		m_Poisson->BiCGStab_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
			DiagVals, DiagCols, DiagRowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	else if (m_PoissonSolverType == POISSONTYPE::GS) {
		// std::cout << "Poisson : GS" << std::endl;
		m_Poisson->GS_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
				std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
	}
	
	ApplyBC_P_2D(ps);

	std::ofstream outF;
	std::string fname("P_ASCII.plt");
	if (m_iter == 1) {
		outF.open(fname.c_str(), std::ios::out);

		outF << "TITLE = VEL" << std::endl;
		outF << "VARIABLES = \"X\", \"Y\", \"FW\", \"FE\", \"FS\", \"FN\", \"F\",\"RW\", \"RE\", \"RS\", \"RN\", \"H\", \"LS\", \"PS\", \"HatDIV\", \"RHS\", \"P_W\", \"P_S\",\"KAPPA\"" << std::endl;
		outF.close();
	}

	outF.open(fname.c_str(), std::ios::app);
	
	outF << std::string("ZONE T=\"") << m_iter
		<< std::string("\", I=") << kNx + 6 << std::string(", J=") << kNy + 6
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;

	// for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	// 	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
				for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
			<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FEVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(FWVec[idx(i, j)] + FEVec[idx(i, j)] + FSVec[idx(i, j)] + FNVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoW[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoE[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoS[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoN[idx(i, j)]) << std::string(",")
			<< static_cast<double>(H[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(div[idx(i, j)]) << std::string(",")
			<< static_cast<double>(rhs[idx(i, j)]) << std::string(",")
			<< static_cast<double>((ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx) << std::string(",")
			<< static_cast<double>((ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;

	outF.close();
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
	std::vector<double> dudX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dudY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dvdY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		dldY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> U_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		V_PGrid((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	double dudXW = 0.0, dudXE = 0.0, dudYS = 0.0, dudYN = 0.0;
	double dvdXW = 0.0, dvdXE = 0.0, dvdYS = 0.0, dvdYN = 0.0;

	std::vector<double> nxVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		nyVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> AWW((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		AWM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		ASS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		ASM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> iRhoEffWVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		iRhoEffSVec((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> aWEff((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		aSEff((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		U_PGrid[idx(i, j)] = 0.5 * (u[idx(i, j)] + u[idx(i + 1, j)]);
		V_PGrid[idx(i, j)] = 0.5 * (v[idx(i, j)] + v[idx(i, j + 1)]);
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsE = ls[idx(i + 1, j)];
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		lsN = ls[idx(i, j + 1)];
		
		// At P grid 
		dudX[idx(i, j)] = (U_PGrid[idx(i + 1, j)] - U_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dudY[idx(i, j)] = (U_PGrid[idx(i, j + 1)] - U_PGrid[idx(i, j - 1)]) / (2.0 * kDy);
		dvdX[idx(i, j)] = (V_PGrid[idx(i + 1, j)] - V_PGrid[idx(i - 1, j)]) / (2.0 * kDx);
		dvdY[idx(i, j)] = (V_PGrid[idx(i, j + 1)] - V_PGrid[idx(i, j - 1)]) / (2.0 * kDy);

		// Jump occurs when computing dudX
		if (lsW * lsE <= 0) {
			// interface lies between u[i - 1, j] and u[i + 1, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsE));
			// |(lsW)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsE)|
			// |(lsW)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsE)|
			uEff = U_PGrid[idx(i + 1, j)] * theta + U_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dudXW = (uEff - U_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dudXE = (U_PGrid[idx(i + 1, j)] - uEff) / (2.0 * (1.0 - theta) * kDx);
			dudX[idx(i, j)] = dudXE * theta + dudXW * (1.0 - theta);

			vEff = V_PGrid[idx(i + 1, j)] * theta + V_PGrid[idx(i - 1, j)] * (1.0 - theta);
			dvdXW = (vEff - V_PGrid[idx(i - 1, j)]) / (2.0 * theta * kDx);
			dvdXE = (V_PGrid[idx(i + 1, j)] - vEff) / (2.0 * (1.0 - theta) * kDx);
			dvdX[idx(i, j)] = dvdXE * theta + dvdXW * (1.0 - theta);
		}

		if (lsS * lsN <= 0) {
			// interface lies between ls[i, j - 1] and u[i, j + 1]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsN));
			// |(lsS)| ===  inside(+) === |(interface)| ===    outside(-)    === |(lsN)|
			// |(lsS)| === theta  * d === |(interface)| === (1 - theta) * d  === |(lsN)|
			uEff = U_PGrid[idx(i, j + 1)] * theta + U_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dudYS = (uEff - U_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dudYN = (U_PGrid[idx(i, j + 1)] - uEff) / (2.0 * (1.0 - theta) * kDy);
			dudY[idx(i, j)] = dudYN * theta + dudYS * (1.0 - theta);

			vEff = V_PGrid[idx(i, j + 1)] * theta + V_PGrid[idx(i, j - 1)] * (1.0 - theta);
			dvdYS = (vEff - V_PGrid[idx(i, j - 1)]) / (2.0 * theta * kDy);
			dvdYN = (V_PGrid[idx(i, j + 1)] - vEff) / (2.0 * (1.0 - theta) * kDy);
			dvdY[idx(i, j)] = dvdYN * theta + dvdYS * (1.0 - theta);
		}
		
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
	
	/*
	// http://ctr.stanford.edu/Summer/SP08/3_1_Moureau.pdf
	// Moreau, V., and O. Desjardins. 
	// "A second-order ghost-fluid method for the primary atomization
	//	 of liquid fuel in air-blast type injectors."
	// Proceedings of the Summer Program. Vol. 143. 2008.
	*/
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		lsW = ls[idx(i - 1, j)];
		lsM = ls[idx(i, j)];
		aW = 0.0;
		aM = 0.0;
		
		nXW = dldX[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps * eps);
		nYW = dldY[idx(i - 1, j)] / (std::sqrt(std::pow(dldX[idx(i - 1, j)], 2.0) + std::pow(dldY[idx(i - 1, j)], 2.0)) + eps * eps);
		nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps * eps);
		nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps * eps);
		nxVec[idx(i, j)] = nXM;
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
		u[idx(i, j)] *= iRhoEff / kDx;

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
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			u[idx(i, j)] = -aEff * iRhoEff / kDx;
		}
		else if (lsW <= 0 && lsM > 0) {
			// interface lies between ls[i - 1, j] and ls[i, j]
			theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
			// |(lsW)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aW = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i - 1, j)] * nXW + dudY[idx(i - 1, j)] * nYW) * nXW
				+ (dvdX[idx(i - 1, j)] * nXW + dvdY[idx(i - 1, j)] * nYW) * nYW)
				+ m_dt * kSigma * m_kappa[idx(i - 1, j)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aW * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			u[idx(i, j)] = aEff * iRhoEff / kDx;
		}

		iRhoEffWVec[idx(i, j)] = iRhoEff;
		aWEff[idx(i, j)] = u[idx(i, j)];
		*/
		u[idx(i, j)] += us[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
		lsM = ls[idx(i, j)];
		lsS = ls[idx(i, j - 1)];
		aS = 0.0;
		aM = 0.0;

		nXS = dldX[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps * eps);
		nYS = dldY[idx(i, j - 1)] / (std::sqrt(std::pow(dldX[idx(i, j - 1)], 2.0) + std::pow(dldY[idx(i, j - 1)], 2.0)) + eps * eps);
		nXM = dldX[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps * eps);
		nYM = dldY[idx(i, j)] / (std::sqrt(std::pow(dldX[idx(i, j)], 2.0) + std::pow(dldY[idx(i, j)], 2.0)) + eps * eps);
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
		v[idx(i, j)] *= iRhoEff / kDy;

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
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
			v[idx(i, j)] = -aEff * iRhoEff / kDy;
		}
		else if (lsS <= 0 && lsM > 0) {
			// interface lies between ls[i, j - 1] and ls[i, j]
			theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
			// |(lsS)| === outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
			// |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
			aS = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j - 1)] * nXS + dudY[idx(i, j - 1)] * nYS) * nXS
				+ (dvdX[idx(i, j - 1)] * nXS + dvdY[idx(i, j - 1)] * nYS) * nYS)
				+ m_dt * kSigma * m_kappa[idx(i, j - 1)];
			aM = 2.0 * m_dt * (kMuH - kMuL)
				* ((dudX[idx(i, j)] * nXM + dudY[idx(i, j)] * nYM) * nXM
				+ (dvdX[idx(i, j)] * nXM + dvdY[idx(i, j)] * nYM) * nYM)
				+ m_dt * kSigma * m_kappa[idx(i, j)];
			aEff = aM * theta + aS * (1.0 - theta);
			iRhoEff = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
			v[idx(i, j)] = aEff * iRhoEff / kDy;
		}

		aSEff[idx(i, j)] = v[idx(i, j)];
		iRhoEffSVec[idx(i, j)] = iRhoEff;
		*/
		v[idx(i, j)] += vs[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy;
	}
	
	/*
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
		<< std::string("\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", SOLUTIONTIME=") << m_iter * 0.1
		<< std::string(", STRANDID=") << m_iter + 1
		<< std::endl;
	
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
			outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
			<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
			<< static_cast<double>((u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5) << std::string(",")
			<< static_cast<double>((v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5) << std::string(",")
			<< static_cast<double>((us[idx(i, j)] + us[idx(i + 1, j)]) * 0.5) << std::string(",")
			<< static_cast<double>((vs[idx(i, j)] + vs[idx(i, j + 1)]) * 0.5) << std::string(",")
			<< static_cast<double>(ls[idx(i, j)]) << std::string(",")
			<< static_cast<double>(ps[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoEffWVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aWEff[idx(i, j)]) << std::string(",")
			<< static_cast<double>(-(ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx) << std::string(",")
			<< static_cast<double>(-iRhoEffWVec[idx(i, j)] * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx) << std::string(",")
			<< static_cast<double>(u[idx(i, j)] - us[idx(i, j)]) << std::string(",")
			<< static_cast<double>(iRhoEffSVec[idx(i, j)]) << std::string(",")
			<< static_cast<double>(aSEff[idx(i, j)]) << std::string(",")
			<< static_cast<double>(-(ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy) << std::string(",")
			<< static_cast<double>(-iRhoEffSVec[idx(i, j)] * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy) << std::string(",")
			<< static_cast<double>(v[idx(i, j)] - vs[idx(i, j)]) << std::string(",")
			<< static_cast<double>(u[idx(i + 1, j)] - u[idx(i, j)]) / kDx << std::string(",")
			<< static_cast<double>(v[idx(i, j + 1)] - v[idx(i, j)]) / kDy << std::string(",")
			<< static_cast<double>((u[idx(i + 1, j)] - u[idx(i, j)]) / kDx + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy) << std::string(",")
			<< static_cast<double>(m_kappa[idx(i, j)]) << std::endl;

	outF.close();
	*/
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
	Vefl = std::max(kMuH / kRhoH, kMuL / kRhoL) * (2.0 / (kDx * kDx) + 2.0 / (kDy * kDy));
	Gefl = std::max(std::sqrt(std::fabs(kG) / kDy), std::sqrt(std::fabs(kG) / kDx));
	
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

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
inline int MACSolver2D::sign(const double& val) {
	return (val > 0) - (val < 0);
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
				// << static_cast<double>(m_kappa[idx(i, j)]) << std::string(",")
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
			resDiv(kNx * kNy, 0.0),
			resLS(kNx * kNy, 0.0),
			resPs(kNx * kNy, 0.0);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)	 {
			resX[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseX + (i + 0.5 - kNumBCGrid) * kDx;
			resY[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseY + (j + 0.5 - kNumBCGrid) * kDy;
			resU[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
			resV[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
			resLS[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ls[idx(i, j)];
			resPs[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ps[idx(i, j)];
			resDiv[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
			 + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
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
