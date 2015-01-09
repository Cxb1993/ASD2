#include "levelset2d.h"

LevelSetSolver2D::LevelSetSolver2D(int nx, int ny, int num_bc_grid, double dx, double dy) :
	kEps(std::min(dx, dy)),
	kThickness(2.0 * kEps),
	kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid),
	kDx(dx), kDy(dy),
	m_adt(kEps * 0.5), m_atime(0.0)  {
}

int LevelSetSolver2D::Solve_LevelSet_2D(
	std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v) {

	// Level set is a P grid, interpolate u and v to P grid
	std::vector<double> tmpU((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> tmpV((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		tmpU[idx(i, j)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
		tmpV[idx(i, j)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
	}
		
	// TVD RK2
	std::vector<double> tmpLS1 = HJWENO5_LS_2D(ls, tmpU, tmpV);
	std::vector<double> tmpLS2 = HJWENO5_LS_2D(tmpLS1, tmpU, tmpV);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		ls[idx(i, j)] = (tmpLS1[idx(i, j)] + tmpLS2[idx(i, j)]) * 0.5;
	}

	return 0;
}

std::vector<double> LevelSetSolver2D::HJWENO5_LS_2D(std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v) {

	std::vector<double> LS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LXM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYP((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> LYM((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> vecF_LSX(kNx + 2 * kNumBCGrid, 0.0), vecF_LSY(kNy + 2 * kNumBCGrid);

	// U : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_LSX[i] = ls[idx(i, j)];
		}

		UnitHJWENO5(vecF_LSX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j)] = FXP[i];
			LXM[idx(i, j)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_LSX.begin(), vecF_LSX.end(), 0.0);
	}

	// U : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_LSY[i] = ls[idx(i, j)];
		}

		UnitHJWENO5(vecF_LSY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; i++) {
			LYP[idx(i, j)] = FYP[i];
			LYM[idx(i, j)] = FYM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_LSY.begin(), vecF_LSY.end(), 0.0);
	}

// combine together with Local Lax-Friedrichs Scheme
double alphaX = 0.0, alphaY = 0.0;

for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		alphaX = u[idx(i, j)];
		alphaY = v[idx(i, j)];

		ls[idx(i, j)]
			= 0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
			+ 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
			- alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
			- alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)]));
	}

return ls;
}

int LevelSetSolver2D::Reinit_Sussman_2D(std::vector<double>& ls) {
	// using 2nd order HJ ENO, based on Sussman's work
	// Sussman, Mark, et al. "An improved level set method for incompressible two-phase flows."
	// Computers & Fluids 27.5 (1998): 663-680.

	// HJENO 2nd order + RK2
	std::vector<double> dLS1 = HJENO_Abs_2D(ls);
	std::vector<double> ls2((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		ls2[idx(i, j)] = ls[idx(i, j)] + m_adt * dLS1[idx(i, j)];
	}
	
	std::vector<double> dLS2 = HJENO_Abs_2D(ls2);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		ls[idx(i, j)] = ls[idx(i, j)] + 0.5 * m_adt * (dLS1[idx(i, j)] + dLS2[idx(i, j)]);
	}

	return 0;
}

std::vector<double>
	LevelSetSolver2D::HJENO_Abs_2D(std::vector<double>& ls) {

	int k1 = 0, kEff = 0;
	double a = 0.0, b = 0.0, c = 0.0;

	// Apply Boundary Condition First
	ApplyBC_P_2D(ls);

	// Divided Difference Table
	std::vector<double> Q0X_k1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		Q1X((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		Q2X((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<double> Q0Y_k1((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		Q1Y((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		Q2Y((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	std::vector<int> signedLS = GetSignedLSNormalized(ls);

	std::vector<double> absdLS((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> L((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<int> kMinVecX((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0),
		kMinVecY((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);

	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++)
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++) {
		Q1X[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
		Q1Y[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		a = (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / kDx;
		b = (ls[idx(i, j)] - 2.0 * ls[idx(i + 1, j)] + ls[idx(i + 2, j)]) / kDx;
		if (fabs(a) <= fabs(b)) {
			c = a;
			k1 = i - 1;
		}
		else {
			c = b;
			k1 = i;
		}

		kMinVecX[idx(i, j)] = k1;
		Q2X[idx(i, j)] = Q1X[idx(i, j)] - 0.5 * kDx * c * (2.0 * (k1 - i) + 1);
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		a = (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / kDx;
		b = (ls[idx(i, j)] - 2.0 * ls[idx(i, j + 1)] + ls[idx(i, j + 2)]) / kDx;
		if (fabs(a) <= fabs(b)) {
			c = a;
			k1 = j - 1;
		}
		else {
			c = b;
			k1 = j;
		}

		kMinVecY[idx(i, j)] = k1;
		Q2Y[idx(i, j)] = Q1Y[idx(i, j)] - 0.5 * kDy * c * (2.0 * (k1 - i) + 1);
	}

	double dPX = 0.0, dMX = 0.0, dPY = 0.0, dMY = 0.0;
	double dLSdX = 0.0, dLSdY = 0.0;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		dPX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
		dMX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
		if (kMinVecX[idx(i, j)] == i) {
			dPX = Q2X[idx(i, j)];
		}
		else {
			dMX = Q2X[idx(i, j)];
		}

		if (dPX * signedLS[idx(i, j)] < 0 && dMX * signedLS[idx(i, j)] < -dPX * signedLS[idx(i, j)]) {
			dLSdX = dPX;
		}
		else if (dPX * signedLS[idx(i, j)] < 0 && dPX * signedLS[idx(i, j)] < -dMX * signedLS[idx(i, j)]) {
			dLSdX = dMX;
		}
		else {
			dLSdX = (dPX + dPX) * 0.5;
		}

		dPY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDx;
		dMY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDx;
		if (kMinVecX[idx(i, j)] == j) {
			dPY = Q2Y[idx(i, j)];
		}
		else {
			dMY = Q2Y[idx(i, j)];
		}

		if (dPY * signedLS[idx(i, j)] < 0 && dMY * signedLS[idx(i, j)] < -dPY * signedLS[idx(i, j)]) {
			dLSdY = dPY;
		}
		else if (dPY * signedLS[idx(i, j)] < 0 && dPY * signedLS[idx(i, j)] < -dMY * signedLS[idx(i, j)]) {
			dLSdY = dMY;
		}
		else {
			dLSdY = (dPY + dPY) * 0.5;
		}

		absdLS[idx(i, j)] = std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY);
	}

	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++)
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++) {
		L[idx(i, j)] = m_signedInitLS[idx(i, j)] * (1.0 - absdLS[idx(i, j)]);
	}

	return L;
}

int LevelSetSolver2D::Reinit_Original_2D(std::vector<double>& ls) {

	return 0;
}

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
int LevelSetSolver2D::sign(const double& val) {
	return (val > 0) - (val < 0);
}

std::vector<int> LevelSetSolver2D::GetSignedLSNormalized(const std::vector<double> &v) {
	std::vector<int> r(v.size());
	// std::transform(v.begin(), v.end(), r.begin(), (static_cast<int>((const int&)sign)));

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		r[idx(i, j)] = sign(v[idx(i, j)]);
	}

	return r;
}

int LevelSetSolver2D::UnitHJWENO5(
	const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n) {
	int stat = 0;
	const double kWeightEps = 1e-6;

	double alphaPlus[3], alphaMinus[3];
	double weightPlus[3], weightMinus[3];

	std::vector<double> IS0(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS1(n + 2 * kNumBCGrid, 0.0);
	std::vector<double> IS2(n + 2 * kNumBCGrid, 0.0);

	// \dfrac{\Delta^+ F}{\Delta x}
	std::vector<double> DFPlus = std::vector<double>(n + 2 * kNumBCGrid, 0.0);

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


int LevelSetSolver2D::SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int LevelSetSolver2D::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int LevelSetSolver2D::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	}

	m_BC->SetBC_P_2D(BC_W, BC_E, BC_S, BC_N);

	return 0;
}

int LevelSetSolver2D::ApplyBC_U_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_U_2D(arr);
	
	return 0;
}

int LevelSetSolver2D::ApplyBC_V_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_V_2D(arr);

	return 0;
}

int LevelSetSolver2D::ApplyBC_P_2D(std::vector<double>& arr) {
	m_BC->ApplyBC_P_2D(arr);
	
	return 0;
}

void LevelSetSolver2D::SetBCConstantUW(double BC_ConstantW) {
	return m_BC->SetBCConstantUW(BC_ConstantW);
}

void LevelSetSolver2D::SetBCConstantUE(double BC_ConstantE) {
	return m_BC->SetBCConstantUE(BC_ConstantE);
}

void LevelSetSolver2D::SetBCConstantUS(double BC_ConstantS) {
	return m_BC->SetBCConstantUS(BC_ConstantS);
}

void LevelSetSolver2D::SetBCConstantUN(double BC_ConstantN) {
	return m_BC->SetBCConstantUN(BC_ConstantN);
}

void LevelSetSolver2D::SetBCConstantVW(double BC_ConstantW) {
	return m_BC->SetBCConstantVW(BC_ConstantW);
}

void LevelSetSolver2D::SetBCConstantVE(double BC_ConstantE) {
	return m_BC->SetBCConstantVE(BC_ConstantE);
}

void LevelSetSolver2D::SetBCConstantVS(double BC_ConstantS) {
	return m_BC->SetBCConstantVS(BC_ConstantS);
}

void LevelSetSolver2D::SetBCConstantVN(double BC_ConstantN) {
	return m_BC->SetBCConstantVN(BC_ConstantN);
}

void LevelSetSolver2D::SetBCConstantPW(double BC_ConstantW) {
	return m_BC->SetBCConstantPW(BC_ConstantW);
}

void LevelSetSolver2D::SetBCConstantPE(double BC_ConstantE) {
	return m_BC->SetBCConstantPE(BC_ConstantE);
}

void LevelSetSolver2D::SetBCConstantPS(double BC_ConstantS) {
	return m_BC->SetBCConstantPS(BC_ConstantS);
}

void LevelSetSolver2D::SetBCConstantPN(double BC_ConstantN) {
	return m_BC->SetBCConstantPN(BC_ConstantN);
}

double LevelSetSolver2D::minmod(double a, double b) {
	double sign = 0.0;

	if (a * b > 0)
		sign = fabs(a) / a;
	else if (a * b == 0 && a != b)
		sign = 0;
	else
		return 0;

	return sign * std::min(fabs(a), fabs(b));
}

inline int LevelSetSolver2D::idx(int i, int j) {
	return i + (kNx + 2 * kNumBCGrid) * j;
}