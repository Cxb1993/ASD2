#include "levelset3d.h"

LevelSetSolver3D::LevelSetSolver3D(int nx, int ny, int nz, int num_bc_grid, double baseX, double baseY, double baseZ, double dx, double dy, double dz) :
	kEps(std::min(std::min(dx, dy), dz)),
	kThickness(2.0 * kEps),
	kNx(nx), kNy(ny), kNz(nz), kNumBCGrid(num_bc_grid),  
	kArrSize(kArrSize * (kNz + 2 * kNumBCGrid)),
	kBaseX(baseX), kBaseY(baseY), kBaseZ(baseZ),
	kDx(dx), kDy(dy), kDz(dz),
	kAdt(kEps * 0.5), m_atime(0.0), kMaxATime(1.5 * kThickness),
	kENOSpatialOrder(2)  {
}

LevelSetSolver3D::LevelSetSolver3D(int nx, int ny, int nz, int num_bc_grid, double baseX, double baseY, double baseZ, double dx, double dy, double dz, double maxTime) :
	kEps(std::min(std::min(dx, dy), dz)),
	kThickness(2.0 * kEps),
	kNx(nx), kNy(ny), kNz(nz), kNumBCGrid(num_bc_grid),
	kArrSize(kArrSize * (kNz + 2 * kNumBCGrid)),
	kDx(dx), kDy(dy), kDz(dz),
	kBaseX(baseX), kBaseY(baseY), kBaseZ(baseZ),
	kAdt(kEps * 0.5), m_atime(0.0), kMaxATime(maxTime),
	kENOSpatialOrder(2)  {
}

std::vector<double> LevelSetSolver3D::UpdateKappa(const std::vector<double>& ls) {
	std::vector<double> dLSdX(kArrSize, 0.0), dLSdY(kArrSize, 0.0), dLSdZ(kArrSize, 0.0);
	std::vector<double> LSSize(kArrSize, 0.0);
	std::vector<double> kappa(kArrSize, 0.0);

	const double eps = 1.0e-100;
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; j < kNz + kNumBCGrid; k++) {
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
	for (int k = kNumBCGrid; j < kNz + kNumBCGrid; k++) {
		kappa[idx(i, j, k)]
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
		kappa[idx(i, j, k)] = std::fabs(std::min(std::fabs(kappa[idx(i, j, k)]), 1.0 / kEps));

		assert(kappa[idx(i, j, k)] == kappa[idx(i, j, k)]);
	}

	return kappa;
}

int LevelSetSolver3D::Solve_LevelSet_3D(
	std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, double dt) {

	// Level set is a P grid, interpolate u and v to P grid
	std::vector<double> tmpU(kArrSize, 0.0), tmpV(kArrSize, 0.0), tmpW(kArrSize, 0.0);

	for (int i = kNumBCGrid - 1; i < kNx + kNumBCGrid + 1; i++)
	for (int j = kNumBCGrid - 1; j < kNy + kNumBCGrid + 1; j++)
	for (int k = kNumBCGrid - 1; k < kNz + kNumBCGrid + 1; k++) {
		tmpU[idx(i, j, k)] = (u[idx(i, j, k)] + u[idx(i + 1, j, k)]) * 0.5;
		tmpV[idx(i, j, k)] = (v[idx(i, j, k)] + v[idx(i, j + 1, k)]) * 0.5;
		tmpW[idx(i, j, k)] = (w[idx(i, j, k)] + w[idx(i, j, k + 1)]) * 0.5;
	}
	
	// TVD RK3
	std::vector<double> DLS1 = HJWENO5_LS_3D(ls, tmpU, tmpV, tmpW);
	
	std::vector<double> LS1(kArrSize, 0.0);
	std::vector<double> LS2(kArrSize, 0.0);
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		LS1[idx(i, j, k)] = ls[idx(i, j, k)] - dt * DLS1[idx(i, j, k)];
	}
	ApplyBC_P_3D(LS1);
	std::vector<double> DLS2 = HJWENO5_LS_3D(LS1, tmpU, tmpV, tmpW);
	
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		LS2[idx(i, j, k)] = 0.75 * ls[idx(i, j, k)] + 0.25 * LS1[idx(i, j, k)] - 0.25 * dt * DLS2[idx(i, j, k)];
	}
	ApplyBC_P_3D(LS2);
	std::vector<double> DLS3 = HJWENO5_LS_3D(LS2, tmpU, tmpV, tmpW);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		ls[idx(i, j, k)] = 1.0 / 3.0 * ls[idx(i, j, k)] + 2.0 / 3.0 * LS2[idx(i, j, k)] - 2.0 / 3.0 * dt * DLS3[idx(i, j, k)];
	}
	ApplyBC_P_3D(ls);
	
	return 0;
}

std::vector<double> LevelSetSolver3D::HJWENO5_LS_3D(std::vector<double>& ls,
	const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w) {

	std::vector<double> L(kArrSize, 0.0);
	std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
	std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
	std::vector<double> LZP(kArrSize, 0.0), LZM(kArrSize, 0.0);
	
	std::vector<double> vecF_LSX(kNx + 2 * kNumBCGrid, 0.0), vecF_LSY(kNy + 2 * kNumBCGrid, 0.0), vecF_LSZ(kNy + 2 * kNumBCGrid, 0.0);

	ApplyBC_P_3D(ls);

	// U : X direction
	std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; j < kNz + kNumBCGrid; k++) {
		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			vecF_LSX[i] = ls[idx(i, j, k)] * u[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_LSX, FXP, FXM, kDx, kNx);

		for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
			LXP[idx(i, j, k)] = FXP[i];
			LXM[idx(i, j, k)] = FXM[i];
		}

		// set all vector elements to zero keeping its size
		std::fill(FXP.begin(), FXP.end(), 0.0);
		std::fill(FXM.begin(), FXM.end(), 0.0);
		std::fill(vecF_LSX.begin(), vecF_LSX.end(), 0.0);
	}

	// U : Y direction	
	std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			vecF_LSY[j] = ls[idx(i, j, k)] * v[idx(i, j, k)];
		}

		UnitHJWENO5(vecF_LSY, FYP, FYM, kDy, kNy);

		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
			LYP[idx(i, j, k)] = FYP[j];
			LYM[idx(i, j, k)] = FYM[j];
		}

		// set all vector elements to zero keeping its size
		std::fill(FYP.begin(), FYP.end(), 0.0);
		std::fill(FYM.begin(), FYM.end(), 0.0);
		std::fill(vecF_LSY.begin(), vecF_LSY.end(), 0.0);
	}

	// Z direction	
	std::vector<double> FZP(kNz + 2 * kNumBCGrid, 0.0), FZM(kNz + 2 * kNumBCGrid, 0.0);
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			vecF_LSZ[j] = ls[idx(i, j, k)] * w[idx(i, j, k)];
		}
		
		UnitHJWENO5(vecF_LSZ, FZP, FZM, kDz, kNz);

		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++) {
			LZP[idx(i, j, k)] = FYP[k];
			LZM[idx(i, j, k)] = FYM[k];
		}

		// set all vector elements to zero keeping its size
		std::fill(FZP.begin(), FZP.end(), 0.0);
		std::fill(FZM.begin(), FZM.end(), 0.0);
		std::fill(vecF_LSZ.begin(), vecF_LSZ.end(), 0.0);
	}

	// combine together with Local Lax-Friedrichs Scheme
	double alphaX = 0.0, alphaY = 0.0, alphaZ = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		alphaX = u[idx(i, j, k)];
		alphaY = v[idx(i, j, k)];
		alphaY = w[idx(i, j, k)];

		L[idx(i, j, k)]
			= 0.5 * (LXP[idx(i, j, k)] + LXM[idx(i, j, k)])
			+ 0.5 * (LYP[idx(i, j, k)] + LYM[idx(i, j, k)])
			+ 0.5 * (LZP[idx(i, j, k)] + LZM[idx(i, j, k)])
			- alphaX * (0.5 * (LXP[idx(i, j, k)] - LXM[idx(i, j, k)]))
			- alphaY * (0.5 * (LYP[idx(i, j, k)] - LYM[idx(i, j, k)]))
			- alphaZ * (0.5 * (LZP[idx(i, j, k)] - LZM[idx(i, j, k)]));
	}

	return L;
}

int LevelSetSolver3D::UnitHJWENO5(
	const std::vector<double>& F, std::vector<double>& FP, std::vector<double>& FM, const double d, const int n) {
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

int LevelSetSolver3D::Reinit_Original_3D(std::vector<double>& ls) {
	m_atime = 0.0;

	std::vector<double> lsInit(kArrSize, 0.0);
	lsInit = ls;

	std::vector<double> absdLS1(kArrSize, 0.0);
	std::vector<double> absdLS2(kArrSize, 0.0);
	std::vector<double> LS1(kArrSize, 0.0);
	std::vector<double> LS2(kArrSize, 0.0);

	while (m_atime < kMaxATime) {
		m_atime += kAdt;
		
		// Apply Boundary Condition First
		ApplyBC_P_3D(ls);

		absdLS1 = ENO_DerivAbsLS_3D(ls, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			LS1[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt 
				* lsInit[idx(i, j, k)] / std::sqrt(std::pow(lsInit[idx(i, j, k)], 2.0) + kEps * kEps)
					* (absdLS1[idx(i, j, k)] - 1.0);
		}

		absdLS2 = ENO_DerivAbsLS_3D(LS1, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			LS2[idx(i, j, k)] = LS1[idx(i, j, k)]
				- kAdt
				* lsInit[idx(i, j, k)] / std::sqrt(std::pow(lsInit[idx(i, j, k)], 2.0) + kEps * kEps)
				* (absdLS2[idx(i, j, k)] - 1.0);
		}

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			ls[idx(i, j, k)] = 0.5 * (LS1[idx(i, j, k)] + LS2[idx(i, j, k)]);
		}
	}

	ApplyBC_P_3D(ls);

	return 0;
}

int LevelSetSolver3D::Reinit_Sussman_3D(std::vector<double>& ls) {
	// using 2nd order HJ ENO, based on Sussman's work
	// Sussman, Mark, et al. "An improved level set method for incompressible two-phase flows."
	// Computers & Fluids 27.5 (1998): 663-680.

	// HJENO 2nd order + RK2
	m_atime = 0.0;

	std::vector<double> lsInit(kArrSize, 0.0);
	lsInit = ls;

	std::vector<double> absdLS1(kArrSize, 0.0);
	std::vector<double> absdLS2(kArrSize, 0.0);
	std::vector<double> L1(kArrSize, 0.0);
	std::vector<double> L2(kArrSize, 0.0);
	std::vector<double> tmpLS(kArrSize, 0.0);
	std::vector<double> sussmanConstraint(kArrSize, 0.0);
	
	while (m_atime < kMaxATime) {
		m_atime += kAdt;
		
		// Apply Boundary Condition First
		ApplyBC_P_3D(ls);

		absdLS1 = ENO_DerivAbsLS_3D(ls, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			L1[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt * sign(lsInit[idx(i, j, k)]) * (absdLS1[idx(i, j, k)] - 1.0);
		}

		absdLS2 = ENO_DerivAbsLS_3D(L1, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			L2[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt * sign(lsInit[idx(i, j, k)]) * (absdLS2[idx(i, j, k)] - 1.0);
		}

		// to apply reinitialization near band level set, use temporary level set variable.
		// no need reinitialization over whole level set domain
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			tmpLS[idx(i, j, k)] = 0.5 * (L1[idx(i, j, k)] + L2[idx(i, j, k)]);
		}

		sussmanConstraint = GetSussmanReinitConstraint(tmpLS, lsInit);
		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			ls[idx(i, j, k)] += sussmanConstraint[idx(i, j, k)];
		}
	}

	ApplyBC_P_3D(ls);
	return 0;
}

std::vector<double> LevelSetSolver3D::GetSussmanReinitConstraint(const std::vector<double>& ls, const std::vector<double>& lsInit) {
	double x[3][3][3], y[3][3][3], z[3][3][3], lsFraction[3][3][3], lsInitFraction[3][3][3],
		Hp[3][3][3], L[3][3][3], dLSAbs[3][3][3], f[3][3][3];
	std::vector<double> sussmanConstraint(kArrSize, 0.0);

	// Sussman's reinitialization
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		for (int li = 0; li < 3; li++)
		for (int lj = 0; lj < 3; lj++)
		for (int lk = 0; lk < 3; lk++) {
			// i - 1/2, i, i + 1/2
			x[li][lj][lk] = kBaseX + (i - kNumBCGrid - 0.5 * (li - 1)) * kDx;
			// j - 1/2, i, j + 1/2
			y[li][lj][lk] = kBaseY + (j - kNumBCGrid - 0.5 * (lj - 1)) * kDy;
			// k - 1/2, i, k + 1/2
			z[li][lj][lk] = kBaseZ + (k - kNumBCGrid - 0.5 * (lk - 1)) * kDz;
		}

		// ls
		lsFraction[0][0][0] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i - 1, j, k - 1)]
			+ ls[idx(i - 1, j - 1, k - 1)] + ls[idx(i, j - 1, k - 1)]);
		lsFraction[0][1][0] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i - 1, j, k - 1)]);
		lsFraction[0][2][0] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i - 1, j + 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i - 1, j, k - 1)]
			+ ls[idx(i, j + 1, k - 1)] + ls[idx(i - 1, j + 1, k - 1)]);
		lsFraction[1][0][0] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j - 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i, j - 1, k - 1)]);
		lsFraction[1][1][0] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j, k - 1)]);
		lsFraction[1][2][0] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i, j + 1, k - 1)]);
		lsFraction[2][0][0] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j - 1, k)] + ls[idx(i + 1, j - 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i + 1, j, k - 1)]
			+ ls[idx(i, j - 1, k - 1)] + ls[idx(i + 1, j - 1, k - 1)]);
		lsFraction[2][1][0] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i + 1, j, k - 1)]);
		lsFraction[2][2][0] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i + 1, j + 1, k)]
			+ ls[idx(i, j, k - 1)] + ls[idx(i + 1, j, k - 1)]
			+ ls[idx(i, j + 1, k - 1)] + ls[idx(i + 1, j + 1, k - 1)]);

		lsFraction[0][0][1] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]);
		lsFraction[0][1][1] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]);
		lsFraction[0][2][1] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i - 1, j + 1, k)]);
		lsFraction[1][0][1] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j - 1, k)]);
		lsFraction[1][1][1] = ls[idx(i, j, k)];
		lsFraction[1][2][1] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]);
		lsFraction[2][0][1] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j - 1, k)] + ls[idx(i + 1, j - 1, k)]);
		lsFraction[2][1][1] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]);
		lsFraction[2][2][1] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i + 1, j + 1, k)]);

		lsFraction[0][0][2] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i - 1, j - 1, k)] + ls[idx(i, j - 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i - 1, j, k + 1)]
			+ ls[idx(i - 1, j - 1, k + 1)] + ls[idx(i, j - 1, k + 1)]);
		lsFraction[0][1][2] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i - 1, j, k + 1)]);
		lsFraction[0][2][2] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i - 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i - 1, j + 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i - 1, j, k + 1)]
			+ ls[idx(i, j + 1, k + 1)] + ls[idx(i - 1, j + 1, k + 1)]);
		lsFraction[1][0][2] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j - 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i, j - 1, k + 1)]);
		lsFraction[1][1][2] = 0.5 * (ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]);
		lsFraction[1][2][2] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i, j + 1, k + 1)]);
		lsFraction[2][0][2] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j - 1, k)] + ls[idx(i + 1, j - 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i + 1, j, k + 1)]
			+ ls[idx(i, j - 1, k + 1)] + ls[idx(i + 1, j - 1, k + 1)]);
		lsFraction[2][1][2] = 0.25 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i + 1, j, k + 1)]);
		lsFraction[2][2][2] = 0.125 * (ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]
			+ ls[idx(i, j + 1, k)] + ls[idx(i + 1, j + 1, k)]
			+ ls[idx(i, j, k + 1)] + ls[idx(i + 1, j, k + 1)]
			+ ls[idx(i, j + 1, k + 1)] + ls[idx(i + 1, j + 1, k + 1)]);

		// lsOrg
		lsInitFraction[0][0][0] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i - 1, j, k - 1)]
			+ lsInit[idx(i - 1, j - 1, k - 1)] + lsInit[idx(i, j - 1, k - 1)]);
		lsInitFraction[0][1][0] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i - 1, j, k - 1)]);
		lsInitFraction[0][2][0] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i - 1, j + 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i - 1, j, k - 1)]
			+ lsInit[idx(i, j + 1, k - 1)] + lsInit[idx(i - 1, j + 1, k - 1)]);
		lsInitFraction[1][0][0] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i, j - 1, k - 1)]);
		lsInitFraction[1][1][0] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j, k - 1)]);
		lsInitFraction[1][2][0] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j + 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i, j + 1, k - 1)]);
		lsInitFraction[2][0][0] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i + 1, j, k - 1)]
			+ lsInit[idx(i, j - 1, k - 1)] + lsInit[idx(i + 1, j - 1, k - 1)]);
		lsInitFraction[2][1][0] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i + 1, j, k - 1)]);
		lsInitFraction[2][2][0] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i, j, k - 1)] + lsInit[idx(i + 1, j, k - 1)]
			+ lsInit[idx(i, j + 1, k - 1)] + lsInit[idx(i + 1, j + 1, k - 1)]);

		lsInitFraction[0][0][1] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)]);
		lsInitFraction[0][1][1] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]);
		lsInitFraction[0][2][1] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i - 1, j + 1, k)]);
		lsInitFraction[1][0][1] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j - 1, k)]);
		lsInitFraction[1][1][1] = lsInit[idx(i, j, k)];
		lsInitFraction[1][2][1] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j + 1, k)]);
		lsInitFraction[2][0][1] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)]);
		lsInitFraction[2][1][1] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]);
		lsInitFraction[2][2][1] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)]);

		lsInitFraction[0][0][2] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i - 1, j, k + 1)]
			+ lsInit[idx(i - 1, j - 1, k + 1)] + lsInit[idx(i, j - 1, k + 1)]);
		lsInitFraction[0][1][2] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i - 1, j, k + 1)]);
		lsInitFraction[0][2][2] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i - 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i - 1, j + 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i - 1, j, k + 1)]
			+ lsInit[idx(i, j + 1, k + 1)] + lsInit[idx(i - 1, j + 1, k + 1)]);
		lsInitFraction[1][0][2] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i, j - 1, k + 1)]);
		lsInitFraction[1][1][2] = 0.5 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j, k + 1)]);
		lsInitFraction[1][2][2] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i, j + 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i, j + 1, k + 1)]);
		lsInitFraction[2][0][2] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i + 1, j, k + 1)]
			+ lsInit[idx(i, j - 1, k + 1)] + lsInit[idx(i + 1, j - 1, k + 1)]);
		lsInitFraction[2][1][2] = 0.25 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i + 1, j, k + 1)]);
		lsInitFraction[2][2][2] = 0.125 * (lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)]
			+ lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i, j, k + 1)] + lsInit[idx(i + 1, j, k + 1)]
			+ lsInit[idx(i, j + 1, k + 1)] + lsInit[idx(i + 1, j + 1, k + 1)]);

		// |\nabla \phi |
		dLSAbs[0][0][0]
			= std::sqrt(std::pow((lsInitFraction[1][0][0] -
			0.25 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j - 1, k)]
			+ lsInit[idx(i - 1, j, k - 1)] + lsInit[idx(i - 1, j - 1, k - 1)])) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][1][0] -
			0.25 * (lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i - 1, j - 1, k - 1)] + lsInit[idx(i, j - 1, k - 1)])) / kDy, 2.0));
		dLSAbs[0][1][0]
			= std::sqrt(std::pow((lsInitFraction[1][1][0] -
			0.5 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j, k - 1)])) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][2][0] - lsInitFraction[0][0][0]) / kDy, 2.0));
		dLSAbs[0][2][0]
			= std::sqrt(std::pow((lsInitFraction[1][2][0] -
			0.25 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j + 1, k)]
			+ lsInit[idx(i - 1, j, k - 1)] + lsInit[idx(i - 1, j + 1, k - 1)])) / kDx, 2.0)
			+ std::pow((0.25 * (lsInit[idx(i - 1, j + 1, k)] + lsInit[idx(i, j + 1, k)]
			+ lsInit[idx(i - 1, j + 1, k - 1)] + lsInit[idx(i, j + 1, k - 1)])
			- lsInitFraction[0][1][0]) / kDy, 2.0));
		dLSAbs[1][0][0]
			= std::sqrt(std::pow((lsInitFraction[2][0][0] - lsInitFraction[0][0][0]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][1][0] -
			0.5 * (lsInit[idx(i, j - 1, k)] + lsInit[idx(i, j - 1, k - 1)])) / kDy, 2.0));
		dLSAbs[1][1][0]
			= std::sqrt(std::pow((lsInitFraction[2][1][0] - lsInitFraction[0][1][0]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][2][0] - lsInitFraction[1][0][0]) / kDy, 2.0));
		dLSAbs[1][2][0]
			= std::sqrt(std::pow((lsInitFraction[2][2][0] - lsInitFraction[0][2][0]) / kDx, 2.0)
			+ std::pow((0.5 * (lsInit[idx(i, j + 1, k)] + lsInit[idx(i, j + 1, k - 1)])
			- lsInitFraction[1][1][0]) / kDy, 2.0));
		dLSAbs[2][0][0]
			= std::sqrt(std::pow((0.25 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i + 1, j, k - 1)] + lsInit[idx(i + 1, j - 1, k - 1)])
			- lsInitFraction[1][0][0]) / kDx, 2.0)
			+ std::pow((lsInitFraction[2][1][0] -
			0.25 * (lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i, j - 1, k - 1)] + lsInit[idx(i + 1, j - 1, k - 1)])) / kDy, 2.0));
		dLSAbs[2][1][0]
			= std::sqrt(std::pow((0.5 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j, k - 1)])
			- lsInitFraction[1][1][0]) / kDx, 2.0) +
			std::pow((lsInitFraction[2][2][0] - lsInitFraction[2][0][0]) / kDy, 2.0));
		dLSAbs[2][2][0]
			= std::sqrt(std::pow((0.25 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i + 1, j, k - 1)] + lsInit[idx(i + 1, j + 1, k - 1)])
			- lsInitFraction[1][2][0]) / kDx, 2.0)
			+ std::pow((0.25 * (lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i, j + 1, k - 1)] + lsInit[idx(i + 1, j + 1, k - 1)])
			- lsInitFraction[2][1][0]) / kDy, 2.0));

		dLSAbs[0][0][1]
			= std::sqrt(std::pow((lsInitFraction[1][0][1] -
			0.5 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j - 1, k)])) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][1][1] -
			0.5 * (lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)])) / kDy, 2.0));
		dLSAbs[0][1][1]
			= std::sqrt(std::pow((lsInitFraction[1][1][1] - lsInit[idx(i - 1, j, k)]) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][2][1] - lsInitFraction[0][0][1]) / kDy, 2.0));
		dLSAbs[0][2][1]
			= std::sqrt(std::pow((lsInitFraction[1][2][1] -
			0.5 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j + 1, k)])) / kDx, 2.0)
			+ std::pow((0.5 * (lsInit[idx(i - 1, j + 1, k)] + lsInit[idx(i, j + 1, k)])
			- lsInitFraction[0][1][1]) / kDy, 2.0));
		dLSAbs[1][0][1]
			= std::sqrt(std::pow((lsInitFraction[2][0][1] - lsInitFraction[0][0][1]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][1][1] - lsInit[idx(i, j - 1, k)]) / kDy, 2.0));
		dLSAbs[1][1][1]
			= std::sqrt(std::pow((lsInitFraction[2][1][1] - lsInitFraction[0][1][1]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][2][1] - lsInitFraction[1][0][1]) / kDy, 2.0));
		dLSAbs[1][2][1]
			= std::sqrt(std::pow((lsInitFraction[2][2][1] - lsInitFraction[0][2][1]) / kDx, 2.0)
			+ std::pow((lsInit[idx(i, j + 1, k)] - lsInitFraction[1][1][1]) / kDy, 2.0));
		dLSAbs[2][0][1]
			= std::sqrt(std::pow((0.5 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j - 1, k)])
			- lsInitFraction[1][0][1]) / kDx, 2.0)
			+ std::pow((lsInitFraction[2][1][1] -
			0.5 * (lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)])) / kDy, 2.0));
		dLSAbs[2][1][1]
			= std::sqrt(std::pow((lsInit[idx(i + 1, j, k)] - lsInitFraction[1][1][1]) / kDx, 2.0) +
			std::pow((lsInitFraction[2][2][1] - lsInitFraction[2][0][1]) / kDy, 2.0));
		dLSAbs[2][2][1]
			= std::sqrt(std::pow((0.5 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j + 1, k)])
			- lsInitFraction[1][2][1]) / kDx, 2.0)
			+ std::pow((0.5 * (lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)])
			- lsInitFraction[2][1][1]) / kDy, 2.0));

		dLSAbs[0][0][2]
			= std::sqrt(std::pow((lsInitFraction[1][0][2] -
			0.25 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j - 1, k)]
			+ lsInit[idx(i - 1, j, k + 1)] + lsInit[idx(i - 1, j - 1, k + 1)])) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][1][2] -
			0.25 * (lsInit[idx(i - 1, j - 1, k)] + lsInit[idx(i, j - 1, k)]
			+ lsInit[idx(i - 1, j - 1, k + 1)] + lsInit[idx(i, j - 1, k + 1)])) / kDy, 2.0));
		dLSAbs[0][1][2]
			= std::sqrt(std::pow((lsInitFraction[1][1][2] -
			0.5 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j, k + 1)])) / kDx, 2.0)
			+ std::pow((lsInitFraction[0][2][2] - lsInitFraction[0][0][2]) / kDy, 2.0));
		dLSAbs[0][2][2]
			= std::sqrt(std::pow((lsInitFraction[1][2][2] -
			0.25 * (lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 1, j + 1, k)]
			+ lsInit[idx(i - 1, j, k + 1)] + lsInit[idx(i - 1, j + 1, k + 1)])) / kDx, 2.0)
			+ std::pow((0.25 * (lsInit[idx(i - 1, j + 1, k)] + lsInit[idx(i, j + 1, k)]
			+ lsInit[idx(i - 1, j + 1, k + 1)] + lsInit[idx(i, j + 1, k + 1)])
			- lsInitFraction[0][1][2]) / kDy, 2.0));
		dLSAbs[1][0][2]
			= std::sqrt(std::pow((lsInitFraction[2][0][2] - lsInitFraction[0][0][2]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][1][2] -
			0.5 * (lsInit[idx(i, j - 1, k)] + lsInit[idx(i, j - 1, k + 1)])) / kDy, 2.0));
		dLSAbs[1][1][2]
			= std::sqrt(std::pow((lsInitFraction[2][1][2] - lsInitFraction[0][1][2]) / kDx, 2.0)
			+ std::pow((lsInitFraction[1][2][2] - lsInitFraction[1][0][2]) / kDy, 2.0));
		dLSAbs[1][2][2]
			= std::sqrt(std::pow((lsInitFraction[2][2][2] - lsInitFraction[0][2][2]) / kDx, 2.0)
			+ std::pow((0.5 * (lsInit[idx(i, j + 1, k)] + lsInit[idx(i, j + 1, k + 1)])
			- lsInitFraction[1][1][2]) / kDy, 2.0));
		dLSAbs[2][0][2]
			= std::sqrt(std::pow((0.25 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i + 1, j, k + 1)] + lsInit[idx(i + 1, j - 1, k + 1)])
			- lsInitFraction[1][0][2]) / kDx, 2.0)
			+ std::pow((lsInitFraction[2][1][2] -
			0.25 * (lsInit[idx(i, j - 1, k)] + lsInit[idx(i + 1, j - 1, k)]
			+ lsInit[idx(i, j - 1, k + 1)] + lsInit[idx(i + 1, j - 1, k + 1)])) / kDy, 2.0));
		dLSAbs[2][1][2]
			= std::sqrt(std::pow((0.5 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j, k + 1)])
			- lsInitFraction[1][1][2]) / kDx, 2.0) +
			std::pow((lsInitFraction[2][2][2] - lsInitFraction[2][0][2]) / kDy, 2.0));
		dLSAbs[2][2][2]
			= std::sqrt(std::pow((0.25 * (lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i + 1, j, k + 1)] + lsInit[idx(i + 1, j + 1, k + 1)])
			- lsInitFraction[1][2][2]) / kDx, 2.0)
			+ std::pow((0.25 * (lsInit[idx(i, j + 1, k)] + lsInit[idx(i + 1, j + 1, k)]
			+ lsInit[idx(i, j + 1, k + 1)] + lsInit[idx(i + 1, j + 1, k + 1)])
			- lsInitFraction[2][1][2]) / kDy, 2.0));

		// Hp, L
		for (int li = 0; li < 3; li++)
		for (int lj = 0; lj < 3; lj++)
		for (int lk = 0; lk < 3; lk++)  {
			if (std::fabs(lsFraction[li][lj][lk]) > kThickness)
				Hp[li][lj][lk] = 0.0;
			else
				Hp[li][lj][lk] = 0.5 / kEps + 0.5 / kEps * cos(M_PI * lsFraction[li][lj][lk] / kEps);
			L[li][lj][lk] = (lsFraction[li][lj][lk] - lsInitFraction[li][lj][lk]) / kAdt;

			f[li][lj][lk] = Hp[li][lj][lk] * dLSAbs[li][lj][lk];
		}

		// Simpson's Rule with 27 point stencil
		double sum1 = 0.0, sum2 = 0.0;
		for (int li = 0; li < 3; li++)
		for (int lj = 0; lj < 3; lj++)
		for (int lk = 0; lk < 3; lk++)  {
			sum1 += Hp[li][lj][lk] * L[li][lj][lk];
			sum2 += Hp[li][lj][lk] * f[li][lj][lk];
		}
		sum1 += 35.0 * Hp[1][1][1] * L[1][1][1];
		sum2 += 35.0 * Hp[1][1][1] * f[1][1][1];

		sussmanConstraint[idx(i, j, k)] -= kAdt * f[1][1][1] * sum1 / sum2;
	}

	return sussmanConstraint;
}

int LevelSetSolver3D::FirstTimeOnlyReinit_Sussman_3D(std::vector<double>& ls) {
	// using 2nd order HJ ENO, based on Sussman's work
	// Sussman, Mark, et al. "An improved level set method for incompressible two-phase flows."
	// Computers & Fluids 27.5 (1998): 663-680.

	// HJENO 2nd order + RK2
	m_atime = 0.0;
	std::vector<double> lsInit(kArrSize, 0.0);
	lsInit = ls;

	std::vector<double> absdLS1(kArrSize, 0.0);
	std::vector<double> absdLS2(kArrSize, 0.0);
	std::vector<double> L1(kArrSize, 0.0);
	std::vector<double> L2(kArrSize, 0.0);

	absdLS1 = ENO_DerivAbsLS_3D(lsInit, lsInit);

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		if (absdLS1[idx(i, j, k)] > 0.0)
			lsInit[idx(i, j, k)] = lsInit[idx(i, j, k)] / absdLS1[idx(i, j, k)];
	}

	while (m_atime < std::max(std::max(kDx * kNx, kDy * kNy), kDz * kNz)) {
		m_atime += kAdt;
		
		// Apply Boundary Condition First
		ApplyBC_P_3D(ls);

		absdLS1 = ENO_DerivAbsLS_3D(ls, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			L1[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt * sign(lsInit[idx(i, j, k)]) * (absdLS1[idx(i, j, k)] - 1.0);
		}

		absdLS2 = ENO_DerivAbsLS_3D(L1, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			L2[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt * sign(lsInit[idx(i, j, k)]) * (absdLS2[idx(i, j, k)] - 1.0);
		}

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			ls[idx(i, j, k)] = 0.5 * (L1[idx(i, j, k)] + L2[idx(i, j, k)]);
		}

	}

	ApplyBC_P_3D(ls);

	return 0;
}

int LevelSetSolver3D::Reinit_MinRK2_3D(std::vector<double>& ls) {
	m_atime = 0.0;

	std::vector<double> lsInit(kArrSize, 0.0);
	lsInit = ls;


	double lambda = 0.0;
	double x[3][3], y[3][3], lsFraction[3][3], lsInitFraction[3][3], Hp[3][3], L[3][3], dLSAbs[3][3], f[3][3];
	std::vector<double> absdLS1(kArrSize, 0.0);
	std::vector<double> absdLS2(kArrSize, 0.0);
	std::vector<double> LS1(kArrSize, 0.0);
	std::vector<double> LS2(kArrSize, 0.0);

	while (m_atime < kMaxATime) {
		m_atime += kAdt;
		
		// Apply Boundary Condition First
		ApplyBC_P_3D(ls);

		absdLS1 = Subcell_DerivAbsLS_3D(ls, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			LS1[idx(i, j, k)] = ls[idx(i, j, k)]
				- kAdt
				* sign(lsInit[idx(i, j, k)])
				* (absdLS1[idx(i, j, k)] - 1.0);
		}

		absdLS2 = Subcell_DerivAbsLS_3D(LS1, lsInit);

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			LS2[idx(i, j, k)] = LS1[idx(i, j, k)]
				- kAdt
				* sign(lsInit[idx(i, j, k)])
				* (absdLS2[idx(i, j, k)] - 1.0);
		}

		for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
			ls[idx(i, j, k)] = 0.5 * (LS1[idx(i, j, k)] + LS2[idx(i, j, k)]);
		}
	}

	ApplyBC_P_3D(ls);
	return 0;
}

std::vector<double>
LevelSetSolver3D::ENO_DerivAbsLS_3D(const std::vector<double>& ls, const std::vector<double>& lsInit) {
	// Use Newton Polynomial
	double x = 0.0, y = 0.0;

	// Divided Difference Table (a.k.a. DDT)
	std::vector<double>
		DDTX((2 * kENOSpatialOrder + 1) * (kENOSpatialOrder + 1), 0.0),
		DDTY((2 * kENOSpatialOrder + 1) * (kENOSpatialOrder + 1), 0.0);
	
	std::vector<double>
		dPX(kArrSize, 0.0),
		dMX(kArrSize, 0.0),
		dPY(kArrSize, 0.0),
		dMY(kArrSize, 0.0),
		dPZ(kArrSize, 0.0),
		dMZ(kArrSize, 0.0);

	std::vector<double>
		dLSdX(kArrSize, 0.0),
		dLSdY(kArrSize, 0.0),
		dLSdZ(kArrSize, 0.0);

	std::vector<double> absdLS(kArrSize, 0.0);

	//	const int kDDTRows = 2 * kENOSpatialOrder + 1;
	// general divided difference table, but not used in this method
	/*
	Ex) 3rd order ENO, total array size is 28. and { } indicates index of element.
	----------------------------------------------------------------
	|D0_(i+3) {6}	|{13}			|{20}			|{27}			|
	|D0_(i+2)		|D1_(i+2, i+3)	|				|				|
	|D0_(i+1)		|D1_(i+1, i+2)	|D2_(i+1,i+3)	|				|
	|D0_(i)			|D1_(i, i+1)	|D2_(i,i+2)		|D3_(i,i+3)		|
	|D0_(i-1)		|D1_(i-1, i)	|D2_(i-1,i+1)	|D3_(i-1,i+2)	|
	|D0_(i-2)		|D1_(i-2, i-1)	|D2_(i-2,i)		|D3_(i-2,i+1)	|
	|D0_(i-3) {0}	|D1_(i-3, i-2)	|D2_(i-3,i-1)r	|D3_(i-3,i)		|
	----------------------------------------------------------------

	Ex) 2nd order ENO, total array size is 28. and { } indicates index of element.
	-------------------------------------------------
	|D0_(i+2)		|				|				|
	|D0_(i+1)		|D1_(i+1, i+2)	|				|
	|D0_(i)			|D1_(i, i+1)	|D2_(i,i+2)		|
	|D0_(i-1)		|D1_(i-1, i)	|D2_(i-1,i+1)	|
	|D0_(i-2)		|D1_(i-2, i-1)	|D2_(i-2,i)		|
	-------------------------------------------------
	*/
	/*
	// init 0th divided difference table
	for (int d = 0; d < 2 * kENOSpatialOrder + 1; d++)
	DDTX[d] = ls[idx(i - kENOSpatialOrder + d, j)];

	// init 1st ~ kENOSpatialOrder(th) divided difference table
	for (int r = 1; r <= kENOSpatialOrder; r++) {
	for (int d = 0; d < kDDTRows - r; d++)
	DDTX[kDDTRows * r + d]
	= (DDTX[kDDTRows * (r - 1) + d + 1] - DDTX[kDDTRows * (r - 1) + d])
	/ (kDx * r);
	}
	*/

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPX[idx(i, j, k)] = (ls[idx(i + 1, j, k)] - ls[idx(i, j, k)]) / kDx
			- kDx * 0.5 * MinMod(
			(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i + 1, j, k)] + ls[idx(i + 2, j, k)]) / (kDx * kDx));
		dMX[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i - 1, j, k)]) / kDx
			+ kDx * 0.5 * MinMod(
			(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
			(ls[idx(i - 2, j, k)] - 2.0 * ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]) / (kDx * kDx));
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPY[idx(i, j, k)] = (ls[idx(i, j + 1, k)] - ls[idx(i, j, k)]) / kDy
			- kDy * 0.5 * MinMod(
			(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j + 1, k)] + ls[idx(i, j + 2, k)]) / (kDy * kDy));
		dMY[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i, j - 1, k)]) / kDy
			+ kDy * 0.5 * MinMod(
			(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
			(ls[idx(i, j - 2, k)] - 2.0 * ls[idx(i, j - 1, k)] + ls[idx(i, j, k)]) / (kDy * kDy));
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPZ[idx(i, j, k)] = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k)]) / kDz
			- kDz * 0.5 * MinMod(
			(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j, k + 1)] + ls[idx(i, j, k + 2)]) / (kDz * kDz));
		dMZ[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i, j, k - 1)]) / kDz
			+ kDz * 0.5 * MinMod(
			(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
			(ls[idx(i, j, k - 2)] - 2.0 * ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]) / (kDz * kDz));
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		if (lsInit[idx(i, j, k)] >= 0)
			absdLS[idx(i, j, k)] = std::sqrt(
			std::max(std::pow(std::min(dPX[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::max(dMX[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::min(dPY[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::max(dMY[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::min(dPZ[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::max(dMZ[idx(i, j, k)], 0.0), 2.0)));
		else
			absdLS[idx(i, j, k)] = std::sqrt(
			std::max(std::pow(std::max(dPX[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::min(dMX[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::max(dPY[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::min(dMY[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::max(dPZ[idx(i, j, k)], 0.0), 2.0),
				std::pow(std::min(dMZ[idx(i, j, k)], 0.0), 2.0)));
	}

	return absdLS;
}

std::vector<double>
LevelSetSolver3D::Subcell_DerivAbsLS_3D(const std::vector<double>& ls, const std::vector<double>& lsInit) {
	// Use Newton Polynomial
	double x = 0.0, y = 0.0;

	// Divided Difference Table (a.k.a. DDT)
	std::vector<double>
		DDTX((2 * kENOSpatialOrder + 1) * (kENOSpatialOrder + 1), 0.0),
		DDTY((2 * kENOSpatialOrder + 1) * (kENOSpatialOrder + 1), 0.0);

	std::vector<double>
		dPX(kArrSize, 0.0),
		dMX(kArrSize, 0.0),
		dPY(kArrSize, 0.0),
		dMY(kArrSize, 0.0),
		dPZ(kArrSize, 0.0),
		dMZ(kArrSize, 0.0);

	std::vector<double>
		dLSdX(kArrSize, 0.0),
		dLSdY(kArrSize, 0.0),
		dLSdZ(kArrSize, 0.0);

	std::vector<double> absdLS(kArrSize, 0.0);

	// new dx or dy using subcell resolution
	double newDxP = 0.0, newDxM = 0.0, newDyP = 0.0, newDyM = 0.0, newDzP = 0.0, newDzM = 0.0;
	// minmod of undivided differences of lsInitial
	double uDDX = 0.0, uDDY = 0.0, uDDZ = 0.0;
	// discriminant of the qudaratic polynomial
	double D = 0.0;
	const double kSmallValue = 1.0e-10;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPX[idx(i, j, k)] = (ls[idx(i + 1, j, k)] - ls[idx(i, j, k)]) / kDx
			- kDx * 0.5 * MinMod(
			(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i + 1, j, k)] + ls[idx(i + 2, j, k)]) / (kDx * kDx));
		dMX[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i - 1, j, k)]) / kDx
			+ kDx * 0.5 * MinMod(
			(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
			(ls[idx(i - 2, j, k)] - 2.0 * ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]) / (kDx * kDx));

		if (ls[idx(i, j, k)] * ls[idx(i + 1, j, k)] < 0) {
			uDDX = MinMod(
				lsInit[idx(i - 1, j, k)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i + 1, j, k)] + lsInit[idx(i + 2, j, k)]);
			// std::fabs(uDDX) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDX) <= kSmallValue : linear interpolation
			if (std::fabs(uDDX) > kSmallValue) {
				// this coefficient is correct? I'm not sure
				// A second order accurate level set method on non - graded adaptive cartesian grids
				// On reinitializing level set functions
				// Two papers .. which one is correct?
				D = std::pow(uDDX * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i + 1, j, k)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i + 1, j, k)];
				newDxP = kDx * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i + 1, j, k)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i + 1, j, k)]) * std::sqrt(D))
					/ uDDX);
			}
			else
				newDxP = kDx * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i + 1, j, k)]);

			dPX[idx(i, j, k)] = -ls[idx(i, j, k)] / newDxP
				- newDxP * 0.5 *  MinMod(
				(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
				(ls[idx(i, j, k)] - 2.0 * ls[idx(i + 1, j, k)] + ls[idx(i + 2, j, k)]) / (kDx * kDx));
		}

		if (ls[idx(i, j, k)] * ls[idx(i - 1, j, k)] < 0) {
			uDDX = MinMod(
				lsInit[idx(i - 1, j, k)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i + 1, j, k)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i - 1, j, k)] + lsInit[idx(i - 2, j, k)]);
			// std::fabs(uDDX) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDX) <= kSmallValue : linear interpolation
			if (std::fabs(uDDX) > kSmallValue) {
				D = std::pow(uDDX * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i - 1, j, k)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i - 1, j, k)];
				newDxM = kDx * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i - 1, j, k)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i - 1, j, k)]) * std::sqrt(D))
					/ uDDX);
			}
			else
				newDxM = kDx * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i - 1, j, k)]);

			dMX[idx(i, j, k)] = ls[idx(i, j, k)] / newDxM
				+ newDxM * 0.5 *  MinMod(
				(ls[idx(i - 1, j, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i + 1, j, k)]) / (kDx * kDx),
				(ls[idx(i - 2, j, k)] - 2.0 * ls[idx(i - 1, j, k)] + ls[idx(i, j, k)]) / (kDx * kDx));
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPY[idx(i, j, k)] = (ls[idx(i, j + 1, k)] - ls[idx(i, j, k)]) / kDy
			- kDy * 0.5 * MinMod(
			(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j + 1, k)] + ls[idx(i, j + 2, k)]) / (kDy * kDy));
		dMY[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i, j - 1, k)]) / kDy
			+ kDy * 0.5 * MinMod(
			(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
			(ls[idx(i, j - 2, k)] - 2.0 * ls[idx(i, j - 1, k)] + ls[idx(i, j, k)]) / (kDy * kDy));

		if (ls[idx(i, j, k)] * ls[idx(i, j + 1, k)] < 0) {
			uDDY = MinMod(
				lsInit[idx(i, j - 1, k)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i, j + 1, k)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i, j + 1, k)] + lsInit[idx(i, j + 2, k)]);
			// std::fabs(uDDY) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDY) <= kSmallValue : linear interpolation
			if (std::fabs(uDDY) > kSmallValue) {
				D = std::pow(uDDY * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i, j + 1, k)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i, j + 1, k)];
				newDyP = kDy * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i, j + 1, k)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i, j + 1, k)]) * std::sqrt(D))
					/ uDDY);
			}
			else
				newDyP = kDy * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i, j + 1, k)]);

			dPY[idx(i, j, k)] = -ls[idx(i, j, k)] / newDyP
				- newDyP * 0.5 *  MinMod(
				(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
				(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j + 1, k)] + ls[idx(i, j + 2, k)]) / (kDy * kDy));
		}
		
		if (ls[idx(i, j, k)] * ls[idx(i, j - 1, k)] < 0) {
			uDDY = MinMod(
				lsInit[idx(i, j - 1, k)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i, j + 1, k)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i, j - 1, k)] + lsInit[idx(i, j - 2, k)]);
			// std::fabs(uDDY) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDY) <= kSmallValue : linear interpolation
			if (std::fabs(uDDY) > kSmallValue) {
				D = std::pow(uDDY * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i, j - 1, k)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i, j - 1, k)];
				newDyM = kDy * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i, j - 1, k)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i, j - 1, k)]) * std::sqrt(D))
					/ uDDY);
			}
			else
				newDyM = kDy * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i, j - 1, k)]);

			dMY[idx(i, j, k)] = ls[idx(i, j, k)] / newDyM
				+ newDyM * 0.5 *  MinMod(
				(ls[idx(i, j - 1, k)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j + 1, k)]) / (kDy * kDy),
				(ls[idx(i, j - 2, k)] - 2.0 * ls[idx(i, j - 1, k)] + ls[idx(i, j, k)]) / (kDy * kDy));
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		// 2nd order ENO Method
		dPZ[idx(i, j, k)] = (ls[idx(i, j, k + 1)] - ls[idx(i, j, k)]) / kDz
			- kDz * 0.5 * MinMod(
			(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
			(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j, k + 1)] + ls[idx(i, j, k + 2)]) / (kDz * kDz));
		dMZ[idx(i, j, k)] = (ls[idx(i, j, k)] - ls[idx(i, j, k - 1)]) / kDz
			+ kDz * 0.5 * MinMod(
			(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
			(ls[idx(i, j, k - 2)] - 2.0 * ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]) / (kDz * kDz));

		if (ls[idx(i, j, k)] * ls[idx(i, j, k + 1)] < 0) {
			uDDZ = MinMod(
				lsInit[idx(i, j, k - 1)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i, j, k + 1)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i, j, k + 1)] + lsInit[idx(i, j, k + 2)]);
			// std::fabs(uDDZ) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDZ) <= kSmallValue : linear interpolation
			if (std::fabs(uDDZ) > kSmallValue) {
				D = std::pow(uDDZ * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i, j, k + 1)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i, j, k + 1)];
				newDzP = kDz * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i, j, k + 1)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i, j, k + 1)]) * std::sqrt(D))
					/ uDDZ);
			}
			else
				newDzP = kDy * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i, j, k + 1)]);

			dPZ[idx(i, j, k)] = -ls[idx(i, j, k)] / newDzP
				- newDzP * 0.5 *  MinMod(
				(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
				(ls[idx(i, j, k)] - 2.0 * ls[idx(i, j, k + 1)] + ls[idx(i, j, k + 2)]) / (kDz * kDz));
		}

		if (ls[idx(i, j, k)] * ls[idx(i, j, k - 1)] < 0) {
			uDDZ = MinMod(
				lsInit[idx(i, j, k - 1)] - 2.0 * lsInit[idx(i, j, k)] + lsInit[idx(i, j, k + 1)],
				lsInit[idx(i, j, k)] - 2.0 * lsInit[idx(i, j, k - 1)] + lsInit[idx(i, j, k - 2)]);
			// std::fabs(uDDZ) > kSmallValue : quadaratic interpolation
			// std::fabs(uDDZ) <= kSmallValue : linear interpolation
			if (std::fabs(uDDZ) > kSmallValue) {
				D = std::pow(uDDZ * 0.5 - lsInit[idx(i, j, k)] - lsInit[idx(i, j, k - 1)], 2.0)
					- 4.0 * lsInit[idx(i, j, k)] * lsInit[idx(i, j, k - 1)];
				newDzM = kDz * (0.5 +
					(lsInit[idx(i, j, k)] - lsInit[idx(i, j, k - 1)]
					- sign(lsInit[idx(i, j, k)] - lsInit[idx(i, j, k - 1)]) * std::sqrt(D))
					/ uDDZ);
			}
			else
				newDzM = kDz * lsInit[idx(i, j, k)] / (lsInit[idx(i, j, k)] - lsInit[idx(i, j, k - 1)]);

			dMZ[idx(i, j, k)] = ls[idx(i, j, k)] / newDzM
				+ newDzM * 0.5 *  MinMod(
				(ls[idx(i, j, k - 1)] - 2.0 * ls[idx(i, j, k)] + ls[idx(i, j, k + 1)]) / (kDz * kDz),
				(ls[idx(i, j, k - 2)] - 2.0 * ls[idx(i, j, k - 1)] + ls[idx(i, j, k)]) / (kDz * kDz));
		}
	}

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++) {
		if (lsInit[idx(i, j, k)] >= 0)
			absdLS[idx(i, j, k)] = std::sqrt(
			std::max(std::pow(std::min(dPX[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::max(dMX[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::min(dPY[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::max(dMY[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::min(dPZ[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::max(dMZ[idx(i, j, k)], 0.0), 2.0)));
		else
			absdLS[idx(i, j, k)] = std::sqrt(
			std::max(std::pow(std::max(dPX[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::min(dMX[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::max(dPY[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::min(dMY[idx(i, j, k)], 0.0), 2.0))
			+ std::max(std::pow(std::max(dPZ[idx(i, j, k)], 0.0), 2.0),
			std::pow(std::min(dMZ[idx(i, j, k)], 0.0), 2.0)));
	}

	return absdLS;
}

std::vector<int> LevelSetSolver3D::GetSignedLSNormalized(const std::vector<double> &v) {
	std::vector<int> r(v.size());
	// std::transform(v.begin(), v.end(), r.begin(), (static_cast<int>((const int&)sign)));
	
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
	for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)	{
		if (std::fabs(v[idx(i, j, k)]) <= kEps)
			r[idx(i, j, k)]
				= static_cast<int>(
				2.0 * ((0.5 + cos(M_PI * v[idx(i, j, k)] / kEps)) / kEps - 0.5));
		else if (std::fabs(v[idx(i, j, k)]) > kEps && v[idx(i, j, k)] > 0)
			r[idx(i, j, k)] = 1;
		else if (std::fabs(v[idx(i, j, k)]) > kEps && v[idx(i, j, k)] < 0)
			r[idx(i, j, k)] = -1;
	}

	return r;
}

int LevelSetSolver3D::SetBC_U_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_U_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int LevelSetSolver3D::SetBC_V_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNz, kNy, kNumBCGrid);
	}

	m_BC->SetBC_V_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int LevelSetSolver3D::SetBC_W_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNz, kNy, kNumBCGrid);
	}

	m_BC->SetBC_W_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int LevelSetSolver3D::SetBC_P_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {
	if (!m_BC) {
		m_BC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	}

	m_BC->SetBC_P_3D(BC_W, BC_E, BC_S, BC_N, BC_B, BC_T);

	return 0;
}

int LevelSetSolver3D::ApplyBC_U_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_U_3D(arr);
	
	return 0;
}

int LevelSetSolver3D::ApplyBC_V_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_V_3D(arr);

	return 0;
}

int LevelSetSolver3D::ApplyBC_W_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_W_3D(arr);

	return 0;
}

int LevelSetSolver3D::ApplyBC_P_3D(std::vector<double>& arr) {
	m_BC->ApplyBC_P_3D(arr);
	
	return 0;
}

// Prescribed bounadry condition (Dirichlet) for U grid
void LevelSetSolver3D::SetBCConstantUW(double BC_ConstantW) {
	return m_BC->SetBCConstantUW(BC_ConstantW);
}

void LevelSetSolver3D::SetBCConstantUE(double BC_ConstantE) {
	return m_BC->SetBCConstantUE(BC_ConstantE);
}

void LevelSetSolver3D::SetBCConstantUS(double BC_ConstantS) {
	return m_BC->SetBCConstantUS(BC_ConstantS);
}

void LevelSetSolver3D::SetBCConstantUN(double BC_ConstantN) {
	return m_BC->SetBCConstantUN(BC_ConstantN);
}

void LevelSetSolver3D::SetBCConstantUB(double BC_ConstantB) {
	return m_BC->SetBCConstantUB(BC_ConstantB);
}

void LevelSetSolver3D::SetBCConstantUT(double BC_ConstantT) {
	return m_BC->SetBCConstantUT(BC_ConstantT);
}

// Prescribed bounadry condition (Dirichlet) for V grid
void LevelSetSolver3D::SetBCConstantVW(double BC_ConstantW) {
	return m_BC->SetBCConstantVW(BC_ConstantW);
}

void LevelSetSolver3D::SetBCConstantVE(double BC_ConstantE) {
	return m_BC->SetBCConstantVE(BC_ConstantE);
}

void LevelSetSolver3D::SetBCConstantVS(double BC_ConstantS) {
	return m_BC->SetBCConstantVS(BC_ConstantS);
}

void LevelSetSolver3D::SetBCConstantVN(double BC_ConstantN) {
	return m_BC->SetBCConstantVN(BC_ConstantN);
}

void LevelSetSolver3D::SetBCConstantVB(double BC_ConstantB) {
	return m_BC->SetBCConstantVB(BC_ConstantB);
}

void LevelSetSolver3D::SetBCConstantVT(double BC_ConstantT) {
	return m_BC->SetBCConstantVT(BC_ConstantT);
}

// Prescribed bounadry condition (Dirichlet) for W grid
void LevelSetSolver3D::SetBCConstantWW(double BC_ConstantW) {
	return m_BC->SetBCConstantWW(BC_ConstantW);
}

void LevelSetSolver3D::SetBCConstantWE(double BC_ConstantE) {
	return m_BC->SetBCConstantWE(BC_ConstantE);
}

void LevelSetSolver3D::SetBCConstantWS(double BC_ConstantS) {
	return m_BC->SetBCConstantWS(BC_ConstantS);
}

void LevelSetSolver3D::SetBCConstantWN(double BC_ConstantN) {
	return m_BC->SetBCConstantWN(BC_ConstantN);
}

void LevelSetSolver3D::SetBCConstantWB(double BC_ConstantB) {
	return m_BC->SetBCConstantWB(BC_ConstantB);
}

void LevelSetSolver3D::SetBCConstantWT(double BC_ConstantT) {
	return m_BC->SetBCConstantWT(BC_ConstantT);
}

// Prescribed bounadry condition (Dirichlet) for P grid
void LevelSetSolver3D::SetBCConstantPW(double BC_ConstantW) {
	return m_BC->SetBCConstantPW(BC_ConstantW);
}

void LevelSetSolver3D::SetBCConstantPE(double BC_ConstantE) {
	return m_BC->SetBCConstantPE(BC_ConstantE);
}

void LevelSetSolver3D::SetBCConstantPS(double BC_ConstantS) {
	return m_BC->SetBCConstantPS(BC_ConstantS);
}

void LevelSetSolver3D::SetBCConstantPN(double BC_ConstantN) {
	return m_BC->SetBCConstantPN(BC_ConstantN);
}

void LevelSetSolver3D::SetBCConstantPB(double BC_ConstantB) {
	return m_BC->SetBCConstantPB(BC_ConstantB);
}

void LevelSetSolver3D::SetBCConstantPT(double BC_ConstantT) {
	return m_BC->SetBCConstantPT(BC_ConstantT);
}

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
inline int LevelSetSolver3D::sign(const double& val) {
	return (val > 0) - (val < 0);
}

double LevelSetSolver3D::MinAbs(const double& val1, const double& val2) {
	if (std::fabs(val1) < std::fabs(val2))
		return val1;
	else
		return val2;

	// impossible
	return std::numeric_limits<double>::quiet_NaN();
}

double LevelSetSolver3D::MinMod(const double& val1, const double& val2) {
	if ((std::fabs(val1) <= std::fabs(val2)) && (val1 * val2 > 0)) 
		return val1;
	else if ((std::fabs(val1) > std::fabs(val2)) && (val1 * val2 > 0))
		return val2;
	else if (val1 * val2 <= 0)
		return 0.0;

	// impossible
	return std::numeric_limits<double>::quiet_NaN();
}

inline int LevelSetSolver3D::idx(int i, int j, int k) {
	return (k + (kNz + 2 * kNumBCGrid) * j + (kNz + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid) * i);
}