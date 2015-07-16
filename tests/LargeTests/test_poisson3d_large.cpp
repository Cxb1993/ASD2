#include "test_poisson3d_large.h"

int test_poisson3D_CG() {
	const int kNx = 64, kNy = 64, kNz = 64, kNumBCGrid = 3;
	const double kBaseX = 0.0, kBaseY = 0.0, kBaseZ = 0.0, kLenX = 1.0, kLenY = 1.0, kLenZ = 1.0;
	const double kDx = kLenX / kNx, kDy = kLenY / kNy, kDz = kLenZ / kNz;

	std::shared_ptr<PoissonSolver3D> PSolver = std::make_shared<PoissonSolver3D>(kNx, kNy, kNz, kNumBCGrid);
	std::shared_ptr<BoundaryCondition3D> PBC = std::make_shared<BoundaryCondition3D>(kNx, kNy, kNz, kNumBCGrid);
	std::vector<double> u((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid) * (kNz + 2 * kNumBCGrid), 0.0);
	std::vector<double> b((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid) * (kNz + 2 * kNumBCGrid), 0.0);
	PBC->SetBC_P_3D("dirichlet", "dirichlet", "dirichlet", "dirichlet", "dirichlet", "dirichlet");
	PBC->SetBCConstantPW(0.0);
	PBC->SetBCConstantPE(0.0);
	PBC->SetBCConstantPS(0.0);
	PBC->SetBCConstantPN(0.0);
	PBC->SetBCConstantPB(0.0);
	PBC->SetBCConstantPT(1.0);
	
	double x = 0.0, y = 0.0, z = 0.0;

	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		x = kBaseX + (i - kNumBCGrid) * kDx;
		y = kBaseY + (j - kNumBCGrid) * kDy;
		z = kBaseZ + (k - kNumBCGrid) * kDz;
		b[idx3_3D(kNy, kNz, i, j, k)] = 0.0;
	}

	std::vector<double> AVals, DiagVals, MVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy * kNz, Ancols = kNx * kNy * kNz;
	MKL_INT size = kNx * kNy * kNz;
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
	// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
	// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
	
	// Set default values, if a current pointer is in interior, it will not be changed.
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		AValsDic.clear();
		AColsDic.clear();
		tmpRowIdx = 0;
		tmpMRowIdx = 0;
		// Add starting rowIdx
		ARowIdx.push_back(rowIdx);
		DiagRowIdx.push_back(MRowIdx);

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["B"] = -1.0 / (kDz * kDz);
		AColsDic["B"] = (i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - 1 - kNumBCGrid);
		AValsDic["S"] = -1.0 / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + kNx * (j - 1 - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid);
		AValsDic["W"] = -1.0 / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid);
		AValsDic["C"] = (1.0 + 1.0) / (kDx * kDx) + (1.0 + 1.0) / (kDy * kDy) + (1.0 + 1.0) / (kDz * kDz);
		AColsDic["C"] = (i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid);
		AValsDic["E"] = -1.0 / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid);
		AValsDic["N"] = -1.0 / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + kNx * (j + 1 - kNumBCGrid) + kNx * kNy * (k - kNumBCGrid);
		AValsDic["T"] = -1.0 / (kDz * kDz);
		AColsDic["T"] = (i - kNumBCGrid) + kNx * (j - kNumBCGrid) + kNx * kNy * (k + 1 - kNumBCGrid);

		if (i == kNumBCGrid && PBC->m_BC_PW == BC3D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid && PBC->m_BC_PW == BC3D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPW);
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC3D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC3D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPE);
		}

		if (j == kNumBCGrid && PBC->m_BC_PS == BC3D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid && PBC->m_BC_PS == BC3D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPS);
		}

		if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC3D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC3D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPN);
		}

		if (k == kNumBCGrid && PBC->m_BC_PB == BC3D::NEUMANN) {
			AColsDic["B"] = -1;
			AValsDic["B"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDz * kDz);
		}
		else if (k == kNumBCGrid && PBC->m_BC_PB == BC3D::DIRICHLET) {
			AColsDic["B"] = -1;
			AValsDic["B"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDz * kDz) * (PBC->m_BC_DirichletConstantPB);
		}

		if (k == kNumBCGrid + kNz - 1 && PBC->m_BC_PT == BC3D::NEUMANN) {
			AColsDic["T"] = -1;
			AValsDic["T"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDz* kDz);
		}
		else if (k == kNumBCGrid + kNz - 1 && PBC->m_BC_PT == BC3D::DIRICHLET) {
			AColsDic["T"] = -1;
			AValsDic["T"] = 0.0;
			b[idx3_3D(kNy, kNz, i, j, k)] -= -1.0 / (kDz * kDz) * (PBC->m_BC_DirichletConstantPT);
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
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);

	PSolver->CG_2FUniformP_3D(u, b, AVals, ACols, ARowIdx, DiagVals, DiagCols, DiagRowIdx, PBC, 10000);
	
	std::string fname_p_base("TestPoisson3D_CG");
	std::ofstream outF;
	std::string fname_p(fname_p_base + "_ASCII.plt");

	outF.open(fname_p.c_str(), std::ios::out);

	outF << "TITLE = VEL" << std::endl;
	outF << "VARIABLES = \"X\", \"Y\", \"Z\", \"P\"" << std::endl;
	outF.close();

	outF.open(fname_p.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"1\", I=") << kNx << std::string(", J=") << kNy << std::string(", K=") << kNz
		<< std::string(", DATAPACKING=POINT")
		<< std::string(", SOLUTIONTIME=") << 0.0
		<< std::string(", STRANDID=") << 1
		<< std::endl;
	
	for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
		<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
		<< kBaseZ + static_cast<double>(k + 0.5 - kNumBCGrid) * kDz << std::string(",")
		<< static_cast<double>(u[idx3_3D(kNy, kNz, i, j, k)]) << std::endl;

	outF.close();

	return 0;
}