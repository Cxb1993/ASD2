#include "test_poisson2d_large.h"

int test_poisson2D_GuassSeidel() {
	// http://farside.ph.utexas.edu/teaching/329/lectures/node71.html
	const int kNx = 64, kNy = 64, kNumBCGrid = 3;
	const double kBaseX = 0.0, kBaseY = 0.0, kLenX = 1.0, kLenY = 1.0;
	const double kDx = kLenX / kNx, kDy = kLenY / kNy;
	
	std::shared_ptr<PoissonSolver2D> PSolver = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);
	std::shared_ptr<BoundaryCondition2D> PBC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	std::vector<double> u((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> b((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	PBC->SetBC_P_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	PBC->SetBCConstantPW(0.0);
	PBC->SetBCConstantPE(0.0);
	PBC->SetBCConstantPS(0.0);
	PBC->SetBCConstantPN(1.0);

	double x = 0.0, y = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		x = kBaseX + (i - kNumBCGrid) * kDx;
		y = kBaseY + (j - kNumBCGrid) * kDy;
		b[i + j * (kNx + 2 * kNumBCGrid)] = 0.0;
	}
	std::vector<double> AVals, DiagVals, MVals;
	std::vector<MKL_INT> ACols, ARowIdx, MCols, MRowIdx;
	MKL_INT rowIdx = 0, tmpRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
	// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
	// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		AValsDic.clear();
		AColsDic.clear();
		tmpRowIdx = 0;
		// Add starting rowIdx
		ARowIdx.push_back(rowIdx);

		// Set default values, if a current pointer is in interior, it will not be changed.
		AValsDic["S"] = -1.0 / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
		AValsDic["W"] = -1.0 / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["C"] = (1.0 + 1.0) / (kDx * kDx) + (1.0 + 1.0) / (kDy * kDy);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["E"] = -1.0 / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["N"] = -1.0 / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;

		if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] += 1.0 / (kDx * kDx);
			b[i + j * (kNx + 2 * kNumBCGrid)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPW);
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] += 1.0 / (kDy * kDy);
			b[i + j * (kNx + 2 * kNumBCGrid)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPE);
		}

		if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] += 1.0 / (kDy * kDy);
			b[i + j * (kNx + 2 * kNumBCGrid)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPS);
		}

		if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] += 1.0 / (kDx * kDx);
			b[i + j * (kNx + 2 * kNumBCGrid)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPN);
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
	}
	ARowIdx.push_back(rowIdx);
	
	PSolver->GS_2FUniform_2D(u, b, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, PBC, 10000);
	
	std::string fname_p_base("TestPoisson2D_GS");
	std::ofstream outF;
	std::string fname_p(fname_p_base + "_ASCII.plt");

	outF.open(fname_p.c_str(), std::ios::out);

	outF << "TITLE = VEL" << std::endl;
	outF << "VARIABLES = \"X\", \"Y\", \"P\"" << std::endl;
	outF.close();

	outF.open(fname_p.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"1\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", DATAPACKING=POINT")
		<< std::string(", SOLUTIONTIME=") << 0.1
		<< std::string(", STRANDID=") << 2
		<< std::endl;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
		outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
		<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
		<< static_cast<double>(u[i + j * (kNx + 2 * kNumBCGrid)]) << std::endl;

	outF.close();

	return 0;
}

int test_poisson2D_CG() {
	const int kNx = 64, kNy = 64, kNumBCGrid = 3;
	const double kBaseX = 0.0, kBaseY = 0.0, kLenX = 1.0, kLenY = 1.0;
	const double kDx = kLenX / kNx, kDy = kLenY / kNy;

	std::shared_ptr<PoissonSolver2D> PSolver = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);
	std::shared_ptr<BoundaryCondition2D> PBC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	std::vector<double> u((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> b((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	PBC->SetBC_P_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	PBC->SetBCConstantPW(0.0);
	PBC->SetBCConstantPE(0.0);
	PBC->SetBCConstantPS(0.0);
	PBC->SetBCConstantPN(1.0);
	
	double x = 0.0, y = 0.0;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
		x = kBaseX + (i - kNumBCGrid) * kDx;
		y = kBaseY + (j - kNumBCGrid) * kDy;
		b[idx3_2D(kNy, i, j)] = 0.0;
	}

	std::vector<double> AVals, DiagVals, MVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
	// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
	// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
	
	// Set default values, if a current pointer is in interior, it will not be changed.
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
		AValsDic["S"] = -1.0 / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
		AValsDic["W"] = -1.0 / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["C"] = (1.0 + 1.0) / (kDx * kDx) + (1.0 + 1.0) / (kDy * kDy);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["E"] = -1.0 / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["N"] = -1.0 / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;

		if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDx * kDx);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPW);
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDy * kDy);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPE);
		}

		if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDy * kDy);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPS);
		}

		if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDx * kDx);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPN);
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

		rowIdx += tmpRowIdx;
		MRowIdx += tmpMRowIdx;
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);

	PSolver->CG_2FUniformP_2D(u, b, AVals, ACols, ARowIdx, DiagVals, DiagCols, DiagRowIdx, PBC, size, 10000);

	std::string fname_p_base("TestPoisson2D_CG");
	std::ofstream outF;
	std::string fname_p(fname_p_base + "_ASCII.plt");

	outF.open(fname_p.c_str(), std::ios::out);

	outF << "TITLE = VEL" << std::endl;
	outF << "VARIABLES = \"X\", \"Y\", \"P\"" << std::endl;
	outF.close();

	outF.open(fname_p.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"1\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", DATAPACKING=POINT")
		<< std::string(", SOLUTIONTIME=") << 0.0
		<< std::string(", STRANDID=") << 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
		<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
		<< static_cast<double>(u[idx3_2D(kNy, i, j)]) << std::endl;

	outF.close();

	return 0;
}

int test_poisson2D_BiCGStab() {
	const int kNx = 64, kNy = 64, kNumBCGrid = 3;
	const double kBaseX = 0.0, kBaseY = 0.0, kLenX = 1.0, kLenY = 1.0;
	const double kDx = kLenX / kNx, kDy = kLenY / kNy;

	std::shared_ptr<PoissonSolver2D> PSolver = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);
	std::shared_ptr<BoundaryCondition2D> PBC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
	std::vector<double> u((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	std::vector<double> b((kNx + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid), 0.0);
	PBC->SetBC_P_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	PBC->SetBCConstantPW(0.0);
	PBC->SetBCConstantPE(0.0);
	PBC->SetBCConstantPS(0.0);
	PBC->SetBCConstantPN(1.0);

	double x = 0.0, y = 0.0;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
		x = kBaseX + (i - kNumBCGrid) * kDx;
		y = kBaseY + (j - kNumBCGrid) * kDy;
		b[i + j * (kNx + 2 * kNumBCGrid)] = 0.0;
	}

	std::vector<double> AVals, DiagVals;
	std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
	MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
	MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
	MKL_INT size = kNx * kNy;
	std::map<std::string, double> AValsDic;
	std::map<std::string, MKL_INT> AColsDic;

	// AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
	// Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
	// At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
	// Set default values, if a current pointer is in interior, it will not be changed.
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
		AValsDic["S"] = -1.0 / (kDy * kDy);
		AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
		AValsDic["W"] = -1.0 / (kDx * kDx);
		AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["C"] = (1.0 + 1.0) / (kDx * kDx) + (1.0 + 1.0) / (kDy * kDy);
		AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["E"] = -1.0 / (kDx * kDx);
		AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
		AValsDic["N"] = -1.0 / (kDy * kDy);
		AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;

		if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::NEUMANN) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid && PBC->m_BC_PW == BC2D::DIRICHLET) {
			AColsDic["W"] = -1;
			AValsDic["W"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDx * kDx);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPW);
		}

		// East boundary
		if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::NEUMANN) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDx * kDx);
		}
		else if (i == kNumBCGrid + kNx - 1 && PBC->m_BC_PE == BC2D::DIRICHLET) {
			AColsDic["E"] = -1;
			AValsDic["E"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDy * kDy);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDx * kDx) * (PBC->m_BC_DirichletConstantPE);
		}

		if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::NEUMANN) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid && PBC->m_BC_PS == BC2D::DIRICHLET) {
			AColsDic["S"] = -1;
			AValsDic["S"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDy * kDy);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPS);
		}

		if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::NEUMANN) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			AValsDic["C"] -= 1.0 / (kDy * kDy);
		}
		else if (j == kNumBCGrid + kNy - 1 && PBC->m_BC_PN == BC2D::DIRICHLET) {
			AColsDic["N"] = -1;
			AValsDic["N"] = 0.0;
			// AValsDic["C"] += 1.0 / (kDx * kDx);
			b[idx3_2D(kNy, i, j)] -= -1.0 / (kDy * kDy) * (PBC->m_BC_DirichletConstantPN);
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

		rowIdx += tmpRowIdx;
		MRowIdx += tmpMRowIdx;
	}
	ARowIdx.push_back(rowIdx);
	DiagRowIdx.push_back(MRowIdx);

	PSolver->BiCGStab_2FUniform_2D(u, b, AVals, ACols, ARowIdx, DiagVals, DiagCols, DiagRowIdx, PBC, size, 10000);

	std::string fname_p_base("TestPoisson2D_BiCGSTAB");
	std::ofstream outF;
	std::string fname_p(fname_p_base + "_ASCII.plt");

	outF.open(fname_p.c_str(), std::ios::out);

	outF << "TITLE = VEL" << std::endl;
	outF << "VARIABLES = \"X\", \"Y\", \"P\"" << std::endl;
	outF.close();

	outF.open(fname_p.c_str(), std::ios::app);

	outF << std::string("ZONE T=\"1\", I=") << kNx << std::string(", J=") << kNy
		<< std::string(", DATAPACKING=POINT")
		<< std::string(", SOLUTIONTIME=") << 0.0
		<< std::string(", STRANDID=") << 1
		<< std::endl;

	for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
		<< kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
		<< static_cast<double>(u[idx3_2D(kNy, i, j)]) << std::endl;

	outF.close();

	return 0;
}
