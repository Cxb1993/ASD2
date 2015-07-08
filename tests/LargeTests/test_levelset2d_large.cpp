#include "test_levelset2d_large.h"

int LevelSetTest_2D_Simple() {
	// Zalesak Problem
	// int nx = 200, ny = 200;
	// double dx = 0.5, dy = 0.5, dt = 0.5;
	int nx = 100, ny = 100;
	double lenX = 10.0, lenY = 10.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.2;
	double dt = cfl * dx;
	int iter = 0, maxIter = 50;
	double curTime = 0.0, maxTime = 628.1;
	double x = 0.0, y = 0.0;
	double baseX = -5.0, baseY = -5.0;
	int orgMass = 0, LSMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> div((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy);
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 2.0, d = 0.0;
	
	for (int i = NBC3; i < nx + NBC3; i++)
	for (int j = NBC3; j < ny + NBC3; j++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];

		err[idx3_2D(ny, i, j)] = 0.0;
	}

	for (int i = 0; i < nx + 2 * NBC3; i++)
	for (int j = 0; j < ny + 2 * NBC3; j++) {
		u[idx3_2D(ny, i, j)] = 1.0;
		v[idx3_2D(ny, i, j)] = 0.0;

	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	LSolver->ApplyBC_P_2D(ls);
	
	std::ostringstream outfname_stream1;
	outfname_stream1 << "2D_LS_SimpleVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "2D_LS_SimpleDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Solve_LevelSet_2D(ls, u, v, dt);
		LSolver->Reinit_Original_2D(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	double errNormInfty = -std::numeric_limits<double>::max();
	double errNorm1 = 0.0, errNorm2 = 0.0, errMass = 0.0;

	LSMass = 0;
	orgMass = 0;

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		err[idx3_2D(ny, i, j)] = std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		errNormInfty = std::max(std::fabs(err[idx3_2D(ny, i, j)]), errNormInfty);
		errNorm1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		errNorm2 += std::pow(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), 2.0);

		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNorm2 = std::sqrt(errNorm2) / static_cast<double>(nx * ny);

	std::cout << "Simple - err(L-Infty):" << errNormInfty << " err(L1):" << errNorm1 << " err(L2):" << errNorm2 << std::endl;
	std::cout << "Simple - err(Mass loss) : " << errMass << "%" << std::endl;

	std::cout << "LevelSet Simple " << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2D_ReinitOnly_Original() {
	int nx = 128, ny = 128;
	double lenX = 4.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 2;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0;
	double baseX = -2.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2D_ReinitOnlyOriginalVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2D_ReinitOnlyOriginalDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_Original_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;

	std::cout << "Reinit Only(2D, Original) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "LevelSet Reinit Only 2D Original " << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();

	return 0;
}

int LevelSetTest_2D_ReinitOnly_Sussman() {
	int nx = 128, ny = 128;
	double lenX = 4.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 2;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0;
	double baseX = -2.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2D_ReinitOnlySussmanVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2D_ReinitOnlySussmanDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;
	
	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;
		
		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;
	}

	std::cout << "Reinit Only(2D, Sussman) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "LevelSet Reinit Only 2D Sussman " << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();

	return 0;
}

int LevelSetTest_2DAxisym_ReinitOnly_Sussman() {
	int nx = 128, ny = 128;
	double lenX = 2.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 2;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0;
	double baseX = 0.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "dirichlet", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "dirichlet", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2DAxisym_ReinitOnlySussmanVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2DAxisym_ReinitOnlySussmanDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
			for (int j = NBC3; j < ny + NBC3; j++) {
				if (ls[idx3_2D(ny, i, j)] > 0.0)
					LSMass++;
			}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;
	}

	std::cout << "Reinit Only(Axisym, Sussman) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "LevelSet Reinit Only 2DAxisym, Sussman" << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();

	return 0;
}

int LevelSetTest_2D_ReinitOnly_MinRK2() {
	int nx = 128, ny = 128;
	double lenX = 10.0, lenY = 10.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 1;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0, r = 1.0, a = 0.7;
	double baseX = -5.0, baseY = -5.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> kappa((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2D_ReinitOnlyMinRK2Vel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2D_ReinitOnlyMinRK2Div_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	double errNearInterfaceL1 = 0.0, errWholeDomainL1 = 0.0;
	double errNearInterfaceLInfty = -std::numeric_limits<double>::max(), errWholeDomainLInfty = -std::numeric_limits<double>::max();

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;

		if (std::fabs(ls[idx3_2D(ny, i, j)]) < 1.2 * std::min(dx, dy)) {
			errNearInterfaceLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errNearInterfaceL1);
			errNearInterfaceL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}

		if (ls[idx3_2D(ny, i, j)] > -0.8) {
			errWholeDomainLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errWholeDomainL1);
			errWholeDomainL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNearInterfaceL1 *= std::min(dx, dy) * std::min(dx, dy);
	errWholeDomainL1 *= std::min(dx, dy) * std::min(dx, dy);

	std::cout << "Reinit Only(2D, MinRK2) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) Near Interface - err(L1) : " << errNearInterfaceL1 << " err(LInfty) : " << errNearInterfaceLInfty << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) WholeDomain - err(L1) : " << errWholeDomainL1 << " err(LInfty) : " << errWholeDomainLInfty << std::endl;
	std::cout << "LevelSet Reinit Only 2D SigneddistanceFunc, MinRK2" << std::endl;

	OutRes_2D(iter + 1, curTime * 1.1, fname_vel_base, fname_div_base,
		u, v, err, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2DAxisym_ReinitOnly_MinRK2() {
	int nx = 64, ny = 128;
	double lenX = 2.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 5;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0, r = 1.0, a = 0.7;
	double baseX = 0.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> kappa((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2DAxisym_ReinitOnlyMinRK2Vel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2DAxisym_ReinitOnlyMinRK2Div_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	double errNearInterfaceL1 = 0.0, errWholeDomainL1 = 0.0;
	double errNearInterfaceLInfty = -std::numeric_limits<double>::max(), errWholeDomainLInfty = -std::numeric_limits<double>::max();

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;

		if (std::fabs(ls[idx3_2D(ny, i, j)]) < 1.2 * std::min(dx, dy)) {
			errNearInterfaceLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errNearInterfaceL1);
			errNearInterfaceL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}

		if (ls[idx3_2D(ny, i, j)] > -0.8) {
			errWholeDomainLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errWholeDomainL1);
			errWholeDomainL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNearInterfaceL1 *= std::min(dx, dy) * std::min(dx, dy);
	errWholeDomainL1 *= std::min(dx, dy) * std::min(dx, dy);

	std::cout << "Reinit Only(2D Axisym, MinRK2) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "Reinit Only(2D Axisym, MinRK2) Near Interface - err(L1) : " << errNearInterfaceL1 << " err(LInfty) : " << errNearInterfaceLInfty << std::endl;
	std::cout << "Reinit Only(2D Axisym, MinRK2) WholeDomain - err(L1) : " << errWholeDomainL1 << " err(LInfty) : " << errWholeDomainLInfty << std::endl;
	std::cout << "LevelSet Reinit Only 2D Axisym SigneddistanceFunc, MinRK2" << std::endl;

	OutRes_2D(iter + 1, curTime * 1.1, fname_vel_base, fname_div_base,
		u, v, err, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2D_ReinitOnly_MinRK2_2() {
	int nx = 128, ny = 128;
	double lenX = 4.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.45;
	double dt = cfl * dx;
	int iter = 0, maxIter = 1;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0, r = 1.0, a = 0.7;
	double baseX = -2.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0, errMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;

		ls[idx3_2D(ny, i, j)] = (std::pow(x - 1.0, 2.0) + std::pow(y - 1.0, 2.0) + 0.1)
			* (std::sqrt(x * x + y * y) - 1.0);

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = std::sqrt(x * x + y * y) - 1.0;
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2D_ReinitOnlyMinRK2_2Vel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2D_ReinitOnlyMinRK2_2Div_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}
		std::cout << iter << std::endl;
	}

	LSMass = 0;
	orgMass = 0;
	double errNearInterfaceL1 = 0.0, errWholeDomainL1 = 0.0;
	double errNearInterfaceLInfty = -std::numeric_limits<double>::max(), errWholeDomainLInfty = -std::numeric_limits<double>::max();
	
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;

		if (std::fabs(ls[idx3_2D(ny, i, j)]) < 1.2 * std::min(dx, dy)) {
			errNearInterfaceLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errNearInterfaceL1);
			errNearInterfaceL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}

		if (ls[idx3_2D(ny, i, j)] > -0.8) {
			errWholeDomainLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errWholeDomainL1);
			errWholeDomainL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}
	}
	
	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNearInterfaceL1 *= std::min(dx, dy) * std::min(dx, dy);	
	errWholeDomainL1 *= std::min(dx, dy) * std::min(dx, dy);
	
	std::cout << "Reinit Only(2D, MinRK2) - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) Near Interface - err(L1) : " << errNearInterfaceL1 << " err(LInfty) : " << errNearInterfaceLInfty << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) WholeDomain - err(L1) : " << errWholeDomainL1 << " err(LInfty) : " << errWholeDomainLInfty << std::endl;
	std::cout << "LevelSet Reinit Only 2D NonsigneddistanceFunc, MinRK2" << std::endl;

	OutRes_2D(iter + 1, curTime * 1.1, fname_vel_base, fname_div_base,
		u, v, err, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2D_ReinitOnlyWithKinks() {
	int nx = 128, ny = 128;
	double lenX = 4.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 1;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0, r = 1.0, a = 0.7;
	double baseX = -2.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 80);
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		
		if (((a - x) / std::sqrt(std::pow(a - x, 2.0) + y * y) >= a / r) &&
			((a + x) / std::sqrt(std::pow(a + x, 2.0) + y * y) >= a / r))
			lsInit[idx3_2D(ny, i, j)] = std::min(
			std::sqrt(x * x + std::pow(y + std::sqrt(r * r - a * a), 2.0)), 
			std::sqrt(x * x + std::pow(y - std::sqrt(r * r - a * a), 2.0)));
		else
			lsInit[idx3_2D(ny, i, j)] = std::min(
			std::sqrt(std::pow(x + a, 2.0) + y * y) - r,
			std::sqrt(std::pow(x - a, 2.0) + y * y) - r);

		ls[idx3_2D(ny, i, j)] = std::min(
			std::sqrt(std::pow(x + a, 2.0) + y * y) - r,
			std::sqrt(std::pow(x - a, 2.0) + y * y) - r);

		if (lsInit[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		err[idx3_2D(ny, i, j)] = 0.0;
	}
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;

		ls[idx3_2D(ny, i, j)] *= std::pow(x - 1, 2.0) + std::pow(y - 1.0, 2.0) + 0.1;
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "testLevelSet2D_ReinitOnlyWithKinksVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testLevelSet2D_ReinitOnlyWithKinksDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;

		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		MSolver->UpdateKappa(ls);

		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}

		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	LSMass = 0,	orgMass = 0; 
	double errMass = 0.0;
	double errNearInterfaceL1 = 0.0, errWholeDomainL1 = 0.0;
	double errNearInterfaceLInfty = -std::numeric_limits<double>::max(), errWholeDomainLInfty = -std::numeric_limits<double>::max();

	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;

		if (std::fabs(ls[idx3_2D(ny, i, j)]) < 1.2 * std::min(dx, dy)) {
			errNearInterfaceLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errNearInterfaceL1);
			errNearInterfaceL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}

		if (std::fabs(ls[idx3_2D(ny, i, j)]) >= 0.1) {
			errWholeDomainLInfty = std::max(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), errWholeDomainL1);
			errWholeDomainL1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		}
	}
	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNearInterfaceLInfty /= nx * ny;
	errNearInterfaceL1 /= nx * ny;
	errWholeDomainLInfty /= nx * ny;
	errWholeDomainL1 /= nx * ny;
	
	std::cout << "Reinit Only - err(Mass loss) : " << errMass << "%" << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) Near Interface - err(L1) : " << errNearInterfaceL1 << " err(LInfty) : " << errNearInterfaceLInfty << std::endl;
	std::cout << "Reinit Only(2D, MinRK2) WholeDomain - err(L1) : " << errWholeDomainL1 << " err(LInfty) : " << errWholeDomainLInfty << std::endl;
	std::cout << "LevelSet Reinit Only 2D with Kinks, MinRK2" << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2D_FirstTimeReinit() {
	int nx = 128, ny = 128;
	double lenX = 4.0, lenY = 4.0;
	double dx = lenX / nx, dy = lenY / ny, cfl = 0.3;
	double dt = cfl * dx;
	int iter = 0, maxIter = 10;
	double curTime = 0.0, maxTime = 10.0;
	double x = 0.0, y = 0.0;
	double baseX = -2.0, baseY = -2.0;
	int orgMass = 0, LSMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy, 2 * std::max(nx, ny));
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 1.0, d = 0.0;
	
	for (int j = NBC3; j < ny + NBC3; j++)
	for (int i = NBC3; i < nx + NBC3; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		d = std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0)) - radius;

		if (d < 0)
			ls[idx3_2D(ny, i, j)] = 1;
		else if (d > 0)
			ls[idx3_2D(ny, i, j)] = -1;
		else
			ls[idx3_2D(ny, i, j)] = 0.0;
		//	ls[idx3_2D(ny, i, j)] = -d;

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];

		err[idx3_2D(ny, i, j)] = 0.0;
	}

	for (int j = 0; j < ny + 2 * NBC3; j++)
	for (int i = 0; i < nx + 2 * NBC3; i++) {
		u[idx3_2D(ny, i, j)] = 0.0;
		v[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->SetBC_U_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_V_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_U_2D(u);
	LSolver->ApplyBC_V_2D(v);
	LSolver->ApplyBC_P_2D(ls);

	LSolver->FirstTimeOnlyReinit_Sussman_2D(ls);


	double rhoH = 1.0, rhoL = 1.0, muH = 1.0, muL = 1.0, gConstant = 0.0, sigma = 0.0, L = 1.0, U = 1.0;
	TIMEORDERENUM timeOrder = TIMEORDERENUM::EULER;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	double maxtime = 2.0;
	int maxiter = 20, niterskip = 1, num_bc_grid = NBC3;
	bool writeVTK = false;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);
	MSolver->AllocateVariables();
	MSolver->UpdateKappa(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "2D_LS_FirstTimeReinitVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "2D_LS_FirstTimeReinitDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, MSolver->m_kappa, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	LSMass = 0;
	orgMass = 0;

	std::cout << "LevelSet Reinit Only " << std::endl;

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}

int LevelSetTest_2D_Sussman621_ReinitSussman() {
	// Zalesak Problem
	// int nx = 200, ny = 200;
	// double dx = 0.5, dy = 0.5, dt = 0.5;
	int nx = 100, ny = 100;
	double dx = 1.0, dy = 1.0, dt = 1.0;
	int iter = 0, maxIter = 630;
	double curTime = 0.0, maxTime = 628.1;
	double x = 0.0, y = 0.0;
	const double baseX = 0.0, baseY = 0.0;
	int orgMass = 0, LSMass = 0;

	std::vector<double> lsInit((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> div((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, baseX, baseY, dx, dy);
	
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	const double radius = 15.0, slotW = 15.0;
	double d = 0.0;
	
	const double centerX = 50.0, centerY = 75.0;

	for (int i = 0; i < nx + 2 * NBC3; i++)
	for (int j = 0; j < ny + 2 * NBC3; j++) {
		x = baseX + (i + 0.5 - NBC3) * dx;
		y = baseY + (j + 0.5 - NBC3) * dy;
		
		ls[idx3_2D(ny, i, j)] = -2.0;
		u[idx3_2D(ny, i, j)] = (M_PI / 314.0) * (50.0 - y);
		v[idx3_2D(ny, i, j)] = (M_PI / 314.0) * (x - 50.0);

		lsInit[idx3_2D(ny, i, j)] = ls[idx3_2D(ny, i, j)];

		err[idx3_2D(ny, i, j)] = 0.0;
	}

	for (int i = NBC3; i < nx + NBC3; i++)
	for (int j = NBC3; j < ny + NBC3; j++) {
		// positive : inside, negative : outside
		x = baseX + (i - NBC3) * dx;
		y = baseY + (j - NBC3) * dy;
		d = std::sqrt(std::pow(x - centerX, 2.0) + std::pow(y - centerY, 2.0)) - radius;

		ls[idx3_2D(ny, i, j)] = -2.0;
		if (d < 0)
			ls[idx3_2D(ny, i, j)] = 2.0;
		else if (d == 0)
			ls[idx3_2D(ny, i, j)] = 0.0;
		else
			ls[idx3_2D(ny, i, j)] = -2.0;

		if (x >= (50 - 0.25 * slotW) && x <= (50 + 0.25 * slotW)) {
			if (y <= 75 + 7.5)
				ls[idx3_2D(ny, i, j)] = -2.0;
		}

		if (ls[idx3_2D(ny, i, j)] > 0)
			orgMass++;
	}

	LSolver->Reinit_Original_2D(ls);
	LSolver->ApplyBC_P_2D(ls);

	std::ostringstream outfname_stream1;
	outfname_stream1 << "2D_LS_Sussman621_SussmanReinitVel_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel_base = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "2D_LS_Sussman621_SussmanReinitDiv_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div_base = outfname_stream2.str();
	PLTTYPE m_PLTType = PLTTYPE::BOTH;

	OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);
	
	while (curTime < maxTime && iter < maxIter) {
		iter++;
		curTime += dt;
		
		LSolver->Solve_LevelSet_2D(ls, u, v, dt);
		LSolver->ApplyBC_P_2D(ls);
		LSolver->Reinit_Original_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		
		OutRes_2D(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);
		
		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (ls[idx3_2D(ny, i, j)] > 0.0)
				LSMass++;
		}
			
		massVec.push_back(LSMass * dx * dy);
		std::cout << iter << std::endl;
	}

	double errNormInfty = -std::numeric_limits<double>::max();
	double errNorm1 = 0.0, errNorm2 = 0.0, errMass = 0.0;

	LSMass = 0;
	orgMass = 0;

	for (int i = NBC3; i < nx + NBC3; i++)
	for (int j = NBC3; j < ny + NBC3; j++) {
		err[idx3_2D(ny, i, j)] = std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		errNormInfty = std::max(std::fabs(err[idx3_2D(ny, i, j)]), errNormInfty);
		errNorm1 += std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]);
		errNorm2 += std::pow(std::fabs(ls[idx3_2D(ny, i, j)] - lsInit[idx3_2D(ny, i, j)]), 2.0);

		if (lsInit[idx3_2D(ny, i, j)] > 0.0)
			orgMass++;

		if (ls[idx3_2D(ny, i, j)] > 0.0)
			LSMass++;
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNorm2 = std::sqrt(errNorm2) / static_cast<double>(nx * ny);

	std::cout << "Sussman 6-2-1 - err(L-Infty):" << errNormInfty << " err(L1):" << errNorm1 << " err(L2):" << errNorm2 << std::endl;
	std::cout << "Sussman 6-2-1 - err(Mass loss) : " << errMass << "%" << std::endl;

	std::cout << "LevelSet Zalesak Disc : Sussman 6-2-1 (2D) " << std::endl;
	
	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH)
		OutResClose();
	return 0;
}
