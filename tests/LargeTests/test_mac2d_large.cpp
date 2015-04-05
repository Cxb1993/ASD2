#include "test_mac2d_large.h"

int MAC2DTest_CavityFlow() {
	// Set initial level set
	const double Re = 100.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.0, muI = 0.01, sigma = 0.0;
	const double densityRatio = 1.0, viscosityRatio = 1.0;
	const double gConstant = 0.0;
	// # of cells
	const int nx = 128, ny = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, lenX = 1.0, lenY = 1.0, cfl = 0.1;
	double x = 0.0, y = 0.0, d = 0.0;

	const int maxtime = 10.0, maxiter = 10000, niterskip = 50, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_CavityVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_CavityDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::RK3;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr,
		L, U, sigma, densityRatio, viscosityRatio, rhoI, muI,
		nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip,
		num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("dirichlet", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(1.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);

	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3(ny, i, j)] = 0.0;
		MSolver->m_v[idx3(ny, i, j)] = 0.0;
		MSolver->m_p[idx3(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3(ny, i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);

	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;

		ls[idx3(ny, i, j)] = 1.0;
	}
	LSolver->ApplyBC_P_2D(ls);

	LSolver->m_signedInitLS = LSolver->GetSignedLSNormalized(ls);
	// LSolver->Reinit_Sussman_2D(ls);

	LSolver->ApplyBC_P_2D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dx, dy) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid)), vhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	std::vector<double> div((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		for (int i = 0; i < nx + 2 * num_bc_grid; i++)
		for (int j = 0; j < ny + 2 * num_bc_grid; j++)
			lsB[idx3(ny, i, j)] = ls[idx3(ny, i, j)];

		// lsB = ls;
		// LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		// LSolver->Reinit_Sussman_2D(ls);
		// LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part

		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, lsB);
		
		// Get intermediate velocity with RK method
		MSolver->GetIntermediateVel(LSolver, lsB, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, lsB, ls, uhat, vhat, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;
		MSolver->m_iter++;
	
		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			double rnorm2 = 0.0;
			std::cout << "Cavity : " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << rnorm2 << std::endl;
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				uhat, vhat, MSolver->m_ps, ls);
			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_SmallAirBubbleRising() {
	// Set initial level set
	const double Re = 1.0, We = 1.0, Fr = 1.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.226, muI = 1.78e-5;
	const double rhoO = 1000, muO = 1.137e-3, sigma = 0.0728;
	const double gConstant = -9.81;
	// # of cells
	const int nx = 80, ny = 120;
	// related to initialize level set
	const double baseX = -0.01, baseY = -0.01, lenX = 0.02, lenY = 0.03, cfl = 0.1;
	double radius = 1.0 / 300.0, x = 0.0, y = 0.0, d = 0.0;
	// const double baseX = -1.0, baseY = -1.0, lenX = 2.0, lenY = 3.0, cfl = 0.5;
	// double radius = 1.0 / 3.0, x = 0.0, y = 0.0, d = 0.0;

	const int maxtime = 2.0, maxiter = 5, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_BubbleRisingVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_BubbleRisingDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::RK1;
	int stat = 0;
	
	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoI, rhoO, muI, muO, gConstant,
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
	
	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3(ny, i, j)] = 0.0;
		MSolver->m_v[idx3(ny, i, j)] = 0.0;
		MSolver->m_p[idx3(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3(ny, i, j)] = 0.0;
	}
	
	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	
	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;
		
		d = std::sqrt(x * x + y * y) - radius;
		
		ls[idx3(ny, i, j)] = -d;
	}
	LSolver->ApplyBC_P_2D(ls);

	LSolver->m_signedInitLS = LSolver->GetSignedLSNormalized(ls);
	LSolver->Reinit_Sussman_2D(ls);

	LSolver->ApplyBC_P_2D(ls);
	
	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dx, dy) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid)), vhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	std::vector<double> div((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	
	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		for (int i = 0; i < nx + 2 * num_bc_grid; i++)
		for (int j = 0; j < ny + 2 * num_bc_grid; j++)
			lsB[idx3(ny, i, j)] = ls[idx3(ny, i, j)];
		
		// lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		
		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, lsB);
		
		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, lsB, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);
		
		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, lsB, ls, uhat, vhat, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		
		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls);
		
		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;
		MSolver->m_iter++;

		
		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			double rnorm2 = 0.0;
			std::cout << "Bubble : " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << rnorm2 << std::endl;
			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	uhat, vhat, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_TaylorInstability() {
	// Set initial level set
	const double L = 1.0, U = 1.0;
	const double rhoI = 1, muI = 1.78e-5;
	const double Re = 1.0, We = 1.0, Fr = 1.0;
	const double rhoO = 3.0, muO = 1.137e-3, sigma = 0.0728;
	const double densityRatio = 3.0, viscosityRatio = 1.0;
	const double gConstant = -2.0;
	// # of cells
	const int nx = 128, ny = 512;
	// related to initialize level set
	const double d = 0.5;
	const double baseX = -d * 0.5, baseY = -d * 2.0, lenX = d, lenY = 4.0 * d, cfl = 0.1;
	double radius = 1.0 / 300.0, x = 0.0, y = 0.0;
	
	const int maxtime = 2.0, maxiter = 20, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("TaylorVel_Re_" + std::to_string(Re));
	const std::string fname_div("TaylorDiv_Re_" + std::to_string(Re));
	const TimeOrderEnum timeOrder = TimeOrderEnum::RK1;
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr,
		L, U, sigma, densityRatio, viscosityRatio, rhoI, muI,
		nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("periodic", "periodic", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("periodic", "periodic", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("periodic", "periodic", "dirichlet", "dirichlet");
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetBCConstantPS(0.0);
	MSolver->SetBCConstantPN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);

	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3(ny, i, j)] = 0.0;
		MSolver->m_v[idx3(ny, i, j)] = 0.0;
		MSolver->m_p[idx3(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3(ny, i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);

	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;

		if (x > -0.1 * d * cos(2.0 * M_PI * x / d))
			ls[idx3(ny, i, j)] = 2.0;
		else if (x < -0.1 * d * cos(2.0 * M_PI * x / d))
			ls[idx3(ny, i, j)] = -2.0;
		else
			ls[idx3(ny, i, j)] = 0.0;
	}
	LSolver->ApplyBC_P_2D(ls);

	LSolver->m_signedInitLS = LSolver->GetSignedLSNormalized(ls);
	LSolver->Reinit_Sussman_2D(ls);

	LSolver->ApplyBC_P_2D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dx, dy) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid)), vhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	std::vector<double> div((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		for (int i = 0; i < nx + 2 * num_bc_grid; i++)
			for (int j = 0; j < ny + 2 * num_bc_grid; j++)
				lsB[idx3(ny, i, j)] = ls[idx3(ny, i, j)];
		// lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part

		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, lsB);
		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, lsB, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, lsB, ls, uhat, vhat, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;
		MSolver->m_iter++;

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			double rnorm2 = 0.0;
			std::cout << "Taylor : " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << rnorm2 << std::endl;
			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	uhat, vhat, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}