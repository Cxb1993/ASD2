#include "test_mac2d_large.h"

int MAC2DTest_CavityFlow() {
	// Set initial level set
	const double Re = 100.0, We = 0.0, FrX = 0.0, FrY = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.0, muI = 0.01, sigma = 0.0;
	const double densityRatio = 1.0, viscosityRatio = 1.0;
	const double gConstant = 0.0;
	// # of cells
	const int nx = 128, ny = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, lenX = 1.0, lenY = 1.0, cfl = 0.3;
	double x = 0.0, y = 0.0, d = 0.0;

	const double maxtime = 10.0;
	const int maxiter = 10000, niterskip = 100, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_CavityVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_CavityDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::RK3;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, FrX, FrY,
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
	lsB = ls;
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
		
		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, ls);
		
		// Get intermediate velocity with RK method
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(div);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, uhat, vhat, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v, uhat, vhat, MSolver->m_ps, ls, lsB);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;
		MSolver->m_iter++;
	
		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Cavity : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}
	// http://stackoverflow.com/a/12836048/743078
	std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

	std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
	std::time_t t = s.count();
	std::size_t fractional_seconds = ms.count() % 1000;

	std::cout << "(Final) Cavity : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_StationaryBubble() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.226, muI = 1.78e-5;
	// const double rhoI = 1000, muI = 1.137e-3;
	const double rhoO = 1000, muO = 1.137e-3, sigma = 0.0728;
	// const double rhoO = 1000, muO = 1.137e-3, sigma = 0.0;
	const double gConstantX = 0.0, gConstantY = 0.0;
	// # of cells
	const int nx = 128, ny = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, lenX = 0.04, lenY = 0.04, cfl = 0.1;
	double radius = 0.01, x = 0.0, y = 0.0, d = 0.0;

	const int maxiter = 10, niterskip = 1, num_bc_grid = 3;
	const double maxtime = 0.06;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_StationaryBubbleVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_StationaryBubbleDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::EULER;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoI, rhoO, muI, muO, gConstantX, gConstantY,
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

		d = std::sqrt(std::pow(x - 0.02, 2.0) + std::pow(y - 0.02, 2.0)) - radius;

		ls[idx3(ny, i, j)] = -d;
	}
	LSolver->ApplyBC_P_2D(ls);

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
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Sussman_2D(ls);
		
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		// MSolver->UpdateKappa(ls);

		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Stationary Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	uhat, vhat, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}
	// http://stackoverflow.com/a/12836048/743078
	std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

	std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
	std::time_t t = s.count();
	std::size_t fractional_seconds = ms.count() % 1000;

	std::cout << "(Final) Small Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_SmallAirBubbleRising() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.226, muI = 1.78e-5;
	const double rhoO = 1000, muO = 1.137e-3, sigma = 0.0728;
	const double gConstantX = 0.0, gConstantY = 9.81;
	// # of cells
	const int nx = 80, ny = 120;
	// related to initialize level set
	const double baseX = -0.01, baseY = -0.01, lenX = 0.02, lenY = 0.03, cfl = 0.1;
	double radius = 1.0 / 300.0, x = 0.0, y = 0.0, d = 0.0;
	// const double baseX = -1.0, baseY = -1.0, lenX = 2.0, lenY = 3.0, cfl = 0.5;
	// double radius = 1.0 / 3.0, x = 0.0, y = 0.0, d = 0.0;

	const double maxtime = 2.0;
	const int maxiter = 10, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_SmallBubbleRisingVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_SmallBubbleRisingDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::EULER;
	int stat = 0;
	
	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoI, rhoO, muI, muO, gConstantX, gConstantY,
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
	LSolver->Reinit_Original_2D(ls);

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
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Original_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		
		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, ls);
		
		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);
		
		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		
		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB);
		
		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;
		
		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());
			
			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;
			
			std::cout << "Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	uhat, vhat, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}
	// http://stackoverflow.com/a/12836048/743078
	std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

	std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
	std::time_t t = s.count();
	std::size_t fractional_seconds = ms.count() % 1000;

	std::cout << "(Final) Small Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_LargeAirBubbleRising() {
	// Set initial level set
	const double Re = 1.0, We = 1.0, Fr = 1.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.226, muI = 1.78e-5;
	const double rhoO = 1000, muO = 1.137e-3, sigma = 0.0728;
	const double gConstantX = 0.0, gConstantY = 9.81;
	// # of cells
	const int nx = 80, ny = 120;
	// related to initialize level set
	const double baseX = -1.0, baseY = -1.0, lenX = 2.0, lenY = 3.0, cfl = 0.5;
	double radius = 1.0 / 3.0, x = 0.0, y = 0.0, d = 0.0;

	const double maxtime = 2.0;
	const int maxiter = 5, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("testMAC2D_LargeBubbleRisingVel_Re_" + std::to_string(Re));
	const std::string fname_div("testMAC2D_LargeBubbleRisingDiv_Re_" + std::to_string(Re));
	const int iterskip = 1;
	const TimeOrderEnum timeOrder = TimeOrderEnum::RK2;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(rhoI, rhoO, muI, muO, gConstantX, gConstantY,
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
	
	LSolver->Reinit_Original_2D(ls);
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
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Original_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part

		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	uhat, vhat, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}
	// http://stackoverflow.com/a/12836048/743078
	std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

	std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
	std::time_t t = s.count();
	std::size_t fractional_seconds = ms.count() % 1000;

	std::cout << "(Final) Large Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DTest_TaylorInstability() {
	// Set initial level set
	const double gConstantX = 0.0, gConstantY = 200;
	// Length Scale = d, Time Scale = \sqrt(d / g), Vel Scale = Length / Time
	const double L = 1.0, U = L / (std::sqrt(L / gConstantY));
	// Froude Number : u0 / (sqrt(g0 * l0))
	const double Re = 1000.0, We = 0.0, FrX = 0.0, FrY = U / std::sqrt(gConstantY * L);
	const double rhoI = 1.0, muI = rhoI / Re * std::sqrt(L * L * L * 9.8);
	const double rhoO = 3.0, muO = rhoO / Re * std::sqrt(L * L * L * 9.8),
		 sigma = 0.0;
	const double densityRatio = 3.0, viscosityRatio = muO / muI;
	
	// # of cells
	const int nx = 128, ny = 512;
	// related to initialize level set
	const double d = 0.5;
	const double baseX = -0.5, baseY = -2.0, lenX = 1.0, lenY = 4.0, cfl = 0.2;
	double radius = 1.0 / 300.0, x = 0.0, y = 0.0;
	
	const double maxtime = 3.0;
	const int maxiter = 10, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("TaylorVel_Re_" + std::to_string(Re));
	const std::string fname_div("TaylorDiv_Re_" + std::to_string(Re));
	const TimeOrderEnum timeOrder = TimeOrderEnum::EULER;
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, FrX, FrY, 
		L, U, sigma, densityRatio, viscosityRatio, rhoI, muI,
		nx, ny, baseX, baseY, lenX, lenY,
		timeOrder, cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("neumann", "neumann", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("neumann", "neumann", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
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

	double distance = 0.0;
	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside, negative : outside
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;
		
		ls[idx3(ny, i, j)] = -y - 0.1 * cos(2.0 * M_PI * x);
	}
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
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);

		// Update F and apply time integration
		MSolver->UpdateJumpCond(MSolver->m_u, MSolver->m_v, ls);
		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		
		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		
		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Taylor : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
		}
		std::fill(uhat.begin(), uhat.end(), 0.0);
		std::fill(vhat.begin(), vhat.end(), 0.0);
	}
	// http://stackoverflow.com/a/12836048/743078
	std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
	std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

	std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
	std::time_t t = s.count();
	std::size_t fractional_seconds = ms.count() % 1000;

	std::cout << "(Final) Taylor : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}