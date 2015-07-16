#include "test_mac2d_large.h"

int MAC2DTest_CavityFlow() {
	// Set initial level set
	const double Re = 100.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoH = 1.0, muH = 0.01, sigma = 0.0;
	const double densityRatio = 1.0, viscosityRatio = 1.0;
	const double gConstant = 0.0;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
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
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_CavityVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_CavityDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK3;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr, GAxis,
		L, U, sigma, densityRatio, viscosityRatio, rhoH, muH,
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
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
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

		ls[idx3_2D(ny, i, j)] = 1.0;
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
		H = MSolver->UpdateHeavisideFunc(ls);
		
		// Get intermediate velocity with RK method
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, H);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(div);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, uhat, vhat, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v, uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

int MAC2DTest_NonSurfaceTension() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.0, muL = 0.1, sigma = 0.0;
	const double rhoH = 1.0, muH = 0.1;
	const double gConstant = 0.0;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	// # of cells
	const int nx = 128, ny = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, lenX = 1.0, lenY = 1.0, cfl = 0.3;
	double radius = 0.25, x = 0.0, y = 0.0, d = 0.0;

	const int maxiter = 10000, niterskip = 100, num_bc_grid = 3;
	const double maxtime = 4.0;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_NonSurfaceTensionVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_NonSurfaceTensionDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK3;
	int stat = 0;

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
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0), HSmooth((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

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
		// positive : inside & gas, negative : outside & liquid 
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(x - 0.5, 2.0) + std::pow(y - 0.5, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_2D(ny, i, j)] = d;
	}

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
	}

	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);

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
		LSolver->ApplyBC_P_2D(ls);
		LSolver->Reinit_MinRK2_2D(ls);

		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		LSolver->ApplyBC_P_2D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

			std::cout << "Non Surface Tension : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

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

	std::cout << "(Final) Non Surface Tension : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
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
	const double rhoL = 1.226, muL = 1.78e-5;
	// const double rhoH = 1000, muH = 1.137e-3;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	// const double rhoL = 1000, muL = 1.137e-3, sigma = 0.0;
	const double gConstant = 0.0; 
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	// # of cells
	const int nx = 128, ny = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, lenX = 0.04, lenY = 0.04, cfl = 0.3;
	double radius = 0.01, x = 0.0, y = 0.0, d = 0.0;

	const int maxiter = 20, niterskip = 1, num_bc_grid = 3;
	const double maxtime = 5.0;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_StationaryBubbleVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_StationaryBubbleDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK3;
	int stat = 0;

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

	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0), HSmooth((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
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
		// positive : inside & gas, negative : outside & liquid 
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(x - 0.02, 2.0) + std::pow(y - 0.02, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_2D(ny, i, j)] = d;
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);
	
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
		LSolver->ApplyBC_P_2D(ls);
		LSolver->Reinit_MinRK2_2D(ls);
		
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		LSolver->ApplyBC_P_2D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

	std::cout << "(Final) Stationary Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
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
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	const double gConstant = 9.81;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	// # of cells
	const int nx = 80, ny = 120;
	// related to initialize level set
	const double baseX = -0.01, baseY = -0.01, lenX = 0.02, lenY = 0.03, cfl = 0.3;
	double radius = 1.0 / 300.0, x = 0.0, y = 0.0, d = 0.0;
	
	const double maxtime = 2.0;
	const int maxiter = 1000, niterskip = 50, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_SmallBubbleRisingVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_SmallBubbleRisingDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK3;
	int stat = 0;	
	
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
	
	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative
	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0), HSmooth((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
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
		
		ls[idx3_2D(ny, i, j)] = d;
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);

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
		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		
		// Update F and apply time integration
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);
		
		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		
		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);
		
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
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	const double gConstant = 9.81;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
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
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_LargeBubbleRisingVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_LargeBubbleRisingDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK2;
	int stat = 0;

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

	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative
	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0), HSmooth((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
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

		ls[idx3_2D(ny, i, j)] = d;
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

		// Update F and apply time integration
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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
	const double gConstant = 200;
	GAXISENUM2D GAxis = GAXISENUM2D::Y;
	// Length Scale = d, Time Scale = \sqrt(d / g), Vel Scale = Length / Time
	const double L = 1.0, U = L / (std::sqrt(L / gConstant));
	// Froude Number : u0 / (sqrt(g0 * l0))
	const double Re = 1000.0, We = 0.0, Fr = U / std::sqrt(gConstant * L);
	const double rhoH = 3.0, muH = rhoH / Re * std::sqrt(L * L * L * 9.8);
	const double rhoL = 1.0, muL = rhoL / Re * std::sqrt(L * L * L * 9.8), sigma = 0.0;
	const double densityRatio = rhoL / rhoH, viscosityRatio = muL / muH;
	
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
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2DAxisym_TaylorVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2DAxisym_TaylorDiv_Re_"
		<< std::to_string(Re)
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str();
	std::string fname_div = outfname_stream2.str();
	const TIMEORDERENUM timeOrder = TIMEORDERENUM::RK3;
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr, GAxis,
		L, U, sigma, densityRatio, viscosityRatio, rhoH, muH,
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
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, baseX, baseY, dx, dy);
	// \phi^n
	std::vector<double> lsB((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative
	std::vector<double> H((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0), HSmooth((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid), 0.0);

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(ny, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(ny, i, j)] = 0.0;
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
		
		ls[idx3_2D(ny, i, j)] = -y - 0.1 * cos(2.0 * M_PI * x);
	}
	LSolver->Reinit_MinRK2_2D(ls);
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
		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Update F and apply time integration
		// Get intermediate velocity
		MSolver->GetIntermediateVel(LSolver, ls, MSolver->m_u, MSolver->m_v, uhat, vhat, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);
		
		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		
		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

int MAC2DAxisymTest_NonSurfaceTension() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1.226, muH = 1.78e-5, sigma = 0.0;
	const double gConstant = 0.0;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;
	// # of cells
	const int nr = 64, nz = 128;
	// related to initialize level set
	const double baseR = 0.0, baseZ = 0.0, lenR = 0.5, lenZ = 1.0, cfl = 0.5;

	const int maxiter = 10, niterskip = 1, num_bc_grid = 3;
	const double maxtime = 1.0;
	const bool writeVTK = false;
	// length of each cell
	const double dr = lenR / nr, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2DAxisym_NonSurfaceTensionVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2DAxisym_NonSurfaceTensionDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	double ambientPressure = 1.0;
	int stat = 0;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nr, nz, baseR, baseZ, lenR, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("axisym", "wall", "inlet", "outlet");
	MSolver->SetBC_V_2D("axisym", "wall", "inlet", "outlet");
	MSolver->SetBC_P_2D("neumann", "wall", "neumann", "pressure");
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(1.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetAmbientPressure(ambientPressure);
	MSolver->SetPLTType(PLTTYPE::BOTH);

	MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nr, nz, num_bc_grid, baseR, baseZ, dr, dz);
	// \phi^n
	std::vector<double> lsB((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0),
		HSmooth((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	LSolver->SetBC_P_2D("axisym", "lsfree", "lsfree", "lsfree");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	double radius = 0.2, r = 0.0, z = 0.0, d = 0.0;
	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		// positive : inside & gas, negative : outside & liquid 
		r = baseR + (i + 0.5 - num_bc_grid) * dr;
		z = baseZ + (j + 0.5 - num_bc_grid) * dz;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(r, 2.0) + std::pow(z - lenZ * 0.5, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_2D(nz, i, j)] = d;
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);
	
	LSolver->ApplyBC_P_2D(ls);

	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(nz, i, j)] = ambientPressure;
		MSolver->m_ps[idx3_2D(nz, i, j)] = ambientPressure;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dr, dz) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		vhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> rhsU((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		rhsV((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> div((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->ApplyBC_P_2D(ls);
		LSolver->Reinit_MinRK2_2D(ls);

		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		rhsU = MSolver->UpdateFU(ls, MSolver->m_u, MSolver->m_v, HSmooth);
		uhat = MSolver->GetUHat(ls, MSolver->m_u, rhsU, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);

		rhsV = MSolver->UpdateFV(ls, MSolver->m_u, uhat, MSolver->m_v, HSmooth);
		vhat = MSolver->GetVHat(ls, MSolver->m_v, rhsV, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence4Poisson(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		LSolver->ApplyBC_P_2D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);

		// MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_curTime += MSolver->m_dt;

		std::ofstream outF;
		std::string fname("RR_ASCII.plt");
		if (MSolver->m_iter == 1) {
			outF.open(fname.c_str(), std::ios::out);

			outF << "TITLE = VEL" << std::endl;
			outF << "VARIABLES = \"R\", \"Z\",\"RHSU\", \"UHAT\", \"RHSV\", \"VHAT\", \"HATDIV\", \"HatRealDiv\"" << std::endl;
			outF.close();
		}
		if (MSolver->m_iter % MSolver->kNIterSkip == 0) {
			std::vector<double> realDiv((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));

			realDiv = MSolver->GetDivergence(uhat, vhat);

			outF.open(fname.c_str(), std::ios::app);

			outF << std::string("ZONE T=\"") << MSolver->m_iter
				<< std::string("\", I=") << nr + 6 << std::string(", J=") << nz + 6
				<< std::string(", SOLUTIONTIME=") << MSolver->m_iter * 0.1
				<< std::string(", STRANDID=") << MSolver->m_iter + 1
				<< std::endl;

			// for (int j = kNumBCGrid; j < kNz + kNumBCGrid; j++)
			// for (int i = kNumBCGrid; i < kNr + kNumBCGrid; i++)
			for (int j = 0; j < nz + 2 * num_bc_grid; j++)
			for (int i = 0; i < nr + 2 * num_bc_grid; i++)
				outF << baseR + static_cast<double>(i + 0.5 - num_bc_grid) * dr << std::string(",")
					<< baseZ + static_cast<double>(j + 0.5 - num_bc_grid) * dz << std::string(",")
					<< static_cast<double>(rhsU[idx3_2D(nz, i, j)]) << std::string(",")
					<< static_cast<double>(uhat[idx3_2D(nz, i, j)]) << std::string(",")
					<< static_cast<double>(rhsV[idx3_2D(nz, i, j)]) << std::string(",")
					<< static_cast<double>(vhat[idx3_2D(nz, i, j)]) << std::string(",")
					<< static_cast<double>(div[idx3_2D(nz, i, j)])	<< std::string(",")
					<< static_cast<double>(realDiv[idx3_2D(nz, i, j)]) << std::endl;
			outF.close();
		}

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Non Surface Tension : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

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

	std::cout << "(Final) Non Surface Tension : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DAxisymTest_StationaryBubble() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	const double gConstant = 0.0;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;
	// # of cells
	const int nr = 64, nz = 128;
	// related to initialize level set
	const double baseR = 0.0, baseZ = 0.0, lenR = 0.02, lenZ = 0.04, cfl = 0.5;
	
	const int maxiter = 2000, niterskip = 50, num_bc_grid = 3;
	const double maxtime = 1.0;
	const bool writeVTK = false;
	// length of each cell
	const double dr = lenR / nr, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2DAxisym_StationaryBubbleVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2DAxisym_StationaryBubbleDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nr, nz, baseR, baseZ, lenR, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("axisym", "dirichlet", "neumann", "neumann");
	MSolver->SetBC_V_2D("axisym", "dirichlet", "neumann", "neumann");
	MSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);

	MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	// MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nr, nz, num_bc_grid, baseR, baseZ, dr, dz);
	// \phi^n
	std::vector<double> lsB((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0),
		HSmooth((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(nz, i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("axisym", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	double radius = 0.01, r = 0.0, z = 0.0, d = 0.0;
	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		// positive : inside & gas, negative : outside & liquid 
		r = baseR + (i + 0.5 - num_bc_grid) * dr;
		z = baseZ + (j + 0.5 - num_bc_grid) * dz;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(r, 2.0) + std::pow(z - lenZ * 0.5, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_2D(nz, i, j)] = d;
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);

	LSolver->ApplyBC_P_2D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dr, dz) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		vhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> rhsU((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		rhsV((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> div((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->ApplyBC_P_2D(ls);
		LSolver->Reinit_MinRK2_2D(ls);

		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		rhsU = MSolver->UpdateFU(ls, MSolver->m_u, MSolver->m_v, HSmooth);
		uhat = MSolver->GetUHat(ls, MSolver->m_u, rhsU, HSmooth);

		rhsV = MSolver->UpdateFV(ls, MSolver->m_u, uhat, MSolver->m_v, HSmooth);
		vhat = MSolver->GetVHat(ls, MSolver->m_v, rhsV, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence4Poisson(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		LSolver->ApplyBC_P_2D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

	std::cout << "(Final) Stationary Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC2DAxisymTest_SmallAirBubbleRising() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	const double gConstant = 9.81;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;
	// # of cells
	const int nr = 64, nz = 128;
	// related to initialize level set
	const double baseR = 0.0, baseZ = -0.01, lenR = 0.01, lenZ = 0.03, cfl = 0.5;
	
	const double maxtime = 2.0;
	const int maxiter = 1000, niterskip = 50, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dr = lenR / nr, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2DAxisym_SmallBubbleRisingVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2DAxisym_SmallBubbleRisingDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nr, nz, baseR, baseZ, lenR, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("axisym", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_V_2D("axisym", "dirichlet", "dirichlet", "dirichlet");
	MSolver->SetBC_P_2D("axisym", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BOTH);

	MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	// MSolver->SetPoissonSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::GS);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nr, nz, num_bc_grid, baseR, baseZ, dr, dz);
	// \phi^n
	std::vector<double> lsB((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative
	std::vector<double> H((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0), HSmooth((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(nz, i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("axisym", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	
	double radius = 1.0 / 300.0, r = 0.0, z = 0.0, d = 0.0;
	for (int j = 0; j < nz + 2 * num_bc_grid; j++)
	for (int i = 0; i < nr + 2 * num_bc_grid; i++) {
		// positive : inside, negative : outside
		r = baseR + (i + 0.5 - num_bc_grid) * dr;
		z = baseZ + (j + 0.5 - num_bc_grid) * dz;

		d = std::sqrt(r * r + z * z) - radius;

		ls[idx3_2D(nz, i, j)] = d;
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);

	LSolver->ApplyBC_P_2D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dr, dz) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		vhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> rhsU((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		rhsV((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> div((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part

		// Update F and apply time integration
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		rhsU = MSolver->UpdateFU(ls, MSolver->m_u, MSolver->m_v, HSmooth);
		uhat = MSolver->GetUHat(ls, MSolver->m_u, rhsU, HSmooth);

		rhsV = MSolver->UpdateFV(ls, MSolver->m_u, uhat, MSolver->m_v, HSmooth);
		vhat = MSolver->GetVHat(ls, MSolver->m_v, rhsV, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence4Poisson(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

int MAC2DAxisymTest_WaterDropletCollison1() {
	// Set initial level set
	const double L = 300 * 10e-6, U = 1.4;
	const double rhoL = 1.226, muL = 1.78e-5;
	const double rhoH = 1000.0, muH = 1.137e-3, sigma = 0.0728;
	const double gConstant = 0.0;
	const double Re = rhoH / muH * U * L, We = rhoH / sigma * L * U * U, Fr = 0.0;
	GAXISENUM2DAXISYM GAxis = GAXISENUM2DAXISYM::Z;

	// # of cells
	const int nr = 120, nz = 320;
	// related to initialize level set
	const double lenR = L * 3.0, lenZ = L * 8.0, cfl = 0.1, baseR = 0.0, baseZ = -0.5 * lenZ;
	const double densityRatio = rhoL / rhoH, viscosityRatio = muL / muH;

	const double maxtime = 2.0;
	const int maxiter = 10, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dr = lenR / nr, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC2D_WaterDropletCollision01Vel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC2D_WaterDropletCollision01Div_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nr))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver2DAxisym> MSolver;
	MSolver = std::make_unique<MACSolver2DAxisym>(Re, We, Fr, GAxis,
		L, U, sigma, densityRatio, viscosityRatio, rhoH, muH,
		nr, nz, baseR, baseZ, lenR, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_2D("axisym", "dirichlet", "neumann", "neumann");
	MSolver->SetBC_V_2D("axisym", "dirichlet", "neumann", "neumann");
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

	MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	const int poissonMaxIter = 20000;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nr, nz, num_bc_grid, baseR, baseZ, dr, dz);
	// \phi^n
	std::vector<double> lsB((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// \phi^{n + 1}
	std::vector<double> ls((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);
	// inside value must be positive levelset, otherwise, negative
	std::vector<double> H((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0),
		HSmooth((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid), 0.0);

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nr + 2 * num_bc_grid; i++)
	for (int j = 0; j < nz + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_v[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_p[idx3_2D(nz, i, j)] = 0.0;
		MSolver->m_ps[idx3_2D(nz, i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	MSolver->ApplyBC_P_2D(MSolver->m_ps);
	MSolver->ApplyBC_P_2D(MSolver->m_p);

	LSolver->SetBC_P_2D("axisym", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	
	double radius = 150 * 10e-6, r = 0.0, z = 0.0, d1 = 0.0, d2 = 0.0;
	for (int j = 3; j < nz + num_bc_grid; j++)
	for (int i = 3; i < nr + num_bc_grid; i++) {
		// positive : inside, negative : outside
		r = baseR + (i + 0.5 - num_bc_grid) * dr;
		z = baseZ + (j + 0.5 - num_bc_grid) * dz;
		
		d1 = std::sqrt(r * r + (z - L * 1.5) * (z - L * 1.5)) - radius;
		d2 = std::sqrt(r * r + (z + L * 1.5) * (z + L * 1.5)) - radius;

		if (d1 < 0)
			MSolver->m_v[idx3_2D(nz, i, j)] = -U * 0.5;

		if (d2 < 0)
			MSolver->m_v[idx3_2D(nz, i, j)] = U * 0.5;

		ls[idx3_2D(nz, i, j)] = std::min(d1, d2);
	}
	LSolver->ApplyBC_P_2D(ls);
	LSolver->Reinit_MinRK2_2D(ls);
	LSolver->ApplyBC_P_2D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dr, dz) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> uhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		vhat((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> rhsU((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid)),
		rhsV((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	std::vector<double> div((nr + 2 * num_bc_grid) * (nz + 2 * num_bc_grid));
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_dt);
		LSolver->Reinit_MinRK2_2D(ls);
		LSolver->ApplyBC_P_2D(ls);
		// Solve Momentum Part

		// Update F and apply time integration
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		// Get intermediate velocity
		rhsU = MSolver->UpdateFU(ls, MSolver->m_u, MSolver->m_v, HSmooth);
		uhat = MSolver->GetUHat(ls, MSolver->m_u, rhsU, HSmooth);

		rhsV = MSolver->UpdateFV(ls, MSolver->m_u, uhat, MSolver->m_v, HSmooth);
		vhat = MSolver->GetVHat(ls, MSolver->m_v, rhsV, HSmooth);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence4Poisson(uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, H, poissonMaxIter);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps, ls, lsB, H);

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

			std::cout << "Droplet Collision : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

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

	std::cout << "(Final) Droplet Collision : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}
