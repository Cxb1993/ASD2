#include "test_mac3d_large.h"

int MAC3DTest_NonSurfaceTension() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.226, muL = 1.78e-5;
	// const double rhoH = 1000, muH = 1.137e-3;
	// const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	const double rhoH = 1.226, muH = 1.78e-5, sigma = 0.0;
	const double gConstant = 0.0;
	GAXISENUM3D GAxis = GAXISENUM3D::Y;
	// # of cells
	const int nx = 128, ny = 128, nz = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, baseZ = 0.0, lenX = 1.0, lenY = 1.0, lenZ = 1.0, cfl = 0.5;
	double radius = 0.2, x = 0.0, y = 0.0, z = 0.0, d = 0.0;

	const int maxiter = 2, niterskip = 1, num_bc_grid = 3;
	const int64_t arrSize = (nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid) * (nz + 2 * num_bc_grid);
	const double maxtime = 0.06;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC3D_NonSurfaceTensionVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC3D_NonSurfaceTensionDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	double ambientPressure = 1.0;
	int stat = 0;

	std::unique_ptr<MACSolver3D> MSolver;
	MSolver = std::make_unique<MACSolver3D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, nz, baseX, baseY, baseY, lenX, lenY, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_3D("wall", "wall", "wall", "wall", "inlet", "outlet");
	MSolver->SetBC_V_3D("wall", "wall", "wall", "wall", "inlet", "outlet");
	MSolver->SetBC_W_3D("wall", "wall", "wall", "wall", "inlet", "outlet");
	MSolver->SetBC_P_3D("wall", "wall", "wall", "wall", "neumann", "pressure");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantUB(0.0);
	MSolver->SetBCConstantUT(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetBCConstantVB(0.0);
	MSolver->SetBCConstantVT(0.0);
	MSolver->SetBCConstantWW(0.0);
	MSolver->SetBCConstantWE(0.0);
	MSolver->SetBCConstantWS(0.0);
	MSolver->SetBCConstantWN(0.0);
	MSolver->SetBCConstantWB(1.0);
	MSolver->SetBCConstantWT(0.0);
	
	MSolver->SetPLTType(PLTTYPE::BINARY);

	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	MSolver->SetImplicitSolver(POISSONTYPE::CG);
	// MSolver->SetPoissonSolver(POISSONTYPE::BICGSTAB);
	const int poissonMaxIter = 2500;

	std::shared_ptr<LevelSetSolver3D> LSolver;
	LSolver = std::make_shared<LevelSetSolver3D>(nx, ny, nz, num_bc_grid, baseX, baseY, baseZ, dx, dy, dz);
	// \phi^n
	std::vector<double> lsB(arrSize, 0.0);
	// \phi^{n + 1}
	std::vector<double> ls(arrSize, 0.0);
	// inside value must be positive levelset, otherwise, negative

	std::vector<double> H(arrSize, 0.0),
		HSmooth(arrSize, 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	MSolver->ApplyBC_U_3D(MSolver->m_u);
	MSolver->ApplyBC_V_3D(MSolver->m_v);
	MSolver->ApplyBC_W_3D(MSolver->m_w);
	MSolver->ApplyBC_P_3D(MSolver->m_ps);
	MSolver->ApplyBC_P_3D(MSolver->m_p);
	LSolver->SetBC_P_3D("neumann", "neumann", "neumann", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	LSolver->SetBCConstantPB(0.0);
	LSolver->SetBCConstantPT(0.0);

	for (int k = 0; k < nz + 2 * num_bc_grid; k++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside & gas, negative : outside & liquid 
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;
		z = baseZ + (k + 0.5 - num_bc_grid) * dz;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(x - lenX * 0.5, 2.0) + std::pow(y - lenY * 0.5, 2.0) + std::pow(z - lenZ * 0.5, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_3D(ny, nz, i, j, k)] = d;

		MSolver->m_u[idx3_3D(ny, nz, i, j, k)] = 0.0;
		MSolver->m_v[idx3_3D(ny, nz, i, j, k)] = 0.0;
		MSolver->m_w[idx3_3D(ny, nz, i, j, k)] = 0.0;
		MSolver->m_p[idx3_3D(ny, nz, i, j, k)] = ambientPressure;
		MSolver->m_ps[idx3_3D(ny, nz, i, j, k)] = ambientPressure;
	}

	LSolver->ApplyBC_P_3D(ls);
	LSolver->Reinit_MinRK2_3D(ls);

	LSolver->ApplyBC_P_3D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(std::min(dx, dy), dz) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;
	MSolver->UpdateAmbientPressure(MSolver->m_dt * ambientPressure);

	std::vector<double> rhsU(arrSize), rhsV(arrSize), rhsW(arrSize);
	std::vector<double> uhat(arrSize), vhat(arrSize), what(arrSize);
	std::vector<double> div(arrSize);
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_3D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_dt);
		LSolver->ApplyBC_P_3D(ls);
		LSolver->Reinit_MinRK2_3D(ls);

		LSolver->ApplyBC_P_3D(ls);

		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		rhsU = MSolver->GetRHSU(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);
		rhsV = MSolver->GetRHSV(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);
		rhsW = MSolver->GetRHSW(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);

		uhat = MSolver->GetUHat(ls, rhsU, H, poissonMaxIter);
		vhat = MSolver->GetVHat(ls, rhsV, H, poissonMaxIter);
		what = MSolver->GetWHat(ls, rhsW, H, poissonMaxIter);

		MSolver->ApplyBC_U_3D(uhat);
		MSolver->ApplyBC_V_3D(vhat);
		MSolver->ApplyBC_W_3D(what);

		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat, what);
		MSolver->ApplyBC_P_3D(MSolver->m_ps);
		LSolver->ApplyBC_P_3D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, MSolver->m_w, H, poissonMaxIter);
		MSolver->ApplyBC_P_3D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v, MSolver->m_w,
			uhat, vhat, what, MSolver->m_ps, ls, lsB, H);

		MSolver->ApplyBC_U_3D(MSolver->m_u);
		MSolver->ApplyBC_V_3D(MSolver->m_v);
		MSolver->ApplyBC_W_3D(MSolver->m_w);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v, MSolver->m_w);
		MSolver->m_curTime += MSolver->m_dt;
		MSolver->UpdateAmbientPressure(MSolver->m_dt * ambientPressure);

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Non SurfaceTension : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

			// MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
			// 	MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);
			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				uhat, vhat, what, MSolver->m_ps, ls);
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
		MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}

int MAC3DTest_StationaryBubble() {
	// Set initial level set
	const double Re = 0.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoL = 1.226, muL = 1.78e-5;
	// const double rhoH = 1000, muH = 1.137e-3;
	const double rhoH = 1000, muH = 1.137e-3, sigma = 0.0728;
	// const double rhoL = 1000, muL = 1.137e-3, sigma = 0.0;
	const double ambientPressure = 1.0;
	const double gConstant = 0.0; 
	GAXISENUM3D GAxis = GAXISENUM3D::Y;
	// # of cells
	const int nx = 128, ny = 128, nz = 128;
	// related to initialize level set
	const double baseX = 0.0, baseY = 0.0, baseZ = 0.0, lenX = 0.04, lenY = 0.04, lenZ = 0.04, cfl = 0.5;
	double radius = 0.01, x = 0.0, y = 0.0, z = 0.0, d = 0.0;

	const int maxiter = 5, niterskip = 1, num_bc_grid = 3;
	const int64_t arrSize = (nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid) * (nz + 2 * num_bc_grid);
	const double maxtime = 0.06;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny, dz = lenZ / nz;
	std::ostringstream outfname_stream1;
	outfname_stream1 << "testMAC3D_StationaryBubbleBubbleRisingVel_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	const std::string fname_vel = outfname_stream1.str();

	std::ostringstream outfname_stream2;
	outfname_stream2 << "testMAC3D_StationaryBubbleBubbleRisingDiv_Re_"
		<< std::to_string(Re) << "_"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nx))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << ny))->str()
		<< "x"
		<< static_cast<std::ostringstream*>(&(std::ostringstream() << nz))->str();
	const std::string fname_div = outfname_stream2.str();
	const int iterskip = 1;
	int stat = 0;

	std::unique_ptr<MACSolver3D> MSolver;
	MSolver = std::make_unique<MACSolver3D>(rhoH, rhoL, muH, muL, gConstant, GAxis,
		L, U, sigma, nx, ny, nz, baseX, baseY, baseY, lenX, lenY, lenZ,
		cfl, maxtime, maxiter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_3D("wall", "wall", "wall", "wall", "wall", "wall");
	MSolver->SetBC_V_3D("wall", "wall", "wall", "wall", "wall", "wall");
	MSolver->SetBC_W_3D("wall", "wall", "wall", "wall", "wall", "wall");
	MSolver->SetBC_P_3D("wall", "wall", "wall", "wall", "wall", "wall");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantUB(0.0);
	MSolver->SetBCConstantUT(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetBCConstantVB(0.0);
	MSolver->SetBCConstantVT(0.0);
	MSolver->SetBCConstantWW(0.0);
	MSolver->SetBCConstantWE(0.0);
	MSolver->SetBCConstantWS(0.0);
	MSolver->SetBCConstantWN(0.0);
	MSolver->SetBCConstantWB(0.0);
	MSolver->SetBCConstantWT(0.0);
	MSolver->UpdateAmbientPressure(MSolver->m_dt * ambientPressure);
	MSolver->SetPLTType(PLTTYPE::BINARY);

	MSolver->SetPoissonSolver(POISSONTYPE::CG);
	MSolver->SetImplicitSolver(POISSONTYPE::CG);
	const int poissonMaxIter = 2500;

	std::shared_ptr<LevelSetSolver3D> LSolver;
	LSolver = std::make_shared<LevelSetSolver3D>(nx, ny, nz, num_bc_grid, baseX, baseY, baseZ, dx, dy, dz);
	// \phi^n
	std::vector<double> lsB(arrSize, 0.0);
	// \phi^{n + 1}
	std::vector<double> ls(arrSize, 0.0);
	// inside value must be positive levelset, otherwise, negative
	
	std::vector<double> H(arrSize, 0.0), 
		HSmooth(arrSize, 0.0);
	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	MSolver->ApplyBC_U_3D(MSolver->m_u);
	MSolver->ApplyBC_V_3D(MSolver->m_v);
	MSolver->ApplyBC_W_3D(MSolver->m_w);
	MSolver->ApplyBC_P_3D(MSolver->m_ps);
	MSolver->ApplyBC_P_3D(MSolver->m_p);
	LSolver->SetBC_P_3D("neumann", "neumann", "neumann", "neumann", "neumann", "neumann");
	LSolver->SetBCConstantPW(0.0);
	LSolver->SetBCConstantPE(0.0);
	LSolver->SetBCConstantPS(0.0);
	LSolver->SetBCConstantPN(0.0);
	LSolver->SetBCConstantPB(0.0);
	LSolver->SetBCConstantPT(0.0);
	
	for (int k = 0; k < nz + 2 * num_bc_grid; k++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++)
	for (int i = 0; i < nx + 2 * num_bc_grid; i++) {
		// positive : inside & gas, negative : outside & liquid 
		x = baseX + (i + 0.5 - num_bc_grid) * dx;
		y = baseY + (j + 0.5 - num_bc_grid) * dy;
		z = baseZ + (k + 0.5 - num_bc_grid) * dz;

		// d - inside : -, outside : +
		d = std::sqrt(std::pow(x - 0.02, 2.0) + std::pow(y - 0.02, 2.0) + std::pow(z - 0.02, 2.0)) - radius;

		// ls - inside : -(gas), outside : +(liquid)
		ls[idx3_3D(ny, nz, i, j, k)] = d;
	}
	LSolver->ApplyBC_P_3D(ls);
	LSolver->Reinit_MinRK2_3D(ls);
	
	LSolver->ApplyBC_P_3D(ls);

	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dx, dy) / U;
	std::cout << " dt : " << MSolver->m_dt << std::endl;

	std::vector<double> rhsU(arrSize), rhsV(arrSize), rhsW(arrSize);
	std::vector<double> uhat(arrSize), vhat(arrSize), what(arrSize);
	std::vector<double> div(arrSize);
	MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
		MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);

	while (MSolver->m_curTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		MSolver->m_iter++;

		lsB = ls;
		LSolver->Solve_LevelSet_3D(ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_dt);
		LSolver->ApplyBC_P_3D(ls);
		LSolver->Reinit_MinRK2_3D(ls);
		
		LSolver->ApplyBC_P_3D(ls);
		
		// Solve Momentum Part
		MSolver->UpdateKappa(ls);
		H = MSolver->UpdateHeavisideFunc(ls);
		HSmooth = MSolver->UpdateSmoothHeavisideFunc(ls);

		rhsU = MSolver->GetRHSU(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);
		rhsV = MSolver->GetRHSV(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);
		rhsW = MSolver->GetRHSW(LSolver, ls, MSolver->m_u, MSolver->m_v, MSolver->m_w, H);

		uhat = MSolver->GetUHat(ls, rhsU, H, poissonMaxIter);
		vhat = MSolver->GetVHat(ls, rhsV, H, poissonMaxIter);
		what = MSolver->GetWHat(ls, rhsW, H, poissonMaxIter);

		MSolver->ApplyBC_U_3D(uhat);
		MSolver->ApplyBC_V_3D(vhat);
		MSolver->ApplyBC_W_3D(what);
		
		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat, what);
		MSolver->ApplyBC_P_3D(MSolver->m_ps);
		LSolver->ApplyBC_P_3D(ls);

		// Solve Poisson equation
		// m_phi = pressure * dt
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, lsB, MSolver->m_u, MSolver->m_v, MSolver->m_w, H, poissonMaxIter);
		MSolver->ApplyBC_P_3D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v, MSolver->m_w,
			uhat, vhat, what, MSolver->m_ps, ls, lsB, H);

		MSolver->ApplyBC_U_3D(MSolver->m_u);
		MSolver->ApplyBC_V_3D(MSolver->m_v);
		MSolver->ApplyBC_W_3D(MSolver->m_w);

		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v, MSolver->m_w);
		MSolver->m_curTime += MSolver->m_dt;

		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			std::chrono::high_resolution_clock::time_point p = std::chrono::high_resolution_clock::now();
			std::chrono::milliseconds ms = std::chrono::duration_cast<std::chrono::milliseconds>(p.time_since_epoch());

			std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(ms);
			std::time_t t = s.count();
			std::size_t fractional_seconds = ms.count() % 1000;

			std::cout << "Stationary Bubble : " << std::ctime(&t) << " " << MSolver->m_iter << " " << MSolver->m_curTime << " " << MSolver->m_dt << " " << std::endl;

			MSolver->OutRes(MSolver->m_iter, MSolver->m_curTime, fname_vel, fname_div,
				MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);
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
		MSolver->m_u, MSolver->m_v, MSolver->m_w, MSolver->m_ps, ls);
	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}
