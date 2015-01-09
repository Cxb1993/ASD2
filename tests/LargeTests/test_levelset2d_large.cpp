#include "test_levelset2d_large.h"

int LevelSetTest_2D_Sussman621_ReinitSussman() {
	// Zalesak Problem
	// int nx = 200, ny = 200;
	// double dx = 0.5, dy = 0.5, dt = 0.5;
	int nx = 100, ny = 100;
	double dx = 1.0, dy = 1.0, dt = 1.0;
	int iter = 0;
	double curTime = 0.0;
	double x = 0.0, y = 0.0;

	int orgMass = 0, LSMass = 0;

	std::vector<double> lsOrg((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> ls((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> err((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> u((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> v((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> div((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> phi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);
	std::vector<double> massVec;

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, NBC3, dx, dy);
	
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
	LSolver->ApplyBC_P_2D(ls);

	double radius = 15.0, slotW = 15.0, d = 0.0;
	double baseX = 0.0, baseY = 0.0;
	double centerX = 50.0, centerY = 75.0;
	for (int i = NBC3; i < nx + NBC3; i++)
	for (int j = NBC3; j < ny + NBC3; j++) {
		// positive : inside, negative : outside
		x = baseX + (i - NBC3) * dx;
		y = baseY + (j - NBC3) * dy;
		d = std::sqrt(std::pow(x - centerX, 2.0) + std::pow(y - centerY, 2.0)) - radius;

		ls[idx3(nx, i, j)] = -2.0;

		if (x >= (50 - 0.25 * slotW) && x <= (50 + 0.25 * slotW)) {
			if (y <= 75 + 7.5)
				ls[idx3(nx, i, j)] = -2.0;
		}

		if (ls[idx3(nx, i, j)] > 0)
			orgMass++;

		u[idx3(nx, i, j)] = (M_PI / 314) * (50.0 - y);
		v[idx3(nx, i, j)] = (M_PI / 314) * (x - 50.0);

		lsOrg[idx3(nx, i, j)] = ls[idx3(nx, i, j)];

		err[idx3(nx, i, j)] = 0.0;
	}

	LSolver->m_signedInitLS = LSolver->GetSignedLSNormalized(ls);
	LSolver->Reinit_Sussman_2D(ls);

	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");
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

	OutRes(iter, curTime, fname_vel_base, fname_div_base,
		u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);

	while (curTime < 628.0) {
		iter++;
		curTime += dt;

		LSolver->Solve_LevelSet_2D(ls, u, v);
		LSolver->Reinit_Sussman_2D(ls);
		LSolver->ApplyBC_P_2D(ls);

		OutRes(iter, curTime, fname_vel_base, fname_div_base,
			u, v, phi, div, ls, nx, ny, dx, dy, baseX, baseY, m_PLTType);
		
		LSMass = 0;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++) {
			if (phi[idx3(nx, i, j)] > 0.0)
				LSMass++;
		}
			
		massVec.push_back(LSMass * dx * dy);
	}

	double errNormInfty = -std::numeric_limits<double>::max();
	double errNorm1 = 0.0, errNorm2 = 0.0, errMass = 0.0;

	LSMass = 0;
	orgMass = 0;

	for (int i = NBC3; i < nx + NBC3; i++)
	for (int j = NBC3; j < ny + NBC3; j++) {
		err[idx3(nx, i, j)] = std::fabs(ls[idx3(nx, i, j)] - lsOrg[idx3(nx, i, j)]);
		errNormInfty = std::max(std::fabs(err[idx3(nx, i, j)]), errNormInfty);
		errNorm1 += std::fabs(ls[idx3(nx, i, j)] - lsOrg[idx3(nx, i, j)]);
		errNorm2 += std::pow(std::fabs(ls[idx3(nx, i, j)] - lsOrg[idx3(nx, i, j)]), 2.0);

		if (lsOrg[idx3(nx, i, j)] < 0.0)
			orgMass++;

		if (ls[idx3(nx, i, j)] < 0.0)
			LSMass++;
	}

	errMass = std::fabs(static_cast<double>(LSMass - orgMass)) / orgMass * 100;
	errNorm2 = std::sqrt(errNorm2) / static_cast<double>(nx * ny);

	std::cout << "Sussman 6-2-1 - err(L-Infty):" << errNormInfty << " err(L1):" << errNorm1 << " err(L2):" << errNorm2 << std::endl;
	std::cout << "Sussman 6-2-1 - err(Mass loss) : " << errMass << "%" << std::endl;

	std::cout << "LevelSet Zalesak Disc : Sussman 6-2-1 (2D) " << std::endl;

	return 0;
}
