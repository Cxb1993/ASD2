#include "mac2d.h"
#include "levelset2d.h"

int GNx = 0, GNy = 0, GNumBCGrid = 0;
int idx(int i, int j);

int main(int argc, char *argv[]) {
	// Set initial level set
	const double Re = 100.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rhoI = 1.0, muI = rhoI * L * U / Re, sigma = 0.0;
	// for 1-phase flow, set ratio to 1.
	const double densityRatio = 1, viscosityRatio = 1;
	// # of cells
	// const int nx = 128, ny = 128, nz = 8;
	const int nx = 128, ny = 128;
	const double lenX = L, lenY = L, cfl = 0.1;
	const int maxtime = 2.0, niter = 30, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const std::string fname_vel("LSVel_Re_" + std::to_string(Re));
	const std::string fname_div("LSDiv_Re_" + std::to_string(Re));
	int iterskip = 1;
	int stat = 0;

	GNx = nx;
	GNy = ny;
	GNumBCGrid = num_bc_grid;

	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr,
		L, U, sigma, densityRatio, viscosityRatio, rhoI, muI,
		nx, ny, lenX, lenY, cfl,
		maxtime, niter, niterskip, num_bc_grid, writeVTK);
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
	MSolver->SetPLTType(PLTTYPE::BINARY);
	MSolver->SetPoissonSolver(POISSONTYPE::MKL);

	std::shared_ptr<LevelSetSolver2D> LSolver;
	LSolver = std::make_shared<LevelSetSolver2D>(nx, ny, num_bc_grid, dx, dy);
	std::vector<double> ls((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	// inside value must be positive levelset, otherwise, negative

	// init level set
	for (int i = num_bc_grid; i < nx + num_bc_grid; i++)
	for (int j = num_bc_grid; j < ny + num_bc_grid; j++) {
		ls[idx(i, j)] = 1.0;
	}
	LSolver->Reinit_Sussman_2D(ls);

	LSolver->m_signedInitLS = LSolver->GetSignedLSNormalized(ls);
	LSolver->SetBC_P_2D("neumann", "neumann", "neumann", "neumann");

	// init velocity and pseudo-pressure
	MSolver->AllocateVariables();

	for (int i = 0; i < nx + 2 * num_bc_grid; i++)
	for (int j = 0; j < ny + 2 * num_bc_grid; j++) {
		MSolver->m_u[idx(i, j)] = 1.0;
		MSolver->m_v[idx(i, j)] = 0.0;
		MSolver->m_w[idx(i, j)] = 0.0;
		MSolver->m_ps[idx(i, j)] = 0.0;
	}

	MSolver->ApplyBC_U_2D(MSolver->m_u);
	MSolver->ApplyBC_V_2D(MSolver->m_v);
	
	// prevent dt == 0.0
	MSolver->m_dt = cfl * std::min(dx, dy) / U;
	MSolver->OutRes(0, 0.0, fname_vel, fname_div, MSolver->m_u, MSolver->m_v, MSolver->m_ps);

	std::vector<double> FU((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid)), FV((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	std::vector<double> uhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid)), vhat((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	std::vector<double> div((nx + 2 * num_bc_grid) * (ny + 2 * num_bc_grid));
	
	while (MSolver->m_totTime < MSolver->kMaxTime && MSolver->m_iter < MSolver->kMaxIter) {
		// Solver Level set part first
		// Have to use \phi^{n+1} for rho, mu, kappa
		LSolver->Solve_LevelSet_2D(ls, MSolver->m_u, MSolver->m_v);
		LSolver->Reinit_Sussman_2D(ls);

		// Solve Momentum Part
		// Update Rho, Mu, Kappa
		if (We != 0.0)
			MSolver->UpdateKappa(ls);
		MSolver->ApplyBC_P_2D(MSolver->m_kappa);

		// Update F and apply time discretization
		FU = MSolver->UpdateFU(LSolver, ls,
			MSolver->m_u, MSolver->m_v);
		FV = MSolver->UpdateFV(LSolver, ls,
			MSolver->m_u, MSolver->m_v);
		
		// Get intermediate velocity
		uhat = MSolver->GetUhat(MSolver->m_u, FU);
		vhat = MSolver->GetVhat(MSolver->m_v, FV);

		MSolver->ApplyBC_U_2D(uhat);
		MSolver->ApplyBC_V_2D(vhat);
		// MSolver->OutRes(MSolver->m_iter, MSolver->m_totTime, fname_vel, fname_div, uhat, vhat, what, MSolver->m_phi);
		// From intermediate velocity, get divergence
		div = MSolver->GetDivergence(uhat, vhat);

		MSolver->ApplyBC_P_2D(MSolver->m_ps);
		MSolver->ApplyBC_P_2D(div);
		// Solve Poisson equation
		// m_phi = pressure * dt / rho
		stat = MSolver->SolvePoisson(MSolver->m_ps, div, ls, uhat, vhat);
		MSolver->ApplyBC_P_2D(MSolver->m_ps);

		stat = MSolver->UpdateVel(MSolver->m_u, MSolver->m_v,
			uhat, vhat, MSolver->m_ps);

		MSolver->ApplyBC_U_2D(MSolver->m_u);
		MSolver->ApplyBC_V_2D(MSolver->m_v);
		
		MSolver->m_dt = MSolver->UpdateDt(MSolver->m_u, MSolver->m_v);
		MSolver->m_totTime += MSolver->m_dt;
		MSolver->m_iter++;
		if ((MSolver->m_iter % MSolver->kNIterSkip) == 0) {
			MSolver->OutRes(MSolver->m_iter, MSolver->m_totTime, fname_vel, fname_div, MSolver->m_u, MSolver->m_v, MSolver->m_ps);
		}
	}

	MSolver->OutResClose();

	MSolver.reset();
	LSolver.reset();

	return 0;
}


inline int idx(int i, int j) {
	return i + (GNx + 2 * GNumBCGrid) * j;
}