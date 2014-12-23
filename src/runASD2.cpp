#include "../include/mac2d.h"

int main(int argc, char *argv[]) {
	// Set initial level set
	const double Re = 100.0, We = 0.0, Fr = 0.0;
	const double L = 1.0, U = 1.0;
	const double rho1 = 1.0, mu1 = rho1 * L * U / Re, sigma = 0.0;
	// for 1-phase flow, set ratio to 1.
	const double densityRatio = 1, viscosityRatio = 1;
	// # of cells
	// const int nx = 128, ny = 128, nz = 8;
	const int nx = 128, ny = 128, nz = 128;
	const double lenX = L, lenY = L, cfl = 0.1;
	const int maxtime = 2.0, niter = 30, niterskip = 1, num_bc_grid = 3;
	const bool writeVTK = false;
	// length of each cell
	const double dx = lenX / nx, dy = lenY / ny;
	const double dz = dx;
	const double lenZ = dz * nz;
	const std::string fname_vel("cavityVel_Re_" + std::to_string(Re));
	const std::string fname_div("cavityDiv_Re_" + std::to_string(Re));
	int iterskip = 1;
	std::unique_ptr<MACSolver2D> MSolver;
	MSolver = std::make_unique<MACSolver2D>(Re, We, Fr,
		L, U, sigma, densityRatio, viscosityRatio, rho1, mu1,
		nx, ny, nz, lenX, lenY, lenZ, cfl,
		maxtime, niter, niterskip, num_bc_grid, writeVTK);
	MSolver->SetBC_U_3D("dirichlet", "dirichlet", "dirichlet", "dirichlet", "periodic", "periodic");
	MSolver->SetBC_V_3D("dirichlet", "dirichlet", "dirichlet", "dirichlet", "periodic", "periodic");
	MSolver->SetBC_W_3D("dirichlet", "dirichlet", "dirichlet", "dirichlet", "periodic", "periodic");
	MSolver->SetBC_P_3D("neumann", "neumann", "neumann", "neumann", "neumann", "neumann");
	MSolver->SetBCConstantUW(0.0);
	MSolver->SetBCConstantUE(0.0);
	MSolver->SetBCConstantUS(0.0);
	MSolver->SetBCConstantUN(0.0);
	MSolver->SetBCConstantVW(0.0);
	MSolver->SetBCConstantVE(0.0);
	MSolver->SetBCConstantVS(0.0);
	MSolver->SetBCConstantVN(0.0);
	MSolver->SetPLTType(PLTTYPE::BINARY);
	MSolver->SetPoissonSolver(POISSONTYPE::CG);
}