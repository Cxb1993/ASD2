#ifndef __MAC_H_
#define __MAC_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <array>
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>

#include "common.h"
#include "data.h"
#include "bc2d.h"
#include "poisson2d.h"

class MACSolver2D {
public:
	const int kNx, kNy, kNz, kNumBCGrid;
	const double kLenX, kLenY, kLenZ, kDx, kDy, kDz;

	const double kRe, kWe, kFr;
	const double kLScale, kUScale, kSigma, kG, kMuScale, kRhoScale;
	// *I : inside fluid, *O : outside fluid
	const double kRhoI, kRhoO, kRhoRatio;
	const double kMuI, kMuO, kMuRatio;

	const int kMaxIter, kNIterSkip;
	const double kCFL, kMaxTime;
	const bool kWriteVTK;
	POISSONTYPE m_PoissonSolverType;
	PLTTYPE m_PLTType;

	const double kEps_div = 1.0e-6;

	std::shared_ptr<BoundaryCondition> m_BC;
	std::shared_ptr<PoissonSolver> m_Poisson;
	int m_iter;
	double m_dt, m_totTime;
	std::string m_outFname;

	// velocity
	std::vector<double> m_u, m_v, m_w;
	// density, curvature, viscosity
	std::vector<double> m_rho, m_kappa, m_mu;

	std::vector<double> m_ps, m_p;

	MACSolver2D();
	MACSolver2D::MACSolver2D(double Re, double We, double Fr,
		double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoI, double mu1,
		int nx, int ny, int nz, double lenX, double lenY, double lenZ, double cfl,
		int maxtime, int maxiter, int niterskip, int num_bc_grid, bool writeVTK);
	MACSolver2D::MACSolver2D(double rhoI, double rhoO, double muI, double muO, double gConstant,
		double L, double U, double sigma,
		int nx, int ny, int nz, double lenX, double lenY, double lenZ, double cfl,
		int maxtime, int maxiter, int niterskip, int num_bc_grid, bool writeVTK);
	~MACSolver2D();

	int AllocateVariables();

	// Related to Level Set Related
	int UpdateKappa(const std::vector<double> ls);
	
	// Convection Term
	std::vector<double> AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v);
	std::vector<double> AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v);
	int UnitHJWENO5(
		const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n);

	// Viscosity Term
	std::vector<double> AddViscosityFU(const std::vector<double>& u, const std::vector<double>& v, const std::vecotr<double>& ls);
	// Surface Term

	// Gravity Term

	// Poisson 
	int SetPoissonSolver(POISSONTYPE type);

	// update velocity using projection
	int UpdateVel(std::vector<double> u, std::vector<double> v,
		std::vector<double> us, std::vector<double> vs, std::vector<double> ps);

	// BC
	int SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int SetBC_W_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int ApplyBC_U_2D(std::vector<double> arr);
	int ApplyBC_V_2D(std::vector<double> arr);
	int ApplyBC_W_2D(std::vector<double> arr);
	int ApplyBC_P_2D(std::vector<double> arr);
	void SetBCConstantUW(double BC_ConstantW);
	void SetBCConstantUE(double BC_ConstantE);
	void SetBCConstantUS(double BC_ConstantS);
	void SetBCConstantUN(double BC_ConstantN);
	void SetBCConstantVW(double BC_ConstantW);
	void SetBCConstantVE(double BC_ConstantE);
	void SetBCConstantVS(double BC_ConstantS);
	void SetBCConstantVN(double BC_ConstantN);
	void SetBCConstantPW(double BC_ConstantW);
	void SetBCConstantPE(double BC_ConstantE);
	void SetBCConstantPS(double BC_ConstantS);
	void SetBCConstantPN(double BC_ConstantN);

	void TDMA(double* a, double* b, double* c, double* d, int n);
	int SetPLTType(PLTTYPE type);
	int OutRes(int iter, double curTime, const std::string fname_vel_base, const std::string fname_div_base,
		double ***u, double ***v, double ***w, double ***phi);
	int OutResClose();

	int idx(int i, int j);
};

#endif __MAC_H_