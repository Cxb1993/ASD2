#ifndef __MAC3D_H_
#define __MAC3D_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <cfloat>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <istream>
#include <fstream>
#include <string>

#include <algorithm>
#include <chrono>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "common.h"
#include "bc3d.h"
#include "levelset3d.h"
#include "poisson3d.h"

// Tecplot  IO
#include "TECIO.h"
#include "TECXXX.h"

class MACSolver3D {
public:
	const int kNx, kNy, kNz, kNumBCGrid;
	const int64_t kArrSize;
	const double kBaseX, kBaseY, kBaseZ, kLenX, kLenY, kLenZ, kDx, kDy, kDz;

	const double kRe, kWe, kFr;
	const double kLScale, kUScale, kSigma, kG, kMuScale, kRhoScale;
	// *H : Heavy fluid(such as liquid), *L : light fluid (such as gas)
	const double kRhoH, kRhoL, kRhoRatio;
	const double kMuH, kMuL, kMuRatio;
	
	const int kMaxIter, kNIterSkip;
	const TIMEORDERENUM kTimeOrder;
	const double kCFL, kMaxTime;
	const bool kWriteVTK;
	POISSONTYPE m_PoissonSolverType;
	PLTTYPE m_PLTType;
	GAXISENUM3D kGAxis;

	const double kEps_div = 1.0e-6;

	std::shared_ptr<BoundaryCondition3D> m_BC;
	std::shared_ptr<PoissonSolver3D> m_Poisson;
	int m_iter;
	double m_dt, m_curTime;
	std::string m_outFname;

	// velocity
	std::vector<double> m_u, m_v, m_w;
	// density, curvature, viscosity
	std::vector<double> m_rho, m_kappa, m_mu;

	std::vector<double> m_ps, m_p;

	// normal vector & tangent vector
	std::vector<double> m_nx, m_ny, m_nz;
	std::vector<double> m_t1x, m_t1y, m_t1z;
	std::vector<double> m_t2x, m_t2y, m_t2z;

	MACSolver3D();
	MACSolver3D(double Re, double We, double Fr, GAXISENUM3D GAxis,
		double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoI, double mu1,
		int nx, int ny, int nz, double baseX, double baseY, double baseZ, double lenX, double lenY, double lenZ,
		TIMEORDERENUM RKOrder, double cfl, double maxtime, int maxiter, int niterskip,
		int num_bc_grid, bool writeVTK);
	MACSolver3D(double rhoI, double rhoO, double muI, double muO,
		double gConstant, GAXISENUM3D GAxis, 
		double L, double U, double sigma,
		int nx, int ny, int nz, double baseX, double baseY, double baseZ, double lenX, double lenY, double lenZ,
		TIMEORDERENUM RKOrder, double cfl, double maxtime, int maxiter, int niterskip,
		int num_bc_grid, bool writeVTK);
	~MACSolver3D();

	int AllocateVariables();

	// Related to Level Set Related
	std::vector<double> UpdateSmoothHeavisideFunc(const std::vector<double>& ls);
	std::vector<double> UpdateHeavisideFunc(const std::vector<double>& ls);
	int UpdateNTK(const std::shared_ptr<LevelSetSolver3D>& LSolver, const std::vector<double>& ls,
		std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> normalVec,
		std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t1Vec,
		std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t2Vec,
		std::vector<double>& kappa);
	int UpdateKappa(const std::vector<double>& ls);

	std::vector<double> UpdateFU(const std::shared_ptr<LevelSetSolver3D>& LSolver,
		const std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		const std::vector<double>& H);
	std::vector<double> UpdateFV(const std::shared_ptr<LevelSetSolver3D>& LSolver,
		const std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		const std::vector<double>& H);
	std::vector<double> UpdateFW(const std::shared_ptr<LevelSetSolver3D>& LSolver,
		const std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		const std::vector<double>& H);

	// Convection Term
	std::vector<double> AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);
	std::vector<double> AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);
	std::vector<double> AddConvectionFW(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);
	int UnitHJWENO5(
		const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n);

	// Viscosity Term
	std::vector<double> AddViscosityFU(
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
			const std::vector<double>& ls, const std::vector<double>& H);
	std::vector<double> AddViscosityFV(
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
			const std::vector<double>& ls, const std::vector<double>& H);
	std::vector<double> AddViscosityFW(
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
			const std::vector<double>& ls, const std::vector<double>& H);
	
	// Gravity Term
	std::vector<double> AddGravityFU();
	std::vector<double> AddGravityFV();
	std::vector<double> AddGravityFW();
	
	// Intermediate Velocity
	int GetIntermediateVel(const std::shared_ptr<LevelSetSolver3D>& LSolver, const std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		std::vector<double>& uhat, std::vector<double>& vhat, std::vector<double>& what,
		const std::vector<double>& H);
	
	// Poisson 
	int SetPoissonSolver(POISSONTYPE type);
	int SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
		const std::vector<double>& ls, const std::vector<double>& lsB,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, const std::vector<double>& H, const int maxiter);

	// update velocity using projection
	std::vector<double> GetDivergence(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);
	int UpdateVel(std::vector<double>& u, std::vector<double>& v, std::vector<double>& w,
		const std::vector<double>& us, const std::vector<double>& vs, const std::vector<double>& ws,
		const std::vector<double>& ps, const std::vector<double>& ls, const std::vector<double>& lsB,
		const std::vector<double>& H);
	double UpdateDt(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);

	// BC
	int SetBC_U_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_V_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_W_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_P_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int ApplyBC_U_3D(std::vector<double>& arr);
	int ApplyBC_V_3D(std::vector<double>& arr);
	int ApplyBC_W_3D(std::vector<double>& arr);
	int ApplyBC_P_3D(std::vector<double>& arr);
	void SetBCConstantUW(double BC_ConstantW);
	void SetBCConstantUE(double BC_ConstantE);
	void SetBCConstantUS(double BC_ConstantS);
	void SetBCConstantUN(double BC_ConstantN);
	void SetBCConstantUB(double BC_ConstantS);
	void SetBCConstantUT(double BC_ConstantN);
	void SetBCConstantVW(double BC_ConstantW);
	void SetBCConstantVE(double BC_ConstantE);
	void SetBCConstantVS(double BC_ConstantS);
	void SetBCConstantVN(double BC_ConstantN);
	void SetBCConstantVB(double BC_ConstantS);
	void SetBCConstantVT(double BC_ConstantN);
	void SetBCConstantWW(double BC_ConstantW);
	void SetBCConstantWE(double BC_ConstantE);
	void SetBCConstantWS(double BC_ConstantS);
	void SetBCConstantWN(double BC_ConstantN);
	void SetBCConstantWB(double BC_ConstantS);
	void SetBCConstantWT(double BC_ConstantN);
	void SetBCConstantPW(double BC_ConstantW);
	void SetBCConstantPE(double BC_ConstantE);
	void SetBCConstantPS(double BC_ConstantS);
	void SetBCConstantPN(double BC_ConstantN);
	void SetBCConstantPB(double BC_ConstantS);
	void SetBCConstantPT(double BC_ConstantN);

	int sign(const double& val);
	int SetPLTType(PLTTYPE type);
	int OutRes(const int iter, const double curTime, const std::string fname_vel_base, const std::string fname_div_base,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
		const std::vector<double>& ps, const std::vector<double>& ls);
	int OutResClose();

	int idx(int i, int j, int k);
};

#endif __MAC3D_H_
