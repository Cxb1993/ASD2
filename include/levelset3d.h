#ifndef __LEVELSET3D_H__
#define __LEVELSET3D_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include <algorithm>

#include "common.h"
#include "bc3d.h"

class LevelSetSolver3D {
public:
	LevelSetSolver3D(int nx, int ny, int nz, int num_bc_grid,
		double baseX, double baseY, double baseZ, double dx, double dy, double dz);
	LevelSetSolver3D(int nx, int ny, int nz, int num_bc_grid,
		double baseX, double baseY, double baseZ, double dx, double dy, double dz, double maxTime);
	
	std::vector<double> UpdateHeavisideFuncDeriv(const std::vector<double>& ls);

	int Solve_LevelSet_3D(std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w, double dt);

	int Reinit_Original_3D(std::vector<double>& ls);
	int Reinit_Sussman_3D(std::vector<double>& ls);
	int Reinit_Sussman_FirstTime_3D(std::vector<double>& ls);
	int FirstTimeOnlyReinit_Sussman_3D(std::vector<double>& ls);
	int Reinit_MinRK2_3D(std::vector<double>& ls);
	
	std::vector<double> GetSussmanReinitConstraint(const std::vector<double>& ls,
		const std::vector<double>& lsInit, const std::vector<double>& heavisideDeriv);

	int sign(const double& val);
	double MinAbs(const double& val1, const double& val2);
	double MinMod(const double& val1, const double& val2);
	std::vector<double> GetSmoothedSignFunc(const std::vector<double> &ls);
	std::vector<int> GetSignedLSNormalized(const std::vector<double>& ls);

	std::vector<double> HJWENO5_LS_3D(std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w);
	int UnitHJWENO5(
		const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n);
	std::vector<double>
		ENO_DerivAbsLS_3D(const std::vector<double>& ls, const std::vector<double>& lsInit);
	std::vector<double>
		Subcell_DerivAbsLS_3D(const std::vector<double>& ls, const std::vector<double>& lsInit);

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

	int idx(int i, int j, int k);
	
	const double kEps, kThickness, kMaxATime;
	const int kNx, kNy, kNz;
	const int kNumBCGrid, kENOSpatialOrder;
	const int64_t kArrSize;
	const double kBaseX, kBaseY, kBaseZ;
	const double kDx, kDy, kDz;
	
	// artificial dt for reinitialization only 
	double m_atime;
	const double kAdt;

	std::shared_ptr<BoundaryCondition3D> m_BC;
};

#endif __LEVELSET3D_H__