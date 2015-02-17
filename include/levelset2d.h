#ifndef __LEVELSET2D_H__
#define __LEVELSET2D_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <cassert>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

#include <algorithm>

#include "common.h"
#include "bc2d.h"

class LevelSetSolver2D {
public:
	LevelSetSolver2D(int nx, int ny, int num_bc_grid, double dx, double dy);
	
	std::vector<double> HJWENO5_LS_2D(std::vector<double>& ls,
		const std::vector<double>& u, const std::vector<double>& v);
	
	int Solve_LevelSet_2D(std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v, double dt);

	int Reinit_Original_2D(std::vector<double>& ls);
	int Reinit_Sussman_2D(std::vector<double>& ls);
	std::vector<double> HJENO_ReinitABS_2D(std::vector<double>& ls);

	int sign(const double& val);
	std::vector<int> GetSignedLSNormalized(const std::vector<double>& ls);
	double NewtonDivdedDifference_X(int order, const std::vector<double> val, int i, int j, int lI, int rI);
	double NewtonDivdedDifference_Y(int order, const std::vector<double> val, int i, int j, int lJ, int rJ);

	int UnitHJWENO5(
		const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n);

	// BC
	int SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
	int ApplyBC_U_2D(std::vector<double>& arr);
	int ApplyBC_V_2D(std::vector<double>& arr);
	int ApplyBC_P_2D(std::vector<double>& arr);
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

	double minmod(double a, double b);
	int idx(int i, int j);
	
	const double kEps, kThickness, kMaxATime;
	const int kNx, kNy, kNumBCGrid;
	const double kDx, kDy;
	std::vector<int> m_signedInitLS;

	// artificial dt for reinitialization only 
	double m_atime;
	const double kAdt;

	std::shared_ptr<BoundaryCondition2D> m_BC;
};

#endif __LEVELSET2D_H__