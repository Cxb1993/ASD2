#ifndef __BC3D_H_
#define __BC3D_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>

#include <array>
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "common.h"

class BoundaryCondition3D {
public:
	BoundaryCondition3D(int nx, int ny, int nz, int num_bc_grid);

	const int64_t kNx, kNy, kNz;
	const int kNumBCGrid;

	// Boundary Condition Variables
	
	BC3D m_BC_UW, m_BC_UE, m_BC_US, m_BC_UN, m_BC_UB, m_BC_UT;
	BC3D m_BC_VW, m_BC_VE, m_BC_VS, m_BC_VN, m_BC_VB, m_BC_VT;
	BC3D m_BC_WW, m_BC_WE, m_BC_WS, m_BC_WN, m_BC_WB, m_BC_WT;
	BC3D m_BC_PW, m_BC_PE, m_BC_PS, m_BC_PN, m_BC_PB, m_BC_PT;
	double m_BC_DirichletConstantUW, m_BC_DirichletConstantUE, m_BC_DirichletConstantUS, m_BC_DirichletConstantUN, m_BC_DirichletConstantUB, m_BC_DirichletConstantUT;
	double m_BC_DirichletConstantVW, m_BC_DirichletConstantVE, m_BC_DirichletConstantVS, m_BC_DirichletConstantVN, m_BC_DirichletConstantVB, m_BC_DirichletConstantVT;
	double m_BC_DirichletConstantWW, m_BC_DirichletConstantWE, m_BC_DirichletConstantWS, m_BC_DirichletConstantWN, m_BC_DirichletConstantWB, m_BC_DirichletConstantWT;
	double m_BC_DirichletConstantPW, m_BC_DirichletConstantPE, m_BC_DirichletConstantPS, m_BC_DirichletConstantPN, m_BC_DirichletConstantPB, m_BC_DirichletConstantPT;

	inline int idx(int i, int j, int k);
	// BC
	/*
	For staggered grid
	*_U_* for U grid
	*_V_* for V grid
	*_P_* for P grid
	*/
	int SetBC_U_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_V_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_W_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_P_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int ApplyBC_U_3D(std::vector<double>& arr);
	int ApplyBC_V_3D(std::vector<double>& arr);
	int ApplyBC_W_3D(std::vector<double>& arr);
	int ApplyBC_P_3D(std::vector<double>& arr);

	void BC_UW(std::vector<double>& arr);
	void BC_UE(std::vector<double>& arr);
	void BC_US(std::vector<double>& arr);
	void BC_UN(std::vector<double>& arr);
	void BC_UB(std::vector<double>& arr);
	void BC_UT(std::vector<double>& arr);
	void BC_VW(std::vector<double>& arr);
	void BC_VE(std::vector<double>& arr);
	void BC_VS(std::vector<double>& arr);
	void BC_VN(std::vector<double>& arr);
	void BC_VB(std::vector<double>& arr);
	void BC_VT(std::vector<double>& arr);
	void BC_WW(std::vector<double>& arr);
	void BC_WE(std::vector<double>& arr);
	void BC_WS(std::vector<double>& arr);
	void BC_WN(std::vector<double>& arr);
	void BC_WB(std::vector<double>& arr);
	void BC_WT(std::vector<double>& arr);
	void BC_PW(std::vector<double>& arr);
	void BC_PE(std::vector<double>& arr);
	void BC_PS(std::vector<double>& arr);
	void BC_PN(std::vector<double>& arr);
	void BC_PB(std::vector<double>& arr);
	void BC_PT(std::vector<double>& arr);

	// Prescribed bounadry condition (3D Periodic)
	void BC_PeriodicUW(std::vector<double>& arr);
	void BC_PeriodicUE(std::vector<double>& arr);
	void BC_PeriodicUS(std::vector<double>& arr);
	void BC_PeriodicUN(std::vector<double>& arr);
	void BC_PeriodicUB(std::vector<double>& arr);
	void BC_PeriodicUT(std::vector<double>& arr);
	void BC_PeriodicVW(std::vector<double>& arr);
	void BC_PeriodicVE(std::vector<double>& arr);
	void BC_PeriodicVS(std::vector<double>& arr);
	void BC_PeriodicVN(std::vector<double>& arr);
	void BC_PeriodicVB(std::vector<double>& arr);
	void BC_PeriodicVT(std::vector<double>& arr);
	void BC_PeriodicWW(std::vector<double>& arr);
	void BC_PeriodicWE(std::vector<double>& arr);
	void BC_PeriodicWS(std::vector<double>& arr);
	void BC_PeriodicWB(std::vector<double>& arr);
	void BC_PeriodicWT(std::vector<double>& arr);
	void BC_PeriodicWN(std::vector<double>& arr);
	void BC_PeriodicPW(std::vector<double>& arr);
	void BC_PeriodicPE(std::vector<double>& arr);
	void BC_PeriodicPS(std::vector<double>& arr);
	void BC_PeriodicPN(std::vector<double>& arr);
	void BC_PeriodicPB(std::vector<double>& arr);
	void BC_PeriodicPT(std::vector<double>& arr);

	// Prescribed bounadry condition (3D Neumann)
	void BC_NeumannUW(std::vector<double>& arr);
	void BC_NeumannUE(std::vector<double>& arr);
	void BC_NeumannUS(std::vector<double>& arr);
	void BC_NeumannUN(std::vector<double>& arr);
	void BC_NeumannUB(std::vector<double>& arr);
	void BC_NeumannUT(std::vector<double>& arr);
	void BC_NeumannVW(std::vector<double>& arr);
	void BC_NeumannVE(std::vector<double>& arr);
	void BC_NeumannVS(std::vector<double>& arr);
	void BC_NeumannVN(std::vector<double>& arr);
	void BC_NeumannVB(std::vector<double>& arr);
	void BC_NeumannVT(std::vector<double>& arr);
	void BC_NeumannWW(std::vector<double>& arr);
	void BC_NeumannWE(std::vector<double>& arr);
	void BC_NeumannWS(std::vector<double>& arr);
	void BC_NeumannWN(std::vector<double>& arr);
	void BC_NeumannWB(std::vector<double>& arr);
	void BC_NeumannWT(std::vector<double>& arr);
	void BC_NeumannPW(std::vector<double>& arr);
	void BC_NeumannPE(std::vector<double>& arr);
	void BC_NeumannPS(std::vector<double>& arr);
	void BC_NeumannPN(std::vector<double>& arr);
	void BC_NeumannPB(std::vector<double>& arr);
	void BC_NeumannPT(std::vector<double>& arr);
	// Prescribed bounadry condition (3D Dirichlet)
	void BC_DirichletUW(std::vector<double>& arr);
	void BC_DirichletUE(std::vector<double>& arr);
	void BC_DirichletUS(std::vector<double>& arr);
	void BC_DirichletUN(std::vector<double>& arr);
	void BC_DirichletUB(std::vector<double>& arr);
	void BC_DirichletUT(std::vector<double>& arr);
	void BC_DirichletVW(std::vector<double>& arr);
	void BC_DirichletVE(std::vector<double>& arr);
	void BC_DirichletVS(std::vector<double>& arr);
	void BC_DirichletVN(std::vector<double>& arr);
	void BC_DirichletVB(std::vector<double>& arr);
	void BC_DirichletVT(std::vector<double>& arr);
	void BC_DirichletWW(std::vector<double>& arr);
	void BC_DirichletWE(std::vector<double>& arr);
	void BC_DirichletWS(std::vector<double>& arr);
	void BC_DirichletWN(std::vector<double>& arr);
	void BC_DirichletWB(std::vector<double>& arr);
	void BC_DirichletWT(std::vector<double>& arr);
	void BC_DirichletPW(std::vector<double>& arr);
	void BC_DirichletPE(std::vector<double>& arr);
	void BC_DirichletPS(std::vector<double>& arr);
	void BC_DirichletPN(std::vector<double>& arr);
	void BC_DirichletPB(std::vector<double>& arr);
	void BC_DirichletPT(std::vector<double>& arr);
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
};

#endif __BC3D_H_