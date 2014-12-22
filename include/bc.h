#ifndef __BC_H_
#define __BC_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <limits>

#include "common.h"

class BoundaryCondition {
public:
	BoundaryCondition(int nr, int nz, int num_bc_grid);
	BoundaryCondition(int nx, int ny, int nz, int num_bc_grid);

	// for 2D axisymmetric
	const int kNr, kNz;
	// for 3D cartesian grid, kNz already declared
	const int kNx, kNy;
	const int kNumBCGrid;

	// Boundary Condition Variables
	
	BC m_BC_UW, m_BC_UE, m_BC_US, m_BC_UN, m_BC_UB, m_BC_UT;
	BC m_BC_VW, m_BC_VE, m_BC_VS, m_BC_VN, m_BC_VB, m_BC_VT;
	BC m_BC_WW, m_BC_WE, m_BC_WS, m_BC_WN, m_BC_WB, m_BC_WT;
	BC m_BC_PW, m_BC_PE, m_BC_PS, m_BC_PN, m_BC_PB, m_BC_PT;
	double m_BC_DirichletConstantUW, m_BC_DirichletConstantUE,
		m_BC_DirichletConstantUS, m_BC_DirichletConstantUN,
		m_BC_DirichletConstantUB, m_BC_DirichletConstantUT;
	double m_BC_DirichletConstantVW, m_BC_DirichletConstantVE,
		m_BC_DirichletConstantVS, m_BC_DirichletConstantVN,
		m_BC_DirichletConstantVB, m_BC_DirichletConstantVT;
	double m_BC_DirichletConstantWW, m_BC_DirichletConstantWE,
		m_BC_DirichletConstantWS, m_BC_DirichletConstantWN,
		m_BC_DirichletConstantWB, m_BC_DirichletConstantWT;
	double m_BC_DirichletConstantPW, m_BC_DirichletConstantPE,
		m_BC_DirichletConstantPS, m_BC_DirichletConstantPN,
		m_BC_DirichletConstantPB, m_BC_DirichletConstantPT;

	// BC
	/*
	For staggered grid
	*_U_* for U grid
	*_V_* for V grid
	*_W_* for W grid
	*_P_* for P grid
	*/
	int SetBC_U_3D(std::string BC_W, std::string BC_E,
		std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_V_3D(std::string BC_W, std::string BC_E,
		std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_W_3D(std::string BC_W, std::string BC_E,
		std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int SetBC_P_3D(std::string BC_W, std::string BC_E,
		std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T);
	int ApplyBC_U_3D(double ***arr);
	int ApplyBC_V_3D(double ***arr);
	int ApplyBC_W_3D(double ***arr);
	int ApplyBC_P_3D(double ***arr);

	void BC_UW(double ***arr);
	void BC_UE(double ***arr);
	void BC_US(double ***arr);
	void BC_UN(double ***arr);
	void BC_UB(double ***arr);
	void BC_UT(double ***arr);
	void BC_VW(double ***arr);
	void BC_VE(double ***arr);
	void BC_VS(double ***arr);
	void BC_VN(double ***arr);
	void BC_VB(double ***arr);
	void BC_VT(double ***arr);
	void BC_WW(double ***arr);
	void BC_WE(double ***arr);
	void BC_WS(double ***arr);
	void BC_WN(double ***arr);
	void BC_WB(double ***arr);
	void BC_WT(double ***arr);
	void BC_PW(double ***arr);
	void BC_PE(double ***arr);
	void BC_PS(double ***arr);
	void BC_PN(double ***arr);
	void BC_PB(double ***arr);
	void BC_PT(double ***arr);

	// Prescribed bounadry condition (2D Periodic)
	void BC_PeriodicUW(double ***arr);
	void BC_PeriodicUE(double ***arr);
	void BC_PeriodicUS(double ***arr);
	void BC_PeriodicUN(double ***arr);
	void BC_PeriodicUB(double ***arr);
	void BC_PeriodicUT(double ***arr);
	void BC_PeriodicVW(double ***arr);
	void BC_PeriodicVE(double ***arr);
	void BC_PeriodicVS(double ***arr);
	void BC_PeriodicVN(double ***arr);
	void BC_PeriodicVB(double ***arr);
	void BC_PeriodicVT(double ***arr);
	void BC_PeriodicWW(double ***arr);
	void BC_PeriodicWE(double ***arr);
	void BC_PeriodicWS(double ***arr);
	void BC_PeriodicWN(double ***arr);
	void BC_PeriodicWB(double ***arr);
	void BC_PeriodicWT(double ***arr);
	void BC_PeriodicPW(double ***arr);
	void BC_PeriodicPE(double ***arr);
	void BC_PeriodicPS(double ***arr);
	void BC_PeriodicPN(double ***arr);
	void BC_PeriodicPB(double ***arr);
	void BC_PeriodicPT(double ***arr);
	// Prescribed bounadry condition (2D Neumann)
	void BC_NeumannUW(double ***arr);
	void BC_NeumannUE(double ***arr);
	void BC_NeumannUS(double ***arr);
	void BC_NeumannUN(double ***arr);
	void BC_NeumannUB(double ***arr);
	void BC_NeumannUT(double ***arr);
	void BC_NeumannVW(double ***arr);
	void BC_NeumannVE(double ***arr);
	void BC_NeumannVS(double ***arr);
	void BC_NeumannVN(double ***arr);
	void BC_NeumannVB(double ***arr);
	void BC_NeumannVT(double ***arr);
	void BC_NeumannWW(double ***arr);
	void BC_NeumannWE(double ***arr);
	void BC_NeumannWS(double ***arr);
	void BC_NeumannWN(double ***arr);
	void BC_NeumannWB(double ***arr);
	void BC_NeumannWT(double ***arr);
	void BC_NeumannPW(double ***arr);
	void BC_NeumannPE(double ***arr);
	void BC_NeumannPS(double ***arr);
	void BC_NeumannPN(double ***arr);
	void BC_NeumannPB(double ***arr);
	void BC_NeumannPT(double ***arr);
	// Prescribed bounadry condition (2D Dirichlet)
	void BC_DirichletUW(double ***arr);
	void BC_DirichletUE(double ***arr);
	void BC_DirichletUS(double ***arr);
	void BC_DirichletUN(double ***arr);
	void BC_DirichletUB(double ***arr);
	void BC_DirichletUT(double ***arr);
	void BC_DirichletVW(double ***arr);
	void BC_DirichletVE(double ***arr);
	void BC_DirichletVS(double ***arr);
	void BC_DirichletVN(double ***arr);
	void BC_DirichletVB(double ***arr);
	void BC_DirichletVT(double ***arr);
	void BC_DirichletWW(double ***arr);
	void BC_DirichletWE(double ***arr);
	void BC_DirichletWS(double ***arr);
	void BC_DirichletWN(double ***arr);
	void BC_DirichletWB(double ***arr);
	void BC_DirichletWT(double ***arr);
	void BC_DirichletPW(double ***arr);
	void BC_DirichletPE(double ***arr);
	void BC_DirichletPS(double ***arr);
	void BC_DirichletPN(double ***arr);
	void BC_DirichletPB(double ***arr);
	void BC_DirichletPT(double ***arr);
	void SetBCConstantUW(double BC_ConstantW);
	void SetBCConstantUE(double BC_ConstantE);
	void SetBCConstantUS(double BC_ConstantS);
	void SetBCConstantUN(double BC_ConstantN);
	void SetBCConstantUB(double BC_ConstantB);
	void SetBCConstantUT(double BC_ConstantT);
	void SetBCConstantVW(double BC_ConstantW);
	void SetBCConstantVE(double BC_ConstantE);
	void SetBCConstantVS(double BC_ConstantS);
	void SetBCConstantVN(double BC_ConstantN);
	void SetBCConstantVB(double BC_ConstantB);
	void SetBCConstantVT(double BC_ConstantT);
	void SetBCConstantWW(double BC_ConstantW);
	void SetBCConstantWE(double BC_ConstantE);
	void SetBCConstantWS(double BC_ConstantS);
	void SetBCConstantWN(double BC_ConstantN);
	void SetBCConstantWB(double BC_ConstantB);
	void SetBCConstantWT(double BC_ConstantT);
	void SetBCConstantPW(double BC_ConstantW);
	void SetBCConstantPE(double BC_ConstantE);
	void SetBCConstantPS(double BC_ConstantS);
	void SetBCConstantPN(double BC_ConstantN);
	void SetBCConstantPB(double BC_ConstantB);
	void SetBCConstantPT(double BC_ConstantT);
};

#endif __BC_H_