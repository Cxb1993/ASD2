#ifndef __BC2D_H_
#define __BC2D_H_

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
#include <vector>

#include "common.h"

class BoundaryCondition2D {
public:
    BoundaryCondition2D(int nx, int ny, int num_bc_grid);
    BoundaryCondition2D(int nx, int ny, double dx, double dy, int num_bc_grid);

    const int kNx, kNy;
    const int kNumBCGrid;
    const double kDx, kDy;
    double m_AmbientPressure;
    // Boundary Condition Variables
    
    BC2D m_BC_UW, m_BC_UE, m_BC_US, m_BC_UN;
    BC2D m_BC_VW, m_BC_VE, m_BC_VS, m_BC_VN;
    BC2D m_BC_PW, m_BC_PE, m_BC_PS, m_BC_PN;
    double m_BC_DirichletConstantUW, m_BC_DirichletConstantUE, m_BC_DirichletConstantUS, m_BC_DirichletConstantUN;
    double m_BC_DirichletConstantVW, m_BC_DirichletConstantVE, m_BC_DirichletConstantVS, m_BC_DirichletConstantVN;
    double m_BC_DirichletConstantPW, m_BC_DirichletConstantPE, m_BC_DirichletConstantPS, m_BC_DirichletConstantPN;

    inline int idx(int i, int j);
    // BC
    /*
    For staggered grid
    *_U_* for U grid
    *_V_* for V grid
    *_P_* for P grid
    */
    int SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
    int SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
    int SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N);
    int ApplyBC_U_2D(std::vector<double>& arr);
    int ApplyBC_V_2D(std::vector<double>& arr);
    int ApplyBC_P_2D(std::vector<double>& arr);

    void BC_UW(std::vector<double>& arr);
    void BC_UE(std::vector<double>& arr);
    void BC_US(std::vector<double>& arr);
    void BC_UN(std::vector<double>& arr);
    void BC_VW(std::vector<double>& arr);
    void BC_VE(std::vector<double>& arr);
    void BC_VS(std::vector<double>& arr);
    void BC_VN(std::vector<double>& arr);
    void BC_PW(std::vector<double>& arr);
    void BC_PE(std::vector<double>& arr);
    void BC_PS(std::vector<double>& arr);
    void BC_PN(std::vector<double>& arr);

    // Prescribed bounadry condition (2D Periodic)
    void BC_PeriodicUW(std::vector<double>& arr);
    void BC_PeriodicUE(std::vector<double>& arr);
    void BC_PeriodicUS(std::vector<double>& arr);
    void BC_PeriodicUN(std::vector<double>& arr);
    void BC_PeriodicVW(std::vector<double>& arr);
    void BC_PeriodicVE(std::vector<double>& arr);
    void BC_PeriodicVS(std::vector<double>& arr);
    void BC_PeriodicVN(std::vector<double>& arr);
    void BC_PeriodicPW(std::vector<double>& arr);
    void BC_PeriodicPE(std::vector<double>& arr);
    void BC_PeriodicPS(std::vector<double>& arr);
    void BC_PeriodicPN(std::vector<double>& arr);

    // Prescribed bounadry condition (2D Neumann)
    void BC_NeumannUW(std::vector<double>& arr);
    void BC_NeumannUE(std::vector<double>& arr);
    void BC_NeumannUS(std::vector<double>& arr);
    void BC_NeumannUN(std::vector<double>& arr);
    void BC_NeumannVW(std::vector<double>& arr);
    void BC_NeumannVE(std::vector<double>& arr);
    void BC_NeumannVS(std::vector<double>& arr);
    void BC_NeumannVN(std::vector<double>& arr);
    void BC_NeumannPW(std::vector<double>& arr);
    void BC_NeumannPE(std::vector<double>& arr);
    void BC_NeumannPS(std::vector<double>& arr);
    void BC_NeumannPN(std::vector<double>& arr);
    // Prescribed bounadry condition (2D Dirichlet)
    void BC_DirichletUW(std::vector<double>& arr);
    void BC_DirichletUE(std::vector<double>& arr);
    void BC_DirichletUS(std::vector<double>& arr);
    void BC_DirichletUN(std::vector<double>& arr);
    void BC_DirichletVW(std::vector<double>& arr);
    void BC_DirichletVE(std::vector<double>& arr);
    void BC_DirichletVS(std::vector<double>& arr);
    void BC_DirichletVN(std::vector<double>& arr);
    void BC_DirichletPW(std::vector<double>& arr);
    void BC_DirichletPE(std::vector<double>& arr);
    void BC_DirichletPS(std::vector<double>& arr);
    void BC_DirichletPN(std::vector<double>& arr);
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
    
    // Free boundary condition for Level Set
    // Extrpolate level set
    void BC_LSFreeBoundaryPW(std::vector<double>& arr);
    void BC_LSFreeBoundaryPE(std::vector<double>& arr);
    void BC_LSFreeBoundaryPS(std::vector<double>& arr);
    void BC_LSFreeBoundaryPN(std::vector<double>& arr);

    void SetAmbientPressure(double ambientPressure);
};

#endif __BC2D_H_