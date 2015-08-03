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
    LevelSetSolver2D(int nx, int ny, int num_bc_grid,
        double baseX, double baseY, double dx, double dy);
    LevelSetSolver2D(int nx, int ny, int num_bc_grid,
        double baseX, double baseY, double dx, double dy, double maxTime);
    
    std::vector<double> UpdateHeavisideFuncDeriv(const std::vector<double>& ls);

    int Solve_LevelSet_2D(std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v, double dt);

    int Reinit_Original_2D(std::vector<double>& ls);
    int Reinit_Sussman_2D(std::vector<double>& ls);
    int Reinit_Sussman_FirstTime_2D(std::vector<double>& ls);
    int FirstTimeOnlyReinit_Sussman_2D(std::vector<double>& ls);
    int Reinit_MinRK2_2D(std::vector<double>& ls);
    
    std::vector<double> GetSussmanReinitConstraint(const std::vector<double>& ls,
        const std::vector<double>& lsInit, const std::vector<double>& heavisideDeriv);

    int sign(const double& val);
    double MinAbs(const double& val1, const double& val2);
    double MinMod(const double& val1, const double& val2);
    std::vector<int> GetSignedLSNormalized(const std::vector<double>& ls);
    std::vector<double> GetSmoothedSignFunc(const std::vector<double> &ls);

    std::vector<double> HJWENO5_LS_2D(std::vector<double>& ls,
        const std::vector<double>& u, const std::vector<double>& v);
    int UnitHJWENO5(
        const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n);

    std::vector<double> HJWENO5_DerivAbsLS_2D(const std::vector<double>& ls, const std::vector<double>& lsInit);
    std::vector<double>
        ENO_DerivAbsLS_2D(const std::vector<double>& ls, const std::vector<double>& lsInit);
    std::vector<double>
        SubcellENO_Min_DerivAbsLS_2D(const std::vector<double>& ls, const std::vector<double>& lsInit);

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

    int idx(int i, int j);
    
    const double kEps, kThickness, kMaxATime;
    const int kNx, kNy, kNumBCGrid, kENOSpatialOrder;
    const int64_t kArrSize;
    const double kBaseX, kBaseY;
    const double kDx, kDy;
    
    // artificial dt for reinitialization only 
    double m_atime;
    const double kAdt;

    std::shared_ptr<BoundaryCondition2D> m_BC;
};

#endif __LEVELSET2D_H__