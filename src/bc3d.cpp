#include "bc3d.h"

BoundaryCondition3D::BoundaryCondition3D(int nx, int ny, int nz, int num_bc_grid)
    : kNx(nx), kNy(ny), kNz(nz), kDx(-1.0), kDy(-1.0), kDz(-1.0), kNumBCGrid(num_bc_grid) {
}

BoundaryCondition3D::BoundaryCondition3D(int nx, int ny, int nz, double dx, double dy, double dz, int num_bc_grid)
    : kNx(nx), kNy(ny), kNz(nz), kDx(dx), kDy(dy), kDz(dz), kNumBCGrid(num_bc_grid) {
}

inline int BoundaryCondition3D::idx(int i, int j, int k) {
    return (static_cast<int64_t>(k) + (kNz + 2 * kNumBCGrid) * static_cast<int64_t>(j) + (kNz + 2 * kNumBCGrid) * (kNy + 2 * kNumBCGrid) * static_cast<int64_t>(i));
}

int BoundaryCondition3D::SetBC_U_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
    std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);
    
    if (BC_W == "periodic")
        m_BC_UW = BC3D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_UW = BC3D::NEUMANN;
    else if (BC_W == "wall" || BC_W == "dirichlet")
        m_BC_UW = BC3D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_UW = BC3D::INLET;
    else if (BC_W == "outlet")
        m_BC_UW = BC3D::OUTLET;
    // not supported
    else
        m_BC_UW = BC3D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_UE = BC3D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_UE = BC3D::NEUMANN;
    else if (BC_E == "wall" || BC_E == "dirichlet")
        m_BC_UE = BC3D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_UE = BC3D::INLET;
    else if (BC_E == "outlet")
        m_BC_UE = BC3D::OUTLET;
    // not supported
    else
        m_BC_UE = BC3D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_US = BC3D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_US = BC3D::NEUMANN;
    else if (BC_S == "wall" || BC_S == "dirichlet")
        m_BC_US = BC3D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_US = BC3D::INLET;
    else if (BC_S == "outlet")
        m_BC_US = BC3D::OUTLET;
    // not supported
    else
        m_BC_US = BC3D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_UN = BC3D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_UN = BC3D::NEUMANN;
    else if (BC_N == "wall" || BC_N == "dirichlet")
        m_BC_UN = BC3D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_UN = BC3D::INLET;
    else if (BC_N == "outlet")
        m_BC_UN = BC3D::OUTLET;
    // not supported
    else
        m_BC_UN = BC3D::CUSTOM;

    if (BC_B == "periodic")
        m_BC_UB = BC3D::PERIODIC;
    else if (BC_B == "neumann")
        m_BC_UB = BC3D::NEUMANN;
    else if (BC_B == "wall" || BC_B == "dirichlet")
        m_BC_UB = BC3D::DIRICHLET;
    else if (BC_B == "inlet")
        m_BC_UB = BC3D::INLET;
    else if (BC_B == "outlet")
        m_BC_UB = BC3D::OUTLET;
    // not supported
    else
        m_BC_UB = BC3D::CUSTOM;

    if (BC_T == "periodic")
        m_BC_UT = BC3D::PERIODIC;
    else if (BC_T == "neumann")
        m_BC_UT = BC3D::NEUMANN;
    else if (BC_T == "wall" || BC_T == "dirichlet")
        m_BC_UT = BC3D::DIRICHLET;
    else if (BC_T == "inlet")
        m_BC_UT = BC3D::INLET;
    else if (BC_T == "outlet")
        m_BC_UT = BC3D::OUTLET;
    // not supported
    else
        m_BC_UT = BC3D::CUSTOM;

    return 0;
}

int BoundaryCondition3D::ApplyBC_U_3D(std::vector<double>& arr) {
    BC_UW(arr);
    BC_UE(arr);
    BC_US(arr);
    BC_UN(arr);
    BC_UB(arr);
    BC_UT(arr);
    
    return 0;
}

int BoundaryCondition3D::SetBC_V_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
    std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);
    
    if (BC_W == "periodic")
        m_BC_VW = BC3D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_VW = BC3D::NEUMANN;
    else if (BC_W == "wall" || BC_W == "dirichlet")
        m_BC_VW = BC3D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_VW = BC3D::INLET;
    else if (BC_W == "outlet")
        m_BC_VW = BC3D::OUTLET;
    // not supported
    else
        m_BC_VW = BC3D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_VE = BC3D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_VE = BC3D::NEUMANN;
    else if (BC_E == "wall" || BC_E == "dirichlet")
        m_BC_VE = BC3D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_VE = BC3D::INLET;
    else if (BC_E == "outlet")
        m_BC_VE = BC3D::OUTLET;
    // not supported
    else
        m_BC_VE = BC3D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_VS = BC3D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_VS = BC3D::NEUMANN;
    else if (BC_S == "wall" || BC_S == "dirichlet")
        m_BC_VS = BC3D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_VS = BC3D::INLET;
    else if (BC_S == "outlet")
        m_BC_VS = BC3D::OUTLET;
    // not supported
    else
        m_BC_VS = BC3D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_VN = BC3D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_VN = BC3D::NEUMANN;
    else if (BC_N == "wall" || BC_N == "dirichlet")
        m_BC_VN = BC3D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_VN = BC3D::INLET;
    else if (BC_N == "outlet")
        m_BC_VN = BC3D::OUTLET;
    // not supported
    else
        m_BC_VN = BC3D::CUSTOM;

    if (BC_B == "periodic")
        m_BC_VB = BC3D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_VB = BC3D::NEUMANN;
    else if (BC_B == "wall" || BC_B == "dirichlet")
        m_BC_VB = BC3D::DIRICHLET;
    else if (BC_B == "inlet")
        m_BC_VB = BC3D::INLET;
    else if (BC_B == "outlet")
        m_BC_VB = BC3D::OUTLET;
    // not supported
    else
        m_BC_VB = BC3D::CUSTOM;

    if (BC_T == "periodic")
        m_BC_VT = BC3D::PERIODIC;
    else if (BC_T == "neumann")
        m_BC_VT = BC3D::NEUMANN;
    else if (BC_T == "wall" || BC_T == "dirichlet")
        m_BC_VT = BC3D::DIRICHLET;
    else if (BC_T == "inlet")
        m_BC_VT = BC3D::INLET;
    else if (BC_T == "outlet")
        m_BC_VT = BC3D::OUTLET;
    // not supported
    else
        m_BC_VT = BC3D::CUSTOM;

    return 0;
}

int BoundaryCondition3D::ApplyBC_V_3D(std::vector<double>& arr) {
    BC_VW(arr);
    BC_VE(arr);
    BC_VS(arr);
    BC_VN(arr);
    BC_VB(arr);
    BC_VT(arr);

    return 0;
}

int BoundaryCondition3D::SetBC_W_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
    std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);

    if (BC_W == "periodic")
        m_BC_WW = BC3D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_WW = BC3D::NEUMANN;
    else if (BC_W == "wall" || BC_W == "dirichlet")
        m_BC_WW = BC3D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_WW = BC3D::INLET;
    else if (BC_W == "outlet")
        m_BC_WW = BC3D::OUTLET;
    // not supported
    else
        m_BC_WW = BC3D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_WE = BC3D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_WE = BC3D::NEUMANN;
    else if (BC_E == "wall" || BC_E == "dirichlet")
        m_BC_WE = BC3D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_WE = BC3D::INLET;
    else if (BC_E == "outlet")
        m_BC_WE = BC3D::OUTLET;
    // not supported
    else
        m_BC_WE = BC3D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_WS = BC3D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_WS = BC3D::NEUMANN;
    else if (BC_S == "wall" || BC_S == "dirichlet")
        m_BC_WS = BC3D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_WS = BC3D::INLET;
    else if (BC_S == "outlet")
        m_BC_WS = BC3D::OUTLET;
    // not supported
    else
        m_BC_WS = BC3D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_WN = BC3D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_WN = BC3D::NEUMANN;
    else if (BC_N == "wall" || BC_N == "dirichlet")
        m_BC_WN = BC3D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_WN = BC3D::INLET;
    else if (BC_N == "outlet")
        m_BC_WN = BC3D::OUTLET;
    // not supported
    else
        m_BC_WN = BC3D::CUSTOM;

    if (BC_B == "periodic")
        m_BC_WB = BC3D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_WB = BC3D::NEUMANN;
    else if (BC_B == "wall" || BC_B == "dirichlet")
        m_BC_WB = BC3D::DIRICHLET;
    else if (BC_B == "inlet")
        m_BC_WB = BC3D::INLET;
    else if (BC_B == "outlet")
        m_BC_WB = BC3D::OUTLET;
    // not supported
    else
        m_BC_WB = BC3D::CUSTOM;

    if (BC_T == "periodic")
        m_BC_WT = BC3D::PERIODIC;
    else if (BC_T == "neumann")
        m_BC_WT = BC3D::NEUMANN;
    else if (BC_T == "wall" || BC_T == "dirichlet")
        m_BC_WT = BC3D::DIRICHLET;
    else if (BC_T == "inlet")
        m_BC_WT = BC3D::INLET;
    else if (BC_T == "outlet")
        m_BC_WT = BC3D::OUTLET;
    // not supported
    else
        m_BC_PT = BC3D::CUSTOM;

    return 0;
}

int BoundaryCondition3D::ApplyBC_W_3D(std::vector<double>& arr) {
    BC_WW(arr);
    BC_WE(arr);
    BC_WS(arr);
    BC_WN(arr);
    BC_WB(arr);
    BC_WT(arr);

    return 0;
}

int BoundaryCondition3D::SetBC_P_3D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
    std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);

    if (BC_W == "periodic")
        m_BC_PW = BC3D::PERIODIC;
    else if (BC_W == "neumann" || BC_W == "wall")
        m_BC_PW = BC3D::NEUMANN;
    else if (BC_W == "dirichlet")
        m_BC_PW = BC3D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_PW = BC3D::INLET;
    else if (BC_W == "outlet")
        m_BC_PW = BC3D::OUTLET;
    else if (BC_W == "pressure")
        m_BC_PW = BC3D::PRESSURE;
    else if (BC_W == "lsfree")
        m_BC_PW = BC3D::LSFREE;
    // not supported
    else
        m_BC_PW = BC3D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_PE = BC3D::PERIODIC;
    else if (BC_E == "neumann" || BC_E == "wall")
        m_BC_PE = BC3D::NEUMANN;
    else if (BC_E == "dirichlet")
        m_BC_PE = BC3D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_PE = BC3D::INLET;
    else if (BC_E == "outlet")
        m_BC_PE = BC3D::OUTLET;
    else if (BC_E == "pressure")
        m_BC_PE = BC3D::PRESSURE;
    else if (BC_E == "lsfree")
        m_BC_PE = BC3D::LSFREE;
    // not supported
    else
        m_BC_PE = BC3D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_PS = BC3D::PERIODIC;
    else if (BC_S == "neumann" || BC_S == "wall")
        m_BC_PS = BC3D::NEUMANN;
    else if (BC_S == "dirichlet")
        m_BC_PS = BC3D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_PS = BC3D::INLET;
    else if (BC_S == "outlet")
        m_BC_PS = BC3D::OUTLET;
    else if (BC_S == "pressure")
        m_BC_PS = BC3D::PRESSURE;
    else if (BC_S == "lsfree")
        m_BC_PS = BC3D::LSFREE;
    // not supported
    else
        m_BC_PS = BC3D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_PN = BC3D::PERIODIC;
    else if (BC_N == "neumann" || BC_N == "wall")
        m_BC_PN = BC3D::NEUMANN;
    else if (BC_N == "dirichlet")
        m_BC_PN = BC3D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_PN = BC3D::INLET;
    else if (BC_N == "outlet")
        m_BC_PN = BC3D::OUTLET;
    else if (BC_N == "pressure")
        m_BC_PN = BC3D::PRESSURE;
    else if (BC_N == "lsfree")
        m_BC_PN = BC3D::LSFREE;
    // not supported
    else
        m_BC_PN = BC3D::CUSTOM;

    if (BC_B == "periodic")
        m_BC_PB = BC3D::PERIODIC;
    else if (BC_B == "neumann" || BC_B == "wall")
        m_BC_PB = BC3D::NEUMANN;
    else if (BC_B == "dirichlet")
        m_BC_PB = BC3D::DIRICHLET;
    else if (BC_B == "inlet")
        m_BC_PB = BC3D::INLET;
    else if (BC_B == "outlet")
        m_BC_PB = BC3D::OUTLET;
    else if (BC_B == "pressure")
        m_BC_PB = BC3D::PRESSURE;
    else if (BC_B == "lsfree")
        m_BC_PB = BC3D::LSFREE;
    // not supported
    else
        m_BC_PB = BC3D::CUSTOM;

    if (BC_T == "periodic")
        m_BC_PT = BC3D::PERIODIC;
    else if (BC_T == "neumann" || BC_T == "wall")
        m_BC_PT = BC3D::NEUMANN;
    else if (BC_T == "dirichlet")
        m_BC_PT = BC3D::DIRICHLET;
    else if (BC_T == "inlet")
        m_BC_PT = BC3D::INLET;
    else if (BC_T == "outlet")
        m_BC_PT = BC3D::OUTLET;
    else if (BC_T == "pressure")
        m_BC_PT = BC3D::PRESSURE;
    else if (BC_T == "lsfree")
        m_BC_PT = BC3D::LSFREE;
    // not supported
    else
        m_BC_PT = BC3D::CUSTOM;

    return 0;
}

int BoundaryCondition3D::ApplyBC_P_3D(std::vector<double>& arr) {
    BC_PW(arr);
    BC_PE(arr);
    BC_PS(arr);
    BC_PN(arr);
    BC_PB(arr);
    BC_PT(arr);

    return 0;
}

void BoundaryCondition3D::BC_UW(std::vector<double>& arr) {
    if (m_BC_UW == BC3D::PERIODIC)
        BC_PeriodicUW(arr);
    else if (m_BC_UW == BC3D::NEUMANN)
        BC_NeumannUW(arr);
    else if (m_BC_UW == BC3D::DIRICHLET)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC3D::INLET)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC3D::OUTLET)
        BC_NeumannUW(arr);
}

void BoundaryCondition3D::BC_UE(std::vector<double>& arr) {
    if (m_BC_UE == BC3D::PERIODIC)
        BC_PeriodicUE(arr);
    else if (m_BC_UE == BC3D::NEUMANN)
        BC_NeumannUE(arr);
    else if (m_BC_UE == BC3D::DIRICHLET)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC3D::INLET)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC3D::OUTLET)
        BC_NeumannUE(arr);
}

void BoundaryCondition3D::BC_US(std::vector<double>& arr) {
    if (m_BC_US == BC3D::PERIODIC)
        BC_PeriodicUS(arr);
    else if (m_BC_US == BC3D::NEUMANN)
        BC_NeumannUS(arr);
    else if (m_BC_US == BC3D::DIRICHLET)
        BC_DirichletUS(arr);
    else if (m_BC_US == BC3D::INLET)
        BC_DirichletUS(arr);
    else if (m_BC_US == BC3D::OUTLET)
        BC_NeumannUS(arr);
}

void BoundaryCondition3D::BC_UN(std::vector<double>& arr) {
    if (m_BC_UN == BC3D::PERIODIC)
        BC_PeriodicUN(arr);
    else if (m_BC_UN == BC3D::NEUMANN)
        BC_NeumannUN(arr);
    else if (m_BC_UN == BC3D::DIRICHLET)
        BC_DirichletUN(arr);
    else if (m_BC_UN == BC3D::INLET)
        BC_DirichletUN(arr);
    else if (m_BC_UN == BC3D::OUTLET)
        BC_NeumannUN(arr);
}

void BoundaryCondition3D::BC_UB(std::vector<double>& arr) {
    if (m_BC_UB == BC3D::PERIODIC)
        BC_PeriodicUB(arr);
    else if (m_BC_UB == BC3D::NEUMANN)
        BC_NeumannUB(arr);
    else if (m_BC_UB == BC3D::DIRICHLET)
        BC_DirichletUB(arr);
    else if (m_BC_UB == BC3D::INLET)
        BC_DirichletUB(arr);
    else if (m_BC_UB == BC3D::OUTLET)
        BC_NeumannUB(arr);
}

void BoundaryCondition3D::BC_UT(std::vector<double>& arr) {
    if (m_BC_UT == BC3D::PERIODIC)
        BC_PeriodicUT(arr);
    else if (m_BC_UT == BC3D::NEUMANN)
        BC_NeumannUT(arr);
    else if (m_BC_UT == BC3D::DIRICHLET)
        BC_DirichletUT(arr);
    else if (m_BC_UT == BC3D::INLET)
        BC_DirichletUT(arr);
    else if (m_BC_UT == BC3D::OUTLET)
        BC_NeumannUT(arr);
}

void BoundaryCondition3D::BC_VW(std::vector<double>& arr) {
    if (m_BC_VW == BC3D::PERIODIC)
        BC_PeriodicVW(arr);
    else if (m_BC_VW == BC3D::NEUMANN)
        BC_NeumannVW(arr);
    else if (m_BC_VW == BC3D::DIRICHLET)
        BC_DirichletVW(arr);
    else if (m_BC_VW == BC3D::INLET)
        BC_DirichletVW(arr);
    else if (m_BC_VW == BC3D::OUTLET)
        BC_NeumannVW(arr);
}

void BoundaryCondition3D::BC_VE(std::vector<double>& arr) {
    if (m_BC_VE == BC3D::PERIODIC)
        BC_PeriodicVE(arr);
    else if (m_BC_VE == BC3D::NEUMANN)
        BC_NeumannVE(arr);
    else if (m_BC_VE == BC3D::DIRICHLET)
        BC_DirichletVE(arr);
    else if (m_BC_VE == BC3D::INLET)
        BC_DirichletVE(arr);
    else if (m_BC_VE == BC3D::OUTLET)
        BC_NeumannVE(arr);
}

void BoundaryCondition3D::BC_VS(std::vector<double>& arr) {
    if (m_BC_VS == BC3D::PERIODIC)
        BC_PeriodicVS(arr);
    else if (m_BC_VS == BC3D::NEUMANN)
        BC_NeumannVS(arr);
    else if (m_BC_VS == BC3D::DIRICHLET)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC3D::INLET)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC3D::OUTLET)
        BC_NeumannVS(arr);
}

void BoundaryCondition3D::BC_VN(std::vector<double>& arr) {
    if (m_BC_VN == BC3D::PERIODIC)
        BC_PeriodicVN(arr);
    else if (m_BC_VN == BC3D::NEUMANN)
        BC_NeumannVN(arr);
    else if (m_BC_VN == BC3D::DIRICHLET)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC3D::INLET)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC3D::OUTLET)
        BC_NeumannVN(arr);
}

void BoundaryCondition3D::BC_VB(std::vector<double>& arr) {
    if (m_BC_VB == BC3D::PERIODIC)
        BC_PeriodicVB(arr);
    else if (m_BC_VB == BC3D::NEUMANN)
        BC_NeumannVB(arr);
    else if (m_BC_VB == BC3D::DIRICHLET)
        BC_DirichletVB(arr);
    else if (m_BC_VB == BC3D::INLET)
        BC_DirichletVB(arr);
    else if (m_BC_VB == BC3D::OUTLET)
        BC_NeumannVB(arr);
}

void BoundaryCondition3D::BC_VT(std::vector<double>& arr) {
    if (m_BC_VT == BC3D::PERIODIC)
        BC_PeriodicVT(arr);
    else if (m_BC_VT == BC3D::NEUMANN)
        BC_NeumannVT(arr);
    else if (m_BC_VT == BC3D::DIRICHLET)
        BC_DirichletVT(arr);
    else if (m_BC_VT == BC3D::INLET)
        BC_DirichletVT(arr);
    else if (m_BC_VT == BC3D::OUTLET)
        BC_NeumannVT(arr);
}

void BoundaryCondition3D::BC_WW(std::vector<double>& arr) {
    if (m_BC_WW == BC3D::PERIODIC)
        BC_PeriodicWW(arr);
    else if (m_BC_WW == BC3D::NEUMANN)
        BC_NeumannWW(arr);
    else if (m_BC_WW == BC3D::DIRICHLET)
        BC_DirichletWW(arr);
    else if (m_BC_WW == BC3D::INLET)
        BC_DirichletWW(arr);
    else if (m_BC_WW == BC3D::OUTLET)
        BC_NeumannWW(arr);
}

void BoundaryCondition3D::BC_WE(std::vector<double>& arr) {
    if (m_BC_WE == BC3D::PERIODIC)
        BC_PeriodicWE(arr);
    else if (m_BC_WE == BC3D::NEUMANN)
        BC_NeumannWE(arr);
    else if (m_BC_WE == BC3D::DIRICHLET)
        BC_DirichletWE(arr);
    else if (m_BC_WE == BC3D::INLET)
        BC_DirichletWE(arr);
    else if (m_BC_WE == BC3D::OUTLET)
        BC_NeumannWE(arr);
}

void BoundaryCondition3D::BC_WS(std::vector<double>& arr) {
    if (m_BC_WS == BC3D::PERIODIC)
        BC_PeriodicWS(arr);
    else if (m_BC_WS == BC3D::NEUMANN)
        BC_NeumannWS(arr);
    else if (m_BC_WS == BC3D::DIRICHLET)
        BC_DirichletWS(arr);
    else if (m_BC_WS == BC3D::INLET)
        BC_DirichletWS(arr);
    else if (m_BC_WS == BC3D::OUTLET)
        BC_NeumannWS(arr);
}

void BoundaryCondition3D::BC_WN(std::vector<double>& arr) {
    if (m_BC_WN == BC3D::PERIODIC)
        BC_PeriodicWN(arr);
    else if (m_BC_WN == BC3D::NEUMANN)
        BC_NeumannWN(arr);
    else if (m_BC_WN == BC3D::DIRICHLET)
        BC_DirichletWN(arr);
    else if (m_BC_WN == BC3D::INLET)
        BC_DirichletWN(arr);
    else if (m_BC_WN == BC3D::OUTLET)
        BC_NeumannWN(arr);
}

void BoundaryCondition3D::BC_WB(std::vector<double>& arr) {
    if (m_BC_WB == BC3D::PERIODIC)
        BC_PeriodicWB(arr);
    else if (m_BC_WB == BC3D::NEUMANN)
        BC_NeumannWB(arr);
    else if (m_BC_WB == BC3D::DIRICHLET)
        BC_DirichletWB(arr);
    else if (m_BC_WB == BC3D::INLET)
        BC_DirichletWB(arr);
    else if (m_BC_WB == BC3D::OUTLET)
        BC_NeumannWB(arr);
}

void BoundaryCondition3D::BC_WT(std::vector<double>& arr) {
    if (m_BC_WT == BC3D::PERIODIC)
        BC_PeriodicWT(arr);
    else if (m_BC_WT == BC3D::NEUMANN)
        BC_NeumannWT(arr);
    else if (m_BC_WT == BC3D::DIRICHLET)
        BC_DirichletWT(arr);
    else if (m_BC_WT == BC3D::INLET)
        BC_DirichletWT(arr);
    else if (m_BC_WT == BC3D::OUTLET)
        BC_NeumannWT(arr);
}

void BoundaryCondition3D::BC_PW(std::vector<double>& arr) {
    if (m_BC_PW == BC3D::PERIODIC)
        BC_PeriodicPW(arr);
    else if (m_BC_PW == BC3D::NEUMANN)
        BC_NeumannPW(arr);
    else if (m_BC_PW == BC3D::DIRICHLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC3D::INLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC3D::OUTLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC3D::PRESSURE)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC3D::LSFREE)
        BC_LSFreeBoundaryPW(arr);
}

void BoundaryCondition3D::BC_PE(std::vector<double>& arr) {
    if (m_BC_PE == BC3D::PERIODIC)
        BC_PeriodicPE(arr);
    else if (m_BC_PE == BC3D::NEUMANN)
        BC_NeumannPE(arr);
    else if (m_BC_PE == BC3D::DIRICHLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC3D::INLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC3D::OUTLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC3D::PRESSURE)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC3D::LSFREE)
        BC_LSFreeBoundaryPE(arr);
}

void BoundaryCondition3D::BC_PS(std::vector<double>& arr) {
    if (m_BC_PS == BC3D::PERIODIC)
        BC_PeriodicPS(arr);
    else if (m_BC_PS == BC3D::NEUMANN)
        BC_NeumannPS(arr);
    else if (m_BC_PS == BC3D::DIRICHLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC3D::INLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC3D::OUTLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC3D::PRESSURE)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC3D::LSFREE)
        BC_LSFreeBoundaryPS(arr);
}

void BoundaryCondition3D::BC_PN(std::vector<double>& arr) {
    if (m_BC_PN == BC3D::PERIODIC)
        BC_PeriodicPN(arr);
    else if (m_BC_PN == BC3D::NEUMANN)
        BC_NeumannPN(arr);
    else if (m_BC_PN == BC3D::DIRICHLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC3D::INLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC3D::OUTLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC3D::PRESSURE)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC3D::LSFREE)
        BC_LSFreeBoundaryPN(arr);
}

void BoundaryCondition3D::BC_PB(std::vector<double>& arr) {
    if (m_BC_PB == BC3D::PERIODIC)
        BC_PeriodicPB(arr);
    else if (m_BC_PB == BC3D::NEUMANN)
        BC_NeumannPB(arr);
    else if (m_BC_PB == BC3D::DIRICHLET)
        BC_DirichletPB(arr);
    else if (m_BC_PB == BC3D::INLET)
        BC_DirichletPB(arr);
    else if (m_BC_PB == BC3D::OUTLET)
        BC_DirichletPB(arr);
    else if (m_BC_PB == BC3D::PRESSURE)
        BC_DirichletPB(arr);
    else if (m_BC_PB == BC3D::LSFREE)
        BC_LSFreeBoundaryPB(arr);
}

void BoundaryCondition3D::BC_PT(std::vector<double>& arr) {
    if (m_BC_PT == BC3D::PERIODIC)
        BC_PeriodicPT(arr);
    else if (m_BC_PT == BC3D::NEUMANN)
        BC_NeumannPT(arr);
    else if (m_BC_PT == BC3D::DIRICHLET)
        BC_DirichletPT(arr);
    else if (m_BC_PT == BC3D::INLET)
        BC_DirichletPT(arr);
    else if (m_BC_PT == BC3D::OUTLET)
        BC_DirichletPT(arr);
    else if (m_BC_PT == BC3D::PRESSURE)
        BC_DirichletPT(arr);
    else if (m_BC_PT == BC3D::LSFREE)
        BC_LSFreeBoundaryPT(arr);
}

void BoundaryCondition3D::BC_PeriodicUW(std::vector<double>& arr) {
    for (int i = 0; i <= kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNx - 1 + i, j, k)];
}

void BoundaryCondition3D::BC_PeriodicUE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + i + kNumBCGrid, j, k)] = arr[idx(i + kNumBCGrid, j, k)];
}

void BoundaryCondition3D::BC_PeriodicUS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_PeriodicUN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + j + kNumBCGrid, k)] = arr[idx(i, j + kNumBCGrid, k)];
}

void BoundaryCondition3D::BC_PeriodicUB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNz + k)];
}

void BoundaryCondition3D::BC_PeriodicUT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + k + kNumBCGrid)] = arr[idx(i, j, k + kNumBCGrid)];
}

void BoundaryCondition3D::BC_NeumannUW(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNumBCGrid, j, k)] = arr[idx(kNumBCGrid + 1, j, k)];

    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j,k )] = arr[idx(kNumBCGrid * 2 - i, j, k)];
}

void BoundaryCondition3D::BC_NeumannUE(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid, j, k)] = arr[idx(kNx + kNumBCGrid - 1, j, k)];

    for (int i = 1; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_NeumannUS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNumBCGrid * 2 - j - 1, k)];
}

void BoundaryCondition3D::BC_NeumannUN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_NeumannUB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNumBCGrid * 2 - k - 1)];
}

void BoundaryCondition3D::BC_NeumannUT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = arr[idx(i, j, kNz + k)];
}

// Prescribed bounadry condition (Dirichlet) for U grid
void BoundaryCondition3D::SetBCConstantUW(double BC_ConstantW) {
    m_BC_DirichletConstantUW = BC_ConstantW;
}

void BoundaryCondition3D::SetBCConstantUE(double BC_ConstantE) {
    m_BC_DirichletConstantUE = BC_ConstantE;
}

void BoundaryCondition3D::SetBCConstantUS(double BC_ConstantS) {
    m_BC_DirichletConstantUS = BC_ConstantS;
}

void BoundaryCondition3D::SetBCConstantUN(double BC_ConstantN) {
    m_BC_DirichletConstantUN = BC_ConstantN;
}

void BoundaryCondition3D::SetBCConstantUB(double BC_ConstantB) {
    m_BC_DirichletConstantUB = BC_ConstantB;
}

void BoundaryCondition3D::SetBCConstantUT(double BC_ConstantT) {
    m_BC_DirichletConstantUT = BC_ConstantT;
}

void BoundaryCondition3D::BC_DirichletUW(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNumBCGrid, j, k)] = m_BC_DirichletConstantUW;

    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(kNumBCGrid * 2 - i, j, k)]
            + 2.0 * m_BC_DirichletConstantUW;
}

void BoundaryCondition3D::BC_DirichletUE(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid, j, k)] = m_BC_DirichletConstantUE;
    
    for (int i = 1; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i, j, k)] = -arr[idx(kNx + i, j, k)]
            + 2.0 * m_BC_DirichletConstantUE;
}

void BoundaryCondition3D::BC_DirichletUS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, kNumBCGrid * 2 - j - 1, k)]
            + 2.0 * m_BC_DirichletConstantUS;
}

void BoundaryCondition3D::BC_DirichletUN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = -arr[idx(i, kNy + j, k)]
            + 2.0 * m_BC_DirichletConstantUN;
}

void BoundaryCondition3D::BC_DirichletUB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, j, kNumBCGrid * 2 - k - 1)]
            + 2.0 * m_BC_DirichletConstantUB;
}

void BoundaryCondition3D::BC_DirichletUT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = -arr[idx(i, j, kNz + k)]
            + 2.0 * m_BC_DirichletConstantUT;
}

void BoundaryCondition3D::BC_PeriodicVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_PeriodicVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + i + kNumBCGrid, j, k)] = arr[idx(i + kNumBCGrid, j, k)];
}

void BoundaryCondition3D::BC_PeriodicVS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j <= kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNy - 1 + j, k)];
}

void BoundaryCondition3D::BC_PeriodicVN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + j + kNumBCGrid, k)] = arr[idx(i, j + kNumBCGrid, k)];
}

void BoundaryCondition3D::BC_PeriodicVB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k <= kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNz + k)];
}

void BoundaryCondition3D::BC_PeriodicVT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + k + kNumBCGrid)] = arr[idx(i, j, k + kNumBCGrid)];
}

void BoundaryCondition3D::BC_NeumannVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNumBCGrid * 2 - i - 1, j, k)];
}

void BoundaryCondition3D::BC_NeumannVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_NeumannVS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNumBCGrid, k)] = arr[idx(i, kNumBCGrid + 1, k)];
    
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNumBCGrid * 2 - j, k)];
}

void BoundaryCondition3D::BC_NeumannVN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid, k)] = arr[idx(i, kNy + kNumBCGrid - 1, k)];

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 1; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_NeumannVB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNumBCGrid * 2 - k - 1)];
}

void BoundaryCondition3D::BC_NeumannVT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = arr[idx(i, j, kNz + k)];
}

// Prescribed bounadry condition (Dirichlet) for V grid
void BoundaryCondition3D::SetBCConstantVW(double BC_ConstantW) {
    m_BC_DirichletConstantVW = BC_ConstantW;
}

void BoundaryCondition3D::SetBCConstantVE(double BC_ConstantE) {
    m_BC_DirichletConstantVE = BC_ConstantE;
}

void BoundaryCondition3D::SetBCConstantVS(double BC_ConstantS) {
    m_BC_DirichletConstantVS = BC_ConstantS;
}

void BoundaryCondition3D::SetBCConstantVN(double BC_ConstantN) {
    m_BC_DirichletConstantVN = BC_ConstantN;
}

void BoundaryCondition3D::SetBCConstantVB(double BC_ConstantB) {
    m_BC_DirichletConstantVB = BC_ConstantB;
}

void BoundaryCondition3D::SetBCConstantVT(double BC_ConstantT) {
    m_BC_DirichletConstantVT = BC_ConstantT;
}

void BoundaryCondition3D::BC_DirichletVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(kNumBCGrid * 2 - i - 1, j, k)]
            + 2.0 * m_BC_DirichletConstantVW;
}

void BoundaryCondition3D::BC_DirichletVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = -arr[idx(kNx + i, j, k)]
            + 2.0 * m_BC_DirichletConstantVE;
}

void BoundaryCondition3D::BC_DirichletVS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNumBCGrid, k)] = m_BC_DirichletConstantVS;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, kNumBCGrid * 2 - j, k)]
            + 2.0 * m_BC_DirichletConstantVS;
}

void BoundaryCondition3D::BC_DirichletVN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid, k)] = m_BC_DirichletConstantVN;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 1; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j, k)] = -arr[idx(i, kNy + j, k)]
            + 2.0 * m_BC_DirichletConstantVN;
}

void BoundaryCondition3D::BC_DirichletVB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, j, kNumBCGrid * 2 - k - 1)]
            + 2.0 * m_BC_DirichletConstantVB;
}

void BoundaryCondition3D::BC_DirichletVT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = -arr[idx(i, j, kNz + k)]
            + 2.0 * m_BC_DirichletConstantVT;
}

void BoundaryCondition3D::BC_PeriodicWW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_PeriodicWE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + i + kNumBCGrid, j, k)] = arr[idx(i + kNumBCGrid, j, k)];
}

void BoundaryCondition3D::BC_PeriodicWS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_PeriodicWN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + j + kNumBCGrid, k)] = arr[idx(i, j + kNumBCGrid, k)];
}

void BoundaryCondition3D::BC_PeriodicWB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k <= kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNz - 1 + k)];
}

void BoundaryCondition3D::BC_PeriodicWT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + k + kNumBCGrid)] = arr[idx(i, j, k + kNumBCGrid)];
}

void BoundaryCondition3D::BC_NeumannWW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNumBCGrid * 2 - i - 1, j, k)];
}

void BoundaryCondition3D::BC_NeumannWE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_NeumannWS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNumBCGrid * 2 - j - 1, k)];
}
    
void BoundaryCondition3D::BC_NeumannWN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = arr[idx(i, kNz + j, k)];
}

void BoundaryCondition3D::BC_NeumannWB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(i, j, kNumBCGrid)] = arr[idx(i, j, kNumBCGrid + 1)];

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNumBCGrid * 2 - k)];
}

void BoundaryCondition3D::BC_NeumannWT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(i, j, kNz + kNumBCGrid)] = arr[idx(i, j, kNz + kNumBCGrid - 1)];

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int k = 1; k < kNumBCGrid; k++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k)] = arr[idx(i, j, kNz + k)];
}

// Prescribed bounadry condition (Dirichlet) for W grid
void BoundaryCondition3D::SetBCConstantWW(double BC_ConstantW) {
    m_BC_DirichletConstantWW = BC_ConstantW;
}

void BoundaryCondition3D::SetBCConstantWE(double BC_ConstantE) {
    m_BC_DirichletConstantWE = BC_ConstantE;
}

void BoundaryCondition3D::SetBCConstantWS(double BC_ConstantS) {
    m_BC_DirichletConstantWS = BC_ConstantS;
}

void BoundaryCondition3D::SetBCConstantWN(double BC_ConstantN) {
    m_BC_DirichletConstantWN = BC_ConstantN;
}

void BoundaryCondition3D::SetBCConstantWB(double BC_ConstantB) {
    m_BC_DirichletConstantWB = BC_ConstantB;
}

void BoundaryCondition3D::SetBCConstantWT(double BC_ConstantT) {
    m_BC_DirichletConstantWT = BC_ConstantT;
}

void BoundaryCondition3D::BC_DirichletWW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(kNumBCGrid * 2 - i - 1, j, k)]
            + 2.0 * m_BC_DirichletConstantVW;
}

void BoundaryCondition3D::BC_DirichletWE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = -arr[idx(kNx + i, j, k)]
            + 2.0 * m_BC_DirichletConstantVE;
}

void BoundaryCondition3D::BC_DirichletWS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, kNumBCGrid * 2 - j - 1, k)]
            + 2.0 * m_BC_DirichletConstantVS;
}

void BoundaryCondition3D::BC_DirichletWN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = -arr[idx(i, kNy + j, k)]
            + 2.0 * m_BC_DirichletConstantVN;
}

void BoundaryCondition3D::BC_DirichletWB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(i, j, kNumBCGrid)] = m_BC_DirichletConstantVB;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, j, kNumBCGrid * 2 - k)]
            + 2.0 * m_BC_DirichletConstantVB;
}

void BoundaryCondition3D::BC_DirichletWT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(i, j, kNy + kNumBCGrid)] = m_BC_DirichletConstantVT;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 1; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k)] = -arr[idx(i, j, kNz + k)]
            + 2.0 * m_BC_DirichletConstantVT;
}

void BoundaryCondition3D::BC_PeriodicPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_PeriodicPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + i, j, k)] = arr[idx(i + kNumBCGrid, j, k)];
}

void BoundaryCondition3D::BC_PeriodicPS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_PeriodicPN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + j + kNumBCGrid, k)] = arr[idx(i, j + kNumBCGrid, k)];
}

void BoundaryCondition3D::BC_PeriodicPB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k <= kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNz + k)];
}

void BoundaryCondition3D::BC_PeriodicPT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + k + kNumBCGrid)] = arr[idx(i, j, k + kNumBCGrid)];
}

void BoundaryCondition3D::BC_NeumannPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(kNumBCGrid * 2 - i - 1, j, k)];
}

void BoundaryCondition3D::BC_NeumannPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = arr[idx(kNx + i, j, k)];
}

void BoundaryCondition3D::BC_NeumannPS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, kNumBCGrid * 2 - j - 1, k)];
}

void BoundaryCondition3D::BC_NeumannPN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = arr[idx(i, kNy + j, k)];
}

void BoundaryCondition3D::BC_NeumannPB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j, kNumBCGrid * 2 - k - 1)];
}

void BoundaryCondition3D::BC_NeumannPT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = arr[idx(i, j, kNz + k)];
}

// Prescribed bounadry condition (Dirichlet) for P grid (center)
void BoundaryCondition3D::SetBCConstantPW(double BC_ConstantW) {
    m_BC_DirichletConstantPW = BC_ConstantW;
}

void BoundaryCondition3D::SetBCConstantPE(double BC_ConstantE) {
    m_BC_DirichletConstantPE = BC_ConstantE;
}

void BoundaryCondition3D::SetBCConstantPS(double BC_ConstantS) {
    m_BC_DirichletConstantPS = BC_ConstantS;
}

void BoundaryCondition3D::SetBCConstantPN(double BC_ConstantN) {
    m_BC_DirichletConstantPN = BC_ConstantN;
}

void BoundaryCondition3D::SetBCConstantPB(double BC_ConstantB) {
    m_BC_DirichletConstantPB = BC_ConstantB;
}

void BoundaryCondition3D::SetBCConstantPT(double BC_ConstantT) {
    m_BC_DirichletConstantPT = BC_ConstantT;
}

void BoundaryCondition3D::SetAmbientPressure(double ambientPressure) {
    m_AmbientPressure = ambientPressure;

    // set ambient pressure as dirichlet condition constant
    if (m_BC_PW == BC3D::INLET || m_BC_PW == BC3D::OUTLET || m_BC_PW == BC3D::PRESSURE)
        SetBCConstantPW(ambientPressure);

    if (m_BC_PE == BC3D::INLET || m_BC_PE == BC3D::OUTLET || m_BC_PE == BC3D::PRESSURE)
        SetBCConstantPE(ambientPressure);

    if (m_BC_PS == BC3D::INLET || m_BC_PS == BC3D::OUTLET || m_BC_PS == BC3D::PRESSURE)
        SetBCConstantPS(ambientPressure);

    if (m_BC_PN == BC3D::INLET || m_BC_PN == BC3D::OUTLET || m_BC_PN == BC3D::PRESSURE)
        SetBCConstantPN(ambientPressure);

    if (m_BC_PB == BC3D::INLET || m_BC_PB == BC3D::OUTLET || m_BC_PB == BC3D::PRESSURE)
        SetBCConstantPS(ambientPressure);

    if (m_BC_PT == BC3D::INLET || m_BC_PT == BC3D::OUTLET || m_BC_PT == BC3D::PRESSURE)
        SetBCConstantPN(ambientPressure);
}

void BoundaryCondition3D::BC_DirichletPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(kNumBCGrid * 2 - i - 1, j, k)]
            + 2.0 * m_BC_DirichletConstantPW;
}

void BoundaryCondition3D::BC_DirichletPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j, k)] = -arr[idx(kNx + i, j, k)]
            + 2.0 * m_BC_DirichletConstantPE;
}

void BoundaryCondition3D::BC_DirichletPS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, kNumBCGrid * 2 - j - 1, k)]
            + 2.0 * m_BC_DirichletConstantPS;
}

void BoundaryCondition3D::BC_DirichletPN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1, k)] = -arr[idx(i, kNy + j, k)]
            + 2.0 * m_BC_DirichletConstantPN;
}

void BoundaryCondition3D::BC_DirichletPB(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, k)] = -arr[idx(i, j, kNumBCGrid * 2 - k - 1)]
            + 2.0 * m_BC_DirichletConstantPB;
}

void BoundaryCondition3D::BC_DirichletPT(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid * 2 - k - 1)] = -arr[idx(i, j, kNz + k)]
            + 2.0 * m_BC_DirichletConstantPT;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPW(std::vector<double>& arr) {
    if (kDx < 0.0)
        std::cout << "BC Initialization Error : kDx not set" << std::endl;

    for (int i = kNumBCGrid - 1; i >= 0; i--)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i + 1, j, k)] + kDx * (arr[idx(i + 1, j, k)] - arr[idx(i + 2, j, k)]) / kDx;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPE(std::vector<double>& arr) {
    if (kDx < 0.0)
        std::cout << "BC Initialization Error : kDx not set" << std::endl;

    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
        arr[idx(kNx + kNumBCGrid + i, j, k)] = arr[idx(kNx + kNumBCGrid + i - 1, j, k)]
        + kDx * (arr[idx(kNx + kNumBCGrid + i - 1, j, k)] - arr[idx(kNx + kNumBCGrid + i - 2, j, k)]) / kDx;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPS(std::vector<double>& arr) {
    if (kDy < 0.0)
        std::cout << "BC Initialization Error : kDy not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid - 1; j >= 0; j--)
    for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
        arr[idx(i, j, k)] = arr[idx(i, j + 1, k)] + kDy * (arr[idx(i, j + 1, k)] - arr[idx(i, j + 2, k)]) / kDy;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPN(std::vector<double>& arr) {
    if (kDy < 0.0)
        std::cout << "BC Initialization Error : kDy not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
    for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
        arr[idx(i, kNy + kNumBCGrid + j, k)] = arr[idx(i, kNy + kNumBCGrid + j - 1, k)]
        + kDy * (arr[idx(i, kNy + kNumBCGrid + j - 1, k)] - arr[idx(i, kNy + kNumBCGrid + j - 2, k)]) / kDy;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPB(std::vector<double>& arr) {
    if (kDz < 0.0)
        std::cout << "BC Initialization Error : kDz not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int k = kNumBCGrid - 1; k >= 0; k--)
        arr[idx(i, j, k)] = arr[idx(i, j, k + 1)] + kDz * (arr[idx(i, j, k + 1)] - arr[idx(i, j, k + 2)]) / kDz;
}

void BoundaryCondition3D::BC_LSFreeBoundaryPT(std::vector<double>& arr) {
    if (kDz < 0.0)
        std::cout << "BC Initialization Error : kDz not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int k = 0; k < kNumBCGrid; k++)
        arr[idx(i, j, kNz + kNumBCGrid + k)] = arr[idx(i, j, kNz + kNumBCGrid + k - 1)]
        + kDz * (arr[idx(i, j, kNz + kNumBCGrid + k - 1)] - arr[idx(i, j, kNz + kNumBCGrid + k - 2)]) / kDz;
}