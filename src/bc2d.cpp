#include "bc2d.h"

BoundaryCondition2D::BoundaryCondition2D(int nx, int ny, int num_bc_grid)
    : kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid), kDx(-1.0), kDy(-1.0) {
}

BoundaryCondition2D::BoundaryCondition2D(int nx, int ny, double dx, double dy, int num_bc_grid)
    : kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid), kDx(dx), kDy(dy) {
}


inline int BoundaryCondition2D::idx(int i, int j) {
    return (j + (kNy + 2 * kNumBCGrid) * (i));
}

int BoundaryCondition2D::SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    
    if (BC_W == "pressure" || BC_E == "pressure" || BC_S == "pressure" || BC_N == "pressure") {
        std::cout << "Boundary Condition Error : constant pressure BC is not allowed at velocity" << std::endl;
        exit(1);
    }

    if (BC_W == "periodic")
        m_BC_UW = BC2D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_UW = BC2D::NEUMANN;
    else if (BC_W == "axisym") {
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_UW = BC2D::AXISYM;
        SetBCConstantUW(0.0);
    }
    else if (BC_W == "wall") {
        m_BC_UW = BC2D::WALL;
        SetBCConstantUW(0.0);
    }
    else if (BC_W == "dirichlet")
        m_BC_UW = BC2D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_UW = BC2D::INLET;
    else if (BC_W == "outlet")
        m_BC_UW = BC2D::OUTLET;
    // not supported
    else
        m_BC_UW = BC2D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_UE = BC2D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_UE = BC2D::NEUMANN;
    else if (BC_E == "axisym") {
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_UE = BC2D::AXISYM;
        SetBCConstantUE(0.0);
    }
    else if (BC_E == "wall") {
        m_BC_UE = BC2D::WALL;
        SetBCConstantUE(0.0);
    }
    else if (BC_E == "dirichlet")
        m_BC_UE = BC2D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_UE = BC2D::INLET;
    else if (BC_E == "outlet")
        m_BC_UE = BC2D::OUTLET;
    // not supported
    else
        m_BC_UE = BC2D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_US = BC2D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_US = BC2D::NEUMANN;
    else if (BC_S == "axisym")
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_US = BC2D::AXISYM;
    else if (BC_S == "wall") {
        m_BC_US = BC2D::WALL;
        SetBCConstantUS(0.0);
    }
    else if (BC_S == "dirichlet")
        m_BC_US = BC2D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_US = BC2D::INLET;
    else if (BC_S == "outlet")
        m_BC_US = BC2D::OUTLET;
    // not supported
    else
        m_BC_US = BC2D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_UN = BC2D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_UN = BC2D::NEUMANN;
    else if (BC_N == "axisym")
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_UN = BC2D::AXISYM;
    else if (BC_N == "wall") {
        m_BC_UN = BC2D::WALL;
        SetBCConstantUN(0.0);
    }
    else if (BC_N == "dirichlet")
        m_BC_UN = BC2D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_UN = BC2D::INLET;
    else if (BC_N == "outlet")
        m_BC_UN = BC2D::OUTLET;
    // not supported
    else
        m_BC_UN = BC2D::CUSTOM;

    return 0;
}

int BoundaryCondition2D::ApplyBC_U_2D(std::vector<double>& arr) {
    BC_UW(arr);
    BC_UE(arr);
    BC_US(arr);
    BC_UN(arr);
    
    return 0;
}

int BoundaryCondition2D::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
    
    if (BC_W == "pressure" || BC_E == "pressure" || BC_S == "pressure" || BC_N == "pressure") {
        std::cout << "Boundary Condition Error : constant pressure BC is not allowed at velocity" << std::endl;
        exit(1);
    }

    if (BC_W == "periodic")
        m_BC_VW = BC2D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_VW = BC2D::NEUMANN;
    else if (BC_W == "axisym")
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_VW = BC2D::AXISYM;
    else if (BC_W == "wall") {
        m_BC_VW = BC2D::WALL;
        SetBCConstantVW(0.0);
    }
    else if (BC_W == "dirichlet")
        m_BC_VW = BC2D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_VW = BC2D::INLET;
    else if (BC_W == "outlet")
        m_BC_VW = BC2D::OUTLET;
    // not supported
    else
        m_BC_VW = BC2D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_VE = BC2D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_VE = BC2D::NEUMANN;
    else if (BC_E == "axisym")
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_VE = BC2D::AXISYM;
    else if (BC_E == "wall") {
        m_BC_VE = BC2D::WALL;
        SetBCConstantVE(0.0);
    }
    else if (BC_E == "dirichlet")
        m_BC_VE = BC2D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_VE = BC2D::INLET;
    else if (BC_E == "outlet")
        m_BC_VE = BC2D::OUTLET;
    // not supported
    else
        m_BC_VE = BC2D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_VS = BC2D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_VS = BC2D::NEUMANN;
    else if (BC_S == "axisym") {
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_VS = BC2D::AXISYM;
        SetBCConstantVS(0.0);
    }
    else if (BC_S == "wall") {
        m_BC_VS = BC2D::WALL;
        SetBCConstantVS(0.0);
    }
    else if (BC_S == "dirichlet")
        m_BC_VS = BC2D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_VS = BC2D::INLET;
    else if (BC_S == "outlet")
        m_BC_VS = BC2D::OUTLET;
    // not supported
    else
        m_BC_VS = BC2D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_VN = BC2D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_VN = BC2D::NEUMANN;
    else if (BC_N == "axisym") {
        // zero normal velocity
        // other variables has zero flux = neumann condition
        m_BC_VN = BC2D::AXISYM;
        SetBCConstantVN(0.0);
    }
    else if (BC_N == "wall") {
        m_BC_VN = BC2D::WALL;
        SetBCConstantVN(0.0);
    }
    else if (BC_N == "dirichlet")
        m_BC_VN = BC2D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_VN = BC2D::INLET;
    else if (BC_N == "outlet")
        m_BC_VN = BC2D::OUTLET;
    // not supported
    else
        m_BC_VN = BC2D::CUSTOM;

    return 0;
}

int BoundaryCondition2D::ApplyBC_V_2D(std::vector<double>& arr) {
    BC_VW(arr);
    BC_VE(arr);
    BC_VS(arr);
    BC_VN(arr);

    return 0;
}

int BoundaryCondition2D::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {

    std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
    std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
    std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
    std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);

    if (BC_W == "periodic")
        m_BC_PW = BC2D::PERIODIC;
    else if (BC_W == "neumann")
        m_BC_PW = BC2D::NEUMANN;
    else if (BC_W == "axisym")
        m_BC_PW = BC2D::AXISYM;
    else if (BC_W == "wall")
        m_BC_PW = BC2D::WALL;
    else if (BC_W == "dirichlet")
        m_BC_PW = BC2D::DIRICHLET;
    else if (BC_W == "inlet")
        m_BC_PW = BC2D::INLET;
    else if (BC_W == "outlet")
        m_BC_PW = BC2D::OUTLET;
    else if (BC_W == "pressure")
        m_BC_PW = BC2D::PRESSURE;
    else if (BC_W == "lsfree")
        m_BC_PW = BC2D::LSFREE;
    // not supported
    else
        m_BC_PW = BC2D::CUSTOM;

    if (BC_E == "periodic")
        m_BC_PE = BC2D::PERIODIC;
    else if (BC_E == "neumann")
        m_BC_PE = BC2D::NEUMANN;
    else if (BC_E == "axisym")
        m_BC_PE = BC2D::AXISYM;
    else if (BC_E == "wall")
        m_BC_PE = BC2D::WALL;
    else if (BC_E == "dirichlet")
        m_BC_PE = BC2D::DIRICHLET;
    else if (BC_E == "inlet")
        m_BC_PE = BC2D::INLET;
    else if (BC_E == "outlet")
        m_BC_PE = BC2D::OUTLET;
    else if (BC_E == "pressure")
        m_BC_PE = BC2D::PRESSURE;
    else if (BC_E == "lsfree")
        m_BC_PE = BC2D::LSFREE;
    // not supported
    else
        m_BC_PE = BC2D::CUSTOM;

    if (BC_S == "periodic")
        m_BC_PS = BC2D::PERIODIC;
    else if (BC_S == "neumann")
        m_BC_PS = BC2D::NEUMANN;
    else if (BC_S == "axisym")
        m_BC_PS = BC2D::AXISYM;
    else if (BC_S == "wall")
        m_BC_PS = BC2D::WALL;
    else if (BC_S == "dirichlet")
        m_BC_PS = BC2D::DIRICHLET;
    else if (BC_S == "inlet")
        m_BC_PS = BC2D::INLET;
    else if (BC_S == "outlet")
        m_BC_PS = BC2D::OUTLET;
    else if (BC_S == "pressure")
        m_BC_PS = BC2D::PRESSURE;
    else if (BC_S == "lsfree")
        m_BC_PS = BC2D::LSFREE;
    // not supported
    else
        m_BC_PS = BC2D::CUSTOM;

    if (BC_N == "periodic")
        m_BC_PN = BC2D::PERIODIC;
    else if (BC_N == "neumann")
        m_BC_PN = BC2D::NEUMANN;
    else if (BC_N == "axisym")
        m_BC_PN = BC2D::AXISYM;
    else if (BC_N == "wall")
        m_BC_PN = BC2D::WALL;
    else if (BC_N == "dirichlet")
        m_BC_PN = BC2D::DIRICHLET;
    else if (BC_N == "inlet")
        m_BC_PN = BC2D::INLET;
    else if (BC_N == "outlet")
        m_BC_PN = BC2D::OUTLET;
    else if (BC_N == "pressure")
        m_BC_PN = BC2D::PRESSURE;
    else if (BC_N == "lsfree")
        m_BC_PN = BC2D::LSFREE;
    // not supported
    else
        m_BC_PN = BC2D::CUSTOM;

    return 0;
}

int BoundaryCondition2D::ApplyBC_P_2D(std::vector<double>& arr) {
    BC_PW(arr);
    BC_PE(arr);
    BC_PS(arr);
    BC_PN(arr);
    
    return 0;
}

void BoundaryCondition2D::BC_UW(std::vector<double>& arr) {
    if (m_BC_UW == BC2D::PERIODIC)
        BC_PeriodicUW(arr);
    else if (m_BC_UW == BC2D::NEUMANN)
        BC_NeumannUW(arr);
    else if (m_BC_UW == BC2D::AXISYM)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC2D::WALL)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC2D::DIRICHLET)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC2D::INLET)
        BC_DirichletUW(arr);
    else if (m_BC_UW == BC2D::OUTLET)
        BC_NeumannUW(arr);
}

void BoundaryCondition2D::BC_UE(std::vector<double>& arr) {
    if (m_BC_UE == BC2D::PERIODIC)
        BC_PeriodicUE(arr);
    else if (m_BC_UE == BC2D::NEUMANN)
        BC_NeumannUE(arr);
    else if (m_BC_UE == BC2D::AXISYM)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC2D::WALL)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC2D::DIRICHLET)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC2D::INLET)
        BC_DirichletUE(arr);
    else if (m_BC_UE == BC2D::OUTLET)
        BC_NeumannUE(arr);
}

void BoundaryCondition2D::BC_US(std::vector<double>& arr) {
    if (m_BC_US == BC2D::PERIODIC)
        BC_PeriodicUS(arr);
    else if (m_BC_US == BC2D::NEUMANN)
        BC_NeumannUS(arr);
    else if (m_BC_US == BC2D::AXISYM)
        BC_NeumannUS(arr);
    else if (m_BC_US == BC2D::WALL)
        BC_DirichletUS(arr);
    else if (m_BC_US == BC2D::DIRICHLET)
        BC_DirichletUS(arr);
    else if (m_BC_US == BC2D::INLET)
        BC_DirichletUS(arr);
    else if (m_BC_US == BC2D::OUTLET)
        BC_NeumannUS(arr);
}

void BoundaryCondition2D::BC_UN(std::vector<double>& arr) {
    if (m_BC_UN == BC2D::PERIODIC)
        BC_PeriodicUN(arr);
    else if (m_BC_UN == BC2D::NEUMANN)
        BC_NeumannUN(arr);
    else if (m_BC_UN == BC2D::AXISYM)
        BC_NeumannUN(arr);
    else if (m_BC_UN == BC2D::WALL)
        BC_DirichletUN(arr);
    else if (m_BC_UN == BC2D::DIRICHLET)
        BC_DirichletUN(arr);
    else if (m_BC_UN == BC2D::INLET)
        BC_DirichletUN(arr);
    else if (m_BC_UN == BC2D::OUTLET)
        BC_NeumannUN(arr);
}

void BoundaryCondition2D::BC_VW(std::vector<double>& arr) {
    if (m_BC_VW == BC2D::PERIODIC)
        BC_PeriodicVW(arr);
    else if (m_BC_VW == BC2D::NEUMANN)
        BC_NeumannVW(arr);
    else if (m_BC_VW == BC2D::AXISYM)
        BC_NeumannVW(arr);
    else if (m_BC_VW == BC2D::WALL)
        BC_DirichletVW(arr);
    else if (m_BC_VW == BC2D::DIRICHLET)
        BC_DirichletVW(arr);
    else if (m_BC_VW == BC2D::INLET)
        BC_DirichletVW(arr);
    else if (m_BC_VW == BC2D::OUTLET)
        BC_NeumannVW(arr);
}

void BoundaryCondition2D::BC_VE(std::vector<double>& arr) {
    if (m_BC_VE == BC2D::PERIODIC)
        BC_PeriodicVE(arr);
    else if (m_BC_VE == BC2D::NEUMANN)
        BC_NeumannVE(arr);
    else if (m_BC_VE == BC2D::AXISYM)
        BC_NeumannVE(arr);
    else if (m_BC_VE == BC2D::WALL)
        BC_DirichletVE(arr);
    else if (m_BC_VE == BC2D::DIRICHLET)
        BC_DirichletVE(arr);
    else if (m_BC_VE == BC2D::INLET)
        BC_DirichletVE(arr);
    else if (m_BC_VE == BC2D::OUTLET)
        BC_NeumannVE(arr);
}

void BoundaryCondition2D::BC_VS(std::vector<double>& arr) {
    if (m_BC_VS == BC2D::PERIODIC)
        BC_PeriodicVS(arr);
    else if (m_BC_VS == BC2D::NEUMANN)
        BC_NeumannVS(arr);
    else if (m_BC_VS == BC2D::AXISYM)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC2D::WALL)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC2D::DIRICHLET)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC2D::INLET)
        BC_DirichletVS(arr);
    else if (m_BC_VS == BC2D::OUTLET)
        BC_NeumannVS(arr);
}

void BoundaryCondition2D::BC_VN(std::vector<double>& arr) {
    if (m_BC_VN == BC2D::PERIODIC)
        BC_PeriodicVN(arr);
    else if (m_BC_VN == BC2D::NEUMANN)
        BC_NeumannVN(arr);
    else if (m_BC_VN == BC2D::AXISYM)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC2D::WALL)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC2D::DIRICHLET)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC2D::INLET)
        BC_DirichletVN(arr);
    else if (m_BC_VN == BC2D::OUTLET)
        BC_NeumannVN(arr);
}

void BoundaryCondition2D::BC_PW(std::vector<double>& arr) {
    if (m_BC_PW == BC2D::PERIODIC)
        BC_PeriodicPW(arr);
    else if (m_BC_PW == BC2D::NEUMANN)
        BC_NeumannPW(arr);
    else if (m_BC_PW == BC2D::AXISYM)
        BC_NeumannPW(arr);
    else if (m_BC_PW == BC2D::WALL)
        BC_NeumannPW(arr);
    else if (m_BC_PW == BC2D::DIRICHLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC2D::INLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC2D::OUTLET)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC2D::PRESSURE)
        BC_DirichletPW(arr);
    else if (m_BC_PW == BC2D::LSFREE)
        BC_LSFreeBoundaryPW(arr);
}

void BoundaryCondition2D::BC_PE(std::vector<double>& arr) {
    if (m_BC_PE == BC2D::PERIODIC)
        BC_PeriodicPE(arr);
    else if (m_BC_PE == BC2D::NEUMANN)
        BC_NeumannPE(arr);
    else if (m_BC_PE == BC2D::AXISYM)
        BC_NeumannPE(arr);
    else if (m_BC_PE == BC2D::WALL)
        BC_NeumannPE(arr);
    else if (m_BC_PE == BC2D::DIRICHLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC2D::INLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC2D::OUTLET)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC2D::PRESSURE)
        BC_DirichletPE(arr);
    else if (m_BC_PE == BC2D::LSFREE)
        BC_LSFreeBoundaryPE(arr);
}

void BoundaryCondition2D::BC_PS(std::vector<double>& arr) {
    if (m_BC_PS == BC2D::PERIODIC)
        BC_PeriodicPS(arr);
    else if (m_BC_PS == BC2D::NEUMANN)
        BC_NeumannPS(arr);
    else if (m_BC_PS == BC2D::AXISYM)
        BC_NeumannPS(arr);
    else if (m_BC_PS == BC2D::WALL)
        BC_NeumannPS(arr);
    else if (m_BC_PS == BC2D::DIRICHLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC2D::INLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC2D::OUTLET)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC2D::PRESSURE)
        BC_DirichletPS(arr);
    else if (m_BC_PS == BC2D::LSFREE)
        BC_LSFreeBoundaryPS(arr);
}

void BoundaryCondition2D::BC_PN(std::vector<double>& arr) {
    if (m_BC_PN == BC2D::PERIODIC)
        BC_PeriodicPN(arr);
    else if (m_BC_PN == BC2D::NEUMANN)
        BC_NeumannPN(arr);
    else if (m_BC_PN == BC2D::AXISYM)
        BC_NeumannPN(arr);
    else if (m_BC_PN == BC2D::WALL)
        BC_NeumannPN(arr);
    else if (m_BC_PN == BC2D::DIRICHLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC2D::INLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC2D::OUTLET)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC2D::PRESSURE)
        BC_DirichletPN(arr);
    else if (m_BC_PN == BC2D::LSFREE)
        BC_LSFreeBoundaryPN(arr);
}

void BoundaryCondition2D::BC_PeriodicUW(std::vector<double>& arr) {
    for (int i = 0; i <= kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNx - 1 + i, j)];
}

void BoundaryCondition2D::BC_PeriodicUE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + i + kNumBCGrid, j)] = arr[idx(i + kNumBCGrid, j)];
}

void BoundaryCondition2D::BC_PeriodicUS(std::vector<double>& arr) {
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNy + j)];
}

void BoundaryCondition2D::BC_PeriodicUN(std::vector<double>& arr) {
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + j + kNumBCGrid)] = arr[idx(i, j + kNumBCGrid)];
}

void BoundaryCondition2D::BC_NeumannUW(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(kNumBCGrid, j)] = arr[idx(kNumBCGrid + 1, j)];

    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j <= kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNumBCGrid * 2 - i, j)];
}

void BoundaryCondition2D::BC_NeumannUE(std::vector<double>& arr) {
    for(int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid, j)] = arr[idx(kNx + kNumBCGrid - 1, j)];

    for (int i = 1; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j <= kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i, j)] = arr[idx(kNx + i, j)];
}

void BoundaryCondition2D::BC_NeumannUS(std::vector<double>& arr) {
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNumBCGrid * 2 - j - 1)];
}

void BoundaryCondition2D::BC_NeumannUN(std::vector<double>& arr) {
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1)] = arr[idx(i, kNy + j)];
}

// Prescribed bounadry condition (Dirichlet) for U grid
void BoundaryCondition2D::SetBCConstantUW(double BC_ConstantW) {
    m_BC_DirichletConstantUW = BC_ConstantW;
}

void BoundaryCondition2D::SetBCConstantUE(double BC_ConstantE) {
    m_BC_DirichletConstantUE = BC_ConstantE;
}

void BoundaryCondition2D::SetBCConstantUS(double BC_ConstantS) {
    m_BC_DirichletConstantUS = BC_ConstantS;
}

void BoundaryCondition2D::SetBCConstantUN(double BC_ConstantN) {
    m_BC_DirichletConstantUN = BC_ConstantN;
}

void BoundaryCondition2D::BC_DirichletUW(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(kNumBCGrid, j)] = m_BC_DirichletConstantUW;

    // to remove corner singularity, equal(=) is added 
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(kNumBCGrid * 2 - i, j)]
            + 2.0 * m_BC_DirichletConstantUW;
}

void BoundaryCondition2D::BC_DirichletUE(std::vector<double>& arr) {
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid, j)] = m_BC_DirichletConstantUE;
    
    for (int i = 1; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i, j)] = -arr[idx(kNx + i, j)]
            + 2.0 * m_BC_DirichletConstantUE;
}

void BoundaryCondition2D::BC_DirichletUS(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(i, kNumBCGrid * 2 - j - 1)]
            + 2.0 * m_BC_DirichletConstantUS;
}

void BoundaryCondition2D::BC_DirichletUN(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1)] = -arr[idx(i, kNy + j)]
            + 2.0 * m_BC_DirichletConstantUN;
}

void BoundaryCondition2D::BC_PeriodicVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNx + i, j)];
}

void BoundaryCondition2D::BC_PeriodicVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + i + kNumBCGrid, j)] = arr[idx(i + kNumBCGrid, j)];
}

void BoundaryCondition2D::BC_PeriodicVS(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j <= kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNy - 1 + j)];
}

void BoundaryCondition2D::BC_PeriodicVN(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + j + kNumBCGrid)] = arr[idx(i, j + kNumBCGrid)];
}

void BoundaryCondition2D::BC_NeumannVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNumBCGrid * 2 - i - 1, j)];
}

void BoundaryCondition2D::BC_NeumannVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j)] = arr[idx(kNx + i, j)];
}

void BoundaryCondition2D::BC_NeumannVS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
        arr[idx(i, kNumBCGrid)] = arr[idx(i, kNumBCGrid + 1)];
    
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNumBCGrid * 2 - j)];
}

void BoundaryCondition2D::BC_NeumannVN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
        arr[idx(i, kNy + kNumBCGrid)] = arr[idx(i, kNy + kNumBCGrid - 1)];

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 1; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j)] = arr[idx(i, kNy + j)];
}

// Prescribed bounadry condition (Dirichlet) for V grid
void BoundaryCondition2D::SetBCConstantVW(double BC_ConstantW) {
    m_BC_DirichletConstantVW = BC_ConstantW;
}

void BoundaryCondition2D::SetBCConstantVE(double BC_ConstantE) {
    m_BC_DirichletConstantVE = BC_ConstantE;
}

void BoundaryCondition2D::SetBCConstantVS(double BC_ConstantS) {
    m_BC_DirichletConstantVS = BC_ConstantS;
}

void BoundaryCondition2D::SetBCConstantVN(double BC_ConstantN) {
    m_BC_DirichletConstantVN = BC_ConstantN;
}

void BoundaryCondition2D::BC_DirichletVW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(kNumBCGrid * 2 - i - 1, j)]
            + 2.0 * m_BC_DirichletConstantVW;
}

void BoundaryCondition2D::BC_DirichletVE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j)] = -arr[idx(kNx + i, j)]
            + 2.0 * m_BC_DirichletConstantVE;
}

void BoundaryCondition2D::BC_DirichletVS(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
        arr[idx(i, kNumBCGrid)] = m_BC_DirichletConstantVS;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(i, kNumBCGrid * 2 - j)]
            + 2.0 * m_BC_DirichletConstantVS;
}

void BoundaryCondition2D::BC_DirichletVN(std::vector<double>& arr) {
    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
        arr[idx(i, kNy + kNumBCGrid)] = m_BC_DirichletConstantVN;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 1; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j)] = -arr[idx(i, kNy + j)]
            + 2.0 * m_BC_DirichletConstantVN;
}

void BoundaryCondition2D::BC_PeriodicPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNx + i, j)];
}

void BoundaryCondition2D::BC_PeriodicPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + i, j)] = arr[idx(i + kNumBCGrid, j)];
}

void BoundaryCondition2D::BC_PeriodicPS(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNy + j)];
}

void BoundaryCondition2D::BC_PeriodicPN(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + j + kNumBCGrid)] = arr[idx(i, j + kNumBCGrid)];
}

void BoundaryCondition2D::BC_NeumannPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(kNumBCGrid * 2 - i - 1, j)];
}

void BoundaryCondition2D::BC_NeumannPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j)] = arr[idx(kNx + i, j)];
}

void BoundaryCondition2D::BC_NeumannPS(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i, kNumBCGrid * 2 - j - 1)];
}

void BoundaryCondition2D::BC_NeumannPN(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1)] = arr[idx(i, kNy + j)];
}

// Prescribed bounadry condition (Dirichlet) for P grid (center)
void BoundaryCondition2D::SetBCConstantPW(double BC_ConstantW) {
    m_BC_DirichletConstantPW = BC_ConstantW;
}

void BoundaryCondition2D::SetBCConstantPE(double BC_ConstantE) {
    m_BC_DirichletConstantPE = BC_ConstantE;
}

void BoundaryCondition2D::SetBCConstantPS(double BC_ConstantS) {
    m_BC_DirichletConstantPS = BC_ConstantS;
}

void BoundaryCondition2D::SetBCConstantPN(double BC_ConstantN) {
    m_BC_DirichletConstantPN = BC_ConstantN;
}

void BoundaryCondition2D::SetAmbientPressure(double ambientPressure) {
    m_AmbientPressure = ambientPressure;

    // set ambient pressure as dirichlet condition constant
    if (m_BC_PW == BC2D::INLET || m_BC_PW == BC2D::OUTLET || m_BC_PW == BC2D::PRESSURE)
        SetBCConstantPW(ambientPressure);

    if (m_BC_PE == BC2D::INLET || m_BC_PE == BC2D::OUTLET || m_BC_PE == BC2D::PRESSURE)
        SetBCConstantPE(ambientPressure);

    if (m_BC_PS == BC2D::INLET || m_BC_PS == BC2D::OUTLET || m_BC_PS == BC2D::PRESSURE)
        SetBCConstantPS(ambientPressure);

    if (m_BC_PN == BC2D::INLET || m_BC_PN == BC2D::OUTLET || m_BC_PN == BC2D::PRESSURE)
        SetBCConstantPN(ambientPressure);
}

void BoundaryCondition2D::BC_DirichletPW(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(kNumBCGrid * 2 - i - 1, j)]
            + 2.0 * m_BC_DirichletConstantPW;
}

void BoundaryCondition2D::BC_DirichletPE(std::vector<double>& arr) {
    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid * 2 - i - 1, j)] = -arr[idx(kNx + i, j)]
            + 2.0 * m_BC_DirichletConstantPE;
}

void BoundaryCondition2D::BC_DirichletPS(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, j)] = -arr[idx(i, kNumBCGrid * 2 - j - 1)]
            + 2.0 * m_BC_DirichletConstantPS;
}

void BoundaryCondition2D::BC_DirichletPN(std::vector<double>& arr) {
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid * 2 - j - 1)] = -arr[idx(i, kNy + j)]
            + 2.0 * m_BC_DirichletConstantPN;
}

void BoundaryCondition2D::BC_LSFreeBoundaryPW(std::vector<double>& arr) {
    if (kDx < 0.0)
        std::cout << "BC Initialization Error : kDx not set" << std::endl;

    for (int i = kNumBCGrid - 1; i >= 0; i--)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(i, j)] = arr[idx(i + 1, j)]
            + kDx * (arr[idx(i + 1, j)] - arr[idx(i + 2, j)]) / kDx;
}

void BoundaryCondition2D::BC_LSFreeBoundaryPE(std::vector<double>& arr) {
    if (kDx < 0.0)
        std::cout << "BC Initialization Error : kDx not set" << std::endl;

    for (int i = 0; i < kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        arr[idx(kNx + kNumBCGrid + i, j)] = arr[idx(kNx + kNumBCGrid + i - 1, j)] 
            + kDx * (arr[idx(kNx + kNumBCGrid + i - 1, j)] - arr[idx(kNx + kNumBCGrid + i - 2, j)]) / kDx;
}

void BoundaryCondition2D::BC_LSFreeBoundaryPS(std::vector<double>& arr) {
    if (kDy < 0.0)
        std::cout << "BC Initialization Error : kDy not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid - 1; j >= 0; j--)
        arr[idx(i, j)] = arr[idx(i, j + 1)]
            + kDy * (arr[idx(i, j + 1)] - arr[idx(i, j + 2)]) / kDy;
}

void BoundaryCondition2D::BC_LSFreeBoundaryPN(std::vector<double>& arr) {
    if (kDy < 0.0)
        std::cout << "BC Initialization Error : kDy not set" << std::endl;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = 0; j < kNumBCGrid; j++)
        arr[idx(i, kNy + kNumBCGrid + j)] = arr[idx(i, kNy + kNumBCGrid + j - 1)]
            + kDy * (arr[idx(i, kNy + kNumBCGrid + j - 1)] - arr[idx(i, kNy + kNumBCGrid + j - 2)]) / kDy;
}

