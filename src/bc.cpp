#include "bc.h"

BoundaryCondition::BoundaryCondition(int nr, int nz, int num_bc_grid) 
	: kNr(nr), kNz(nz), kNx(0), kNy(0), kNumBCGrid(num_bc_grid) {
}

BoundaryCondition::BoundaryCondition(int nx, int ny, int nz, int num_bc_grid)
	: kNr(0), kNx(nx), kNy(ny), kNz(nz), kNumBCGrid(num_bc_grid) {
}

int BoundaryCondition::SetBC_U_3D(std::string BC_W, std::string BC_E,
	std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

	std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
	std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
	std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
	std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
	std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
	std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);

	if (BC_W == "periodic")
		m_BC_UW = BC::PERIODIC;
	else if (BC_W == "neumann")
		m_BC_UW = BC::NEUMANN;
	else if (BC_W == "wall" || BC_W == "dirichlet")
		m_BC_UW = BC::DIRICHLET;
	// not supported
	else
		m_BC_UW = BC::CUSTOM;

	if (BC_E == "periodic")
		m_BC_UE = BC::PERIODIC;
	else if (BC_E == "neumann")
		m_BC_UE = BC::NEUMANN;
	else if (BC_E == "wall" || BC_E == "dirichlet")
		m_BC_UE = BC::DIRICHLET;
	// not supported
	else
		m_BC_UE = BC::CUSTOM;

	if (BC_S == "periodic")
		m_BC_US = BC::PERIODIC;
	else if (BC_S == "neumann")
		m_BC_US = BC::NEUMANN;
	else if (BC_S == "wall" || BC_S == "dirichlet")
		m_BC_US = BC::DIRICHLET;
	// not supported
	else
		m_BC_US = BC::CUSTOM;

	if (BC_N == "periodic")
		m_BC_UN = BC::PERIODIC;
	else if (BC_N == "neumann")
		m_BC_UN = BC::NEUMANN;
	else if (BC_N == "wall" || BC_N == "dirichlet")
		m_BC_UN = BC::DIRICHLET;
	// not supported
	else
		m_BC_UN = BC::CUSTOM;

	if (BC_B == "periodic")
		m_BC_UB = BC::PERIODIC;
	else if (BC_B == "neumann")
		m_BC_UB = BC::NEUMANN;
	else if (BC_B == "wall" || BC_B == "dirichlet")
		m_BC_UB = BC::DIRICHLET;
	// not supported
	else
		m_BC_UB = BC::CUSTOM;

	if (BC_T == "periodic")
		m_BC_UT = BC::PERIODIC;
	else if (BC_T == "neumann")
		m_BC_UT = BC::NEUMANN;
	else if (BC_T == "wall" || BC_T == "dirichlet")
		m_BC_UT = BC::DIRICHLET;
	// not supported
	else
		m_BC_UT = BC::CUSTOM;

	return 0;
}

int BoundaryCondition::ApplyBC_U_3D(double ***arr) {
	BC_UW(arr);
	BC_UE(arr);
	BC_US(arr);
	BC_UN(arr);
	BC_UB(arr);
	BC_UT(arr);

	return 0;
}

int BoundaryCondition::SetBC_V_3D(std::string BC_W, std::string BC_E,
	std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

	std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
	std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
	std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
	std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
	std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
	std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);

	if (BC_W == "periodic")
		m_BC_VW = BC::PERIODIC;
	else if (BC_W == "neumann")
		m_BC_VW = BC::NEUMANN;
	else if (BC_W == "wall" || BC_W == "dirichlet")
		m_BC_VW = BC::DIRICHLET;
	// not supported
	else
		m_BC_VW = BC::CUSTOM;

	if (BC_E == "periodic")
		m_BC_VE = BC::PERIODIC;
	else if (BC_E == "neumann")
		m_BC_VE = BC::NEUMANN;
	else if (BC_E == "wall" || BC_E == "dirichlet")
		m_BC_VE = BC::DIRICHLET;
	// not supported
	else
		m_BC_VE = BC::CUSTOM;

	if (BC_S == "periodic")
		m_BC_VS = BC::PERIODIC;
	else if (BC_S == "neumann")
		m_BC_VS = BC::NEUMANN;
	else if (BC_S == "wall" || BC_S == "dirichlet")
		m_BC_VS = BC::DIRICHLET;
	// not supported
	else
		m_BC_VS = BC::CUSTOM;

	if (BC_N == "periodic")
		m_BC_VN = BC::PERIODIC;
	else if (BC_N == "neumann")
		m_BC_VN = BC::NEUMANN;
	else if (BC_N == "wall" || BC_N == "dirichlet")
		m_BC_VN = BC::DIRICHLET;
	// not supported
	else
		m_BC_VN = BC::CUSTOM;

	if (BC_B == "periodic")
		m_BC_VB = BC::PERIODIC;
	else if (BC_B == "neumann")
		m_BC_VB = BC::NEUMANN;
	else if (BC_B == "wall" || BC_B == "dirichlet")
		m_BC_VB = BC::DIRICHLET;
	// not supported
	else
		m_BC_VB = BC::CUSTOM;

	if (BC_T == "periodic")
		m_BC_VT = BC::PERIODIC;
	else if (BC_T == "neumann")
		m_BC_VT = BC::NEUMANN;
	else if (BC_T == "wall" || BC_T == "dirichlet")
		m_BC_VT = BC::DIRICHLET;
	// not supported
	else
		m_BC_VT = BC::CUSTOM;

	return 0;
}

int BoundaryCondition::ApplyBC_V_3D(double ***arr) {
	BC_VW(arr);
	BC_VE(arr);
	BC_VS(arr);
	BC_VN(arr);
	BC_VB(arr);
	BC_VT(arr);

	return 0;
}

int BoundaryCondition::SetBC_W_3D(std::string BC_W, std::string BC_E,
	std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

	std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
	std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
	std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
	std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
	std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
	std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);

	if (BC_W == "periodic")
		m_BC_WW = BC::PERIODIC;
	else if (BC_W == "neumann")
		m_BC_WW = BC::NEUMANN;
	else if (BC_W == "wall" || BC_W == "dirichlet")
		m_BC_WW = BC::DIRICHLET;
	// not supported
	else
		m_BC_WW = BC::CUSTOM;

	if (BC_E == "periodic")
		m_BC_WE = BC::PERIODIC;
	else if (BC_E == "neumann")
		m_BC_WE = BC::NEUMANN;
	else if (BC_E == "wall" || BC_E == "dirichlet")
		m_BC_WE = BC::DIRICHLET;
	// not supported
	else
		m_BC_WE = BC::CUSTOM;

	if (BC_S == "periodic")
		m_BC_WS = BC::PERIODIC;
	else if (BC_S == "neumann")
		m_BC_WS = BC::NEUMANN;
	else if (BC_S == "wall" || BC_S == "dirichlet")
		m_BC_WS = BC::DIRICHLET;
	// not supported
	else
		m_BC_WS = BC::CUSTOM;

	if (BC_N == "periodic")
		m_BC_WN = BC::PERIODIC;
	else if (BC_N == "neumann")
		m_BC_WN = BC::NEUMANN;
	else if (BC_N == "wall" || BC_N == "dirichlet")
		m_BC_WN = BC::DIRICHLET;
	// not supported
	else
		m_BC_WN = BC::CUSTOM;

	if (BC_B == "periodic")
		m_BC_WB = BC::PERIODIC;
	else if (BC_B == "neumann")
		m_BC_WB = BC::NEUMANN;
	else if (BC_B == "wall" || BC_B == "dirichlet")
		m_BC_WB = BC::DIRICHLET;
	// not supported
	else
		m_BC_WB = BC::CUSTOM;

	if (BC_T == "periodic")
		m_BC_WT = BC::PERIODIC;
	else if (BC_T == "neumann")
		m_BC_WT = BC::NEUMANN;
	else if (BC_T == "wall" || BC_T == "dirichlet")
		m_BC_WT = BC::DIRICHLET;
	// not supported
	else
		m_BC_WT = BC::CUSTOM;

	return 0;
}

int BoundaryCondition::ApplyBC_W_3D(double ***arr) {
	BC_WW(arr);
	BC_WE(arr);
	BC_WS(arr);
	BC_WN(arr);
	BC_WB(arr);
	BC_WT(arr);

	return 0;
}

int BoundaryCondition::SetBC_P_3D(std::string BC_W, std::string BC_E,
	std::string BC_S, std::string BC_N, std::string BC_B, std::string BC_T) {

	std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
	std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
	std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
	std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
	std::transform(BC_B.begin(), BC_B.end(), BC_B.begin(), ::tolower);
	std::transform(BC_T.begin(), BC_T.end(), BC_T.begin(), ::tolower);
	if (BC_W == "periodic")
		m_BC_PW = BC::PERIODIC;
	else if (BC_W == "neumann")
		m_BC_PW = BC::NEUMANN;
	else if (BC_W == "wall" || BC_W == "dirichlet")
		m_BC_PW = BC::DIRICHLET;
	// not supported
	else
		m_BC_PW = BC::CUSTOM;

	if (BC_E == "periodic")
		m_BC_PE = BC::PERIODIC;
	else if (BC_E == "neumann")
		m_BC_PE = BC::NEUMANN;
	else if (BC_E == "wall" || BC_E == "dirichlet")
		m_BC_PE = BC::DIRICHLET;
	// not supported
	else
		m_BC_PE = BC::CUSTOM;

	if (BC_S == "periodic")
		m_BC_PS = BC::PERIODIC;
	else if (BC_S == "neumann")
		m_BC_PS = BC::NEUMANN;
	else if (BC_S == "wall" || BC_S == "dirichlet")
		m_BC_PS = BC::DIRICHLET;
	// not supported
	else
		m_BC_PS = BC::CUSTOM;

	if (BC_N == "periodic")
		m_BC_PN = BC::PERIODIC;
	else if (BC_N == "neumann")
		m_BC_PN = BC::NEUMANN;
	else if (BC_N == "wall" || BC_N == "dirichlet")
		m_BC_PN = BC::DIRICHLET;
	// not supported
	else
		m_BC_PN = BC::CUSTOM;

	if (BC_B == "periodic")
		m_BC_PB = BC::PERIODIC;
	else if (BC_B == "neumann")
		m_BC_PB = BC::NEUMANN;
	else if (BC_B == "wall" || BC_B == "dirichlet")
		m_BC_PB = BC::DIRICHLET;
	// not supported
	else
		m_BC_PB = BC::CUSTOM;

	if (BC_T == "periodic")
		m_BC_PT = BC::PERIODIC;
	else if (BC_T == "neumann")
		m_BC_PT = BC::NEUMANN;
	else if (BC_T == "wall" || BC_T == "dirichlet")
		m_BC_PT = BC::DIRICHLET;
	// not supported
	else
		m_BC_PT = BC::CUSTOM;

	return 0;
}

int BoundaryCondition::ApplyBC_P_3D(double ***arr) {
	BC_PW(arr);
	BC_PE(arr);
	BC_PS(arr);
	BC_PN(arr);
	BC_PB(arr);
	BC_PT(arr);

	return 0;
}

void BoundaryCondition::BC_UW(double ***arr) {
	if (m_BC_UW == BC::PERIODIC)
		this->BC_PeriodicUW(arr);
	else if (m_BC_UW == BC::NEUMANN)
		this->BC_NeumannUW(arr);
	else if (m_BC_UW == BC::DIRICHLET)
		this->BC_DirichletUW(arr);
}

void BoundaryCondition::BC_UE(double ***arr) {
	if (m_BC_UE == BC::PERIODIC)
		this->BC_PeriodicUE(arr);
	else if (m_BC_UE == BC::NEUMANN)
		this->BC_NeumannUE(arr);
	else if (m_BC_UE == BC::DIRICHLET)
		this->BC_DirichletUE(arr);
}

void BoundaryCondition::BC_US(double ***arr) {
	if (m_BC_US == BC::PERIODIC)
		this->BC_PeriodicUS(arr);
	else if (m_BC_US == BC::NEUMANN)
		this->BC_NeumannUS(arr);
	else if (m_BC_US == BC::DIRICHLET)
		this->BC_DirichletUS(arr);
}

void BoundaryCondition::BC_UN(double ***arr) {
	if (m_BC_UN == BC::PERIODIC)
		this->BC_PeriodicUN(arr);
	else if (m_BC_UN == BC::NEUMANN)
		this->BC_NeumannUN(arr);
	else if (m_BC_UN == BC::DIRICHLET)
		this->BC_DirichletUN(arr);
}

void BoundaryCondition::BC_UB(double ***arr) {
	if (m_BC_UB == BC::PERIODIC)
		this->BC_PeriodicUB(arr);
	else if (m_BC_UB == BC::NEUMANN)
		this->BC_NeumannUB(arr);
	else if (m_BC_UB == BC::DIRICHLET)
		this->BC_DirichletUB(arr);
}

void BoundaryCondition::BC_UT(double ***arr) {
	if (m_BC_UT == BC::PERIODIC)
		this->BC_PeriodicUT(arr);
	else if (m_BC_UT == BC::NEUMANN)
		this->BC_NeumannUT(arr);
	else if (m_BC_UT == BC::DIRICHLET)
		this->BC_DirichletUT(arr);
}

void BoundaryCondition::BC_VW(double ***arr) {
	if (m_BC_VW == BC::PERIODIC)
		this->BC_PeriodicVW(arr);
	else if (m_BC_VW == BC::NEUMANN)
		this->BC_NeumannVW(arr);
	else if (m_BC_VW == BC::DIRICHLET)
		this->BC_DirichletVW(arr);
}

void BoundaryCondition::BC_VE(double ***arr) {
	if (m_BC_VE == BC::PERIODIC)
		this->BC_PeriodicVE(arr);
	else if (m_BC_VE == BC::NEUMANN)
		this->BC_NeumannVE(arr);
	else if (m_BC_VE == BC::DIRICHLET)
		this->BC_DirichletVE(arr);
}

void BoundaryCondition::BC_VS(double ***arr) {
	if (m_BC_VS == BC::PERIODIC)
		this->BC_PeriodicVS(arr);
	else if (m_BC_VS == BC::NEUMANN)
		this->BC_NeumannVS(arr);
	else if (m_BC_VS == BC::DIRICHLET)
		this->BC_DirichletVS(arr);
}

void BoundaryCondition::BC_VN(double ***arr) {
	if (m_BC_VN == BC::PERIODIC)
		this->BC_PeriodicVN(arr);
	else if (m_BC_VN == BC::NEUMANN)
		this->BC_NeumannVN(arr);
	else if (m_BC_VN == BC::DIRICHLET)
		this->BC_DirichletVN(arr);
}

void BoundaryCondition::BC_VB(double ***arr) {
	if (m_BC_VB == BC::PERIODIC)
		this->BC_PeriodicVB(arr);
	else if (m_BC_VB == BC::NEUMANN)
		this->BC_NeumannVB(arr);
	else if (m_BC_VB == BC::DIRICHLET)
		this->BC_DirichletVB(arr);
}

void BoundaryCondition::BC_VT(double ***arr) {
	if (m_BC_VT == BC::PERIODIC)
		this->BC_PeriodicVT(arr);
	else if (m_BC_VT == BC::NEUMANN)
		this->BC_NeumannVT(arr);
	else if (m_BC_VT == BC::DIRICHLET)
		this->BC_DirichletVT(arr);
}

void BoundaryCondition::BC_WW(double ***arr) {
	if (m_BC_WW == BC::PERIODIC)
		this->BC_PeriodicWW(arr);
	else if (m_BC_WW == BC::NEUMANN)
		this->BC_NeumannWW(arr);
	else if (m_BC_WW == BC::DIRICHLET)
		this->BC_DirichletWW(arr);
}

void BoundaryCondition::BC_WE(double ***arr) {
	if (m_BC_WE == BC::PERIODIC)
		this->BC_PeriodicWE(arr);
	else if (m_BC_WE == BC::NEUMANN)
		this->BC_NeumannWE(arr);
	else if (m_BC_WE == BC::DIRICHLET)
		this->BC_DirichletWE(arr);
}

void BoundaryCondition::BC_WS(double ***arr) {
	if (m_BC_WS == BC::PERIODIC)
		this->BC_PeriodicWS(arr);
	else if (m_BC_WS == BC::NEUMANN)
		this->BC_NeumannWS(arr);
	else if (m_BC_WS == BC::DIRICHLET)
		this->BC_DirichletWS(arr);
}

void BoundaryCondition::BC_WN(double ***arr) {
	if (m_BC_WN == BC::PERIODIC)
		this->BC_PeriodicWN(arr);
	else if (m_BC_WN == BC::NEUMANN)
		this->BC_NeumannWN(arr);
	else if (m_BC_WN == BC::DIRICHLET)
		this->BC_DirichletWN(arr);
}

void BoundaryCondition::BC_WB(double ***arr) {
	if (m_BC_WB == BC::PERIODIC)
		this->BC_PeriodicWB(arr);
	else if (m_BC_WB == BC::NEUMANN)
		this->BC_NeumannWB(arr);
	else if (m_BC_WB == BC::DIRICHLET)
		this->BC_DirichletWB(arr);
}

void BoundaryCondition::BC_WT(double ***arr) {
	if (m_BC_WT == BC::PERIODIC)
		this->BC_PeriodicWT(arr);
	else if (m_BC_WT == BC::NEUMANN)
		this->BC_NeumannWT(arr);
	else if (m_BC_WT == BC::DIRICHLET)
		this->BC_DirichletWT(arr);
}

void BoundaryCondition::BC_PW(double ***arr) {
	if (m_BC_PW == BC::PERIODIC)
		this->BC_PeriodicPW(arr);
	else if (m_BC_PW == BC::NEUMANN)
		this->BC_NeumannPW(arr);
	else if (m_BC_PW == BC::DIRICHLET)
		this->BC_DirichletPW(arr);
}

void BoundaryCondition::BC_PE(double ***arr) {
	if (m_BC_PE == BC::PERIODIC)
		this->BC_PeriodicPE(arr);
	else if (m_BC_PE == BC::NEUMANN)
		this->BC_NeumannPE(arr);
	else if (m_BC_PE == BC::DIRICHLET)
		this->BC_DirichletPE(arr);
}

void BoundaryCondition::BC_PS(double ***arr) {
	if (m_BC_PS == BC::PERIODIC)
		this->BC_PeriodicPS(arr);
	else if (m_BC_PS == BC::NEUMANN)
		this->BC_NeumannPS(arr);
	else if (m_BC_PS == BC::DIRICHLET)
		this->BC_DirichletPS(arr);
}

void BoundaryCondition::BC_PN(double ***arr) {
	if (m_BC_PN == BC::PERIODIC)
		this->BC_PeriodicPN(arr);
	else if (m_BC_PN == BC::NEUMANN)
		this->BC_NeumannPN(arr);
	else if (m_BC_PN == BC::DIRICHLET)
		this->BC_DirichletPN(arr);
}

void BoundaryCondition::BC_PB(double ***arr) {
	if (m_BC_PB == BC::PERIODIC)
		this->BC_PeriodicPB(arr);
	else if (m_BC_PB == BC::NEUMANN)
		this->BC_NeumannPB(arr);
	else if (m_BC_PB == BC::DIRICHLET)
		this->BC_DirichletPB(arr);
}

void BoundaryCondition::BC_PT(double ***arr) {
	if (m_BC_PT == BC::PERIODIC)
		this->BC_PeriodicPT(arr);
	else if (m_BC_PT == BC::NEUMANN)
		this->BC_NeumannPT(arr);
	else if (m_BC_PT == BC::DIRICHLET)
		this->BC_DirichletPT(arr);
}

void BoundaryCondition::BC_PeriodicUW(double ***arr) {
	for (int i = 0; i <= kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNx - 1 + i][j][k];
}

void BoundaryCondition::BC_PeriodicUE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + i + kNumBCGrid][j][k] = arr[i + kNumBCGrid + 1][j][k];
}

void BoundaryCondition::BC_PeriodicUS(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_PeriodicUN(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + j + kNumBCGrid][k] = arr[i][j + kNumBCGrid][k];
}

void BoundaryCondition::BC_PeriodicUB(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNz + k];
}

void BoundaryCondition::BC_PeriodicUT(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + k + kNumBCGrid] = arr[i][j][k + kNumBCGrid];
}

void BoundaryCondition::BC_NeumannUW(double ***arr) {
	for (int i = 0; i <= kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNumBCGrid * 2 + 1 - i][j][k];
}

void BoundaryCondition::BC_NeumannUE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_NeumannUS(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNumBCGrid * 2 - j - 1][k];
}

void BoundaryCondition::BC_NeumannUN(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_NeumannUB(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNumBCGrid * 2 - k - 1];
}

void BoundaryCondition::BC_NeumannUT(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = arr[i][j][kNz + k];
}

// Prescribed bounadry condition (Dirichlet) for U grid
void BoundaryCondition::SetBCConstantUW(double BC_ConstantW) {
	m_BC_DirichletConstantUW = BC_ConstantW;
}

void BoundaryCondition::SetBCConstantUE(double BC_ConstantE) {
	m_BC_DirichletConstantUE = BC_ConstantE;
}

void BoundaryCondition::SetBCConstantUS(double BC_ConstantS) {
	m_BC_DirichletConstantUS = BC_ConstantS;
}

void BoundaryCondition::SetBCConstantUN(double BC_ConstantN) {
	m_BC_DirichletConstantUN = BC_ConstantN;
}

void BoundaryCondition::SetBCConstantUB(double BC_ConstantB) {
	m_BC_DirichletConstantUB = BC_ConstantB;
}

void BoundaryCondition::SetBCConstantUT(double BC_ConstantT) {
	m_BC_DirichletConstantUT = BC_ConstantT;
}

void BoundaryCondition::BC_DirichletUW(double ***arr) {
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
			arr[kNumBCGrid][j][k] = m_BC_DirichletConstantUW;

	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[kNumBCGrid * 2 - i][j][k]
					+ 2.0 * m_BC_DirichletConstantUW;
}

void BoundaryCondition::BC_DirichletUE(double ***arr) {
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
			arr[kNx + kNumBCGrid][j][k] = m_BC_DirichletConstantUE;

	for (int i = 1; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i][j][k] = -arr[kNx + i][j][k]
				+ 2.0 * m_BC_DirichletConstantUE;
}

void BoundaryCondition::BC_DirichletUS(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][kNumBCGrid * 2 - j - 1][k]
				+ 2.0 * m_BC_DirichletConstantUS;
}

void BoundaryCondition::BC_DirichletUN(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = -arr[i][kNy + j][k]
				+ 2.0 * m_BC_DirichletConstantUN;
}

void BoundaryCondition::BC_DirichletUB(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][j][kNumBCGrid * 2 - k - 1]
				+ 2.0 * m_BC_DirichletConstantUB;
}

void BoundaryCondition::BC_DirichletUT(double ***arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++) 
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = -arr[i][j][kNz + k]
				+ 2.0 * m_BC_DirichletConstantUT;
}

void BoundaryCondition::BC_PeriodicVW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_PeriodicVE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + i + kNumBCGrid][j][k] = arr[i + kNumBCGrid][j][k];
}

void BoundaryCondition::BC_PeriodicVS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j <= kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNy - 1 + j][k];
}

void BoundaryCondition::BC_PeriodicVN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + j + kNumBCGrid][k] = arr[i][j + kNumBCGrid + 1][k];
}

void BoundaryCondition::BC_PeriodicVB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNz + k];
}

void BoundaryCondition::BC_PeriodicVT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + k + kNumBCGrid] = arr[i][j][k + kNumBCGrid];
}

void BoundaryCondition::BC_NeumannVW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNumBCGrid * 2 - i - 1][j][k];
}

void BoundaryCondition::BC_NeumannVE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_NeumannVS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j <= kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNumBCGrid * 2 + 1 - j][k];
}

void BoundaryCondition::BC_NeumannVN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_NeumannVB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNumBCGrid * 2 - k - 1];
}

void BoundaryCondition::BC_NeumannVT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = arr[i][j][kNz + k];
}

// Prescribed bounadry condition (Dirichlet) for V grid
void BoundaryCondition::SetBCConstantVW(double BC_ConstantW) {
	m_BC_DirichletConstantVW = BC_ConstantW;
}

void BoundaryCondition::SetBCConstantVE(double BC_ConstantE) {
	m_BC_DirichletConstantVE = BC_ConstantE;
}

void BoundaryCondition::SetBCConstantVS(double BC_ConstantS) {
	m_BC_DirichletConstantVS = BC_ConstantS;
}

void BoundaryCondition::SetBCConstantVN(double BC_ConstantN) {
	m_BC_DirichletConstantVN = BC_ConstantN;
}

void BoundaryCondition::SetBCConstantVB(double BC_ConstantB) {
	m_BC_DirichletConstantVB = BC_ConstantB;
}

void BoundaryCondition::SetBCConstantVT(double BC_ConstantT) {
	m_BC_DirichletConstantVT = BC_ConstantT;
}

void BoundaryCondition::BC_DirichletVW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[kNumBCGrid * 2 - i - 1][j][k]
				+ 2.0 * m_BC_DirichletConstantVW;
}

void BoundaryCondition::BC_DirichletVE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = -arr[kNx + i][j][k]
				+ 2.0 * m_BC_DirichletConstantVE;
}

void BoundaryCondition::BC_DirichletVS(double ***arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
			arr[i][kNumBCGrid][k] = m_BC_DirichletConstantVS;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][kNumBCGrid * 2 - j][k]
					+ 2.0 * m_BC_DirichletConstantVS;
}

void BoundaryCondition::BC_DirichletVN(double ***arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
		for (int k = 0; k < kNz + 2 * kNumBCGrid; k++)
			arr[i][kNy + kNumBCGrid][k] = m_BC_DirichletConstantVN;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 1; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j][k] = -arr[i][kNy + j][k]
					+ 2.0 * m_BC_DirichletConstantVN;
}

void BoundaryCondition::BC_DirichletVB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][j][kNumBCGrid * 2 - k - 1]
				+ 2.0 * m_BC_DirichletConstantVB;
}

void BoundaryCondition::BC_DirichletVT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = -arr[i][j][kNz + k]
				+ 2.0 * m_BC_DirichletConstantVT;
}

void BoundaryCondition::BC_PeriodicWW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_PeriodicWE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[kNx + i + kNumBCGrid][j][k] = arr[i + kNumBCGrid][j][k];
}

void BoundaryCondition::BC_PeriodicWS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_PeriodicWN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + j + kNumBCGrid][k] = arr[i][j + kNumBCGrid][k];
}

void BoundaryCondition::BC_PeriodicWB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k <= kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNz - 1 + k];
}

void BoundaryCondition::BC_PeriodicWT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + k + kNumBCGrid] = arr[i][j][k + kNumBCGrid + 1];
}

void BoundaryCondition::BC_NeumannWW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNumBCGrid * 2 - i - 1][j][k];
}

void BoundaryCondition::BC_NeumannWE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_NeumannWS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNumBCGrid * 2 - j - 1][k];
}

void BoundaryCondition::BC_NeumannWN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_NeumannWB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k <= kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNumBCGrid * 2 + 1 - k];
}

void BoundaryCondition::BC_NeumannWT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = arr[i][j][kNz + k];
}

// Prescribed bounadry condition (Dirichlet) for V grid
void BoundaryCondition::SetBCConstantWW(double BC_ConstantW) {
	m_BC_DirichletConstantWW = BC_ConstantW;
}

void BoundaryCondition::SetBCConstantWE(double BC_ConstantE) {
	m_BC_DirichletConstantWE = BC_ConstantE;
}

void BoundaryCondition::SetBCConstantWS(double BC_ConstantS) {
	m_BC_DirichletConstantWS = BC_ConstantS;
}

void BoundaryCondition::SetBCConstantWN(double BC_ConstantN) {
	m_BC_DirichletConstantWN = BC_ConstantN;
}

void BoundaryCondition::SetBCConstantWB(double BC_ConstantB) {
	m_BC_DirichletConstantWB = BC_ConstantB;
}

void BoundaryCondition::SetBCConstantWT(double BC_ConstantT) {
	m_BC_DirichletConstantWT = BC_ConstantT;
}

void BoundaryCondition::BC_DirichletWW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[kNumBCGrid * 2 - i - 1][j][k]
				+ 2.0 * m_BC_DirichletConstantWW;
}

void BoundaryCondition::BC_DirichletWE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = -arr[kNx + i][j][k]
				+ 2.0 * m_BC_DirichletConstantWE;
}

void BoundaryCondition::BC_DirichletWS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][kNumBCGrid * 2 - j - 1][k]
				+ 2.0 * m_BC_DirichletConstantWS;
}

void BoundaryCondition::BC_DirichletWN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid + 1; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][kNy + kNumBCGrid * 2 - j - 1][k]
				+ 2.0 * m_BC_DirichletConstantWN;
}

void BoundaryCondition::BC_DirichletWB(double ***arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
			arr[i][j][kNumBCGrid] = m_BC_DirichletConstantWB;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][j][kNumBCGrid * 2 - k]
				+ 2.0 * m_BC_DirichletConstantWB;
}

void BoundaryCondition::BC_DirichletWT(double ***arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
		for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
			arr[i][j][kNz + kNumBCGrid] = m_BC_DirichletConstantWT;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 1; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k] = -arr[i][j][kNz + k]
					+ 2.0 * m_BC_DirichletConstantWT;
}

void BoundaryCondition::BC_PeriodicPW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_PeriodicPE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + i + kNumBCGrid][j][k] = arr[i + kNumBCGrid][j][k];
}

void BoundaryCondition::BC_PeriodicPS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_PeriodicPN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + j + kNumBCGrid][k] = arr[i][j + kNumBCGrid][k];
}

void BoundaryCondition::BC_PeriodicPB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNz + k];
}

void BoundaryCondition::BC_PeriodicPT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + k + kNumBCGrid] = arr[i][j][k + kNumBCGrid];
}

void BoundaryCondition::BC_NeumannPW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[kNumBCGrid * 2 - i - 1][j][k];
}

void BoundaryCondition::BC_NeumannPE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = arr[kNx + i][j][k];
}

void BoundaryCondition::BC_NeumannPS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = arr[i][kNumBCGrid * 2 - j - 1][k];
}

void BoundaryCondition::BC_NeumannPN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = arr[i][kNy + j][k];
}

void BoundaryCondition::BC_NeumannPB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = arr[i][j][kNumBCGrid * 2 - k - 1];
}

void BoundaryCondition::BC_NeumannPT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = arr[i][j][kNz + k];
}

// Prescribed bounadry condition (Dirichlet) for P grid (center)
void BoundaryCondition::SetBCConstantPW(double BC_ConstantW) {
	m_BC_DirichletConstantPW = BC_ConstantW;
}

void BoundaryCondition::SetBCConstantPE(double BC_ConstantE) {
	m_BC_DirichletConstantPE = BC_ConstantE;
}

void BoundaryCondition::SetBCConstantPS(double BC_ConstantS) {
	m_BC_DirichletConstantPS = BC_ConstantS;
}

void BoundaryCondition::SetBCConstantPN(double BC_ConstantN) {
	m_BC_DirichletConstantPN = BC_ConstantN;
}

void BoundaryCondition::SetBCConstantPB(double BC_ConstantB) {
	m_BC_DirichletConstantPB = BC_ConstantB;
}

void BoundaryCondition::SetBCConstantPT(double BC_ConstantT) {
	m_BC_DirichletConstantPT = BC_ConstantT;
}

void BoundaryCondition::BC_DirichletPW(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[kNumBCGrid * 2 - i - 1][j][k]
				+ 2.0 * m_BC_DirichletConstantPW;
}

void BoundaryCondition::BC_DirichletPE(double ***arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[kNx + kNumBCGrid * 2 - i - 1][j][k] = -arr[kNx + i][j][k]
				+ 2.0 * m_BC_DirichletConstantPE;
}

void BoundaryCondition::BC_DirichletPS(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][kNumBCGrid * 2 - j - 1][k]
				+ 2.0 * m_BC_DirichletConstantPS;
}

void BoundaryCondition::BC_DirichletPN(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			for (int k = kNumBCGrid; k < kNz + kNumBCGrid; k++)
				arr[i][kNy + kNumBCGrid * 2 - j - 1][k] = -arr[i][kNy + j][k]
				+ 2.0 * m_BC_DirichletConstantPN;
}

void BoundaryCondition::BC_DirichletPB(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][k] = -arr[i][j][kNumBCGrid * 2 - k - 1]
				+ 2.0 * m_BC_DirichletConstantPB;
}

void BoundaryCondition::BC_DirichletPT(double ***arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			for (int k = 0; k < kNumBCGrid; k++)
				arr[i][j][kNz + kNumBCGrid * 2 - k - 1] = -arr[i][j][kNz + k]
				+ 2.0 * m_BC_DirichletConstantPT;
}
