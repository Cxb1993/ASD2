#include "bc2d.h"

BoundaryCondition2D::BoundaryCondition2D(int nx, int ny, int num_bc_grid)
	: kNx(nx), kNy(ny), kNumBCGrid(num_bc_grid) {
}

int BoundaryCondition2D::SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {

	std::transform(BC_W.begin(), BC_W.end(), BC_W.begin(), ::tolower);
	std::transform(BC_E.begin(), BC_E.end(), BC_E.begin(), ::tolower);
	std::transform(BC_S.begin(), BC_S.end(), BC_S.begin(), ::tolower);
	std::transform(BC_N.begin(), BC_N.end(), BC_N.begin(), ::tolower);
	
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

	return 0;
}


int BoundaryCondition2D::ApplyBC_V_2D(std::vector<double>& arr) {
	BC_VW(arr);
	BC_VE(arr);
	BC_VS(arr);
	BC_VN(arr);

	return 0;
}

int BoundaryCondition2D::SetBC_P_2D(std::string BC_W, std::string BC_E,
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
	if (m_BC_UW == BC::PERIODIC)
		this->BC_PeriodicUW(arr);
	else if (m_BC_UW == BC::NEUMANN)
		this->BC_NeumannUW(arr);
	else if (m_BC_UW == BC::DIRICHLET)
		this->BC_DirichletUW(arr);
}


void BoundaryCondition2D::BC_UE(std::vector<double>& arr) {
	if (m_BC_UE == BC::PERIODIC)
		this->BC_PeriodicUE(arr);
	else if (m_BC_UE == BC::NEUMANN)
		this->BC_NeumannUE(arr);
	else if (m_BC_UE == BC::DIRICHLET)
		this->BC_DirichletUE(arr);
}


void BoundaryCondition2D::BC_US(std::vector<double>& arr) {
	if (m_BC_US == BC::PERIODIC)
		this->BC_PeriodicUS(arr);
	else if (m_BC_US == BC::NEUMANN)
		this->BC_NeumannUS(arr);
	else if (m_BC_US == BC::DIRICHLET)
		this->BC_DirichletUS(arr);
}


void BoundaryCondition2D::BC_UN(std::vector<double>& arr) {
	if (m_BC_UN == BC::PERIODIC)
		this->BC_PeriodicUN(arr);
	else if (m_BC_UN == BC::NEUMANN)
		this->BC_NeumannUN(arr);
	else if (m_BC_UN == BC::DIRICHLET)
		this->BC_DirichletUN(arr);
}


void BoundaryCondition2D::BC_VW(std::vector<double>& arr) {
	if (m_BC_VW == BC::PERIODIC)
		this->BC_PeriodicVW(arr);
	else if (m_BC_VW == BC::NEUMANN)
		this->BC_NeumannVW(arr);
	else if (m_BC_VW == BC::DIRICHLET)
		this->BC_DirichletVW(arr);
}


void BoundaryCondition2D::BC_VE(std::vector<double>& arr) {
	if (m_BC_VE == BC::PERIODIC)
		this->BC_PeriodicVE(arr);
	else if (m_BC_VE == BC::NEUMANN)
		this->BC_NeumannVE(arr);
	else if (m_BC_VE == BC::DIRICHLET)
		this->BC_DirichletVE(arr);
}


void BoundaryCondition2D::BC_VS(std::vector<double>& arr) {
	if (m_BC_VS == BC::PERIODIC)
		this->BC_PeriodicVS(arr);
	else if (m_BC_VS == BC::NEUMANN)
		this->BC_NeumannVS(arr);
	else if (m_BC_VS == BC::DIRICHLET)
		this->BC_DirichletVS(arr);
}


void BoundaryCondition2D::BC_VN(std::vector<double>& arr) {
	if (m_BC_VN == BC::PERIODIC)
		this->BC_PeriodicVN(arr);
	else if (m_BC_VN == BC::NEUMANN)
		this->BC_NeumannVN(arr);
	else if (m_BC_VN == BC::DIRICHLET)
		this->BC_DirichletVN(arr);
}


void BoundaryCondition2D::BC_PW(std::vector<double>& arr) {
	if (m_BC_PW == BC::PERIODIC)
		this->BC_PeriodicPW(arr);
	else if (m_BC_PW == BC::NEUMANN)
		this->BC_NeumannPW(arr);
	else if (m_BC_PW == BC::DIRICHLET)
		this->BC_DirichletPW(arr);
}


void BoundaryCondition2D::BC_PE(std::vector<double>& arr) {
	if (m_BC_PE == BC::PERIODIC)
		this->BC_PeriodicPE(arr);
	else if (m_BC_PE == BC::NEUMANN)
		this->BC_NeumannPE(arr);
	else if (m_BC_PE == BC::DIRICHLET)
		this->BC_DirichletPE(arr);
}


void BoundaryCondition2D::BC_PS(std::vector<double>& arr) {
	if (m_BC_PS == BC::PERIODIC)
		this->BC_PeriodicPS(arr);
	else if (m_BC_PS == BC::NEUMANN)
		this->BC_NeumannPS(arr);
	else if (m_BC_PS == BC::DIRICHLET)
		this->BC_DirichletPS(arr);
}


void BoundaryCondition2D::BC_PN(std::vector<double>& arr) {
	if (m_BC_PN == BC::PERIODIC)
		this->BC_PeriodicPN(arr);
	else if (m_BC_PN == BC::NEUMANN)
		this->BC_NeumannPN(arr);
	else if (m_BC_PN == BC::DIRICHLET)
		this->BC_DirichletPN(arr);
}


void BoundaryCondition2D::BC_PeriodicUW(std::vector<double>& arr) {
	for (int i = 0; i <= kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNx - 1 + i + (kNx + 2 * kNumBCGrid) * j];
}


void BoundaryCondition2D::BC_PeriodicUE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + i + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j] = arr[i + kNumBCGrid + 1 + (kNx + 2 * kNumBCGrid) * j];
}


void BoundaryCondition2D::BC_PeriodicUS(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)];
}


void BoundaryCondition2D::BC_PeriodicUN(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j + kNumBCGrid)] = arr[i + (kNx + 2 * kNumBCGrid) * (j + kNumBCGrid)];
}


void BoundaryCondition2D::BC_NeumannUW(std::vector<double>& arr) {
	for (int i = 0; i <= kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNumBCGrid * 2 + 1 - i + (kNx + 2 * kNumBCGrid) * j];
}


void BoundaryCondition2D::BC_NeumannUE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j] = arr[kNx + i + (kNx + 2 * kNumBCGrid) * j];
}


void BoundaryCondition2D::BC_NeumannUS(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
				arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 - j - 1)];
}


void BoundaryCondition2D::BC_NeumannUN(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j - 1)] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)];
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
		arr[kNumBCGrid + (kNx + 2 * kNumBCGrid) * j] = m_BC_DirichletConstantUW;

	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[kNumBCGrid * 2 - i + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantUW;
}

void BoundaryCondition2D::BC_DirichletUE(std::vector<double>& arr) {
	for (int j = 0; j < kNy + 2 * kNumBCGrid; j++)
		arr[kNx + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j] = m_BC_DirichletConstantUE;

	for (int i = 1; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i + (kNx + 2 * kNumBCGrid) * j] = -arr[kNx + i + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantUE;
}

void BoundaryCondition2D::BC_DirichletUS(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 - j - 1)]
				+ 2.0 * m_BC_DirichletConstantUS;
}

void BoundaryCondition2D::BC_DirichletUN(std::vector<double>& arr) {
	for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j - 1)] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)]
				+ 2.0 * m_BC_DirichletConstantUN;
}

void BoundaryCondition2D::BC_PeriodicVW(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNx + i + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_PeriodicVE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			arr[kNx + i + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j] = arr[i + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_PeriodicVS(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j <= kNumBCGrid; j++)
				arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy - 1 + j)];
}

void BoundaryCondition2D::BC_PeriodicVN(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
				arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j + kNumBCGrid)] = arr[i + (kNx + 2 * kNumBCGrid) * (j + kNumBCGrid + 1)];
}

void BoundaryCondition2D::BC_NeumannVW(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_NeumannVE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j] = arr[kNx + i + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_NeumannVS(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j <= kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 + 1 - j)];
}

void BoundaryCondition2D::BC_NeumannVN(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j - 1)] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)];
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
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantVW;
}

void BoundaryCondition2D::BC_DirichletVE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j] = -arr[kNx + i + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantVE;
}

void BoundaryCondition2D::BC_DirichletVS(std::vector<double>& arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
			arr[i + (kNx + 2 * kNumBCGrid) * kNumBCGrid] = m_BC_DirichletConstantVS;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 - j)]
					+ 2.0 * m_BC_DirichletConstantVS;
}

void BoundaryCondition2D::BC_DirichletVN(std::vector<double>& arr) {
	for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid)] = m_BC_DirichletConstantVN;

	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 1; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j)] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)]
				+ 2.0 * m_BC_DirichletConstantVN;
}

void BoundaryCondition2D::BC_PeriodicPW(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNx + i + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_PeriodicPE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + i + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j] = arr[i + kNumBCGrid + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_PeriodicPS(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)];
}

void BoundaryCondition2D::BC_PeriodicPN(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j + kNumBCGrid)] = arr[i + (kNx + 2 * kNumBCGrid) * (j + kNumBCGrid)];
}

void BoundaryCondition2D::BC_NeumannPW(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_NeumannPE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j] = arr[kNx + i + (kNx + 2 * kNumBCGrid) * j];
}

void BoundaryCondition2D::BC_NeumannPS(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 - j - 1)];
}


void BoundaryCondition2D::BC_NeumannPN(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j - 1)] = arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)];
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

void BoundaryCondition2D::BC_DirichletPW(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantPW;
}

void BoundaryCondition2D::BC_DirichletPE(std::vector<double>& arr) {
	for (int i = 0; i < kNumBCGrid; i++)
		for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
			arr[kNx + kNumBCGrid * 2 - i - 1 + (kNx + 2 * kNumBCGrid) * j] = -arr[kNx + i + (kNx + 2 * kNumBCGrid) * j]
				+ 2.0 * m_BC_DirichletConstantPE;
}

void BoundaryCondition2D::BC_DirichletPS(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * j] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNumBCGrid * 2 - j - 1)]
				+ 2.0 * m_BC_DirichletConstantPS;
}

void BoundaryCondition2D::BC_DirichletPN(std::vector<double>& arr) {
	for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
		for (int j = 0; j < kNumBCGrid; j++)
			arr[i + (kNx + 2 * kNumBCGrid) * (kNy + kNumBCGrid * 2 - j - 1)] = -arr[i + (kNx + 2 * kNumBCGrid) * (kNy + j)]
				+ 2.0 * m_BC_DirichletConstantPN;
}
