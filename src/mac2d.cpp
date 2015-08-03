#include "mac2d.h"

MACSolver2D::MACSolver2D(double Re, double We, double Fr, GAXISENUM2D GAxis,
    double L, double U, double sigma, double densityRatio, double viscosityRatio, double rhoH, double muH,
    int nx, int ny, double baseX, double baseY, double lenX, double lenY,
    TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
    bool writeVTK) :
    kRe(Re), kWe(We), kFr(Fr),
    kLScale(L), kUScale(U), kSigma(sigma),
    kG(kFr * L / (U * U)), kGAxis(GAxis), kRhoScale(We * sigma / (L * U * U)), kMuScale(kRhoScale * L * U / Re),
    kRhoH(rhoH), kRhoL(rhoH * densityRatio), kRhoRatio(densityRatio),
    kMuH(muH), kMuL(muH * viscosityRatio), kMuRatio(viscosityRatio),
    kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
    kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
    kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
    kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
    kArrSize(
    static_cast<int64_t>(kNx + 2 * kNumBCGrid) *
    static_cast<int64_t>(kNy + 2 * kNumBCGrid)) {

    // positive level set : inside
    // negative level set : outside
    m_iter = 0;
    m_curTime = 0.0;
    m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::MACSolver2D(double rhoH, double rhoL, double muH, double muL, double gConstant, GAXISENUM2D GAxis,
    double L, double U, double sigma, int nx, int ny, double baseX, double baseY, double lenX, double lenY,
    TIMEORDERENUM timeOrder, double cfl, double maxtime, int maxIter, int niterskip, int num_bc_grid,
    bool writeVTK) :
    kRhoScale(rhoH), kMuScale(muH), kG(gConstant), kLScale(L), kUScale(U), kSigma(sigma),
    kRhoH(rhoH), kRhoL(rhoL), kRhoRatio(rhoL / rhoH), 
    kMuH(muH), kMuL(muL), kMuRatio(muL / muH),
    kRe(rhoH * L * U / muH), kWe(rhoH * L * U * U / sigma), kFr(U * U / (gConstant * L)), kGAxis(GAxis),
    kNx(nx), kNy(ny), kBaseX(baseX), kBaseY(baseY), kLenX(lenX), kLenY(lenY),
    kDx(lenX / static_cast<double>(nx)), kDy(lenY / static_cast<double>(ny)),
    kTimeOrder(timeOrder), kCFL(cfl), kMaxTime(maxtime), kMaxIter(maxIter), kNIterSkip(niterskip),
    kNumBCGrid(num_bc_grid), kWriteVTK(writeVTK),
    kArrSize(
    static_cast<int64_t>(kNx + 2 * kNumBCGrid) *
    static_cast<int64_t>(kNy + 2 * kNumBCGrid)) {

    // positive level set : inside
    // negative level set : outside
    m_iter = 0;
    m_curTime = 0.0;
    m_dt = std::min(0.005, kCFL * std::min(lenX / static_cast<double>(nx), lenY / static_cast<double>(ny)) / U);
}

MACSolver2D::~MACSolver2D() {
    if (m_BC)
        m_BC.reset();
}

// Deallocated automatically
int MACSolver2D::AllocateVariables() {
    m_u = std::vector<double>(kArrSize, 0);
    m_v = std::vector<double>(kArrSize, 0);
    
    m_rho = std::vector<double>(kArrSize, 0);
    m_kappa = std::vector<double>(kArrSize, 0);
    m_mu = std::vector<double>(kArrSize, 0);
    
    m_ps = std::vector<double>(kArrSize, 0);
    m_p = std::vector<double>(kArrSize, 0);

    m_J11 = std::vector<double>(kArrSize, 0);
    m_J12 = std::vector<double>(kArrSize, 0);
    m_J21 = std::vector<double>(kArrSize, 0);
    m_J22 = std::vector<double>(kArrSize, 0);

    m_nx = std::vector<double>(kArrSize, 0);
    m_ny = std::vector<double>(kArrSize, 0);

    m_t1x = std::vector<double>(kArrSize, 0);
    m_t1y = std::vector<double>(kArrSize, 0);
    m_t1z = std::vector<double>(kArrSize, 0);

    m_t2x = std::vector<double>(kArrSize, 0);
    m_t2y = std::vector<double>(kArrSize, 0);
    m_t2z = std::vector<double>(kArrSize, 0);

    return 0;
}

std::vector<double> MACSolver2D::UpdateHeavisideFunc(const std::vector<double>& ls) {
    std::vector<double> Heaviside(kArrSize, 0.0);
    const double eps = std::min(kDx, kDy) * 1.5;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
        // inside
        if (ls[idx(i, j)] >= 0.0)
            Heaviside[idx(i, j)] = 1.0;
        else
            Heaviside[idx(i, j)] = 0.0;
    }

    return Heaviside;
}

std::vector<double> MACSolver2D::UpdateSmoothHeavisideFunc(const std::vector<double>& ls) {
    std::vector<double> Heaviside(kArrSize, 0.0);
    const double eps = std::min(kDx, kDy) * 1.5;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
        // inside
        if (ls[idx(i, j)] > eps)
            Heaviside[idx(i, j)] = 1.0;
        // outside
        else if (ls[idx(i, j)] < -eps)
            Heaviside[idx(i, j)] = 0.0;
        else
            Heaviside[idx(i, j)] 
                = 0.5 * (1.0 + ls[idx(i, j)] / eps
                    + 1.0 / M_PI * sin(M_PI * ls[idx(i, j)] / eps));
    }
        
    return Heaviside;
}

int MACSolver2D::UpdateNTK(const std::shared_ptr<LevelSetSolver2D>& LSolver, const std::vector<double>& ls,
    std::tuple<std::vector<double>&, std::vector<double>&> normalVec,
    std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t1Vec,
    std::tuple<std::vector<double>&, std::vector<double>&, std::vector<double>&> t2Vec,
    std::vector<double>& kappa) {
    // Lervåg, Karl Yngve, Bernhard Müller, and Svend Tollak Munkejord.
    // "Calculation of the interface curvature and normal vector with the level-set method."
    //  Computers & Fluids 84 (2013) : 218 - 230.
    //
    // Macklin, Paul, and John S. Lowengrub.
    //  "A new ghost cell/level set method for moving boundary problems: application to tumor growth."
    // Journal of scientific computing 35.2-3 (2008): 266-299.
    const double eta = 0.0001;
    const double eps = 1.0e-100;

    std::vector<double> Q(kArrSize, 0.0), absDerivLS(kArrSize, 0.0);

    absDerivLS = LSolver->ENO_DerivAbsLS_2D(ls, ls);

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        Q[idx(i, j)] = std::fabs(1.0 - absDerivLS[idx(i, j)]);
    }
    
    double dLSdX = 0.0, dLSdY = 0.0;
    // determination function
    int Dx = 0.0, Dy = 0.0;
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        // complex level set, such as droplet merging
        if (Q[idx(i, j)] > eta) {
            // Using Dierctional Direction First
            if (Q[idx(i - 1, j)] < eta && Q[idx(i + 1, j)] >= eta)
                Dx = -1;
            else if (Q[idx(i - 1, j)] >= eta && Q[idx(i + 1, j)] < eta)
                Dx = 1;
            else if (Q[idx(i - 1, j)] < eta && Q[idx(i, j)] < eta && Q[idx(i + 1, j)] < eta)
                Dx = 0;
            else if (Q[idx(i - 1, j)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i + 1, j)] >= eta)
                Dx = 0;
            else
                Dx = 2;

            // determination funciton
            if (Q[idx(i, j - 1)] < eta && Q[idx(i, j + 1)] >= eta)
                Dy = -1;
            else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j + 1)] < eta)
                Dy = 1;
            else if (Q[idx(i, j - 1)] < eta && Q[idx(i, j)] < eta && Q[idx(i, j + 1)] < eta)
                Dy = 0;
            else if (Q[idx(i, j - 1)] >= eta && Q[idx(i, j)] >= eta && Q[idx(i, j + 1)] >= eta)
                Dy = 0;
            else
                // undetermined
                Dy = 2;

            if (Dx == -1)
                dLSdX = (ls[idx(i, j)] - ls[idx(i - 1, j)]) / kDx;
            else if (Dx == 1)
                dLSdX = (ls[idx(i + 1, j)] - ls[idx(i, j)]) / kDx;
            else if (Dx == 0)
                dLSdX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);

            if (Dy == -1)
                dLSdY = (ls[idx(i, j)] - ls[idx(i, j - 1)]) / kDy;
            else if (Dy == 1)
                dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j)]) / kDy;
            else if (Dy == 0)
                dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);

            // Curve Fitting Method
            if (Dx == 2) {
                // Breadfirst Search

                // Curve Fitting

                // Normal Vector

                // Tangent Vector

                // Curvature

            }

        }
        // simple level set
        else {
            dLSdX = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
            dLSdY = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
            std::get<0>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + eps) * dLSdX;
            std::get<1>(normalVec)[idx(i, j)] = 1.0 / std::sqrt(dLSdX * dLSdX + dLSdY * dLSdY + eps) * dLSdY;
        }
    }

    return 0;
}

int MACSolver2D::UpdateKappa(const std::vector<double>& ls) {
    std::vector<double> dLSdX(kArrSize, 0.0), dLSdY(kArrSize, 0.0);
    std::vector<double> LSSize(kArrSize, 0.0);
    const double eps = 1.0e-100;
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        dLSdX[idx(i, j)] = (ls[idx(i + 1, j)] - ls[idx(i - 1, j)]) / (2.0 * kDx);
        dLSdY[idx(i, j)] = (ls[idx(i, j + 1)] - ls[idx(i, j - 1)]) / (2.0 * kDy);
        
        LSSize[idx(i, j)] = std::sqrt(dLSdX[idx(i, j)] * dLSdX[idx(i, j)] + dLSdY[idx(i, j)] * dLSdY[idx(i, j)] + eps);

        if (LSSize[idx(i, j)] < eps)
            perror("Div/0 Err in computing kappa");
    }

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        m_kappa[idx(i, j)]
            = (dLSdX[idx(i, j)] * dLSdX[idx(i, j)] * (ls[idx(i, j - 1)] - 2.0 * ls[idx(i, j)] + ls[idx(i, j + 1)]) / (kDy * kDy) // phi^2_x \phi_yy
            // - 2.0 * dLSdX[idx(i, j)] * dLSdY[idx(i, j)] * (ls[idx(i + 1, j + 1)] - ls[idx(i - 1, j + 1)] -ls[idx(i + 1, j - 1)] + ls[idx(i - 1, j - 1)]) / (4.0 * kDx * kDy) //2 \phi_x \phi_y \phi_xy
            - 2.0 * dLSdX[idx(i, j)] * dLSdY[idx(i, j)] * (dLSdX[idx(i, j + 1)] - dLSdX[idx(i, j - 1)]) / (2.0 * kDy) //2 \phi_x \phi_y \phi_xy
            + dLSdY[idx(i, j)] * dLSdY[idx(i, j)] * (ls[idx(i - 1, j)] - 2.0 * ls[idx(i, j)] + ls[idx(i + 1, j)]) / (kDx * kDx)) // phi^2_y \phi_xx
                / std::pow(LSSize[idx(i, j)], 3.0);

        // curvature is limiited so that under-resolved regions do not erroneously contribute large surface tensor forces
        m_kappa[idx(i, j)] = sign(m_kappa[idx(i, j)]) * std::fabs(std::min(std::fabs(m_kappa[idx(i, j)]), 1.0 / std::min(kDx, kDy)));
        
        assert(m_kappa[idx(i, j)] == m_kappa[idx(i, j)]);
    }

    return 0;
}

std::vector<double> MACSolver2D::UpdateFU(const std::shared_ptr<LevelSetSolver2D>& LSolver,
    const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& H) {
    
    std::vector<double> cU(kArrSize, 0.0), vU(kArrSize, 0.0), gU(kArrSize, 0.0);
    std::vector<double> rhsU(kArrSize, 0.0);

    // Convection term
    cU = this->AddConvectionFU(u, v);

    // Viscous term
    vU = this->AddViscosityFU(u, v, ls, H);

    gU = this->AddGravityFU();

    // Get RHS(Right Hand Side)
    // level set
    double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
    double theta = 0.0, iRhoEff = 0.0;
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        rhsU[idx(i, j)] = -cU[idx(i, j)] + vU[idx(i, j)] + gU[idx(i, j)];
    }
    
    return rhsU;
}

std::vector<double> MACSolver2D::UpdateFV(const std::shared_ptr<LevelSetSolver2D>& LSolver,
    const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& H) {

    std::vector<double> cV(kArrSize, 0.0), vV(kArrSize, 0.0), gV(kArrSize, 0.0);
    std::vector<double> rhsV(kArrSize, 0.0);

    // Convection term
    cV = this->AddConvectionFV(u, v);

    // Viscous term
    vV = this->AddViscosityFV(u, v, ls, H);

    gV = this->AddGravityFU();

    // Get RHS(Right Hand Side)
    // level set
    double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
    double theta = 0.0, iRhoEff = 0.0;
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
        rhsV[idx(i, j)] = -cV[idx(i, j)] + vV[idx(i, j)] + gV[idx(i, j)];
    }
    
    return rhsV;
}

std::vector<double> MACSolver2D::AddConvectionFU(const std::vector<double>& u, const std::vector<double>& v) {
    std::vector<double> cU(kArrSize, 0.0);
    std::vector<double> tmpV(kArrSize, 0.0);
    std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
    std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
    
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        tmpV[idx(i, j)] = (v[idx(i, j)] + v[idx(i - 1, j)] + v[idx(i - 1, j + 1)] + v[idx(i, j + 1)]) * 0.25;

    std::vector<double> vecF_UX(kNx + 2 * kNumBCGrid, 0.0), vecF_UY(kNy + 2 * kNumBCGrid);

    // U : X direction
    std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
            vecF_UX[i] = u[idx(i, j)] * u[idx(i, j)];
        }

        UnitHJWENO5(vecF_UX, FXP, FXM, kDx, kNx);

        for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
            LXP[idx(i, j)] = FXP[i];
            LXM[idx(i, j)] = FXM[i];
        }

        // set all vector elements to zero keeping its size
        std::fill(FXP.begin(), FXP.end(), 0.0);
        std::fill(FXM.begin(), FXM.end(), 0.0);
        std::fill(vecF_UX.begin(), vecF_UX.end(), 0.0);
    }
    
    // U : Y direction  
    std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
        for (int j = 0; j < kNy + 2 *kNumBCGrid; j++) {
            vecF_UY[j] = u[idx(i, j)] * tmpV[idx(i, j)];
        }

        UnitHJWENO5(vecF_UY, FYP, FYM, kDy, kNy);

        for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
            LYP[idx(i, j)] = FYP[j];
            LYM[idx(i, j)] = FYM[j];
        }

        // set all vector elements to zero keeping its size
        std::fill(FYP.begin(), FYP.end(), 0.0);
        std::fill(FYM.begin(), FYM.end(), 0.0);
        std::fill(vecF_UY.begin(), vecF_UY.end(), 0.0);
    }
    
    // combine together with Local Lax-Friedrichs Scheme
    double alphaX = 0.0, alphaY = 0.0, rhoEff = 0.0;

    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        alphaX = u[idx(i, j)];
        alphaY = tmpV[idx(i, j)];
        
        cU[idx(i, j)]
                = (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
                + 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
                - alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
                - alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])));
        
        if (isnan(cU[idx(i, j)]) || isinf(cU[idx(i, j)])) {
            std::cout << (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
                + 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])) << " " << -alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
                - alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])) << std::endl;
        }
        assert(cU[idx(i, j)] == cU[idx(i, j)]);
        if (std::isnan(cU[idx(i, j)]) || std::isinf(cU[idx(i, j)])) {
            std::cout << "U-convection term nan/inf error : " << i << " " << j << " " << cU[idx(i, j)] << std::endl;
            exit(1);
        }
    }

    return cU;
}

std::vector<double> MACSolver2D::AddConvectionFV(const std::vector<double>& u, const std::vector<double>& v) {
    std::vector<double> cV(kArrSize, 0.0);
    std::vector<double> tmpU(kArrSize, 0.0);
    std::vector<double> LXP(kArrSize, 0.0), LXM(kArrSize, 0.0);
    std::vector<double> LYP(kArrSize, 0.0), LYM(kArrSize, 0.0);
    
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++)
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        tmpU[idx(i, j)] = (u[idx(i, j)] + u[idx(i, j - 1)] + u[idx(i + 1, j - 1)] + u[idx(i + 1, j)]) * 0.25;

    std::vector<double> vecF_VX(kNx + 2 * kNumBCGrid, 0.0), vecF_VY(kNy + 2 * kNumBCGrid);
    
    // V : X direction
    std::vector<double> FXP(kNx + 2 * kNumBCGrid, 0.0), FXM(kNx + 2 * kNumBCGrid, 0.0);
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
            vecF_VX[i] = v[idx(i, j)] * tmpU[idx(i, j)];
        }

        UnitHJWENO5(vecF_VX, FXP, FXM, kDx, kNx);

        for (int i = 0; i < kNx + 2 * kNumBCGrid; i++) {
            LXP[idx(i, j)] = FXP[i];
            LXM[idx(i, j)] = FXM[i];
        }

        // set all vector elements to zero keeping its size
        std::fill(FXP.begin(), FXP.end(), 0.0);
        std::fill(FXM.begin(), FXM.end(), 0.0);
        std::fill(vecF_VX.begin(), vecF_VX.end(), 0.0);
    }

    // V : Y direction  
    std::vector<double> FYP(kNy + 2 * kNumBCGrid, 0.0), FYM(kNy + 2 * kNumBCGrid, 0.0);
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
        for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
            vecF_VY[j] = v[idx(i, j)] * v[idx(i, j)];
        }

        UnitHJWENO5(vecF_VY, FYP, FYM, kDy, kNy);

        for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
            LYP[idx(i, j)] = FYP[j];
            LYM[idx(i, j)] = FYM[j];
        }

        // set all vector elements to zero keeping its size
        std::fill(FYP.begin(), FYP.end(), 0.0);
        std::fill(FYM.begin(), FYM.end(), 0.0);
        std::fill(vecF_VY.begin(), vecF_VY.end(), 0.0);
    }

    // combine together with Local Lax-Friedrichs Scheme
    double alphaX = 0.0, alphaY = 0.0;

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
        alphaX = tmpU[idx(i, j)];
        alphaY = v[idx(i, j)];
        
        cV[idx(i, j)]
            = (0.5 * (LXP[idx(i, j)] + LXM[idx(i, j)])
            + 0.5 * (LYP[idx(i, j)] + LYM[idx(i, j)])
            - alphaX * (0.5 * (LXP[idx(i, j)] - LXM[idx(i, j)]))
            - alphaY * (0.5 * (LYP[idx(i, j)] - LYM[idx(i, j)])));

        assert(cV[idx(i, j)] == cV[idx(i, j)]);
        if (std::isnan(cV[idx(i, j)]) || std::isinf(cV[idx(i, j)])) {
            std::cout << "V-convection term nan/inf error : " << i << " " << j << " " << cV[idx(i, j)] << std::endl;
            exit(1);
        }
    }

    return cV;
}

int MACSolver2D::UnitHJWENO5(
    const std::vector<double> &F, std::vector<double> &FP, std::vector<double> &FM, const double d, const int n) {
    int stat = 0;
    const double kWeightEps = 1e-6;

    double alphaPlus[3], alphaMinus[3];
    double weightPlus[3], weightMinus[3];

    std::vector<double> IS0(n + 2 * kNumBCGrid, 0.0);
    std::vector<double> IS1(n + 2 * kNumBCGrid, 0.0);
    std::vector<double> IS2(n + 2 * kNumBCGrid, 0.0);
    
    // \dfrac{\Delta^+ F}{\Delta x}
    std::vector<double> DFPlus = std::vector<double>(n + 2  * kNumBCGrid, 0.0);

    // Compute Smoothness for phi^{-}
    for (int i = 0; i < n + 2 * kNumBCGrid - 1; i++)
        DFPlus[i] = F[i + 1] - F[i];

    for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
        IS0[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i - 3] - 2.0 * DFPlus[i - 2] + DFPlus[i - 1]) / d, 2.0)
            + 1.0 / 4.0
            * std::pow((DFPlus[i - 3] - 4.0 * DFPlus[i - 2] + 3.0 * DFPlus[i - 1]) / d, 2.0);
        IS1[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i - 2] - 2.0 * DFPlus[i - 1] + DFPlus[i]) / d, 2.0)
            + 1.0 / 4.0 * std::pow((DFPlus[i - 2] - DFPlus[i]) / d, 2.0);
        IS2[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i - 1] - 2.0 * DFPlus[i] + DFPlus[i + 1]) / d, 2.0)
            + 1.0 / 4.0
            * std::pow((3.0 * DFPlus[i - 1] - 4.0 * DFPlus[i] + DFPlus[i + 1]) / d, 2.0);
    }

    // phi^{-}
    for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
        alphaMinus[0] = 0.1 / std::pow((kWeightEps + IS0[i]), 2.0);
        alphaMinus[1] = 0.6 / std::pow((kWeightEps + IS1[i]), 2.0);
        alphaMinus[2] = 0.3 / std::pow((kWeightEps + IS2[i]), 2.0);

        weightMinus[0] = alphaMinus[0]
            / (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);
        weightMinus[1] = alphaMinus[1]
            / (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);
        weightMinus[2] = alphaMinus[2]
            / (alphaMinus[0] + alphaMinus[1] + alphaMinus[2]);

        FM[i] =
            weightMinus[0]
            * (1.0 / 3.0 * DFPlus[i - 3] - 7.0 / 6.0 * DFPlus[i - 2] + 11.0 / 6.0 * DFPlus[i - 1]) / d +
            weightMinus[1]
            * (-1.0 / 6.0 * DFPlus[i - 2] + 5.0 / 6.0 * DFPlus[i - 1] + 1.0 / 3.0 * DFPlus[i]) / d +
            weightMinus[2]
            * (1.0 / 3.0 * DFPlus[i - 1] + 5.0 / 6.0 * DFPlus[i] - 1.0 / 6.0 * DFPlus[i + 1]) / d;

        assert(FM[i] == FM[i]);
    }

    // Compute Smoothness for phi^{+}
    for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
        IS0[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i + 2] - 2.0 * DFPlus[i + 1] + DFPlus[i]) / d, 2.0)
            + 1.0 / 4.0
            * std::pow((DFPlus[i + 2] - 4.0 * DFPlus[i + 1] + 3.0 * DFPlus[i]) / d, 2.0);
        IS1[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i + 1] - 2.0 * DFPlus[i] + DFPlus[i - 1]) / d, 2.0)
            + 1.0 / 4.0 * std::pow((DFPlus[i + 1] - DFPlus[i - 1]) / d, 2.0);
        IS2[i]
            = 13.0 / 12.0
            * std::pow((DFPlus[i] - 2.0 * DFPlus[i - 1] + DFPlus[i - 2]) / d, 2.0)
            + 1.0 / 4.0
            * std::pow((3.0 * DFPlus[i] - 4.0 * DFPlus[i - 1] + DFPlus[i - 2]) / d, 2.0);
    }

    // phi^{+}
    for (int i = kNumBCGrid; i < n + kNumBCGrid; i++) {
        alphaPlus[0] = 0.1 / std::pow((kWeightEps + IS0[i]), 2.0);
        alphaPlus[1] = 0.6 / std::pow((kWeightEps + IS1[i]), 2.0);
        alphaPlus[2] = 0.3 / std::pow((kWeightEps + IS2[i]), 2.0);

        weightPlus[0] = alphaPlus[0]
            / (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);
        weightPlus[1] = alphaPlus[1]
            / (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);
        weightPlus[2] = alphaPlus[2]
            / (alphaPlus[0] + alphaPlus[1] + alphaPlus[2]);

        FP[i] =
            weightPlus[0]
            * (1.0 / 3.0 * DFPlus[i + 2] - 7.0 / 6.0 * DFPlus[i + 1] + 11.0 / 6.0 * DFPlus[i]) / d +
            weightPlus[1]
            * (-1.0 / 6.0 * DFPlus[i + 1] + 5.0 / 6.0 * DFPlus[i] + 1.0 / 3.0 * DFPlus[i - 1]) / d +
            weightPlus[2]
            * (1.0 / 3.0 * DFPlus[i] + 5.0 / 6.0 * DFPlus[i - 1] - 1.0 / 6.0 * DFPlus[i - 2]) / d;

        assert(FP[i] == FP[i]);
    }

    return 0;
}

std::vector<double> MACSolver2D::AddViscosityFU(const std::vector<double>& u, const std::vector<double>& v, 
    const std::vector<double>& ls, const std::vector<double>& H) {
    // This is incompressible viscous flow, which means velocity is CONTINUOUS!
    std::vector<double> dU(kArrSize, 0.0);
    std::vector<double> mu(kArrSize, 0.0);

    std::vector<double> visXVec(kArrSize, 0.0),
        visYVec(kArrSize, 0.0),
        UxxVec(kArrSize, 0.0),
        UyyVec(kArrSize, 0.0),
        VxyVec(kArrSize, 0.0),
        VxyNVec(kArrSize, 0.0),
        VxySVec(kArrSize, 0.0);

    std::vector<double> muU_W(kArrSize, 0.0),
        muU_E(kArrSize, 0.0),
        muU_S(kArrSize, 0.0),
        muU_N(kArrSize, 0.0),
        muV_W(kArrSize, 0.0),
        muV_E(kArrSize, 0.0),
        resJ11W(kArrSize, 0.0), 
        resJ11E(kArrSize, 0.0),
        resJ12S(kArrSize, 0.0), 
        resJ12N(kArrSize, 0.0),
        resJ21W(kArrSize, 0.0),
        resJ21E(kArrSize, 0.0),
        checkDiv(kArrSize, 0.0);
    
    if (kRe <= 0.0) {
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            dU[idx(i, j)] = 0.0;
        
        return dU;
    }

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

    // subcell
    double theta = 0.0;

    // need for updating rho
    /*
    -------------------------
    |           |           |
    |           |           |
    |           |           |
    ---------lsUNHalf---------
    |           |           |
    |   lsW   lsUM   lsM    |
    |           |           |
    ---------lsUSHalf--------
    |           |           |
    |           |           |
    |           |           |
    -------------------------
    */
    double lsW = 0.0, lsM = 0.0, lsUNHalf = 0.0, lsUSHalf = 0.0;
    /*
    -------------------------
    |           |           |
    |         lsUN          |
    |           |           |
    ----lsVW_N-------lsVN----
    |           |           |
  lsUW        lsUM  (i,j)  lsUE
    |           |           |
    ----lsVW---------lsVM----
    |           |           |
    |         lsUS          |
    |           |           |
    -------------------------
    */
    // u_xx, u_yy
    double lsUW = 0.0, lsUE = 0.0, lsUS = 0.0, lsUN = 0.0, lsUM = 0.0;
    // v_xy
    double lsVW_N = 0.0, lsVW = 0.0, lsVM = 0.0, lsVN = 0.0;
    double uW = 0.0, uE = 0.0, uS = 0.0, uN = 0.0, uM = 0.0;
    double vW = 0.0, vW_N = 0.0, vM = 0.0, vN = 0.0;
    double muU_X_W = 0.0, muU_X_E = 0.0, muU_Y_S = 0.0, muU_Y_N = 0.0;
    double muV_X_S = 0.0, muV_X_N = 0.0;
    double rhoU_X_W = 0.0, rhoU_X_E = 0.0, rhoU_Y_S = 0.0, rhoU_Y_N = 0.0;
    double rhoV_X_S = 0.0, rhoV_X_N = 0.0;
    double visX = 0.0, visY = 0.0;

    // jump condition
    double JUW = 0.0, JUE = 0.0, JUS = 0.0, JUN = 0.0, JUM = 0.0;
    double JVW = 0.0, JVW_N = 0.0, JVM = 0.0, JVN = 0.0;
    // effective Jump condition, effective u(uEff), effective mu (muEff), and effective rho (rhoEff)
    // J is a jump condition and defined at P grid
    double JEff = 0.0, JO = 0.0, uEff = 0.0, vEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
    const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        visX = 0.0, visY = 0.0;
        muU_X_W = 0.0; muU_X_E = 0.0; muU_Y_S = 0.0; muU_Y_N = 0.0;
        muV_X_S = 0.0; muV_X_N = 0.0;
        
        lsUW = 0.5 * (ls[idx(i - 2, j)] + ls[idx(i - 1, j)]);
        lsUE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
        lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
        lsUS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
        lsUN = 0.5 * (ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
        lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
        lsVW_N = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)]);
        lsVM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
        lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
        uM = u[idx(i, j)];
        uW = u[idx(i - 1, j)];
        uE = u[idx(i + 1, j)];
        uS = u[idx(i, j - 1)];
        uN = u[idx(i, j + 1)];
        vW = v[idx(i - 1, j)];
        vW_N = v[idx(i - 1, j + 1)];
        vM = v[idx(i, j)];
        vN = v[idx(i, j + 1)];
        
        lsW = ls[idx(i - 1, j)];
        lsM = ls[idx(i, j)];
        lsUNHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j + 1)] + ls[idx(i, j + 1)]);
        lsUSHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
        
        // V part
        /*
        -------------------------------------
        |           |           |           |
        |         lsUN          |           |
        |           |           |           |
        ----lsVW_N-------lsVN----------------
        |           |           |           |
        lsUW      lsUM  (i,j)  lsUE         |
        |           |           |           |
        ----lsVW---------lsVM----------------
        |           |           |           |
        |         lsUS        lsUE_S        |
        |           |           |           |
        -------------------------------------
        */

        /*
        -------------------------
        |           |           |
        |           |           |
        |           |           |
        ---------lsUNHalf---------
        |           |           |
        |   lsW   lsUM   lsM    |
        |           |           |
        ---------lsUSHalf--------
        |           |           |
        |           |           |
        |           |           |
        -------------------------
        */

        double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;
        if (lsW >= 0 && lsM >= 0) {
            rhoEffWE = kRhoH;
            iRhoEffWE = 1.0 / kRhoH;
        }
        else if (lsW <= 0 && lsM <= 0) {
            rhoEffWE = kRhoL;
            iRhoEffWE = 1.0 / kRhoL;
        }
        else if (lsW > 0 && lsM <= 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
            // |(lsW)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsM)|
            // |(lsW)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
            // rhoEffWE = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
            rhoEffWE = kRhoH * theta + kRhoL * (1.0 - theta);
            iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
        }
        else if (lsW <= 0 && lsM > 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
            // |(lsW)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
            // |(lsW)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
            // rhoEffWE = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
            rhoEffWE = kRhoL * theta + kRhoH * (1.0 - theta);
            iRhoEffWE = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
        }

        if (lsUSHalf >= 0 && lsUNHalf >= 0) {
            rhoEffSN = kRhoH;
            iRhoEffSN = 1.0 / kRhoH;
        }
        else if (lsUSHalf <= 0 && lsUNHalf <= 0) {
            rhoEffSN = kRhoL;
            iRhoEffSN = 1.0 / kRhoL;
        }
        else if (lsUSHalf > 0 && lsUNHalf <= 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
            // |(lsUSHalf)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsUNHalf)|
            // |(lsUSHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
            // rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
            rhoEffSN = kRhoH * theta + kRhoL * (1.0 - theta);
            iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
        }
        else if (lsUSHalf <= 0 && lsUNHalf > 0) {
            // interface lies between lsUSHalf and lsUNHalf
            theta = std::fabs(lsUSHalf) / (std::fabs(lsUSHalf) + std::fabs(lsUNHalf));
            // |(lsUSHalf)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsUNHalf)|
            // |(lsUSHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsUNHalf)|
            // rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
            rhoEffSN = kRhoL * theta + kRhoH * (1.0 - theta);
            iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
        }

        muU_X_W = mu[idx(i - 1, j)] * (uM - uW) / kDx;
        muU_X_E = mu[idx(i, j)] * (uE - uM) / kDx;
        muU_Y_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
            * (uM - uS) / kDy;
        muU_Y_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
            * (uN - uM) / kDy;
        muV_X_S = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
            * (vM - vW) / kDx;
        muV_X_N = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j + 1)] + mu[idx(i, j + 1)])
            * (vN - vW_N) / kDx;

        visX = 2.0 * (muU_X_E - muU_X_W) / kDx;
        visY = (muU_Y_N - muU_Y_S) / kDy + (muV_X_N - muV_X_S) / kDy;
        // checkDiv[idx(i, j)] = (muU_X_E - muU_X_W) / kDx + (muV_X_N - muV_X_S) / kDy;
        // dU[idx(i, j)] = visX / rhoEffWE + visY / rhoEffSN;
        dU[idx(i, j)] = visX * iRhoEffWE + visY * iRhoEffSN;
        // dU[idx(i, j)] = visX + visY;

        if (std::isnan(dU[idx(i, j)]) || std::isinf(dU[idx(i, j)])) {
            std::cout << "U-viscosity term nan/inf error : " << i << " " << j << " " << dU[idx(i, j)] << std::endl;
            exit(1);
        }
        assert(dU[idx(i, j)] == dU[idx(i, j)]);
    }

    return dU;
}

std::vector<double> MACSolver2D::AddViscosityFV(const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& ls, const std::vector<double>& H) {
    // This is incompressible viscous flow, which means velocity is CONTINUOUS!
    std::vector<double> dV(kArrSize, 0.0);
    std::vector<double> mu(kArrSize, 0.0);

    std::vector<double> visXVec(kArrSize, 0.0),
        visYVec(kArrSize, 0.0),
        visX2Vec(kArrSize, 0.0),
        visY2Vec(kArrSize, 0.0),
        UxyVec(kArrSize, 0.0),
        UxyWVec(kArrSize, 0.0),
        UxyEVec(kArrSize, 0.0),
        VxxVec(kArrSize, 0.0),
        VyyVec(kArrSize, 0.0),
        dVConst(kArrSize, 0.0);
    std::vector<double> iRhoHVec(kArrSize, 0.0),
        iRhoVVec(kArrSize, 0.0);

    std::vector<double> muV_W(kArrSize, 0.0),
        muV_E(kArrSize, 0.0),
        muV_S(kArrSize, 0.0),
        muV_N(kArrSize, 0.0),
        muU_S(kArrSize, 0.0),
        muU_N(kArrSize, 0.0),
        checkDiv(kArrSize, 0.0);

    if (kRe <= 0.0) {
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            dV[idx(i, j)] = 0.0;

        return dV;
    }

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        mu[idx(i, j)] = kMuL + (kMuH - kMuL) * H[idx(i, j)];

    // level set
    double theta = 0.0;
    /*
    -------------------------------------
    |           |           |           |
    |           |    lsM    |           |
    |           |           |           |
    ---------lsVWHalf----lsVEHalf--------
    |           |           |           |
    |           |    lsS    |           |
    |           |           |           |
    -------------------------------------
    */
    double lsVWHalf = 0.0, lsVEHalf = 0.0, lsM = 0.0, lsS = 0.0;
    /*
    ----------------lsVN-----------------
    |           |           |           |
    |         lsUM  (i,j)  lsUE         |
    |           |           |           |
    -----lsVW-------lsVM--------lsVE-----
    |           |           |           |
    |         lsUS        lsUE_S        |
    |           |           |           |
    ----------------lsVS-----------------
    */
    // need for v_xx, v_yy
    double lsVW = 0.0, lsVE = 0.0, lsVS = 0.0, lsVN = 0.0, lsVM = 0.0;
    // need for u_xy
    double lsUM = 0.0, lsUS = 0.0, lsUE = 0.0, lsUE_S = 0.0;
    // jump condition
    double JVW = 0.0, JVE = 0.0, JVS = 0.0, JVN = 0.0, JVM = 0.0;
    double JUM = 0.0, JUS = 0.0, JUE = 0.0, JUE_S = 0.0;
    double uS = 0.0, uM = 0.0, uE = 0.0, uE_S = 0.0;
    double vW = 0.0, vE = 0.0, vS = 0.0, vN = 0.0, vM = 0.0;
    double muV_X_W = 0.0, muV_X_E = 0.0, muV_Y_S = 0.0, muV_Y_N = 0.0;
    double muU_Y_W = 0.0, muU_Y_E = 0.0;
    double rhoV_X_W = 0.0, rhoV_X_E = 0.0, rhoV_Y_S = 0.0, rhoV_Y_N = 0.0;
    double rhoU_Y_W = 0.0, rhoU_Y_E = 0.0;
    double visX = 0.0, visY = 0.0;
    const double kNuI = kMuH / kRhoH, kNuO = kMuL / kRhoL;

    // effective Jump condition, effective v (vEff), effective mu (muEff), effective rho (rhoEff),
    double JEff = 0.0, JO = 0.0, vEff = 0.0, uEff = 0.0, muEff = 0.0, rhoEff = 0.0, nuEff = 0.0;
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
        muV_X_W = 0.0; muV_X_E = 0.0; muV_Y_S = 0.0; muV_Y_N = 0.0;
        muU_Y_W = 0.0; muU_Y_E = 0.0;
        visX = 0.0, visY = 0.0;
        
        lsVW = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i - 1, j)]);
        lsVE = 0.5 * (ls[idx(i + 1, j - 1)] + ls[idx(i + 1, j)]);
        lsVM = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i, j)]);
        lsVS = 0.5 * (ls[idx(i, j - 2)] + ls[idx(i, j - 1)]);
        lsVN = 0.5 * (ls[idx(i, j)] + ls[idx(i, j + 1)]);
        lsUM = 0.5 * (ls[idx(i - 1, j)] + ls[idx(i, j)]);
        lsUS = 0.5 * (ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
        lsUE = 0.5 * (ls[idx(i, j)] + ls[idx(i + 1, j)]);
        lsUE_S = 0.5 * (ls[idx(i, j - 1)] + ls[idx(i + 1, j - 1)]);
        vW = v[idx(i - 1, j)];
        vE = v[idx(i + 1, j)];
        vS = v[idx(i, j - 1)];
        vN = v[idx(i, j + 1)];
        vM = v[idx(i, j)];
        uM = u[idx(i, j)];
        uE = u[idx(i + 1, j)];
        uS = u[idx(i, j - 1)];
        uE_S = u[idx(i + 1, j - 1)];
        lsVWHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i - 1, j)] + ls[idx(i - 1, j - 1)] + ls[idx(i, j - 1)]);
        lsVEHalf = 0.25 * (ls[idx(i, j)] + ls[idx(i + 1, j)] + ls[idx(i + 1, j - 1)] + ls[idx(i, j - 1)]);
        lsM = ls[idx(i, j)];
        lsS = ls[idx(i, j - 1)];
        
        // U part
        /*
        ----------------lsVN-----------------
        |           |           |           |
        |         lsUM  (i,j)  lsUE         |
        |           |           |           |
        -----lsVW-------lsVM--------lsVE-----
        |           |           |           |
        |         lsUS        lsUE_S        |
        |           |           |           |
        ----------------lsVS-----------------
        */

        /*
        -------------------------------------
        |           |           |           |
        |           |    lsM    |           |
        |           |           |           |
        ---------lsVWHalf----lsVEHalf--------
        |           |           |           |
        |           |    lsS    |           |
        |           |           |           |
        -------------------------------------
        */
        double rhoEffWE = 0.0, rhoEffSN = 0.0, iRhoEffWE = 0.0, iRhoEffSN = 0.0;
        if (lsVWHalf >= 0 && lsVEHalf >= 0) {
            rhoEffWE = kRhoH;
            iRhoEffWE = 1.0 / kRhoH;
        }
        else if (lsVWHalf <= 0 && lsVEHalf <= 0) {
            rhoEffWE = kRhoL;
            iRhoEffWE = 1.0 / kRhoL;
        }
        else if (lsVWHalf > 0 && lsVEHalf <= 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
            // |(lsVWHalf)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsVEHalf)|
            // |(lsVWHalf)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
            // rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
            rhoEffWE = kRhoH * theta + kRhoL * (1.0 - theta);
            iRhoEffWE = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
        }
        else if (lsVWHalf <= 0 && lsVEHalf > 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsVWHalf) / (std::fabs(lsVWHalf) + std::fabs(lsVEHalf));
            // |(lsVWHalf)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsVEHalf)|
            // |(lsVWHalf)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsVEHalf)|
            // rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
            rhoEffWE = kRhoL * theta + kRhoH * (1.0 - theta);
            iRhoEffWE = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
        }

        if (lsS >= 0 && lsM >= 0) {
            rhoEffSN = kRhoH;
            iRhoEffSN = 1.0 / kRhoH;
        }
        else if (lsS <= 0 && lsM <= 0) {
            rhoEffSN = kRhoL;
            iRhoEffSN = 1.0 / kRhoL;
        }
        else if (lsS > 0 && lsM <= 0) {
            // interface lies between u[i - 1, j] and u[i, j]
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
            // |(lsS)| === inside(+)  === |(interface)| ===    outside(-)   === |(lsM)|
            // |(lsS)| === theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
            // rhoEffU = kRhoL * kRhoH / (kRhoL * theta + kRhoH * (1.0 - theta));
            rhoEffSN = kRhoH * theta + kRhoL * (1.0 - theta);
            iRhoEffSN = 1.0 / (kRhoL * kRhoH) / (1.0 / kRhoL * theta + 1.0 / kRhoH * (1.0 - theta));
        }
        else if (lsS <= 0 && lsM > 0) {
            // interface lies between lsS and lsM
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
            // |(lsS)| ===  outside(-) === |(interface)| ===    inside(+)    === |(lsM)|
            // |(lsS)| ===  theta * d  === |(interface)| === (1 - theta) * d === |(lsM)|
            // rhoEffU = kRhoH * kRhoL / (kRhoH * theta + kRhoL * (1.0 - theta));
            rhoEffSN = kRhoL * theta + kRhoH * (1.0 - theta);
            iRhoEffSN = 1.0 / (kRhoH * kRhoL) / (1.0 / kRhoH * theta + 1.0 / kRhoL * (1.0 - theta));
        }
        iRhoEffWE = 1.0 / rhoEffWE;
        iRhoEffSN = 1.0 / rhoEffSN;
        muU_Y_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
            * (uM - uS) / kDy;
        muU_Y_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
            * (uE - uE_S) / kDy;
        muV_X_W = 0.25 * (mu[idx(i, j)] + mu[idx(i - 1, j)] + mu[idx(i - 1, j - 1)] + mu[idx(i, j - 1)])
            * (vM - vW) / kDx;
        muV_X_E = 0.25 * (mu[idx(i, j)] + mu[idx(i + 1, j)] + mu[idx(i + 1, j - 1)] + mu[idx(i, j - 1)])
            * (vE - vM) / kDx;
        muV_Y_S = mu[idx(i, j - 1)] * (vM - vS) / kDy;
        muV_Y_N = mu[idx(i, j)] * (vN - vM) / kDy;

        muV_W[idx(i, j)] = muV_X_W;
        muV_E[idx(i, j)] = muV_X_E;
        muV_S[idx(i, j)] = muV_Y_S;
        muV_N[idx(i, j)] = muV_Y_N;
        muU_S[idx(i, j)] = muU_Y_W;
        muU_N[idx(i, j)] = muU_Y_E;

        visX = (muU_Y_E - muU_Y_W) / kDx + (muV_X_E - muV_X_W) / kDx;
        visY = 2.0 * (muV_Y_N - muV_Y_S) / kDy;
        checkDiv[idx(i, j)] = (muV_Y_N - muV_Y_S) / kDy + (muU_Y_E - muU_Y_W) / kDx;
        
        dV[idx(i, j)] = visX * iRhoEffWE + visY * iRhoEffSN;
        
        iRhoHVec[idx(i, j)] = iRhoEffWE;
        iRhoVVec[idx(i, j)] = iRhoEffSN;
        
        visXVec[idx(i, j)] = visX * iRhoEffWE;
        visYVec[idx(i, j)] = visY * iRhoEffSN;
        visX2Vec[idx(i, j)] = visX;
        visY2Vec[idx(i, j)] = visY;
        assert(dV[idx(i, j)] == dV[idx(i, j)]);
        if (std::isnan(dV[idx(i, j)]) || std::isinf(dV[idx(i, j)])) {
            std::cout << "V-viscosity term nan/inf error : " << i << " " << j << " " << dV[idx(i, j)] << std::endl;
            exit(1);
        }
    }
    
    return dV;  
}

std::vector<double> MACSolver2D::AddGravityFU() {
    std::vector<double> gU(kArrSize, 0.0);
    
    if ((kFr != 0 || kG != 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis == GAXISENUM2D::X) {
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            gU[idx(i, j)] = -kG;
        }
    }
    else {
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            gU[idx(i, j)] = 0.0;
        }
    }

    return gU;
}

std::vector<double> MACSolver2D::AddGravityFV() {
    std::vector<double> gV(kArrSize, 0.0);

    if ((kFr != 0 || kG != 0.0) && !isnan(kFr) && !isinf(kFr) && kGAxis == GAXISENUM2D::Y) {
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            gV[idx(i, j)] = -kG;
        }
    }
    else {
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            gV[idx(i, j)] = 0.0;
        }
    }

    return gV;
}

int MACSolver2D::GetIntermediateVel(const std::shared_ptr<LevelSetSolver2D>& LSolver,
    const std::vector<double>& ls, const std::vector<double>& u, const std::vector<double>& v,
    std::vector<double>& uhat, std::vector<double>& vhat, const std::vector<double>& H) {
        
    // Update rhs
    if (kTimeOrder == TIMEORDERENUM::EULER) {
        std::vector<double> FU1(kArrSize, 0.0);
        std::vector<double> FV1(kArrSize, 0.0);

        FU1 = UpdateFU(LSolver, ls, u, v, H);
        FV1 = UpdateFV(LSolver, ls, u, v, H);
        
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            uhat[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            vhat[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
        }
    }
    else if (kTimeOrder == TIMEORDERENUM::RK2) {
        std::vector<double> FU1(kArrSize, 0.0), FV1(kArrSize, 0.0),
            FU2(kArrSize, 0.0), FV2(kArrSize, 0.0);

        // FU1 & FV1 : L(u^(0))
        // FU2 & FV2 : u^(1)
        FU1 = UpdateFU(LSolver, ls, u, v, H);
        FV1 = UpdateFV(LSolver, ls, u, v, H);
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
        }
        std::fill(FU1.begin(), FU1.end(), 0.0);
        std::fill(FV1.begin(), FV1.end(), 0.0);
        
        // FU1 & FV1 : L(u^(1))
        // FU2 & FV2 : u^(2)
        ApplyBC_U_2D(FU2);
        ApplyBC_V_2D(FV2);
        FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
        FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
        }
        
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            uhat[idx(i, j)] = 0.5 * u[idx(i, j)] + 0.5 * FU2[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            vhat[idx(i, j)] = 0.5 * v[idx(i, j)] + 0.5 * FV2[idx(i, j)];
        }
    }
    else if (kTimeOrder == TIMEORDERENUM::RK3) {
        std::vector<double> FU1(kArrSize, 0.0), FV1(kArrSize, 0.0),
            FU2(kArrSize, 0.0), FV2(kArrSize, 0.0);

        // FU1 & FV1 : L(u^(0))
        // FU2 & FV2 : u^(1)
        FU1 = UpdateFU(LSolver, ls, u, v, H);
        FV1 = UpdateFV(LSolver, ls, u, v, H);
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU2[idx(i, j)] = u[idx(i, j)] + m_dt * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV2[idx(i, j)] = v[idx(i, j)] + m_dt * FV1[idx(i, j)];
        }
        std::fill(FU1.begin(), FU1.end(), 0.0);
        std::fill(FV1.begin(), FV1.end(), 0.0);
        
        // FU1 & FV1 : L(u^(1))
        // FU2 & FV2 : u^(1) + \delta t L (u^(1))
        ApplyBC_U_2D(FU2);
        ApplyBC_V_2D(FV2);
        FU1 = UpdateFU(LSolver, ls, FU2, FV2, H);
        FV1 = UpdateFV(LSolver, ls, FU2, FV2, H);
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU2[idx(i, j)] = FU2[idx(i, j)] + m_dt * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV2[idx(i, j)] = FV2[idx(i, j)] + m_dt * FV1[idx(i, j)];
        }
        std::fill(FU1.begin(), FU1.end(), 0.0);
        std::fill(FV1.begin(), FV1.end(), 0.0);
        
        // FU1 & FV1 : u^(2)
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU1[idx(i, j)] = 0.75 * u[idx(i, j)] + 0.25 * FU2[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV1[idx(i, j)] = 0.75 * v[idx(i, j)] + 0.25 * FV2[idx(i, j)];
        }
        std::fill(FU2.begin(), FU2.end(), 0.0);
        std::fill(FV2.begin(), FV2.end(), 0.0);

        // FU2 & FV2 : L(u^(2))
        // FU1 & FV1 : u^(2) + \delta t L (u^(2))
        ApplyBC_U_2D(FU1);
        ApplyBC_V_2D(FV1);
        FU2 = UpdateFU(LSolver, ls, FU1, FV1, H);
        FV2 = UpdateFV(LSolver, ls, FU1, FV1, H);
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            FU1[idx(i, j)] = FU1[idx(i, j)] + m_dt * FU2[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            FV1[idx(i, j)] = FV1[idx(i, j)] + m_dt * FV2[idx(i, j)];
        }
        std::fill(FU2.begin(), FU2.end(), 0.0);
        std::fill(FV2.begin(), FV2.end(), 0.0);

        // FU1 & FV1 : u^(2) + \delta t L (u^(2))
        // FU2 & FV2 : doesn't need. set to zero.
        for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
            uhat[idx(i, j)] = 1.0 / 3.0 * u[idx(i, j)] + 2.0 / 3.0 * FU1[idx(i, j)];
        }
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
            vhat[idx(i, j)] = 1.0 / 3.0 * v[idx(i, j)] + 2.0 / 3.0 * FV1[idx(i, j)];
        }
    }

    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        if (std::isnan(uhat[idx(i, j)]) || std::isinf(uhat[idx(i, j)])) {
            std::cout << "Uhat term nan/inf error : " << i << " " << j << " " << uhat[idx(i, j)] << std::endl;
            exit(1);
        }
    }

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
        if (std::isnan(vhat[idx(i, j)]) || std::isinf(vhat[idx(i, j)])) {
            std::cout << "Vhat term nan/inf error : " << i << " " << j << " " << vhat[idx(i, j)] << std::endl;
            exit(1);
        }
    }
    
    return 0;
}

int MACSolver2D::SetPoissonSolver(POISSONTYPE type) {
    m_PoissonSolverType = type;
    if (!m_Poisson)
        m_Poisson = std::make_shared<PoissonSolver2D>(kNx, kNy, kNumBCGrid);

    return 0;
}

int MACSolver2D::SolvePoisson(std::vector<double>& ps, const std::vector<double>& div,
    const std::vector<double>& ls, const std::vector<double>& lsB,
    const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& H, const int maxIter) {
    if (!m_Poisson) {
        perror("Solver method for Poisson equations are not set. Please add SetPoissonSolver Method to running code");
    }
    std::vector<double> rhs(kArrSize, 0.0);
    std::vector<double> rhoW(kArrSize, 0.0),
        rhoE(kArrSize, 0.0),
        rhoS(kArrSize, 0.0),
        rhoN(kArrSize, 0.0);
    std::vector<double> iRhoW(kArrSize, 0.0),
        iRhoE(kArrSize, 0.0),
        iRhoS(kArrSize, 0.0),
        iRhoN(kArrSize, 0.0);
    std::vector<double> FWVec(kArrSize, 0.0),
        FEVec(kArrSize, 0.0),
        FSVec(kArrSize, 0.0),
        FNVec(kArrSize, 0.0);
    std::vector<double> aWVec(kArrSize, 0.0),
        aEVec(kArrSize, 0.0),
        aCVec(kArrSize, 0.0),
        aSVec(kArrSize, 0.0),
        aNVec(kArrSize, 0.0);

    /*
    Solver \nabla \cdot ((\nabla p^*) / (\rho)) = rhs

    Liu, Xu-Dong, Ronald P. Fedkiw, and Myungjoo Kang.
    "A boundary condition capturing method for Poisson's equation on irregular domains."
    Journal of Computational Physics 160.1 (2000): 151-178.
    for level set jump condition
    [p^*] - 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N = dt \sigma \kappa
    [p^*] = 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa

    [p^*_x \ rho] = [((2 mu u_x)_x  + (mu(u_y + v_x))_y) / rho]
    [p^*_y \ rho] = [((mu(u_y + v_x))_x  + (2 mu v_y)_y  ) / rho]
    However, for solving poisson equation,
    [p^*_x \ rho] = 0.0
    [p^*_y \ rho] = 0.0
    why?
    */
    double eps = 1.0e-100;
    double lsW = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0, lsM = 0.0;
    double FW = 0.0, FE = 0.0, FS = 0.0, FN = 0.0;
    // derivative
    double duWX = 0.0, duWY = 0.0, duEX = 0.0, duEY = 0.0, duSX = 0.0, duSY = 0.0, duNX = 0.0, duNY = 0.0, duMX = 0.0, duMY = 0.0;
    double dvWX = 0.0, dvWY = 0.0, dvEX = 0.0, dvEY = 0.0, dvSX = 0.0, dvSY = 0.0, dvNX = 0.0, dvNY = 0.0, dvMX = 0.0, dvMY = 0.0;
    double dlWX = 0.0, dlWY = 0.0, dlEX = 0.0, dlEY = 0.0, dlSX = 0.0, dlSY = 0.0, dlNX = 0.0, dlNY = 0.0, dlMX = 0.0, dlMY = 0.0;
    // normal and tangent vector variable (n, t1, t2)
    double nXW = 0.0, nYW = 0.0, nXE = 0.0, nYE = 0.0, nXS = 0.0, nYS = 0.0, nXN = 0.0, nYN = 0.0, nXM = 0.0, nYM = 0.0;
    // jump at grid node (if interface is at grid node, jump occurs and aW, aE, aS, aN, aM describe that condition)
    // aW, aE, aS, aN, aM : at P grid
    double aW = 0.0, aE = 0.0, aS = 0.0, aN = 0.0, aM = 0.0;
    // jump condition [u]_\Gamma = a, [\beta u_n]_\Gamma = b
    double a = 0.0, b = 0.0;
    double thetaH = 0.0, theta = 0.0, uEff = 0.0, vEff = 0.0, aEff = 0.0, kappaEff = 0.0, rhoEff = 0.0, nXEff = 0.0, nYEff = 0.0;
    
    if (kWe != 0.0 && !isnan(kWe) && !isinf(kWe)) {
        UpdateKappa(ls);
        ApplyBC_P_2D(m_kappa);
    }
    // A Matrix is (nx * ny * nz) X (nx * ny * nz) matrix, which is very very huge. hence use sparse blas
    std::vector<double> AVals, DiagVals;
    std::vector<MKL_INT> ACols, ARowIdx, DiagCols, DiagRowIdx;
    MKL_INT rowIdx = 0, MRowIdx = 0, tmpRowIdx = 0, tmpMRowIdx = 0, colIdx = 0;
    MKL_INT Anrows = kNx * kNy, Ancols = kNx * kNy;
    MKL_INT size = kNx * kNy;
    
    double dUdXW = 0.0, dUdXE = 0.0, dUdYS = 0.0, dUdYN = 0.0;
    double dVdXW = 0.0, dVdXE = 0.0, dVdYS = 0.0, dVdYN = 0.0;
    // stored coef for A matrix, Dictionary but it is ordered
    std::map<std::string, double> AValsDic;
    std::map<std::string, MKL_INT> AColsDic;

    for (int i = 0; i < kNx + 2 * kNumBCGrid; i++)
    for (int j = 0; j < kNy + 2 * kNumBCGrid; j++) {
        ps[idx(i, j)] = 0.0;
        rhs[idx(i, j)] = 0.0;
    }
    
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        lsW = ls[idx(i - 1, j)];
        lsE = ls[idx(i + 1, j)];
        lsM = ls[idx(i, j)];
        lsS = ls[idx(i, j - 1)];
        lsN = ls[idx(i, j + 1)];

        FW = 0.0;
        FE = 0.0;
        FS = 0.0;
        FN = 0.0;

        // thetaH = portion of kRhoH
        // theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsW >= 0.0 && lsM >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsW < 0.0 && lsM < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsW >= 0.0 && lsM < 0.0) {
            thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
        }
        else if (lsW < 0.0 && lsM >= 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
        }
        
        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        // coefficient
        iRhoW[idx(i, j)] = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
        
        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // kRhoH - Liquid, kRhoL - Gas
        // - level set : H = 1, + level set : H = 0;
        // Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
        FW = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
        FW *= iRhoW[idx(i, j)] / (kDx * kDx);
        aWVec[idx(i, j)] = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
        
        // thetaH = portion of kRhoH
        // theta = portion of fluid cell adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsM >= 0.0 && lsE >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsM < 0.0 && lsE < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsM >= 0.0 && lsE < 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsM) + std::fabs(lsE));
            theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
        }
        else if (lsM < 0.0 && lsE >= 0.0) {
            thetaH = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
            theta = std::fabs(lsE) / (std::fabs(lsM) + std::fabs(lsE));
        }

        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        iRhoE[idx(i, j)] = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i + 1, j)] * (1.0 - theta);

        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // - level set : H = 1, + level set : H = 0; 
        // Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
        FE = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);
        FE *= iRhoE[idx(i, j)] / (kDx * kDx);
        aEVec[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i + 1, j)] - H[idx(i, j)]);
        // thetaH = portion of kRhoH
        // theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsS >= 0.0 && lsM >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsS < 0.0 && lsM < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsS >= 0.0 && lsM < 0.0) {
            thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
        }
        else if (lsS < 0.0 && lsM >= 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
        }

        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        iRhoS[idx(i, j)] = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
        
        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // kRhoH - Liquid, kRhoL - Gas
        // - level set : H = 1, + level set : H = 0; 
        // Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
        FS = -m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
        FS *= iRhoS[idx(i, j)] / (kDy * kDy);
        
        // thetaH = portion of kRhoH
        // theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsM >= 0.0 && lsN >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsM < 0.0 && lsN < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsM >= 0.0 && lsN < 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsM) + std::fabs(lsN));
            theta = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
        }
        else if (lsM < 0.0 && lsN >= 0.0) {
            thetaH = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
            theta = std::fabs(lsN) / (std::fabs(lsM) + std::fabs(lsN));
        }

        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        iRhoN[idx(i, j)] = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j + 1)] * (1.0 - theta);
        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // kRhoH - Liquid, kRhoL - Gas
        // - level set : H = 1, + level set : H = 0; 
        // Bubble : (-ls)-kRhoL-Gas, (+ls)-kRhoH-Liquid
        FN = m_dt * (-kSigma * kappaEff) * (H[idx(i, j + 1)] - H[idx(i, j)]);
        FN *= iRhoN[idx(i, j)] / (kDy * kDy);
        
        // poisson equation form should be -\beta \nabla p = f
        // iRhoEff has already negative value of rhoEff, then it is just a coefficient.
        // For discretization of pressure gradient, additional term is negative and it goes to RHS
        // Then it is a positive value and don't worrry about the sign
        rhs[idx(i, j)] -= FW + FE + FS + FN;
        FWVec[idx(i, j)] = FW;
        FEVec[idx(i, j)] = FE;
        FSVec[idx(i, j)] = FS;
        FNVec[idx(i, j)] = FN;
        
        assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
        if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
            std::cout << "right hand side of poisson equation nan/inf error : " << i << " " << j << " " 
                << rhs[idx(i, j)] << std::endl;
            exit(1);
        }
    }
    
    // Original value of RHS
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        rhs[idx(i, j)] -= div[idx(i, j)];

    // An order of A matrix coef. is very important, hence reverse j order
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++) {
        // AValsDic and AColsDic contain coef.s of A matrix of each row (a row of A matrix = a each point of 3D grid)
        // Based on boundary condition, a coef.s of A matrix change and are saved in AvalsDic and AColsDic
        // At last, traverse AValsDic and AColsDic, add nonzero value of A matrix
        AValsDic.clear();
        AColsDic.clear();
        tmpRowIdx = 0;
        tmpMRowIdx = 0;
        if ((m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM || m_BC->m_BC_PW == BC2D::WALL) &&
            (m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM || m_BC->m_BC_PE == BC2D::WALL) &&
            (m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM || m_BC->m_BC_PS == BC2D::WALL) &&
            (m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM || m_BC->m_BC_PN == BC2D::WALL) &&
            (i == kNumBCGrid && j == kNumBCGrid))
            continue;
        // Add starting rowIdx
        ARowIdx.push_back(rowIdx);
        DiagRowIdx.push_back(MRowIdx);

        iRhoS[idx(i, j)] *= -1.0;
        iRhoW[idx(i, j)] *= -1.0;
        iRhoE[idx(i, j)] *= -1.0;
        iRhoN[idx(i, j)] *= -1.0;

        // Set default values, if a current pointer is in interior, it will not be changed.
        AValsDic["S"] = iRhoS[idx(i, j)] / (kDy * kDy);
        AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx;
        AValsDic["W"] = iRhoW[idx(i, j)] / (kDx * kDx);
        AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
        AValsDic["C"] = -(iRhoW[idx(i, j)] + iRhoE[idx(i, j)]) / (kDx * kDx)
            - (iRhoS[idx(i, j)] + iRhoN[idx(i, j)]) / (kDy * kDy);
        AColsDic["C"] = (i - kNumBCGrid) + (j - kNumBCGrid) * kNx;
        AValsDic["E"] = iRhoE[idx(i, j)] / (kDx * kDx);
        AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx;
        AValsDic["N"] = iRhoN[idx(i, j)] / (kDy * kDy);
        AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx;
        
        if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::NEUMANN) {
            AColsDic["W"] = -1;
            AValsDic["W"] = 0.0;
            AValsDic["C"] += iRhoW[idx(i, j)] / (kDx * kDx);
        }
        else if (i == kNumBCGrid
            && (m_BC->m_BC_PW == BC2D::DIRICHLET || m_BC->m_BC_PW == BC2D::INLET
            || m_BC->m_BC_PW == BC2D::OUTLET || m_BC->m_BC_PW == BC2D::PRESSURE
            || m_BC->m_BC_PW == BC2D::WALL)) {
            AColsDic["W"] = -1;
            AValsDic["W"] = 0.0;
            rhs[idx(i, j)] -= iRhoW[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPW);
        }
        else if (i == kNumBCGrid && m_BC->m_BC_PW == BC2D::PERIODIC) {
            AValsDic["W"] = iRhoW[idx(kNumBCGrid + kNx - 1, j)];
        }

        // East boundary
        if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC2D::NEUMANN) {
            AColsDic["E"] = -1;
            AValsDic["E"] = 0.0;
            AValsDic["C"] += iRhoE[idx(i, j)] / (kDx * kDx);
        }
        else if (i == kNumBCGrid + kNx - 1
            && (m_BC->m_BC_PE == BC2D::DIRICHLET || m_BC->m_BC_PE == BC2D::INLET
            || m_BC->m_BC_PE == BC2D::OUTLET || m_BC->m_BC_PE == BC2D::PRESSURE
            || m_BC->m_BC_PE == BC2D::WALL)) {
            AColsDic["E"] = -1;
            AValsDic["E"] = 0.0;
            rhs[idx(i, j)] -= iRhoE[idx(i, j)] / (kDx * kDx) * (m_BC->m_BC_DirichletConstantPE);
        }
        else if (i == kNumBCGrid + kNx - 1 && m_BC->m_BC_PE == BC2D::PERIODIC) {
            AValsDic["E"] = iRhoE[idx(kNumBCGrid, j)];
        }

        if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::NEUMANN) {
            AColsDic["S"] = -1;
            AValsDic["S"] = 0.0;
            AValsDic["C"] += iRhoS[idx(i, j)] / (kDy * kDy);
        }
        else if (j == kNumBCGrid
            && (m_BC->m_BC_PS == BC2D::DIRICHLET || m_BC->m_BC_PS == BC2D::INLET
            || m_BC->m_BC_PS == BC2D::OUTLET || m_BC->m_BC_PS == BC2D::PRESSURE
            || m_BC->m_BC_PS == BC2D::WALL)) {
            AColsDic["S"] = -1;
            AValsDic["S"] = 0.0;
            rhs[idx(i, j)] -= iRhoS[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPS);
        }
        else if (j == kNumBCGrid && m_BC->m_BC_PS == BC2D::PERIODIC) {
            AValsDic["S"] = iRhoS[idx(i, kNumBCGrid + kNy - 1)];
        }

        if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC2D::NEUMANN) {
            AColsDic["N"] = -1;
            AValsDic["N"] = 0.0;
            AValsDic["C"] += iRhoN[idx(i, j)] / (kDy * kDy);
        }
        else if (j == kNumBCGrid + kNy - 1
            && (m_BC->m_BC_PN == BC2D::DIRICHLET || m_BC->m_BC_PN == BC2D::INLET
            || m_BC->m_BC_PN == BC2D::OUTLET || m_BC->m_BC_PN == BC2D::PRESSURE
            || m_BC->m_BC_PN == BC2D::WALL)) {
            AColsDic["N"] = -1;
            AValsDic["N"] = 0.0;
            rhs[idx(i, j)] -= iRhoN[idx(i, j)] / (kDy * kDy) * (m_BC->m_BC_DirichletConstantPN);
        }
        else if (j == kNumBCGrid + kNy - 1 && m_BC->m_BC_PN == BC2D::PERIODIC) {
            AValsDic["N"] = iRhoN[idx(i, kNumBCGrid + kNy + 1)];
        }

        // fix ps[idx(kNumBCGrid, kNumBCGrid)] = 0.0 for neumann condition of all boundary
        if ((m_BC->m_BC_PW == BC2D::NEUMANN || m_BC->m_BC_PW == BC2D::AXISYM || m_BC->m_BC_PW == BC2D::WALL) &&
            (m_BC->m_BC_PE == BC2D::NEUMANN || m_BC->m_BC_PE == BC2D::AXISYM || m_BC->m_BC_PE == BC2D::WALL) &&
            (m_BC->m_BC_PS == BC2D::NEUMANN || m_BC->m_BC_PS == BC2D::AXISYM || m_BC->m_BC_PS == BC2D::WALL) &&
            (m_BC->m_BC_PN == BC2D::NEUMANN || m_BC->m_BC_PN == BC2D::AXISYM || m_BC->m_BC_PN == BC2D::WALL)) {
            if (i == kNumBCGrid + 1 && j == kNumBCGrid) {
                AColsDic["W"] = -1;
                AValsDic["W"] = 0.0;
                rhs[idx(i, j)] -= iRhoW[idx(i, j)] / (kDx * kDx) * (-ps[idx(i, j)]);
            }
            else if (i == kNumBCGrid && j == kNumBCGrid + 1) {
                AColsDic["S"] = -1;
                AValsDic["S"] = 0.0;
                rhs[idx(i, j)] -= iRhoS[idx(i, j)] / (kDy * kDy) * (-ps[idx(i, j)]);
            }

            AColsDic["S"] = (i - kNumBCGrid) + (j - 1 - kNumBCGrid) * kNx - 1;
            AColsDic["W"] = (i - 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx - 1;
            AColsDic["E"] = (i + 1 - kNumBCGrid) + (j - kNumBCGrid) * kNx - 1;
            AColsDic["N"] = (i - kNumBCGrid) + (j + 1 - kNumBCGrid) * kNx - 1;

            size = kNx * kNy - 1;
        }

        // add non zero values to AVals and ACols
        // KEEP ORDER OF PUSH_BACK!!

        if (AColsDic["S"] >= 0) {
            tmpRowIdx++;
            AVals.push_back(AValsDic["S"]);
            ACols.push_back(AColsDic["S"]);
        }
        
        if (AColsDic["W"] >= 0) {
            tmpRowIdx++;
            AVals.push_back(AValsDic["W"]);
            ACols.push_back(AColsDic["W"]);
        }

        // Center, it doens't neeed to check
        tmpRowIdx++;
        AVals.push_back(AValsDic["C"]);
        ACols.push_back(AColsDic["C"]);

        if (AColsDic["E"] >= 0) {
            tmpRowIdx++;
            AVals.push_back(AValsDic["E"]);
            ACols.push_back(AColsDic["E"]);
        }

        if (AColsDic["N"] >= 0) {
            tmpRowIdx++;
            AVals.push_back(AValsDic["N"]);
            ACols.push_back(AColsDic["N"]);
        }
        
        tmpMRowIdx++;
        DiagVals.push_back(AValsDic["C"]);
        DiagCols.push_back(AColsDic["C"]);
        
        rowIdx += tmpRowIdx;
        MRowIdx += tmpMRowIdx;
        
        assert(rhs[idx(i, j)] == rhs[idx(i, j)]);
        if (std::isnan(rhs[idx(i, j)]) || std::isinf(rhs[idx(i, j)])) {
            std::cout << "right hand side of poisson equation nan/inf error : " 
                << i << " " << j << " " << rhs[idx(i, j)] << std::endl;
            exit(1);
        }
    }
    ARowIdx.push_back(rowIdx);
    DiagRowIdx.push_back(MRowIdx);
    
    if (m_PoissonSolverType == POISSONTYPE::MKL) {
        // For Level set solver test only, legacy poisson equation solver
        /*
        Solver \nabla \cdot ((\nabla p^*) / (\rho)) = rhs

        Liu, Xu-Dong, Ronald P. Fedkiw, and Myungjoo Kang.
        "A boundary condition capturing method for Poisson's equation on irregular domains."
        Journal of Computational Physics 160.1 (2000): 151-178.
        for level set jump condition
        [p^*] - 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N = dt \sigma \kappa
        [p^*] = 2 dt [mu] (\del u \cdot n, \del v \cdot n) \cdot N + dt \sigma \kappa

        [p^*_x \ rho] = [((2 mu u_x)_x  + (mu(u_y + v_x))_y) / rho]
        [p^*_y \ rho] = [((mu(u_y + v_x))_x  + (2 mu v_y)_y  ) / rho]
        However, for solving poisson equation, 
        [p^*_x \ rho] = 0.0
        [p^*_y \ rho] = 0.0
        why? 
        */

        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
            rhs[idx(i, j)] = -div[idx(i, j)];
        
        m_Poisson->MKL_2FUniform_2D(ps, rhs, kLenX, kLenY, kDx, kDy, m_BC);
        
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
            if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
                std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;

    }
    else if (m_PoissonSolverType == POISSONTYPE::CG) {
        // std::cout << "Poisson : CG" << std::endl;
        m_Poisson->CG_2FUniformP_2D(ps, rhs, AVals, ACols, ARowIdx,
            DiagVals, DiagCols, DiagRowIdx, m_BC, size, maxIter);
    
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
                std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
    }
    else if (m_PoissonSolverType == POISSONTYPE::BICGSTAB) {
        // std::cout << "Poisson : BiCG" << std::endl;
        m_Poisson->BiCGStab_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx,
            DiagVals, DiagCols, DiagRowIdx, m_BC, size, maxIter);

        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
                std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
    }
    else if (m_PoissonSolverType == POISSONTYPE::GS) {
        // std::cout << "Poisson : GS" << std::endl;
        m_Poisson->GS_2FUniform_2D(ps, rhs, AVals, ACols, ARowIdx, kLenX, kLenY, kDx, kDy, m_BC, maxIter);

        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            if (std::isnan(ps[idx(i, j)]) || std::isinf(ps[idx(i, j)]))
                std::cout << "Pseudo-p nan/inf error : " << i << " " << j << " " << ps[idx(i, j)] << std::endl;
    }
    
    ApplyBC_P_2D(ps);

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        m_p[idx(i, j)] = ps[idx(i, j)] / m_dt;
    
    return 0;
}

std::vector<double> MACSolver2D::GetDivergence(const std::vector<double>& u, const std::vector<double>& v) {
    std::vector<double> div(kArrSize, 0.0);

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        div[idx(i, j)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
                        + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
    
    return div;
}

int MACSolver2D::UpdateVel(std::vector<double>& u, std::vector<double>& v,
    const std::vector<double>& us, const std::vector<double>& vs,
    const std::vector<double>& ps, const std::vector<double>& ls, const std::vector<double>& lsB,
    const std::vector<double>& H) {
    
    // velocity update after solving poisson equation
    // ps = p * dt
    double lsW = 0.0, lsM = 0.0, lsE = 0.0, lsS = 0.0, lsN = 0.0;
    double uEff = 0.0, vEff = 0.0, rhoEff = 0.0, theta = 0.0, thetaH = 0.0, iRhoEff = 0.0;
    double nXEff = 0.0, nYEff = 0.0, kappaEff = 0.0;
    double nXW = 0.0, nYW = 0.0, nXS = 0.0, nYS = 0.0, nXM = 0.0, nYM = 0.0;
    double aW = 0.0, aS = 0.0, aM = 0.0, aEff = 0.0;
    const double eps = 1.0e-100;
    
    std::vector<double> AWW(kArrSize, 0.0), AWM(kArrSize, 0.0),
        ASS(kArrSize, 0.0), ASM(kArrSize, 0.0);
    std::vector<double> iRhoEffWVec(kArrSize, 0.0), iRhoEffSVec(kArrSize, 0.0);

    std::vector<double> aWEff(kArrSize, 0.0),
        aSEff(kArrSize, 0.0);
    
    /*
    // http://ctr.stanford.edu/Summer/SP08/3_1_Moureau.pdf
    // Moreau, V., and O. Desjardins. 
    // "A second-order ghost-fluid method for the primary atomization
    //   of liquid fuel in air-blast type injectors."
    // Proceedings of the Summer Program. Vol. 143. 2008.
    */
    for (int i = kNumBCGrid + 1; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        lsW = ls[idx(i - 1, j)];
        lsM = ls[idx(i, j)];
        
        // thetaH = portion of kRhoH
        // theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsW >= 0.0 && lsM >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsW < 0.0 && lsM < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsW >= 0.0 && lsM < 0.0) {
            thetaH = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
        }
        else if (lsW < 0.0 && lsM >= 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsW) + std::fabs(lsM));
            theta = std::fabs(lsW) / (std::fabs(lsW) + std::fabs(lsM));
        }

        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        iRhoEff = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i - 1, j)] * (1.0 - theta);
        
        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // + level set : H = 1, - level set : H = 0; 
        // Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid
        u[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i - 1, j)]);
        u[idx(i, j)] *= iRhoEff / kDx;

        aWEff[idx(i, j)] = u[idx(i, j)];
        iRhoEffWVec[idx(i, j)] = iRhoEff;
        u[idx(i, j)] += us[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i - 1, j)]) / kDx;
    }

    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid + 1; j < kNy + kNumBCGrid; j++) {
        lsM = ls[idx(i, j)];
        lsS = ls[idx(i, j - 1)];
        aS = 0.0;
        aM = 0.0;

        // thetaH = portion of kRhoH
        // theta = portion of fluid adjacent to lsM, such as |lsW| / (|lsW| + |lsM|), |lsE| / (|lsWE| + |lsM|), and so on
        if (lsS >= 0.0 && lsM >= 0.0)  {
            thetaH = 1.0; theta = 0.0;
        }
        else if (lsS < 0.0 && lsM < 0.0) {
            thetaH = 0.0; theta = 0.0;
        }
        else if (lsS >= 0.0 && lsM < 0.0) {
            thetaH = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
        }
        else if (lsS < 0.0 && lsM >= 0.0) {
            thetaH = std::fabs(lsM) / (std::fabs(lsS) + std::fabs(lsM));
            theta = std::fabs(lsS) / (std::fabs(lsS) + std::fabs(lsM));
        }

        rhoEff = kRhoH * thetaH + kRhoL * (1.0 - thetaH);
        iRhoEff = 1.0 / rhoEff;
        kappaEff = m_kappa[idx(i, j)] * theta + m_kappa[idx(i, j - 1)] * (1.0 - theta);
        
        // pressure jump = Liquid - Gas = -sigma * kappa - [\rho] * normal * gvector
        // + level set : H = 1, - level set : H = 0; 
        // Bubble : (+ls)-kRhoL-Gas, (-ls)-kRhoH-Liquid

        v[idx(i, j)] = m_dt * (-kSigma * kappaEff) * (H[idx(i, j)] - H[idx(i, j - 1)]);
        v[idx(i, j)] *= iRhoEff / kDy;

        aSEff[idx(i, j)] = v[idx(i, j)];
        iRhoEffSVec[idx(i, j)] = iRhoEff;
        v[idx(i, j)] += vs[idx(i, j)] - iRhoEff * (ps[idx(i, j)] - ps[idx(i, j - 1)]) / kDy;
    }
    
    ApplyBC_U_2D(u);
    ApplyBC_V_2D(v);
    
    return 0;
}

double MACSolver2D::UpdateDt(const std::vector<double>& u, const std::vector<double>& v) {
    double uAMax = 0.0;
    double vAMax = 0.0;
    double Cefl = 0.0, Vefl = 0.0, Gefl = 0.0, Sefl = 0.0;
    double dt = std::numeric_limits<double>::max();

    // get maximum of absolute value
    for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
    for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++) {
        uAMax = std::max(uAMax, std::fabs((u[idx(i + 1, j)] * u[idx(i, j)]) * 0.5));
        vAMax = std::max(vAMax, std::fabs((v[idx(i, j + 1)] * v[idx(i, j)]) * 0.5));
    }
    
    Cefl = uAMax / kDx + vAMax / kDy;
    Vefl = std::max(kMuH, kMuL) * (2.0 / (kDx * kDx + kDy * kDy));
    Gefl = std::max(std::sqrt(std::fabs(kG) / kDy), std::sqrt(std::fabs(kG) / kDx));
    Sefl = std::sqrt((kSigma / std::min(kDx, kDy))
        / (std::min(kRhoH, kRhoL) * std::pow(std::min(kDx, kDy), 2.0)));
    
    dt = std::min(dt,
        kCFL / (0.5 * (Cefl + Vefl +
        std::sqrt(std::pow(Cefl + Vefl, 2.0) +
        4.0 * Gefl * Gefl + 4.0 * Sefl * Sefl))));

    if (std::isnan(dt) || std::isinf(dt)) {
        std::cout << "dt nan/inf error : Cefl, Vefl, Gefl : " << Cefl << " " << Vefl << " " << Gefl << std::endl;
    }

    return dt;
}

int MACSolver2D::SetBC_U_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
    if (!m_BC) {
        m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
    }

    m_BC->SetBC_U_2D(BC_W, BC_E, BC_S, BC_N);

    return 0;
}

int MACSolver2D::SetBC_V_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
    if (!m_BC) {
        m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
    }

    m_BC->SetBC_V_2D(BC_W, BC_E, BC_S, BC_N);

    return 0;
}

int MACSolver2D::SetBC_P_2D(std::string BC_W, std::string BC_E, std::string BC_S, std::string BC_N) {
    if (!m_BC) {
        m_BC = std::make_shared<BoundaryCondition2D>(kNx, kNy, kNumBCGrid);
    }

    m_BC->SetBC_P_2D(BC_W, BC_E, BC_S, BC_N);

    return 0;
}

int MACSolver2D::ApplyBC_U_2D(std::vector<double>& arr) {
    m_BC->ApplyBC_U_2D(arr);
    return 0;
}

int MACSolver2D::ApplyBC_V_2D(std::vector<double>& arr) {
    m_BC->ApplyBC_V_2D(arr);
    return 0;
}

int MACSolver2D::ApplyBC_P_2D(std::vector<double>& arr) {
    m_BC->ApplyBC_P_2D(arr);
    return 0;
}

void MACSolver2D::SetBCConstantUW(double BC_ConstantW) {
    return m_BC->SetBCConstantUW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantUE(double BC_ConstantE) {
    return m_BC->SetBCConstantUE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantUS(double BC_ConstantS) {
    return m_BC->SetBCConstantUS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantUN(double BC_ConstantN) {
    return m_BC->SetBCConstantUN(BC_ConstantN);
}

void MACSolver2D::SetBCConstantVW(double BC_ConstantW) {
    return m_BC->SetBCConstantVW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantVE(double BC_ConstantE) {
    return m_BC->SetBCConstantVE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantVS(double BC_ConstantS) {
    return m_BC->SetBCConstantVS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantVN(double BC_ConstantN) {
    return m_BC->SetBCConstantVN(BC_ConstantN);
}

void MACSolver2D::SetBCConstantPW(double BC_ConstantW) {
    return m_BC->SetBCConstantPW(BC_ConstantW);
}

void MACSolver2D::SetBCConstantPE(double BC_ConstantE) {
    return m_BC->SetBCConstantPE(BC_ConstantE);
}

void MACSolver2D::SetBCConstantPS(double BC_ConstantS) {
    return m_BC->SetBCConstantPS(BC_ConstantS);
}

void MACSolver2D::SetBCConstantPN(double BC_ConstantN) {
    return m_BC->SetBCConstantPN(BC_ConstantN);
}

// http://stackoverflow.com/questions/11990030/c-sign-function-from-matlab
inline int MACSolver2D::sign(const double& val) {
    return (val > 0) - (val < 0);
}

inline int MACSolver2D::idx(int i, int j) {
    return (j + (kNy + 2 * kNumBCGrid) * (i));
}

int MACSolver2D::SetPLTType(PLTTYPE type) {
    m_PLTType = type;

    return 0;
}

int MACSolver2D::OutRes(const int iter, const double curTime, const std::string fname_vel_base, const std::string fname_div_base,
    const std::vector<double>& u, const std::vector<double>& v,
    const std::vector<double>& ps, const std::vector<double>& ls) {
    if (m_PLTType == PLTTYPE::ASCII || m_PLTType == PLTTYPE::BOTH) {
        std::ofstream outF;
        std::string fname_vel(fname_vel_base + "_ASCII.plt");
        std::string fname_div(fname_div_base + "_ASCII.plt");
        if (m_iter == 0) {
            outF.open(fname_vel.c_str(), std::ios::out);

            outF << "TITLE = VEL" << std::endl;
            outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"PS\" " << std::endl;
            outF.close();

            outF.open(fname_div.c_str(), std::ios::out);
            outF << "TITLE = DIV" << std::endl;
            outF << "VARIABLES = \"X\", \"Y\", \"LS\", \"DIV\", \"PS\" " << std::endl;
            outF.close();
        }

        std::vector<double>
            resU(kArrSize, 0.0),
            resV(kArrSize, 0.0);

        std::vector<double> resDiv = GetDivergence(u, v);
        
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            resU[idx(i, j)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;

        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
            resV[idx(i, j)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
        
        outF.open(fname_vel.c_str(), std::ios::app);
        
        outF << std::string("ZONE T=\"") << iter
            << std::string("\", I=") << kNx << std::string(", J=") << kNy
            << std::string(", SOLUTIONTIME=") << curTime
            << std::string(", STRANDID=") << iter + 1
            << std::endl;

        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
            outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
                << kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
                << static_cast<double>(resU[idx(i, j)]) << std::string(",")
                << static_cast<double>(resV[idx(i, j)]) << std::string(",")
                << static_cast<double>(ls[idx(i, j)]) << std::string(",")
                << static_cast<double>(ps[idx(i, j)]) << std::endl;

        outF.close();
        
        outF.open(fname_div.c_str(), std::ios::app);

        outF << std::string("ZONE T=\"") << iter
            << std::string("\", I=") << kNx << std::string(", J=") << kNy
            << std::string(", SOLUTIONTIME=") << curTime
            << std::string(", STRANDID=") << iter + 1
            << std::endl;
        
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)
        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
            outF << kBaseX + static_cast<double>(i + 0.5 - kNumBCGrid) * kDx << std::string(",")
                << kBaseY + static_cast<double>(j + 0.5 - kNumBCGrid) * kDy << std::string(",")
                << static_cast<double>(ls[idx(i, j)]) << std::string(",")
                << static_cast<double>(resDiv[idx(i, j)]) << std::string(",")
                // << static_cast<double>(m_kappa[idx(i, j)]) << std::string(",")
                << static_cast<double>(ps[idx(i, j)]) << std::endl;

        outF.close();
    }

    if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH) {
        INTEGER4 whichFile = 0, stat = 0;
        std::string fname_vel(fname_vel_base + "_BINARY.szplt");
        std::string fname_div(fname_div_base + "_BINARY.szplt");

        std::vector<double>
            resX(kNx * kNy, 0.0),
            resY(kNx * kNy, 0.0),
            resU(kNx * kNy, 0.0),
            resV(kNx * kNy, 0.0),
            resDiv(kNx * kNy, 0.0),
            resLS(kNx * kNy, 0.0),
            resPs(kNx * kNy, 0.0);

        for (int i = kNumBCGrid; i < kNx + kNumBCGrid; i++)
        for (int j = kNumBCGrid; j < kNy + kNumBCGrid; j++)  {
            resX[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseX + (i + 0.5 - kNumBCGrid) * kDx;
            resY[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = kBaseY + (j + 0.5 - kNumBCGrid) * kDy;
            resU[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i, j)] + u[idx(i + 1, j)]) * 0.5;
            resV[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (v[idx(i, j)] + v[idx(i, j + 1)]) * 0.5;
            resLS[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ls[idx(i, j)];
            resPs[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = ps[idx(i, j)];
            resDiv[(i - kNumBCGrid) + kNx * (j - kNumBCGrid)] = (u[idx(i + 1, j)] - u[idx(i, j)]) / kDx
             + (v[idx(i, j + 1)] - v[idx(i, j)]) / kDy;
        }

        INTEGER4 Debug = 1;
        INTEGER4 VIsDouble = 1; /* 0 = Single precision, 1 = Double precision*/
        INTEGER4 FileType = 0; /* 0 = Full, 1 = Grid only, 2 = Solution only*/
        INTEGER4 FileFormat = 1; /* 0 = plt, 1 = szplt*/

        if (m_iter == 0) {
            /*
            * Open the file and write the tecplot datafile
            * header information
            */
            stat = TECINI142(const_cast<char *>(std::string("VELOCITY").c_str()),  /* Name of the entire dataset.  */
                            const_cast<char *>(std::string("X, Y, U, V, LS, PS").c_str()),  
                            /* Defines the variables for the data file. Each zone must contain each of the variables listed here. 
                            * The order of the variables in the list is used to define the variable number (e.g. X is Var 1).*/
                            const_cast<char *>(fname_vel.c_str()),
                            const_cast<char *>(std::string(".").c_str()),      /* Scratch Directory */
                                &FileFormat, &FileType, &Debug, &VIsDouble);
        }
        else {
            whichFile = 1;
            stat = TECFIL142(&whichFile);
        }

        /* Set the parameters for TecZne */
        INTEGER4 ZoneType = 0; /* sets the zone type to
                               * ordered
                               */
        /* Create an IJ-ordered zone, by using IMax and JMax
        * values that are greater than one, and setting KMax to one.
        */
        INTEGER4 IMax = kNx, JMax = kNy, KMax = 1;

        double   SolTime = curTime;
        INTEGER4 StrandID = iter + 1; /* StaticZone */
        INTEGER4 ParentZn = 0; /* used for surface streams */

        INTEGER4 ICellMax = 0, JCellMax = 0, KCellMax = 0; /* not used */

        INTEGER4 IsBlock = 1; /* Block */

        INTEGER4 NFConns = 0; /* this example does not use face neighbors */
        INTEGER4 FNMode = 0;
        INTEGER4 TotalNumFaceNodes = 1, TotalNumBndryFaces = 1, TotalNumBndryConn = 1;
        INTEGER4 ShrConn = 0;

        /* Create an Ordered Zone */
        stat = TECZNE142((char*) std::to_string(StrandID).c_str(), &ZoneType,
            &IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax,
            &SolTime, &StrandID, &ParentZn, &IsBlock,
            &NFConns, &FNMode, &TotalNumFaceNodes,
            &TotalNumBndryFaces, &TotalNumBndryConn,
            NULL, NULL, NULL, &ShrConn);
            
        INTEGER4 DIsDouble = 1;  /* set DIsDouble to 0, for float
                                 * values.
                                 */

        INTEGER4 ARRSIZEVAL = IMax * JMax * KMax;

        stat = TECDAT142(&ARRSIZEVAL, resX.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resY.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resU.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resV.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);

        if (m_iter == 0) {
            stat = TECINI142(const_cast<char *>(std::string("DIVERGENCE").c_str()),  /* Name of the entire dataset.  */
                const_cast<char *>(std::string("X, Y, LS, DIV, PS").c_str()),
                /* Defines the variables for the data file. Each zone must contain each of the variables listed here.
                * The order of the variables in the list is used to define the variable number (e.g. X is Var 1).*/
                const_cast<char *>(fname_div.c_str()),
                const_cast<char *>(std::string(".").c_str()),      /* Scratch Directory */
                &FileFormat, &FileType, &Debug, &VIsDouble);
        }

        whichFile = 2;
        stat = TECFIL142(&whichFile);

        /* Create an Ordered Zone */
        stat = TECZNE142((char*) std::to_string(StrandID).c_str(), &ZoneType,
            &IMax, &JMax, &KMax, &ICellMax, &JCellMax, &KCellMax,
            &SolTime, &StrandID, &ParentZn, &IsBlock,
            &NFConns, &FNMode, &TotalNumFaceNodes,
            &TotalNumBndryFaces, &TotalNumBndryConn,
            NULL, NULL, NULL, &ShrConn);

        stat = TECDAT142(&ARRSIZEVAL, resX.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resY.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resDiv.data(), &DIsDouble);
        stat = TECDAT142(&ARRSIZEVAL, resPs.data(), &DIsDouble);
    }

    return 0;
}

int MACSolver2D::OutResClose() {
    INTEGER4 stat, whichFile;
    whichFile = 1;
    stat = TECFIL142(&whichFile);

    // close first file (velocity)
    stat = TECEND142();
    stat = TECEND142();

    return 0;
}
