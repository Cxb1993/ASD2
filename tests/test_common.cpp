#include "test_common.h"

inline int idx3(int nx, int i, int j) {
	return i + (nx + 2 * NBC3) * j;
}

int OutRes(int iter, double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& phi, const std::vector<double>& div,
	const std::vector<double>& ls, const int nx, const int ny,
	const double dx, const double dy, const double baseX, const double baseY, PLTTYPE m_PLTType) {

	if (m_PLTType == PLTTYPE::ASCII || m_PLTType == PLTTYPE::BOTH) {
		std::ofstream outF;
		std::string fname_vel(fname_vel_base + "_ASCII.plt");
		std::string fname_div(fname_div_base + "_ASCII.plt");
		if (iter == 0) {
			outF.open(fname_vel.c_str(), std::ios::out);

			outF << "TITLE = VEL" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"LS\", \"PHI\" " << std::endl;
			outF.close();

			outF.open(fname_div.c_str(), std::ios::out);
			outF << "TITLE = DIV" << std::endl;
			outF << "VARIABLES = \"X\", \"Y\", \"LS\", \"DIV\", \"PHI\" " << std::endl;
			outF.close();
		}

		std::vector<double>
			resU((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resV((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);

		std::vector<double> resDiv = div;

		for (int i = NBC3 + 1; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++)
			resU[idx3(nx, i, j)] = (u[idx3(nx, i, j)] + u[idx3(nx, i + 1, j)]) * 0.5;

		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3 + 1; j < ny + NBC3; j++)
			resV[idx3(nx, i, j)] = (v[idx3(nx, i, j)] + v[idx3(nx, i, j + 1)]) * 0.5;

		outF.open(fname_vel.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << nx << std::string(", J=") << ny
			<< std::string(", DATAPACKING=POINT")
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter
			<< std::endl;
		
		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++)
			outF << baseX + static_cast<double>(i + 0.5 - NBC3) * dx << std::string(",")
			<< baseY + static_cast<double>(j + 0.5 - NBC3) * dy << std::string(",")
			<< static_cast<double>(resU[idx3(nx, i, j)]) << std::string(",")
			<< static_cast<double>(resV[idx3(nx, i, j)]) << std::string(",")
			<< static_cast<double>(ls[idx3(nx, i, j)]) << std::string(",")
			<< static_cast<double>(phi[idx3(nx, i, j)]) << std::endl;

		outF.close();

		outF.open(fname_div.c_str(), std::ios::app);

		outF << std::string("ZONE T=\"") << iter
			<< std::string("\", I=") << nx << std::string(", J=") << ny
			<< std::string(", DATAPACKING=POINT")
			<< std::string(", SOLUTIONTIME=") << curTime
			<< std::string(", STRANDID=") << iter
			<< std::endl;
		for (int i = NBC3; i < nx + NBC3; i++)
		for (int j = NBC3; j < ny + NBC3; j++)
			outF << baseX + static_cast<double>(i + 0.5 - NBC3) * dx << std::string(",")
			<< baseY + static_cast<double>(j + 0.5 - NBC3) * dy << std::string(",")
			<< static_cast<double>(ls[idx3(nx, i, j)]) << std::string(",")
			<< static_cast<double>(resDiv[idx3(nx, i, j)]) << std::string(",")
			<< static_cast<double>(phi[idx3(nx, i, j)]) << std::endl;

		outF.close();
	}

	if (m_PLTType == PLTTYPE::BINARY || m_PLTType == PLTTYPE::BOTH) {
		INTEGER4 whichFile = 0, stat = 0;
		std::string fname_vel(fname_vel_base + "_BINARY.szplt");
		std::string fname_div(fname_div_base + "_BINARY.szplt");

		std::vector<double>
			resX((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resY((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resU((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resV((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resLS((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0),
			resPhi((nx + 2 * NBC3) * (ny + 2 * NBC3), 0.0);

		std::vector<double> resDiv = div;

		for (int j = NBC3; j < ny + NBC3; j++)
		for (int i = NBC3; i < nx + NBC3; i++) {
			resX[idx3(nx, i, j)] = baseX + (i + 0.5 - NBC3) * nx;
			resY[idx3(nx, i, j)] = baseY + (j + 0.5 - NBC3) * ny;
			resU[idx3(nx, i, j)] = (u[idx3(nx, i, j)] + u[idx3(nx, i + 1, j)]) * 0.5;
			resV[idx3(nx, i, j)] = (v[idx3(nx, i, j)] + v[idx3(nx, i, j + 1)]) * 0.5;
			resLS[idx3(nx, i, j)] = ls[idx3(nx, i, j)];
			resPhi[idx3(nx, i, j)] = phi[idx3(nx, i, j)];
		}

		INTEGER4 Debug = 1;
		INTEGER4 VIsDouble = 1; /* 0 = Single precision, 1 = Double precision*/
		INTEGER4 FileType = 0; /* 0 = Full, 1 = Grid only, 2 = Solution only*/
		INTEGER4 FileFormat = 1; /* 0 = plt, 1 = szplt*/

		if (iter == 0) {
			/*
			* Open the file and write the tecplot datafile
			* header information
			*/
			stat = TECINI142(const_cast<char *>(std::string("VELOCITY").c_str()),  /* Name of the entire dataset.  */
				const_cast<char *>(std::string("X, Y, U, V, LS, PHI").c_str()),
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
		INTEGER4 IMax = nx, JMax = ny, KMax = 0;

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

		INTEGER4 DIsDouble = 1;  // set DIsDouble to 0, for float values.

		INTEGER4 ARRSIZEVAL = IMax * JMax * KMax;

		stat = TECDAT142(&ARRSIZEVAL, resX.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resY.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resU.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resV.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resLS.data(), &DIsDouble);
		stat = TECDAT142(&ARRSIZEVAL, resPhi.data(), &DIsDouble);

		if (iter == 0) {
			stat = TECINI142(const_cast<char *>(std::string("DIVERGENCE").c_str()),  /* Name of the entire dataset.  */
				const_cast<char *>(std::string("X, Y, LS, DIV, PHI").c_str()),
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
		stat = TECDAT142(&ARRSIZEVAL, resPhi.data(), &DIsDouble);
	}

	return 0;
}

int OutResClose() {
	INTEGER4 stat, whichFile;
	whichFile = 1;
	stat = TECFIL142(&whichFile);

	// close first file (velocity)
	stat = TECEND142();
	stat = TECEND142();

	return 0;
}
