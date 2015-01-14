#ifndef __TEST_COMMON_H__
#define __TEST_COMMON_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// Tecplot  IO
#include "TECIO.h"
#include "TECXXX.h"

#include "common.h"

#define NBC3 (3)

int idx3(int ny, int i, int j);
int OutRes(int iter, double curTime, const std::string fname_vel_base, const std::string fname_div_base,
	const std::vector<double>& u, const std::vector<double>& v,
	const std::vector<double>& phi, const std::vector<double>& div, 
	const std::vector<double>& ls, const int nx, const int ny,
	const double dx, const double dy, const double baseX, const double baseY, PLTTYPE m_PLTType);
int OutResClose();

#endif __TEST_COMMON_H__