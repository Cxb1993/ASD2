#ifndef __TEST_POISSON2DLARGE_H__
#define __TEST_POISSON2DLARGE_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>
#include <istream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../test_common.h"
#include "../include/common.h"
#include "../include/common.h"
#include "../include/poisson2d.h"

int test_poisson2D_GuassSeidel();
int test_poisson2D_CG();
int test_poisson2D_BiCGStab();

#endif __TEST_POISSON2DLARGE_H__