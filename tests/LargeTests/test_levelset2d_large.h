#ifndef __TEST_LEVEL2DLARGE_H__
#define __TEST_LEVEL2DLARGE_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "../test_common.h"
#include "../include/common.h"
#include "../include/levelset2d.h"

int LevelSetTest_2D_Simple();
int LevelSetTest_2D_ReinitOnly();
int LevelSetTest_2D_Sussman621_ReinitSussman();

#endif __TEST_LEVEL2DLARGE_H__