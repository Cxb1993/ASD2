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
#include "../include/mac2d.h"
#include "../include/mac2daxisym.h"

int ENOTest_1D_Burgers();

int LevelSetTest_2D_Simple();
int LevelSetTest_2D_ReinitOnly_Original();
int LevelSetTest_2D_ReinitOnly_Sussman();
int LevelSetTest_2D_ReinitOnly_MinRK2();
int LevelSetTest_2D_ReinitOnly_MinRK2_2();
int LevelSetTest_2D_ReinitOnlyWithKinks();
int LevelSetTest_2D_Sussman621_ReinitSussman();
int LevelSetTest_2D_FirstTimeReinit();

int LevelSetTest_2DAxisym_ReinitOnly_Sussman();

#endif __TEST_LEVEL2DLARGE_H__
