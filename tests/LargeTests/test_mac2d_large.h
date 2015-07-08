#ifndef __TEST_MAC2DLARGE_H__
#define __TEST_MAC2DLARGE_H__

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
#include "../include/mac2d.h"
#include "../include/mac2daxisym.h"
#include "../include/levelset2d.h"

int MAC2DTest_CavityFlow();
int MAC2DTest_NonSurfaceTension();
int MAC2DTest_StationaryBubble();
int MAC2DTest_SmallAirBubbleRising();
int MAC2DTest_LargeAirBubbleRising();
int MAC2DTest_TaylorInstability();

int MAC2DAxisymTest_NonSurfaceTension();
int MAC2DAxisymTest_StationaryBubble();
int MAC2DAxisymTest_SmallAirBubbleRising();
int MAC2DAxisymTest_WaterDropletCollison1();

#endif __TEST_MAC2DLARGE_H__