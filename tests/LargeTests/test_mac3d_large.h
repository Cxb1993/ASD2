#ifndef __TEST_MAC3DLARGE_H__
#define __TEST_MAC3DLARGE_H__

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
#include "../include/mac3d.h"
#include "../include/levelset3d.h"

int MAC3DTest_StationaryBubble();
int MAC3DTest_NonSurfaceTension();

#endif __TEST_MAC3DLARGE_H__