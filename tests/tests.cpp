#include "test_common.h"
#include "LargeTests/test_levelset2d_large.h"
#include "LargeTests/test_mac2d_large.h"
#include "LargeTests/test_poisson2d_large.h"

int main(int argc, char *argv[]) {
	mkl_set_num_threads(4);

	// LevelSetTest_2D_Simple();
	// LevelSetTest_2D_ReinitOnly();
	// LevelSetTest_2D_ReinitOnly2();
	// LevelSetTest_2D_ReinitOnlyWithKinks();
	// LevelSetTest_2D_FirstTimeReinit();
	// LevelSetTest_2D_Sussman621_ReinitSussman();
	// test_poisson_GuassSeidel();
	// test_poisson_CG();
	// test_poisson_BiCGStab();
	// MAC2DTest_CavityFlow();
	// MAC2DTest_StationaryBubble();
	MAC2DTest_SmallAirBubbleRising();
	// MAC2DTest_LargeAirBubbleRising();
	// MAC2DTest_TaylorInstability();
}