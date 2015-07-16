#include "test_common.h"
#include "LargeTests/test_levelset2d_large.h"
#include "LargeTests/test_mac2d_large.h"
#include "LargeTests/test_mac3d_large.h"
#include "LargeTests/test_poisson2d_large.h"
#include "LargeTests/test_poisson3d_large.h"

int main(int argc, char *argv[]) {
	mkl_set_num_threads(4);

	// LevelSetTest_2D_Simple();
	// LevelSetTest_2D_ReinitOnly_Original();
	// LevelSetTest_2D_ReinitOnly_Sussman();
	// LevelSetTest_2D_ReinitOnly_MinRK2();
	// LevelSetTest_2D_ReinitOnly_MinRK2_2();
	// LevelSetTest_2D_ReinitOnlyWithKinks();
	// LevelSetTest_2DAxisym_ReinitOnly_Sussman();
	// LevelSetTest_2DAxisym_ReinitOnly_MinRK2();
	// LevelSetTest_2D_FirstTimeReinit();
	// LevelSetTest_2D_Sussman621_ReinitSussman();

	// test_poisson2D_GuassSeidel(;)
	// test_poisson2D_CG();
	// test_poisson2D_BiCGStab();
	// test_poisson3D_CG();

	// MAC2DTest_CavityFlow();
	// MAC2DTest_NonSurfaceTension();
	// MAC2DTest_StationaryBubble();
	// MAC2DTest_SmallAirBubbleRising();
	// MAC2DTest_LargeAirBubbleRising();
	// MAC2DTest_TaylorInstability();

	MAC2DAxisymTest_NonSurfaceTension();
	// MAC2DAxisymTest_StationaryBubble();
	// MAC2DAxisymTest_SmallAirBubbleRising();
	// MAC2DAxisymTest_WaterDropletCollison1();
	
	// MAC3DTest_NonSurfaceTension();
	// MAC3DTest_StationaryBubble();
}