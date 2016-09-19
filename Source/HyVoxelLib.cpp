#include "HyVoxelPrivatePCH.h"

#include "ProcModeling/Architecture/CodeToMesh.h"
#include "ProcModeling/Spaceship/SpaceshipGenerator.h"
#include "Common/WhiteNoise.h"
#include "Common/HyJobsImpl.h"
#include "Voxel/Layer/BlockDataLayer.h"
#include "Voxel/Voxelization.h"
#include "Voxel/Visualize/Contour.h"
#include "Voxel/HyVoxelImpl.h"

#include "HyVoxelLib.h"

class HyVoxelLibImpl : public IHyVoxelLib
{
	virtual void Init() override
	{
		HyVoxel::CWhiteNoise::initialize();
	}

	virtual HyVoxel::IHyJobMgrInterface* QueryHyJobMgrInterface() override
	{
		return new HyVoxel::JobMgr();
	}

	virtual HyVoxel::IHyArchitectureInterface* QueryHyArchitectureInterface() override
	{
		return new HyVoxel::HyArchitectureInterfaceImpl();
	}

	virtual HyVoxel::IHyVoxelInterface* QueryHyVoxelInterface() override
	{
		return new HyVoxel::HyVoxelInterfaceImpl();
	}

	virtual HyVoxel::IHySpaceshipInterface* QuerySpaceshipInterface() override
	{
		return new HyVoxel::HySpaceshipImpl();
	}

	virtual void Release() { delete this; }

	// TODO: liang
	// to be removed
	virtual void GetHVFuncs(HVFuncs* outFuncs) override
	{
		outFuncs->TestVoxelization = TestVoxelization;
		outFuncs->TestMeshToVoxel = TestMeshToVoxel;
	}
};

#if HY_ENABLE_UNIT_TESTS
TEST(DummyTest, ProveIsHit)
{
	EXPECT_EQ(1, 1);
}

extern "C" __declspec(dllexport) int RunAllUnitTests(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
#endif

extern "C" __declspec(dllexport) IHyVoxelLib* QueryIHyVoxelLib()
{
	return new HyVoxelLibImpl();
}

