#pragma once

#pragma pack(push, 8)
class IHyVoxelInterfaceBase
{
public:
	virtual ~IHyVoxelInterfaceBase() {}
	virtual void Release() = 0;
};
#pragma pack(pop)

#include "HyJobs.h"
#include "HyContainer.h"
#include "HyMath.h"
#include "HyVoxelCell.h"
#include "HyVoxelLayer.h"
#include "HyVoxel.h"
#include "HyArchitecture.h"
#include "HySpaceship.h"

#include <string>
#include <vector>
#include <map>
#include <set>

#pragma pack(push, 8)

namespace HyVoxel
{
	/// Cell Dimensions in decimeters
	const float CELL_WIDTH = 30.0;

	/// How many octree levels will be included in every scene.
	/// Each level is eight times larger than the level before.
	const int CELL_LOD_MAX = 13;

	/// First level of the world octree that is actually visible
	const int CELL_LOD_MIN = 2;

}

// TODO: liang
// to be removed
struct HVFuncs;

class IHyVoxelLib : public IHyVoxelInterfaceBase
{
public:
	virtual void Init() = 0;
	virtual HyVoxel::IHyJobMgrInterface* QueryHyJobMgrInterface() = 0;
	virtual HyVoxel::IHyArchitectureInterface* QueryHyArchitectureInterface() = 0;
	virtual HyVoxel::IHyVoxelInterface* QueryHyVoxelInterface() = 0;
	virtual HyVoxel::IHySpaceshipInterface* QuerySpaceshipInterface() = 0;

	// TODO: liang
	// to be removed
	virtual void GetHVFuncs(HVFuncs* outFuncs) = 0;
};

// TODO: liang
// to be removed
#if 1
typedef void(*TESTVOXELIZATION)(const char* datPath, const HyVoxel::Vector& boxSize, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, bool isCubeContour);
typedef void(*TESTMESHTOVOXEL)(std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, const HyVoxel::Vector& boxSize, bool isCubeContour);

struct HVFuncs
{
	TESTVOXELIZATION			TestVoxelization;
	TESTMESHTOVOXEL				TestMeshToVoxel;
};
#endif

#pragma pack(pop)
