#pragma once

#include "mapindex.h"

#include "HyVoxelConfig.h"
#include "Voxel/Layer/BlockDataLayer.h"

namespace HyVoxel {

struct VisibleCellSet
{
	OutVector<CellId>*		cells = nullptr;
	bool					isValid;
	HyVoxel::Vector			viewPos;

	VisibleCellSet()
	: isValid(false)
	{}

	~VisibleCellSet() { delete cells; }

	inline bool IsOutOfRange(const HyVoxel::Vector& inViewPos)
	{
		const float VIEW_POS_RADIUS = 4 * CELL_WIDTH;

		if (!isValid)
			return true;

		HyVoxel::Vector offset = inViewPos - viewPos;
		if (offset.GetAbsMax() > VIEW_POS_RADIUS)
			return true;

		return false;
	}

	void GenerateCellSetAround(const HyVoxel::Vector& pos);
};

class ClipmapView
{
public:
	BlockDataLayer			blockLayer;

	bool CalculateNextSceneIfNeeded();
	bool ContourLoop(std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals);

protected:
	HyVoxel::Vector					viewPos;
	VisibleCellSet			curScene;

	TVFSet<CellId>			cpuPendingCells;
};

}