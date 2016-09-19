#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

#include <vector>

namespace HyVoxel {

	#pragma pack(push, 8)

	class IContourThreadContext;

	class IVoxelLayer : public IHyVoxelInterfaceBase
	{
	public:
		virtual void GetContourData(
			CellId cell,
			IContourThreadContext* contourThreadContext,
			bool& empty,
			IHyVoxelInterfaceBase* threadContext) = 0;

		virtual void GetCellVoxelMargin1(
			CellId cell,
			CellVoxelMargin1& outData,
			bool& isEmpty) = 0;

		virtual void Clear() {}
		virtual void UpdateLod(const CellId* cellIds, int count) {}

		virtual IHyVoxelInterfaceBase* CreateThreadContext() { return NULL; }
	};

	#pragma pack(pop)

}
