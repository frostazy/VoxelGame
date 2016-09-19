#pragma once

#include "HyVoxelLib.h"

#include "Visualize/Contour.h"
#include "Layer/BlockDataLayer.h"
#include "Material/MaterialLibrary.h"

namespace HyVoxel {

class ContourThreadContextImpl : public IContourThreadContext
{
public:
	ContourThreadContext ctx;

	virtual void Release() override { delete this; }
	virtual void Reset() override
	{
		ctx.voxelData->clear();
	}
};

class HyVoxelInterfaceImpl : public IHyVoxelInterface
{
public:
	virtual void Release() override { delete this; }

	virtual IVoxelLayer* CreateBlockDataLayer() override
	{
		return new BlockDataLayer();
	}

	virtual void Voxelize(IVoxelLayer* blockLayer, const VoxelizationInputMesh& inputMesh, const Matrix& transformMtx, bool isAdd, OutVector<HyVoxel::CellId>*& changedCells) override;

    virtual void VoxelizeHF(const VoxelizationInputMesh& inputMesh, const Matrix& transformMtx, const int nSegmentsX, const int nSegmentsY, const int nSegmentsZ, VoxelHF*& voxelHF) override;

	virtual IContourThreadContext* CreateContourThreadContext() override
	{
		return new ContourThreadContextImpl();
	}

	virtual void Contour(CellId cellId, IContourThreadContext* contourThreadContext, OutputMesh& outMesh) override;

	virtual void GenerateCubeMesh(CellId cellId, const CellVoxelMargin1& voxelData, OutputMesh& outMesh) override;

	virtual void GetVisibleCellsAround(const Vector& pos, OutVector<CellId>*& outCells) override;

	// TODO: liang
	// to be removed
	virtual void LoadMaterial(const char* path) override
	{
		readMaterialDefinitions(path, HyVoxel::gMaterialLibrary);
	}
};

}
