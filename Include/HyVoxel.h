#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

#include <cstdint>
#include <vector>

namespace HyVoxel {

#pragma pack(push, 8)

    class IContourThreadContext : public IHyVoxelInterfaceBase
    {
    public:
        virtual void Reset() = 0;
    };

    struct VoxelHF : public IHyVoxelInterfaceBase
    {
        virtual void GetIndices(int idx, int*& outData, int& outLength) = 0;
        virtual void GetSegments(int& _nSegmentsX, int& _nSegmentsY, int& _nSegmentsZ) = 0;
        virtual int  GetVoxelIndex(int IndexX, int IndexY, int IndexZ) = 0;
        virtual int  GetVoxelNumber() = 0;
        virtual void GetMinMaxPos(HyVoxel::Vector& _minPos, HyVoxel::Vector& _maxPos) = 0;
        virtual void GetVoxDimReciprocal(HyVoxel::Vector& _voxDimReciprocal) = 0;
        virtual void GetNewVerticesPerTriangle(int idx, HyVoxel::Vector*& outData, int& outLength) = 0;
        virtual void GetNewVerticesPerVoxel(int idx, HyVoxel::Vector*& outData, int& outLength) = 0;
    };

    class IHyVoxelInterface : public IHyVoxelInterfaceBase
    {
    public:
        virtual IVoxelLayer* CreateBlockDataLayer() = 0;

        virtual void Voxelize(IVoxelLayer* blockLayer, const VoxelizationInputMesh& inputMesh, const Matrix& transformMtx, bool isAdd, OutVector<HyVoxel::CellId>*& changedCells) = 0;

        virtual void VoxelizeHF(const VoxelizationInputMesh& inputMesh, const Matrix& transformMtx, const int nSegmentsX, const int nSegmentsY, const int nSegmentsZ, VoxelHF*& voxelHF) = 0;

        virtual IContourThreadContext* CreateContourThreadContext() = 0;
        virtual void Contour(CellId cellId, IContourThreadContext* contourThreadContext, OutputMesh& outMesh) = 0;

        virtual void GenerateCubeMesh(CellId cellId, const CellVoxelMargin1& voxelData, OutputMesh& outMesh) = 0;

        virtual void GetVisibleCellsAround(const Vector& pos, OutVector<CellId>*& outCells) = 0;

        // TODO: liang
        // to be removed
        virtual void LoadMaterial(const char* path) = 0;
    };

#pragma pack(pop)
}
