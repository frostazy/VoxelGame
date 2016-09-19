#pragma once

#include "HyVoxelLib.h"

#include "Common/ExternalMutex.h"
#include "HyVoxelConfig.h"

namespace HyVoxel {

class BlockDataLayer : public IVoxelLayer
{
public:
	class ThreadContext : public IHyVoxelInterfaceBase
	{
	public:
		virtual void Release() override { delete this; }

		BlockVoxelData* kernel[3][3][3];

		ThreadContext()
		{
			for (int nz = 0; nz < 3; nz++)
				for (int ny = 0; ny < 3; ny++)
					for (int nx = 0; nx < 3; nx++)
					{
						kernel[nz][ny][nx] = VF_NEW BlockVoxelData();
					}
		}

		~ThreadContext()
		{
			for (int nz = 0; nz < 3; nz++)
				for (int ny = 0; ny < 3; ny++)
					for (int nx = 0; nx < 3; nx++)
					{
						VF_DELETE(kernel[nz][ny][nx]);
					}
		}
	};

	/// This interface allows the CBlockData object to persist blocks in an external storage
	class IBlockIO
	{
	public:
		/// Loads voxel data for the specified cell
		virtual BlockVoxelData* loadCell(CellId cell)
		{
			return NULL;
		}
		virtual ~IBlockIO() {};
		/// Saves voxel data for the specified cell
		virtual void saveCell(CellId cell, BlockVoxelData* data) {}
	};

	BlockDataLayer()
		:blockIO(NULL)
	{}

	virtual void Release() override { delete this; }

	/// Returns the voxel data buffer for the specified cell
	BlockVoxelData* fetchData(
		/// Cell ID for the cell to be retrieved
		CellId cell,
		/// If set to TRUE will create an empty cell buffer if no buffer is found for the cell
		bool create);

	virtual void GetContourData(
		CellId cell,
		IContourThreadContext* contourThreadContext,
		bool& empty,
		IHyVoxelInterfaceBase* threadContext) override;

	virtual void GetCellVoxelMargin1(
		CellId cell,
		CellVoxelMargin1& outData,
		bool& isEmpty) override;

	virtual void Clear() override;

	virtual IHyVoxelInterfaceBase* CreateThreadContext() override
	{
		return new ThreadContext();
	}

	// TODO: liang
	// This method is quite heavy. We'll process the cells one by one in the future.
	// modifiedCells should be on CELL_LOD_MIN
	virtual void UpdateLod(const CellId* cellIds, int count) override;


	void updateBlockLOD(CellId cell);

	/// Sets the IBlockIO interface
	void setBlockIOHandler(IBlockIO* blockIOHandler) { blockIO = blockIOHandler; }

	/// Loads voxel data for the specified cell using the assigned IBlockIO handler
	BlockVoxelData* loadCell(CellId cell);


protected:
	/// Returns a cache voxel data buffer for the specified cell
	virtual BlockVoxelData* fetchCacheData(
		/// Cell ID for the cell to be retrieved
		CellId cell,
		/// If set to true will create an empty cell buffer if no buffer is found for the cell
		bool create);

// TODO: liang
// test-only
public:
// private:
	TVFMap<CellId, BlockVoxelData*> blockCache;
	ExternalMutex::Mutex lock;
	IBlockIO* blockIO;

	/// Copies voxel data from a block and its neighbors
	void blockContourData(CellId cell, ContourVoxelData* data, bool& empty, BlockVoxelData* blockKernel[3][3][3]);

	// TODO: liang
	// verify the logic in this method
	/// Reads the entire Cell contents into a buffer
	void GetContourVoxelData(
		/// ID of the cell
		CellId cell,
		/// A buffer where the voxel data will be copied.
		ContourVoxelData* data,
		/// A flag notifying the entire cell is empty and could be discarded by the caller
		bool& empty,
		/// A thread context. This parameter is ignored since the CBlockData object does not require any threadsafe work buffers
		ThreadContext* threadContext);

};

int getContourIntersections(HyVoxel::Vector* result, ContourVoxelData* data, int x, int y, int z);

int getVoxelIntersections(
	Vector* result,
	BlockDataLayer* blocks,
	BlockVoxelData* data,
	CellId cell,
	int x, int y, int z,
	bool air = true
	);

}