/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

#include "HyVoxelConfig.h"

#include "Common/ExternalMutex.h"

namespace HyVoxel
{
	/// A set of cells
	typedef TVFSet<CellId> Cells;

	/// An index of which Cells are empty and which ones contain information. For offline mode only.
	class MapIndex
	{
	public:
		MapIndex(int aKeyLevel, char* aIndexPath);
		~MapIndex();
		/// Returns whether a cell has data
		bool cellIsEmpty(CellId cell, int level, int x, int y, int z);
	private:
		TVFMap<CellId, Cells*> cache;
		int keyLevel;
		char* indexPath;
		Cells* topCell;
	protected:
		Cells* loadSection(CellId key);
		ExternalMutex::Mutex lock;
	};

	/// A sorting criteria for cells in a scene so closer cells appear first
	bool CellSort(const CellId& a, const CellId& b);

	/// Given a cell's ID computes its eight children
	inline void getChildCells(CellId cell, CellId childCells[2][2][2])
	{
		int level, xc, yc, zc;
		unpackCellId(cell, level, xc, yc, zc);

		// Load children data into array
		for (int qz = 0; qz < 2; qz++)
			for (int qx = 0; qx < 2; qx++)
				for (int qy = 0; qy < 2; qy++)
				{
					// Compute Id of child cell
					CellId childId = packCellId(
										 level - 1,
										 2*xc + qx,
										 2*yc + qy,
										 2*zc + qz);

					childCells[qx][qy][qz] = childId;
				}
	}



}