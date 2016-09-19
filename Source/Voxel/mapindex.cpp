/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#include "HyVoxelPrivatePCH.h"

#include "mapindex.h"
#include <fstream>

using namespace HyVoxel;

// 260 was taken from windef.h
#ifndef MAX_PATH
	#define MAX_PATH  260
#endif

MapIndex::MapIndex(int aKeyLevel, char* aIndexPath)
{
	keyLevel = aKeyLevel;
	indexPath = aIndexPath;
	topCell = NULL;
	if (aIndexPath != NULL)
	{
		topCell = loadSection(0);
	}
}

MapIndex::~MapIndex()
{
}

bool MapIndex::cellIsEmpty(CellId cell, int level, int x, int y, int z)
{
	if (topCell == NULL)
	{
		return false;
	}

	Cells* cells;
	if (level < keyLevel)
	{
		int shifts = keyLevel - level;
		CellId key = packCellId(keyLevel, x >> shifts, y >> shifts, z >> shifts);
		lock.Lock();
		cells = cache[key];
		lock.Unlock();
		if (cells == NULL)
		{
			cells = loadSection(key);
			lock.Lock();
			cache[key] = cells;
			lock.Unlock();
		}
	}
	else
	{
		cells = topCell;
	}
	return cells->find(cell) == cells->end();
}

Cells* MapIndex::loadSection(CellId key)
{
	char filename[MAX_PATH];
	Cells* result = VF_NEW Cells();
	if (topCell == NULL || topCell->find(key) != topCell->end())
	{
		sprintf_s(filename, MAX_PATH, "%s\\%I64d.idx", indexPath, key);
		std::ifstream cellIndexFile(filename, std::ios::in | std::ios::binary);
		CellId cell;
		while (!cellIndexFile.eof())
		{
			cellIndexFile.read((char*)(&cell), sizeof(CellId));
			result->insert(cell);
		}
		cellIndexFile.close();
	}
	return result;
}

bool HyVoxel::CellSort(const CellId& a, const CellId& b)
{
	return a < b;
}
