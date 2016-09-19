#pragma once

#include "mapindex.h"

namespace HyVoxel {


/// Number of voxels in a Cell
const int BLOCK_DATA_SIZE = BLOCK_DIMENSION*BLOCK_DIMENSION*BLOCK_DIMENSION;

/// Actual number of voxels in a Cell once margins are considered
const int BLOCK_CUBE_SIZE = BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE;
/// Maximum level to retain block data
const int MAX_BLOCK_LOD = CELL_LOD_MAX - 2;//

const double VECTOR_MIN_VALUE = -1.0;
const double VECTOR_MAX_VALUE = 2.0;
const double VECTOR_MAX_CODE = 255.0 / (VECTOR_MAX_VALUE - VECTOR_MIN_VALUE);
const double VECTOR_DEFAULT_VALUE = ((int)(floor((VECTOR_MAX_CODE*(0.5 - VECTOR_MIN_VALUE)) + 0.5))) / VECTOR_MAX_CODE + VECTOR_MIN_VALUE;

// TODO: liang
// rename these variables.
extern const int VoxelCorners[8][3];
extern const int VoxelMovementDirections[7][3];
extern const int VoxelEdgeEndpoints[3][3];
extern const int VoxelEdgeLinks[3][4][3];
extern const int VoxelEdges[12][2];
extern const int VoxelEdgeNeighbors[12][2][3];
extern const int VoxelTriangle[12][3][3];
extern const int VoxelNeighbors[12][4][3];

void adjustCellCoordinates(CellId& cell, int& x, int& y, int& z);


inline unsigned char encodeVectorCoord(double coord)
{
	double code = std::min(VECTOR_MAX_VALUE, std::max(VECTOR_MIN_VALUE, coord)) - VECTOR_MIN_VALUE;
	code = floor((VECTOR_MAX_CODE*code) + 0.5);

	return (unsigned char)code;
}

inline double decodeVectorCoord(unsigned char coord)
{
	double value = (double)coord / VECTOR_MAX_CODE + VECTOR_MIN_VALUE;

	return value;
}





}