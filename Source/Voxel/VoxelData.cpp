#include "HyVoxelPrivatePCH.h"

#include "VoxelData.h"

namespace HyVoxel {


	
/*
vertices

	Y		   X
	|		  /
	|		 /
	|  5-----------6
	| /|   /	  /|
	|/ |  /		 / |
	4-----------7  |
	|  |/		|  |
	|  1--------|--2
	| /  		| /
	|/			|/
	0-----------3----------------Z

*/
const int VoxelCorners[8][3] =
{
	{ 0, 0, 0 },
	{ 1, 0, 0 },
	{ 1, 0, 1 },
	{ 0, 0, 1 },
	{ 0, 1, 0 },
	{ 1, 1, 0 },
	{ 1, 1, 1 },
	{ 0, 1, 1 }
};

const int VoxelMovementDirections[7][3] =
{
	{ 1, 0, 0 },
	{ 0, 0, 1 },
	{ 0, 1, 0 },
	{ 1, 0, 1 },
	{ 1, 1, 0 },
	{ 0, 1, 1 },
	{ 1, 1, 1 }
};

/*
Local edges

	Y		   X
	|		  /
	|		 /
	|   -----------
	| /|   /	  /|
	|/ |  /		 / |
	|-----------   |
	|  |/		|  |
	1  /--------|--
	| 0  		| /
	|/			|/
	0----2-----------------------Z

*/

// connected with vertex #0, there are 3 edges
// the coordinates of the other end of each edge are:
const int VoxelEdgeEndpoints[3][3] =
{
	{ 1, 0, 0 },
	{ 0, 1, 0 },
	{ 0, 0, 1 }
};
// connected with vertex #0, there are 3 edges
// each edge has 4 voxels around it
// here is the position of each voxel assuming the original voxel locates in {0, 0, 0}
const int VoxelEdgeLinks[3][4][3] =
{
	{ { 0, 0, 0 },{ 0, 0, -1 },{ 0, -1, -1 },{ 0, -1, 0 } },
	{ { 0, 0, 0 },{ -1, 0, 0 },{ -1, 0, -1 },{ 0, 0, -1 } },
	{ { 0, 0, 0 },{ 0, -1, 0 },{ -1, -1, 0 },{ -1, 0, 0 } }
};


/*
edges

	Y		   X
	|		  /
	|		 /
	|   -------5--- 
	| 4|   /	  /|
	|/ 9  /		 6 |
	 -------7---   10
	|  |/		|  |
	8   ----1---|-- 
	| 0  	   11 2
	|/			|/
	 -----3----- ----------------Z

*/

//	we have 12 edges
//	edge #n has 2 vertices, here we record their index
const int VoxelEdges[12][2] =
{
	{ 0, 1 },
	{ 1, 2 },
	{ 2, 3 },
	{ 3, 0 },
	{ 4, 5 },
	{ 5, 6 },
	{ 6, 7 },
	{ 7, 4 },
	{ 0, 4 },
	{ 1, 5 },
	{ 2, 6 },
	{ 3, 7 }
};

//	we have 12 edges
//	edge #n has 2 vertices, here we record their position coordinates
const int VoxelEdgeNeighbors[12][2][3] =
{
	{ { 0, 0, 0 },{ 1, 0, 0 } },
	{ { 1, 0, 0 },{ 1, 0, 1 } },
	{ { 0, 0, 1 },{ 1, 0, 1 } },
	{ { 0, 0, 0 },{ 0, 0, 1 } },
	{ { 0, 1, 0 },{ 1, 1, 0 } },
	{ { 1, 1, 0 },{ 1, 1, 1 } },
	{ { 0, 1, 1 },{ 1, 1, 1 } },
	{ { 0, 1, 0 },{ 0, 1, 1 } },
	{ { 0, 0, 0 },{ 0, 1, 0 } },
	{ { 1, 0, 0 },{ 1, 1, 0 } },
	{ { 1, 0, 1 },{ 1, 1, 1 } },
	{ { 0, 0, 1 },{ 0, 1, 1 } }
};

/*
edges

	Y		   X
	|		  /
	|		 /
	|   -------5---
	| 4|   /	  /|
	|/ 9  /		 6 |
	-------7---   10
	|  |/		|  |
	8   ----1---|--
	| 0  	   11 2
	|/			|/
	-----3----- ----------------Z

*/

//TODO: liang
//	rename
//	VoxelTriangle ? ? ? !!!WO CAO NI DAYE!
//
//	we have 12 edges
//	edge #n has 3 neighbor voxels, each has offset...
//	assuming the original voxel locates at{ 0, 0, 0 }
const int VoxelTriangle[12][3][3] =
{
	{ { 0,  0, -1 },{ 0, -1, -1 },{ 0, -1,  0 } },
	{ { 1,  0,  0 },{ 1, -1,  0 },{ 0, -1,  0 } },
	{ { 0, -1,  0 },{ 0, -1,  1 },{ 0,  0,  1 } },
	{ { 0, -1,  0 },{ -1, -1,  0 },{ -1,  0,  0 } },
	{ { 0,  1,  0 },{ 0,  1, -1 },{ 0,  0, -1 } },
	{ { 0,  1,  0 },{ 1,  1,  0 },{ 1,  0,  0 } },
	{ { 0,  0,  1 },{ 0,  1,  1 },{ 0,  1,  0 } },
	{ { -1,  0,  0 },{ -1,  1,  0 },{ 0,  1,  0 } },
	{ { -1,  0,  0 },{ -1,  0, -1 },{ 0,  0, -1 } },
	{ { 0,  0, -1 },{ 1,  0, -1 },{ 1,  0,  0 } },
	{ { 1,  0,  0 },{ 1,  0,  1 },{ 0,  0,  1 } },
	{ { 0,  0,  1 },{ -1,  0,  1 },{ -1,  0,  0 } }
};
// we have 12 edges
// edge #n has 4 neighbor voxels, including the original voxel
// assuming the original voxel locates at{ 0, 0, 0 }
const int VoxelNeighbors[12][4][3] =
{
	{ { 0, -1, -1 },{ 0,  0, -1 },{ 0,  0,  0 },{ 0, -1,  0 } },
	{ { 0, -1,  0 },{ 0,  0,  0 },{ 1,  0,  0 },{ 1, -1,  0 } },
	{ { 0, -1,  0 },{ 0,  0,  0 },{ 0,  0,  1 },{ 0, -1,  1 } },
	{ { -1, -1,  0 },{ -1,  0,  0 },{ 0,  0,  0 },{ 0, -1,  0 } },
	{ { 0,  0, -1 },{ 0,  1, -1 },{ 0,  1,  0 },{ 0,  0,  0 } },
	{ { 0,  0,  0 },{ 0,  1,  0 },{ 1,  1,  0 },{ 1,  0,  0 } },
	{ { 0,  0,  0 },{ 0,  1,  0 },{ 0,  1,  1 },{ 0,  0,  1 } },
	{ { -1,  0,  0 },{ -1,  1,  0 },{ 0,  1,  0 },{ 0,  0,  0 } },
	{ { -1,  0, -1 },{ -1,  0,  0 },{ 0,  0,  0 },{ 0,  0, -1 } },
	{ { 0,  0, -1 },{ 0,  0,  0 },{ 1,  0,  0 },{ 1,  0, -1 } },
	{ { 0,  0,  0 },{ 0,  0,  1 },{ 1,  0,  1 },{ 1,  0,  0 } },
	{ { -1,  0,  0 },{ -1,  0,  1 },{ 0,  0,  1 },{ 0,  0,  0 } }
};


void adjustCellCoordinates(CellId& cell, int& x, int& y, int& z)
{
	int level, xc, yc, zc;
	unpackCellId(cell, level, xc, yc, zc);

	while (x < 0 || x >= BLOCK_DIMENSION ||
		y < 0 || y >= BLOCK_DIMENSION ||
		z < 0 || z >= BLOCK_DIMENSION)
	{
		if (x < 0)
		{
			x = BLOCK_DIMENSION + x;
			xc--;
		}
		if (y < 0)
		{
			y = BLOCK_DIMENSION + y;
			yc--;
		}
		if (z < 0)
		{
			z = BLOCK_DIMENSION + z;
			zc--;
		}
		if (x >= BLOCK_DIMENSION)
		{
			x = x - BLOCK_DIMENSION;
			xc++;
		}
		if (y >= BLOCK_DIMENSION)
		{
			y = y - BLOCK_DIMENSION;
			yc++;
		}
		if (z >= BLOCK_DIMENSION)
		{
			z = z - BLOCK_DIMENSION;
			zc++;
		}
	}

	cell = packCellId(level, xc, yc, zc);
}

}