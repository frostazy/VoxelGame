#pragma once

#include "HyVoxelLib.h"

namespace HyVoxel
{

	/// A record used to create intersection lists during the voxelization process.
	struct RayIntersection
	{
		Vector position;
		Vector normal;
		int material;
		int next;
	};

	/// An intersection list produced during the voxelization process.
	struct RayIntersectionList
	{
		int x, y, z;
		int head;
		int tail;
		int count;
	};

	/// Limit for the intersections to be tracked during voxelization
	const unsigned int MAX_INTERSECTIONS = 40 * BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE;

	/// Limit for the intersections to be tracked along one ray
	const unsigned int MAX_TEST_INTERSECTIONS = 0x7FFFF;

}
