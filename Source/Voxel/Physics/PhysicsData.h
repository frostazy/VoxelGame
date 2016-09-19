/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

#include "HyVoxelConfig.h"
#include "mapindex.h"
//#include "CellData.h"
#include "PhysicsMaterials.h"

// contains various common data structures used in the physics/stampmesh pipeline

namespace HyVoxel
{
	namespace Physics
	{
		// records the number of voxels per material class
		struct CMassInfo
		{
			int totalVoxels = 0;
			int* perMatVoxels;
			CMassInfo()
			{
				perMatVoxels = (int*)VF_RAWCALLOC(cMATCLASSES, sizeof(int));
			}
			~CMassInfo()
			{
				VF_FREE(perMatVoxels);
			}
		};

		// container for noiseBrush
		struct CNoiseBrush
		{
			~CNoiseBrush()
			{
				VF_FREE(field);
			}
			unsigned short dimension[2];	// 0 -> cube side length; 1 -> number of worley cells
			unsigned short* field;			// "bit" field defining the mesh
		};

		// container for voxelBuffer
		class CVoxelBuffer
		{
		public:
			void clear()
			{
				voxelRanges.clear();
				voxelIds.clear();
			}
			struct nVoxel
			{
				unsigned char v[3];			// voxel vector
				int type;					// associated fragment type
			};
			struct nVoxelRange
			{
				int cellId = INT_MIN;		// cell padding used to compute voxel ids within the cell
				double pos[3];				// world pos of the cell

				// ranges
				int minx = INT_MAX;
				int miny = INT_MAX;
				int minz = INT_MAX;
				int maxx = INT_MIN;
				int maxy = INT_MIN;
				int maxz = INT_MIN;
			};
			TVFMap<int, nVoxelRange> voxelRanges;	// maps fragment type -> min/max voxel in (x,y,z)
			TVFMap<int, nVoxel> voxelIds;			// maps a unique voxel id -> voxel container
			int material = -1;					// hit position material id
		};

		// TODO: liang
		// uncomment
		#if 0
		// container for fragment render meshes
		struct CMovingSolid
		{
			CellMesh* mesh;
			double origin[3];
			double shift[3];
			double scale[3];
			double transform[16];
			double rot[4];
			clock_t clock;
			double buoyantForce;
			char material;
			CMovingSolid()
			{
				rot[0] = 0;
				rot[1] = 0;
				rot[2] = 0;
				rot[3] = 1;
				buoyantForce = -1;
				material = -1;
			}
		};
		#endif
	};
};
