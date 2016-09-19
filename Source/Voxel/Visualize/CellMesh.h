#pragma once

#include "Common/FastQuadrics.h"

#include "Voxel/Material/MaterialLibrary.h"
#include "Voxel/mapindex.h"

namespace HyVoxel {

// A map containing the instance name (id) and the internal index of the instance
extern TVFMap<VFString, int> instanceNameIdMap;
// A map containing the internal index and the instance name (id) of the instance
extern TVFMap<int, VFString> instanceIdNameMap;


class CellMesh
{
public:
	/// Number of different mediums in the world. Each medium may be subject to a different rendering logic. Medium 0 is reserved to the world geometry. Medium 1 denotes tree crowns and other forms of foliage that will be replaced by billboards, but still must used for casting shadows.
	const static int MEDIUM_MAX = 4;
	/// This medium is used for solid world geometry
	const static int MEDIUM_SOLID = 0;
	/// This medium is used for foliage
	const static int MEDIUM_FOLIAGE = 1;
	/// This medium is use for a placement preview
	const static int MEDIUM_PREVIEW = 2;
	/// This medium is use for glass and transparent solids
	const static int MEDIUM_WATER = 3;


	/// Cell identifier
	CellId cellId;
	/// A mesh containing the polygons that appear in the Cell
	CFastQuadrics fqs[MEDIUM_MAX];

	struct MediumMesh
	{
		void SetFaceCount(int faceCount)
		{
			vertexCount = faceCount * 3;
			vertices.resize(vertexCount);
			normals.resize(vertexCount);
			faceNormals.resize(vertexCount);

			materials.resize(vertexCount);
			lightColors.resize(vertexCount);
		}

		int vertexCount;
		std::vector<HyVoxel::Vector>			vertices;
		std::vector<HyVoxel::Vector>			normals;
		std::vector<HyVoxel::Vector>			faceNormals;

		std::vector<unsigned int>	materials;
		std::vector<unsigned int>	lightColors;
	};
	MediumMesh meshes[MEDIUM_MAX];

	/// A mesh containing the polygons in the Cell seam before simplification
	CFastQuadrics seamMesh[MEDIUM_MAX][2][3];



	void ProcessMesh(CMaterialLibrary& materialLibrary);
};

}