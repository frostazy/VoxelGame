#include "HyVoxelPrivatePCH.h"

#include "CellMesh.h"
#include "Voxel/mapindex.h"
#include "HyVoxelConfig.h"
#include "Common/Vector.h"

using namespace HyVoxel::Algebra;

namespace HyVoxel {

TVFMap<VFString, int> instanceNameIdMap;
TVFMap<int, VFString> instanceIdNameMap;

void CellMesh::ProcessMesh(CMaterialLibrary& materialLibrary)
{
	int level, xc, yc, zc;
	unpackCellId(cellId, level, xc, yc, zc);

	float cellSize = (float)((1 << level)*CELL_WIDTH);
	Vector up = Vector_withValues(0.0f, 1.0f, 0.0f);

	for (int medium = 0; medium < MEDIUM_MAX; ++medium)
	{
		CFastQuadrics& fq = this->fqs[medium];
		if (fq.faceCount == 0)
			continue;

		// TODO: liang
		// Currently, we only support solid.
		if (medium != MEDIUM_SOLID)
			continue;

		MediumMesh& mesh = meshes[medium];
		mesh.SetFaceCount(fq.faceCount);

		for (int vi = 0; vi < fq.vertexCount; ++vi)
		{
			FQ_Vertex& v = fq.vertices[vi];
			int type = INT_MAX;


		}
	}
}


}