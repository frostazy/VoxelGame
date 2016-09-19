
#include "HyVoxelPrivatePCH.h"

#include "ClipmapView.h"
#include "HyVoxelImpl.h"
#include "Voxel/Visualize/Contour.h"
#include "Common/HyContainerImpl.h"

namespace HyVoxel {



// TODO: liang
// The algorithm in this method should be very straightforward,
// but it is a mess right now.
// Revise it when having time.
void VisibleCellSet::GenerateCellSetAround(const HyVoxel::Vector& pos)
{
	HyVoxelInterfaceImpl impl;
	impl.GetVisibleCellsAround(pos, cells);

	viewPos = pos;
	isValid = true;
}

void HyVoxelInterfaceImpl::GetVisibleCellsAround(const Vector& pos, OutVector<CellId>*& outCells)
{
	std::set<CellId> cells;

	const float DECIMETER_TO_CENTIMETER = 10.0f;
	float oneOverCellWidth = 1.0f / (CELL_WIDTH * DECIMETER_TO_CENTIMETER);

	// flip y,z
	int x = pos.x * oneOverCellWidth;
	int y = pos.z * oneOverCellWidth;
	int z = pos.y * oneOverCellWidth;

	const int radiusLOD0 = 3;

	// Iterate over desired levels for scene
	int prevLODExtents[2][3];
	for (int level = CELL_LOD_MIN; level < CELL_LOD_MAX; level++)
	{
		// Find scene coordinates of new cell
		int scale = (1 << level);
		int xc = x / scale;
		int yc = y / scale;
		int zc = z / scale;

		int extents[2][3];
		bool extentsSet = false;
		//int delta = (level == CELL_LOD_MIN)? 6 : 3;
		int delta;
		if (level == CELL_LOD_MIN)
		{
			delta = radiusLOD0;
		}
		else if (level == CELL_LOD_MIN + 1)
		{
			delta = 0;
		}
		else if (level == CELL_LOD_MIN + 2)
		{
			delta = 0;
		}
		else
		{
			delta = 0;
		}
		int padLeftX = (xc - delta) % 2;
		int padLeftY = (yc - delta) % 2;
		int padLeftZ = (zc - delta) % 2;
		int padRightX = (xc + delta) % 2;
		int padRightY = (yc + delta) % 2;
		int padRightZ = (zc + delta) % 2;
		if (level > CELL_LOD_MIN)
		{
			while (xc - delta - padLeftX >= prevLODExtents[0][0])
			{
				padLeftX++;
			}
			while (yc - delta - padLeftY >= prevLODExtents[0][1])
			{
				padLeftY++;
			}
			while (zc - delta - padLeftZ >= prevLODExtents[0][2])
			{
				padLeftZ++;
			}
			while (xc + delta + padRightX <= prevLODExtents[1][0] + 1)
			{
				padRightX++;
			}
			while (yc + delta + padRightY <= prevLODExtents[1][1] + 1)
			{
				padRightY++;
			}
			while (zc + delta + padRightZ <= prevLODExtents[1][2] + 1)
			{
				padRightZ++;
			}
		}
		padLeftX += (xc - delta - padLeftX) % 2;
		padLeftY += (yc - delta - padLeftY) % 2;
		padLeftZ += (zc - delta - padLeftZ) % 2;
		padRightX += (xc + delta - padRightX) % 2;
		padRightY += (yc + delta - padRightY) % 2;
		padRightZ += (zc + delta - padRightZ) % 2;

		for (int x = xc - delta - padLeftX; x < xc + delta + padRightX; x++)
			for (int y = yc - delta - padLeftY; y < yc + delta + padRightY; y++)
				for (int z = zc - delta - padLeftZ; z < zc + delta + padRightZ; z++)
				{
					if (level > CELL_LOD_MIN &&
						x >= prevLODExtents[0][0] && x <= prevLODExtents[1][0] &&
						y >= prevLODExtents[0][1] && y <= prevLODExtents[1][1] &&
						z >= prevLODExtents[0][2] && z <= prevLODExtents[1][2])
					{
						continue;
					}

					int cellLevel = level;
					int cellX = x;
					int cellY = y;
					int cellZ = z;

					// Create cellId
					CellId cell = packCellId(cellLevel, cellX, cellY, cellZ);
					cells.insert(cells.end(), cell);

					// compute LOD extents
					if (extentsSet)
					{
						extents[0][0] = std::min(cellX / 2, extents[0][0]);
						extents[0][1] = std::min(cellY / 2, extents[0][1]);
						extents[0][2] = std::min(cellZ / 2, extents[0][2]);
						extents[1][0] = std::max(cellX / 2, extents[1][0]);
						extents[1][1] = std::max(cellY / 2, extents[1][1]);
						extents[1][2] = std::max(cellZ / 2, extents[1][2]);
					}
					else
					{
						extentsSet = true;
						extents[0][0] = cellX / 2;
						extents[0][1] = cellY / 2;
						extents[0][2] = cellZ / 2;
						extents[1][0] = cellX / 2;
						extents[1][1] = cellY / 2;
						extents[1][2] = cellZ / 2;
					}
				}
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
			{
				prevLODExtents[i][j] = extents[i][j];
			}
	}

	OutVectorImpl<CellId>* oVecImpl = new OutVectorImpl<CellId>();
	oVecImpl->insert(oVecImpl->begin(), cells.begin(), cells.end());

	oVecImpl->SetOutVectorValues(); outCells = oVecImpl;

	// TODO: liang
	// currently, we don't handle the case resolution > 1
	#if 0
	while (resolution > 1)
	{
		Scene* hiResScene = VF_NEW Scene();
		for (Scene::iterator i = scene->begin(); i != scene->end(); ++i)
		{
			int level, xc, yc, zc;
			unpackCellId(*i, level, xc, yc, zc);
			if (level <= CELL_LOD_MIN + 3)
			{
				hiResScene->insert(*i);
			}
			else
			{
				for (int dz = 0; dz < 2; dz++)
					for (int dx = 0; dx < 2; dx++)
						for (int dy = 0; dy < 2; dy++)
						{
							CellId hiResCell = packCellId(level - 1, 2 * xc + dx, 2 * yc + dy, 2 * zc + dz);
							hiResScene->insert(hiResCell);
						}
			}
		}
		VF_DELETE scene;
		scene = hiResScene;
		resolution--;
	}
	return scene;
	#endif

}

bool ClipmapView::CalculateNextSceneIfNeeded()
{
	if (!curScene.IsOutOfRange(viewPos))
		return false;

	curScene.GenerateCellSetAround(HyVoxel::Vector(0.0f));

	cpuPendingCells.insert(curScene.cells->pdata, curScene.cells->pdata + curScene.cells->length);

	return true;
}

bool ClipmapView::ContourLoop(std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals)
{
	ContourThreadContext contourCtx;

	std::vector<CellId> occupiedCells;
	while (cpuPendingCells.size() > 0)
	{
		CellId cellId = *(cpuPendingCells.begin());
		cpuPendingCells.erase(cellId);

		bool isEmpty = true;

		contourCtx.voxelData->clear();
		BlockDataLayer::ThreadContext	threadContext;
		blockLayer.GetContourVoxelData(cellId, contourCtx.voxelData, isEmpty, &threadContext);

		if (!isEmpty)
		{
			occupiedCells.push_back(cellId);

			CellMesh newCellMesh;
			newCellMesh.cellId = cellId;
			GenerateCellMesh(&contourCtx,
				&gMaterialLibrary,
				&newCellMesh,
				true);

			int lod, xc, yc, zc;
			unpackCellId(cellId, lod, xc, yc, zc);

			if (lod > CELL_LOD_MIN)
				continue;

			const float DECIMETER_TO_CENTIMETER = 10.0f;
			float scale = CELL_WIDTH * (1 << lod) * DECIMETER_TO_CENTIMETER;
			HyVoxel::Vector offset(xc*scale, zc*scale, yc*scale);

			CFastQuadrics& fq = newCellMesh.fqs[CellMesh::MEDIUM_SOLID];
			for (int itFace = 0; itFace < fq.faceCount; ++itFace)
			{
				FQ_Face& face = fq.faces[itFace];
				FQ_Vertex& vfV0 = fq.vertices[face[0]];
				FQ_Vertex& vfV1 = fq.vertices[face[1]];
				FQ_Vertex& vfV2 = fq.vertices[face[2]];

				HyVoxel::Vector v0(vfV0.x, vfV0.z, vfV0.y);
				HyVoxel::Vector v1(vfV1.x, vfV1.z, vfV1.y);
				HyVoxel::Vector v2(vfV2.x, vfV2.z, vfV2.y);

				v0 = v0*scale + offset;
				v1 = v1*scale + offset;
				v2 = v2*scale + offset;

				HyVoxel::Vector n = HyVoxel::Vector::CrossProduct(v2 - v0, v1 - v0);
				n.Normalize();

				verts.push_back(v0); normals.push_back(n);
				tris.push_back(verts.size() - 1);
				verts.push_back(v1); normals.push_back(n);
				tris.push_back(verts.size() - 1);
				verts.push_back(v2); normals.push_back(n);
				tris.push_back(verts.size() - 1);
			}
		}
	}

	return true;
}

}

