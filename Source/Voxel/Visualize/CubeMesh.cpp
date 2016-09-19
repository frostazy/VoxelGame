#include "HyVoxelPrivatePCH.h"

#include "Voxel/HyVoxelImpl.h"
#include "Voxel/mapindex.h"

#include "Common/HyContainerImpl.h"

namespace HyVoxel {


struct VoxelVisibility
{
	static const unsigned char X_PLUS	= 0x1 << 0;
	static const unsigned char X_MINUS	= 0x1 << 1;
	static const unsigned char Y_PLUS	= 0x1 << 2;
	static const unsigned char Y_MINUS	= 0x1 << 3;
	static const unsigned char Z_PLUS	= 0x1 << 4;
	static const unsigned char Z_MINUS	= 0x1 << 5;

	static const unsigned char INVISIBLE = 0x3f;

	unsigned char vis[BLOCK_DIMENSION][BLOCK_DIMENSION][BLOCK_DIMENSION];

	VoxelVisibility()
	{
		memset(vis, 0, sizeof(vis));
	}
};

const unsigned char VISIBILITY_FLAGS[] = {
	VoxelVisibility::X_PLUS,
	VoxelVisibility::X_MINUS,
	VoxelVisibility::Y_PLUS,
	VoxelVisibility::Y_MINUS,
	VoxelVisibility::Z_PLUS,
	VoxelVisibility::Z_MINUS,
};

void CalculateVisibility(const CellVoxelMargin1& voxelData, VoxelVisibility& vis)
{
	const int neighbors[6][3] = {
		{	1,	0,	0 },
		{  -1,	0,	0 },
		{	0,	1,	0 },
		{	0, -1,	0 },
		{	0,	0,	1 },
		{	0,	0, -1 },
	};

	for (int xi = 1; xi < BLOCK_DIMENSION + 1; ++xi)
		for (int yi = 1; yi < BLOCK_DIMENSION + 1; ++yi)
			for (int zi = 1; zi < BLOCK_DIMENSION + 1; ++zi)
			{
				unsigned char v = 0;
				for (int iNeig = 0; iNeig < 6; ++iNeig)
				{
					int xn = xi + neighbors[iNeig][0],
						yn = yi + neighbors[iNeig][1],
						zn = zi + neighbors[iNeig][2];

					bool isSolid = voxelData.getType(CellVoxelMargin1::Index(xn, yn, zn)) != 0;
					if (isSolid)
						v |= VISIBILITY_FLAGS[iNeig];
				}

				vis.vis[xi-1][yi-1][zi-1] = v;
			}
}

static void AddVoxelAt(float x, float y, float z, float cubeWidth, unsigned char visFlag, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals)
{
	float scale = cubeWidth;

	// flip y and z
	HyVoxel::Vector posOffset(x, z, y);

	// left handed
	float VERTS[] = {
		0.0f,	0.0f,	0.0f,
		1.0f,	0.0f,	0.0f,
		1.0f,	0.0f,	1.0f,
		0.0f,	0.0f,	1.0f,
		0.0f,	1.0f,	0.0f,
		1.0f,	1.0f,	0.0f,
		1.0f,	1.0f,	1.0f,
		0.0f,	1.0f,	1.0f,
	};

	// flip y and z
	// CCW
	const int FACES[6][2][3] = {
		{ { 2,5,1 }, { 2,6,5 } },	// X+
		{ { 7,0,4 }, { 7,3,0 } },	// X-
		{ { 7,2,3 }, { 7,6,2 } },	// Z+
		{ { 0,5,4 }, { 0,1,5 } },	// Z-
		{ { 6,4,5 }, { 6,7,4 } },	// Y+
		{ { 3,1,0 }, { 3,2,1 } },	// Y-
	};

	int triOffset = (int)verts.size();

	for (int iFace = 0; iFace < 6; ++iFace)
	{
		if (visFlag & VISIBILITY_FLAGS[iFace])
			continue;

		for (int iTri = 0; iTri < 2; ++iTri)
		{
			int i0 = FACES[iFace][iTri][0],
				i1 = FACES[iFace][iTri][1],
				i2 = FACES[iFace][iTri][2];

			HyVoxel::Vector v0(VERTS[i0 * 3 + 0], VERTS[i0 * 3 + 1], VERTS[i0 * 3 + 2]);
			HyVoxel::Vector v1(VERTS[i1 * 3 + 0], VERTS[i1 * 3 + 1], VERTS[i1 * 3 + 2]);
			HyVoxel::Vector v2(VERTS[i2 * 3 + 0], VERTS[i2 * 3 + 1], VERTS[i2 * 3 + 2]);

			v0 = v0 * scale + posOffset;
			v1 = v1 * scale + posOffset;
			v2 = v2 * scale + posOffset;

			HyVoxel::Vector n = HyVoxel::Vector::CrossProduct(v2 - v0, v1 - v0);

			verts.push_back(v0); normals.push_back(n);
			tris.push_back((int)verts.size() - 1);
			verts.push_back(v1); normals.push_back(n);
			tris.push_back((int)verts.size() - 1);
			verts.push_back(v2); normals.push_back(n);
			tris.push_back((int)verts.size() - 1);
		}
	}
}

void HyVoxelInterfaceImpl::GenerateCubeMesh(CellId cellId, const CellVoxelMargin1& voxelData, OutputMesh& outMesh)
{
	OutVectorImpl<HyVoxel::Vector>* overts = new OutVectorImpl<HyVoxel::Vector>();
	OutVectorImpl<HyVoxel::Vector>* onormals = new OutVectorImpl<HyVoxel::Vector>();
	OutVectorImpl<int>* otris = new OutVectorImpl<int>();

	VoxelVisibility vis;
	CalculateVisibility(voxelData, vis);

	int lod, xc, yc, zc;
	unpackCellId(cellId, lod, xc, yc, zc);

	const float oneOnDiv = 1.0f / BLOCK_DIMENSION;
	double scale = (1 << lod) * CELL_WIDTH * 10.0f; // unit: centimeter
	float cubeWidth = scale * oneOnDiv;

	for (int xi = 1; xi < BLOCK_DIMENSION + 1; ++xi)
		for (int yi = 1; yi < BLOCK_DIMENSION + 1; ++yi)
			for (int zi = 1; zi < BLOCK_DIMENSION + 1; ++zi)
			{
				bool isSolid = voxelData.getType(CellVoxelMargin1::Index(xi, yi, zi)) != 0;

				if (!isSolid)
					continue;

				unsigned char visFlag = vis.vis[xi - 1][yi - 1][zi - 1];
				if (visFlag == VoxelVisibility::INVISIBLE)
					continue;

				float xoff = (xc + (xi - 1) * oneOnDiv) * scale;
				float yoff = (yc + (yi - 1) * oneOnDiv) * scale;
				float zoff = (zc + (zi - 1) * oneOnDiv) * scale;

				AddVoxelAt(xoff, yoff, zoff, cubeWidth, visFlag,
					*overts, *otris, *onormals);
			}

	overts->SetOutVectorValues(); outMesh.verts = overts;
	onormals->SetOutVectorValues(); outMesh.normals = onormals;
	otris->SetOutVectorValues(); outMesh.tris = otris;
}

}