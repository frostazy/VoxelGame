#include "HyVoxelPrivatePCH.h"

#include "BlockDataLayer.h"
#include "Common/qef.h"
#include "Common/Vector.h"
#include "HyVoxelConfig.h"
#include "HyVoxelLib.h"
#include "Voxel/VoxelData.h"
#include "Voxel/Visualize/Contour.h"
#include "Voxel/HyVoxelImpl.h"

using namespace HyVoxel;
using namespace HyVoxel::Algebra;

namespace
{

inline void getNeighborCoords(int m, int& nm, int& dm)
{
	if (m < BLOCK_MARGIN)
	{
		nm = 0;
		dm = BLOCK_DIMENSION - BLOCK_MARGIN;
	}
	else if (m >= BLOCK_DIMENSION + BLOCK_MARGIN)
	{
		nm = 2;
		dm = -(BLOCK_MARGIN + BLOCK_DIMENSION);
	}
	else
	{
		nm = 1;
		dm = -BLOCK_MARGIN;
	}
}

bool ContourKernelBoundaryConflict(ContourVoxelData* contourData, BlockVoxelData* blockKernel[3][3][3], int x, int y, int z)
{
	if (contourData == NULL)
	{
		return false;
	}

	int cx, cy, cz, dx, dy, dz;
	getNeighborCoords(x, cx, dx);
	getNeighborCoords(y, cy, dy);
	getNeighborCoords(z, cz, dz);
	BlockVoxelData* kernelData = blockKernel[cz][cy][cx];

	int kernelMat = VOXEL_ONLY_COORDS;
	if (kernelData != NULL)
	{
		BlockVoxelData::Index kernelIndex(x + dx, y + dy, z + dz);
		if (kernelData->hasType(kernelIndex, VOXEL_HAS_MATERIAL))
		{
			kernelMat = kernelData->getMaterial(kernelIndex);
		}
	}

	for (int nz = z; nz <= z + 1; nz++)
	{
		for (int nx = x; nx <= x + 1; nx++)
		{
			for (int ny = y; ny <= y + 1; ny++)
			{
				getNeighborCoords(nx, cx, dx);
				getNeighborCoords(ny, cy, dy);
				getNeighborCoords(nz, cz, dz);
				BlockVoxelData* neighborData = blockKernel[cz][cy][cx];

				int neighborMaterial = VOXEL_ONLY_COORDS;
				if (neighborData != NULL)
				{
					BlockVoxelData::Index neighborIndex(nx + dx, ny + dy, nz + dz);
					if (neighborData->hasType(neighborIndex, VOXEL_HAS_MATERIAL))
					{
						neighborMaterial = neighborData->getMaterial(neighborIndex);
					}
				}

				if ((neighborMaterial != kernelMat) && ((neighborMaterial == VOXEL_ONLY_COORDS) || (kernelMat == VOXEL_ONLY_COORDS)))
				{
					return true;
				}
			}
		}
	}

	return false;
}

int getKernelIntersections(HyVoxel::Vector* result, BlockVoxelData* blockKernel[3][3][3], int x, int y, int z)
{
	int pointCount = 0;

	for (int edge = 0; edge < 12; edge++)
	{
		int nsx = x + VoxelEdgeNeighbors[edge][0][0];
		int nsy = y + VoxelEdgeNeighbors[edge][0][1];
		int nsz = z + VoxelEdgeNeighbors[edge][0][2];

		int nx, ny, nz, dx, dy, dz;
		getNeighborCoords(nsx, nx, dx);
		getNeighborCoords(nsy, ny, dy);
		getNeighborCoords(nsz, nz, dz);
		nsx += dx;
		nsy += dy;
		nsz += dz;
		BlockVoxelData* nsData = blockKernel[nz][ny][nx];

		int ndx = x + VoxelEdgeNeighbors[edge][1][0];
		int ndy = y + VoxelEdgeNeighbors[edge][1][1];
		int ndz = z + VoxelEdgeNeighbors[edge][1][2];

		getNeighborCoords(ndx, nx, dx);
		getNeighborCoords(ndy, ny, dy);
		getNeighborCoords(ndz, nz, dz);
		ndx += dx;
		ndy += dy;
		ndz += dz;
		BlockVoxelData* ndData = blockKernel[nz][ny][nx];

		int nsMaterial = VOXEL_ONLY_COORDS;
		if (nsData != NULL)
		{
			BlockVoxelData::Index nsIndex(nsx, nsy, nsz);
			if (nsData->hasType(nsIndex, VOXEL_HAS_MATERIAL))
			{
				nsMaterial = nsData->getMaterial(nsIndex);
			}
		}

		int ndMaterial = VOXEL_ONLY_COORDS;
		if (ndData != NULL)
		{
			BlockVoxelData::Index ndIndex(ndx, ndy, ndz);
			if (ndData->hasType(ndIndex, VOXEL_HAS_MATERIAL))
			{
				ndMaterial = ndData->getMaterial(ndIndex);
			}
		}

		if ((nsMaterial != ndMaterial) && ((nsMaterial == VOXEL_ONLY_COORDS) || (ndMaterial == VOXEL_ONLY_COORDS)))
		{
			double points[3][3];
			int count = 0;

			for (int neighbor = 0; neighbor < 2; neighbor++)
			{
				int px = x + VoxelTriangle[edge][neighbor][0];
				int py = y + VoxelTriangle[edge][neighbor][1];
				int pz = z + VoxelTriangle[edge][neighbor][2];

				getNeighborCoords(px, nx, dx);
				getNeighborCoords(py, ny, dy);
				getNeighborCoords(pz, nz, dz);
				px += dx;
				py += dy;
				pz += dz;
				BlockVoxelData* pData = blockKernel[nz][ny][nx];

				if (pData != NULL)
				{
					BlockVoxelData::Index pIndex(px, py, pz);
					if (pData->hasType(pIndex, VOXEL_HAS_VECTOR))
					{
						pData->getVector(pIndex, points[neighbor][0], points[neighbor][1], points[neighbor][2]);

						points[neighbor][0] += VoxelTriangle[edge][neighbor][0];
						points[neighbor][1] += VoxelTriangle[edge][neighbor][1];
						points[neighbor][2] += VoxelTriangle[edge][neighbor][2];

						count++;
					}
				}
			}
			if (count == 2)
			{
				HyVoxel::Vector v0 = HyVoxel::Algebra::Vector_withValues((float)points[0][0], (float)points[0][1], (float)points[0][2]);
				HyVoxel::Vector v1 = HyVoxel::Algebra::Vector_withValues((float)points[1][0], (float)points[1][1], (float)points[1][2]);

				result[pointCount++] = v0;
				result[pointCount++] = v1;
			}
		}
	}

	return pointCount;
}

}

namespace HyVoxel {


void BlockDataLayer::GetContourVoxelData(
	/// ID of the cell
	CellId cell,
	/// A buffer where the voxel data will be copied.
	ContourVoxelData* data,
	/// A flag notifying the entire cell is empty and could be discarded by the caller
	bool& empty,
	/// A thread context. This parameter is ignored since the CBlockData object does not require any threadsafe work buffers
	ThreadContext* threadContext)
{
	int level, xc, yc, zc;
	unpackCellId(cell, level, xc, yc, zc);

	// The cell for contouring has a margin overlap with neighbor cells.
	// For this reason it is necessary to also load the 26 neighbor cells.
	// The "blocks" array contains this kernel of cells, where the cell at [1, 1, 1] is the center of the Kernel.
	// Copies must be stored in the thread context to prevent trashing from other threads while working with them
	ThreadContext* tc = threadContext;
	bool anyData = false;
	lock.Lock();
	BlockVoxelData* blocks[3][3][3];
	for (int nz = 0; nz < 3; nz++)
		for (int ny = 0; ny < 3; ny++)
			for (int nx = 0; nx < 3; nx++)
			{
				CellId ncell = packCellId(level, xc + nx - 1, yc + ny - 1, zc + nz - 1);
				BlockVoxelData* buffer = fetchData(ncell, false);
				if (buffer != NULL)
				{
					tc->kernel[nz][ny][nx]->copy(*buffer);
					blocks[nz][ny][nx] = tc->kernel[nz][ny][nx];
					anyData = true;
				}
				else
				{
					blocks[nz][ny][nx] = NULL;
				}
			}
	lock.Unlock();

	// TODO: liang
	// it's safe to release the lock here.

	if (anyData)
		blockContourData(cell, data, empty, blocks);
}

void BlockDataLayer::GetContourData(
	CellId cell,
	IContourThreadContext* contourThreadContext,
	bool& empty,
	IHyVoxelInterfaceBase* threadContext)
{
	this->GetContourVoxelData(cell,
		static_cast<ContourThreadContextImpl*>(contourThreadContext)->ctx.voxelData,
		empty,
		static_cast<ThreadContext*>(threadContext));
}

void BlockDataLayer::GetCellVoxelMargin1(
	CellId cell,
	CellVoxelMargin1& outData,
	bool& isEmpty)
{
	isEmpty = true;
	outData.clear();

	int lod, xc, yc, zc;
	unpackCellId(cell, lod, xc, yc, zc);

	int voxelRange[3][2] = {
		{ BLOCK_DIMENSION - 1, BLOCK_DIMENSION - 1},
		{ 0, BLOCK_DIMENSION - 1 },
		{ 0, 0}
	};

	for (int xcOffset = -1; xcOffset <= 1; ++xcOffset)
		for (int ycOffset = -1; ycOffset <= 1; ++ycOffset)
			for (int zcOffset = -1; zcOffset <= 1; ++zcOffset)
			{
				int xci = xc + xcOffset,
					yci = yc + ycOffset,
					zci = zc + zcOffset;
				CellId srcCellId = packCellId(lod, xci, yci, zci);
				BlockVoxelData* srcCell = fetchData(srcCellId, false);
				if (srcCell == NULL)
					continue;

				for (int xi = voxelRange[xcOffset+1][0]; xi <= voxelRange[xcOffset + 1][1]; ++xi)
					for (int yi = voxelRange[ycOffset + 1][0]; yi <= voxelRange[ycOffset + 1][1]; ++yi)
						for (int zi = voxelRange[zcOffset + 1][0]; zi <= voxelRange[zcOffset + 1][1]; ++zi)
						{
							Voxel voxel;
							srcCell->getVoxel(BlockVoxelData::Index(xi, yi, zi), voxel);
							if (voxel.type != 0)
							{
								int x = xi + BLOCK_DIMENSION * xcOffset + 1,
									y = yi + BLOCK_DIMENSION * ycOffset + 1,
									z = zi + BLOCK_DIMENSION * zcOffset + 1;
								outData.setVoxel(CellVoxelMargin1::Index(x, y, z), voxel);
								isEmpty = false;
							}

						}
			}
}

void BlockDataLayer::Clear()
{
	for (auto& kv : blockCache)
	{
		delete kv.second;
	}
	blockCache.clear();
}

void BlockDataLayer::UpdateLod(const CellId* cellIds, int count)
{
	TVFSet<CellId> updateCells;
	for (int i = 0; i < count; ++i)
	{
		CellId cellId = cellIds[i];

		int level, xc, yc, zc;
		unpackCellId(cellId, level, xc, yc, zc);
		if (level > CELL_LOD_MIN)
			continue;

		for (level = 1; level < CELL_LOD_MAX - CELL_LOD_MIN; level++)
		{
			CellId parent = packCellId(level + CELL_LOD_MIN, xc >> level, yc >> level, zc >> level);
			updateCells.insert(parent);
		}
	}

	for (auto& cellId : updateCells)
	{
		updateBlockLOD(cellId);
	}
}

static bool calculatePlanesQEF(Vector* childTriangles, BlockDataLayer* blockLayer, BlockVoxelData* data, CellId cell, int x, int y, int z,
	int sx, int sy, int sz, double &vx, double &vy, double &vz, MaterialId &mat, bool &vector, bool &bleedX, bool &bleedY, bool &bleedZ)
{
	//set default values
	mat = VOXEL_ONLY_COORDS;
	vector = false;
	vx = VECTOR_DEFAULT_VALUE;
	vy = VECTOR_DEFAULT_VALUE;
	vz = VECTOR_DEFAULT_VALUE;
	bleedX = false;
	bleedY = false;
	bleedZ = false;

	if (data == NULL)
	{
		return true;
	}

	bool hasVector = false;
	bool solid = false;
	TVFMap<MaterialId, int> matCounts;

	//check if there is any vector inside the parent scope
	for (int dz = 0; (dz <= sz) && !(hasVector && solid); dz++)
	{
		for (int dx = 0; (dx <= sx) && !(hasVector && solid); dx++)
		{
			for (int dy = 0; (dy <= sy) && !(hasVector && solid); dy++)
			{
				CellId nsCell = cell;
				int nsx = x + dx;
				int nsy = y + dy;
				int nsz = z + dz;
				adjustCellCoordinates(nsCell, nsx, nsy, nsz);
				BlockVoxelData* nsData = (cell == nsCell) ? data : blockLayer->fetchData(nsCell, false);

				if (nsData != NULL)
				{
					BlockVoxelData::Index nsIndex(nsx, nsy, nsz);

					MaterialId m = VOXEL_ONLY_COORDS;
					if (nsData->hasType(nsIndex, VOXEL_HAS_MATERIAL))
					{
						m = nsData->getMaterial(nsIndex);
						if (m != 0)
						{
							solid = true;
							matCounts[m]++;
						}
					}

					//take material in 0,0,0
					if (dx == 0 && dy == 0 && dz == 0)
					{
						mat = m;
					}

					if (!hasVector)
					{
						hasVector = nsData->hasType(nsIndex, VOXEL_HAS_VECTOR);
					}
				}
			}
		}
	}

	if (mat != 0 && mat != VOXEL_ONLY_COORDS)
	{
		int count = matCounts[mat];
		for (TVFMap<MaterialId, int>::const_iterator matIndex = matCounts.begin(); matIndex != matCounts.end(); ++matIndex)
		{
			int mCount = matIndex->second;
			if (mCount > count)
			{
				mat = matIndex->first;
				count = matIndex->second;
			}
		}
	}

	if (!hasVector)
	{
		return true;
	}

	bool air_material[3] = { false, false, false };
	bool material_air[3] = { false, false, false };
	const int SOLVER_MAX_POINTS = 64;
	const int MAX_Planes = 3 * 3 * 3 * 24;
	int childrenTriangleCount = 0;
	double maxSolution = 6.0;

	//contour to extract the surfaces that affect the scope of the parent voxel
	for (int dz = 0; dz <= sz; dz++)
	{
		for (int dx = 0; dx <= sx; dx++)
		{
			for (int dy = 0; dy <= sy; dy++)
			{
				CellId nsCell = cell;
				int nsx = x + dx;
				int nsy = y + dy;
				int nsz = z + dz;
				adjustCellCoordinates(nsCell, nsx, nsy, nsz);
				BlockVoxelData* nsData = (cell == nsCell) ? data : blockLayer->fetchData(nsCell, false);

				int nsMaterial = VOXEL_ONLY_COORDS;
				if (nsData != NULL)
				{
					BlockVoxelData::Index nsIndex(nsx, nsy, nsz);
					if (nsData->hasType(nsIndex, VOXEL_HAS_MATERIAL))
					{
						MaterialId m = nsData->getMaterial(nsIndex);
						if (!solid || (solid && (m != 0)))
						{
							nsMaterial = m;
						}
					}
				}

				for (int edge = 0; edge < 3; edge++)
				{
					int nx = dx + VoxelEdgeEndpoints[edge][0];
					int ny = dy + VoxelEdgeEndpoints[edge][1];
					int nz = dz + VoxelEdgeEndpoints[edge][2];
					if (nx>sx || ny>sy || nz>sz)
					{
						continue;
					}

					CellId ndCell = cell;
					int ndx = x + nx;
					int ndy = y + ny;
					int ndz = z + nz;

					adjustCellCoordinates(ndCell, ndx, ndy, ndz);
					BlockVoxelData* ndData = NULL;
					if (ndCell == cell)
					{
						ndData = data;
					}
					else if (ndCell == nsCell)
					{
						ndData = nsData;
					}
					else
					{
						ndData = blockLayer->fetchData(ndCell, false);
					}

					int ndMaterial = VOXEL_ONLY_COORDS;
					if (ndData != NULL)
					{
						BlockVoxelData::Index ndIndex(ndx, ndy, ndz);
						if (ndData->hasType(ndIndex, VOXEL_HAS_MATERIAL))
						{
							MaterialId m = ndData->getMaterial(ndIndex);
							if (!solid || (solid && (m != 0)))
							{
								ndMaterial = m;
							}
						}
					}

					if ((nsMaterial != ndMaterial) && ((nsMaterial == VOXEL_ONLY_COORDS) || (ndMaterial == VOXEL_ONLY_COORDS)))
					{
						HyVoxel::Vector triangles[4];
						for (int neighbor = 0; neighbor < 4; neighbor++)
						{
							CellId pCell = cell;
							int px = x + dx + VoxelEdgeLinks[edge][neighbor][0];
							int py = y + dy + VoxelEdgeLinks[edge][neighbor][1];
							int pz = z + dz + VoxelEdgeLinks[edge][neighbor][2];
							adjustCellCoordinates(pCell, px, py, pz);
							BlockVoxelData* pData = NULL;
							if (pCell == cell)
							{
								pData = data;
							}
							else if (pCell == nsCell)
							{
								pData = nsData;
							}
							else if (pCell == ndCell)
							{
								pData = ndData;
							}
							else
							{
								pData = blockLayer->fetchData(pCell, false);
							}

							triangles[neighbor] = HyVoxel::Algebra::Vector_withValues((float)VECTOR_DEFAULT_VALUE, (float)VECTOR_DEFAULT_VALUE, (float)VECTOR_DEFAULT_VALUE);
							if (pData != NULL)
							{
								BlockVoxelData::Index pIndex(px, py, pz);
								if (pData->hasType(pIndex, VOXEL_HAS_VECTOR))
								{
									pData->getVector(pIndex, triangles[neighbor].x, triangles[neighbor].y, triangles[neighbor].z);
								}
							}
							triangles[neighbor].x += dx + VoxelEdgeLinks[edge][neighbor][0];
							triangles[neighbor].y += dy + VoxelEdgeLinks[edge][neighbor][1];
							triangles[neighbor].z += dz + VoxelEdgeLinks[edge][neighbor][2];
						}

						if (nsMaterial == VOXEL_ONLY_COORDS)
						{
							air_material[edge] = true;
							childTriangles[childrenTriangleCount++] = triangles[0];
							childTriangles[childrenTriangleCount++] = triangles[1];
							childTriangles[childrenTriangleCount++] = triangles[2];
							childTriangles[childrenTriangleCount++] = triangles[0];
							childTriangles[childrenTriangleCount++] = triangles[2];
							childTriangles[childrenTriangleCount++] = triangles[3];
						}
						else
						{
							material_air[edge] = true;
							childTriangles[childrenTriangleCount++] = triangles[0];
							childTriangles[childrenTriangleCount++] = triangles[2];
							childTriangles[childrenTriangleCount++] = triangles[1];
							childTriangles[childrenTriangleCount++] = triangles[0];
							childTriangles[childrenTriangleCount++] = triangles[3];
							childTriangles[childrenTriangleCount++] = triangles[2];
						}
					}
				}
			}
		}
	}

	bleedX = air_material[0] && material_air[0];
	bleedY = air_material[1] && material_air[1];
	bleedZ = air_material[2] && material_air[2];

	//calculate the vector for the parent
	if (childrenTriangleCount > 0)
	{
		vector = true;

		double matrix[MAX_Planes][3];
		double vector[MAX_Planes];
		// Mass point for the voxel intersections
		Vector originOffset = Vector_withValues(1.0f, 1.0f, 1.0f);

		//load the plane equations for the QEF solution
		int planeCount = 0;
		for (int idx = 0; idx < childrenTriangleCount; idx += 3)
		{
			HyVoxel::Vector& A = childTriangles[idx];
			HyVoxel::Vector& B = childTriangles[idx + 1];
			HyVoxel::Vector& C = childTriangles[idx + 2];

			HyVoxel::Vector ab = HyVoxel::Algebra::Vector_subtract(B, A);
			HyVoxel::Vector ac = HyVoxel::Algebra::Vector_subtract(C, A);

			HyVoxel::Vector normal = HyVoxel::Algebra::Vector_cross(ab, ac);

			if (HyVoxel::Algebra::Vector_normalize(&normal))
			{
				matrix[planeCount][0] = normal.x;
				matrix[planeCount][1] = normal.y;
				matrix[planeCount][2] = normal.z;

				HyVoxel::Vector point = HyVoxel::Algebra::Vector_subtract(A, originOffset);
				vector[planeCount] = HyVoxel::Algebra::Vector_dot(normal, point);

				planeCount++;
			}
		}

		//decimate
		if (planeCount >= SOLVER_MAX_POINTS)
		{
			int newPlaneCount = 0;
			for (int i = 0; i<planeCount; i++)
			{
				if (i % 2 == 0)
				{
					matrix[newPlaneCount][0] = matrix[i][0];
					matrix[newPlaneCount][1] = matrix[i][1];
					matrix[newPlaneCount][2] = matrix[i][2];
					vector[newPlaneCount] = vector[i];
					newPlaneCount++;
				}
			}
			planeCount = newPlaneCount;
		}

		if (planeCount >= SOLVER_MAX_POINTS)
		{
			//too many planes, return the average of the child vectors
			int innerPoints = 0;
			double sumX = 0;
			double sumY = 0;
			double sumZ = 0;

			for (int idx = 0; idx < childrenTriangleCount; idx += 3)
			{
				HyVoxel::Vector& v = childTriangles[idx];

				sumX += v.x;
				sumY += v.y;
				sumZ += v.z;
				innerPoints++;
			}

			vx = sumX / innerPoints;
			vy = sumY / innerPoints;
			vz = sumZ / innerPoints;

			if (bleedX || bleedY || bleedZ)
			{
				return false;
			}
			else
			{
				return true;
			}
		}

		int finalPlaneCount = planeCount;
		if (finalPlaneCount < 3)
		{
			for (int plane = 0; plane<3 - planeCount; plane++)
			{
				matrix[planeCount + plane][0] = matrix[0][0];
				matrix[planeCount + plane][1] = matrix[0][1];
				matrix[planeCount + plane][2] = matrix[0][2];
				vector[planeCount + plane] = vector[0];
			}

			finalPlaneCount = 3;
		}

		//solve the QEF for the voxel scope
		double px, py, pz;
		HyVoxel::Algebra::QEF::evaluate(matrix, vector, finalPlaneCount, px, py, pz);

		vx = px + originOffset.x;
		vy = py + originOffset.y;
		vz = pz + originOffset.z;

		if (bleedX || bleedY || bleedZ)
		{
			if ((abs(px)>maxSolution) || (abs(py)>maxSolution) || (abs(pz)>maxSolution))
			{
				return false;
			}

			//check the error of the solution only if there is bleeding in the voxel scope
			//the error is calculated as the average of the distances from the QEF solution point to all the planes
			double error = 0;
			HyVoxel::Algebra::QEFMatrix quadrics(0.0);
			for (int idx = 0; idx < planeCount; idx++)
			{
				double plane[4];
				plane[0] = matrix[idx][0];
				plane[1] = matrix[idx][1];
				plane[2] = matrix[idx][2];
				plane[3] = -1 * vector[idx];

				HyVoxel::Algebra::QEFMatrix q(plane);
				quadrics += q;
			}

			error = HyVoxel::Algebra::vertex_error(quadrics, px, py, pz);
			if (error > 0.05)
			{
				return false;
			}
		}
	}
	else
	{
		hasVector = false;
	}

	return true;
}


void BlockDataLayer::updateBlockLOD(CellId cell)
{
	BlockVoxelData* blockData = NULL;
	int level, xc, yc, zc;
	unpackCellId(cell, level, xc, yc, zc);

	// Only parent cells need to be computed
	if (level == CELL_LOD_MIN)
	{
		return;
	}

	const int maxChildTriangleVertex = 12 * 2 * 3; //12 edges * 2 triangles * 3 vertices
	const int maxTriangleVertex = 8 * maxChildTriangleVertex;

	Vector* childrenTriangles = VF_ALLOC(Vector, maxTriangleVertex);
	//for solutions and pending operations in the voxel
	//first byte for operations, second byte for solutions
	int* operations = VF_ALLOC(int, BLOCK_DATA_SIZE);

	memset(operations, 0, BLOCK_DATA_SIZE*sizeof(int));

	int pendingBufferIndex = 0;
	int processingBufferIndex;
	TVFSet<int> pendingOperations[2];

	double lowValue = VECTOR_MIN_VALUE;
	double highValue = VECTOR_MAX_VALUE;

	// The parent cell is broken into eight octants. Each octant has its data defined by a
	// child cell.
	for (int qz = 0; qz < 2; qz++)
	{
		for (int qx = 0; qx < 2; qx++)
		{
			for (int qy = 0; qy < 2; qy++)
			{
				// Compute Id of child cell
				int cellX = 2 * xc + qx;
				int cellY = 2 * yc + qy;
				int cellZ = 2 * zc + qz;

				CellId childId = packCellId(level - 1, cellX, cellY, cellZ);

				// Get a pointer to the memory buffer where child voxels are stored
				BlockVoxelData* childData = fetchData(childId, false);

				// If the child contains no data, go on and process next one.
				if (childData == NULL)
				{
					continue;
				}

				// Get pointer to memory buffer for parent cell
				if (blockData == NULL)
				{
					blockData = fetchData(cell, true);
					blockData->clear();
				}

				// Iterate over every voxel in the cell.
				// Each parent voxel must get its value from the eight corresponding child voxels
				for (int parentZ = 0; parentZ < BLOCK_DIMENSION / 2; parentZ++)
				{
					for (int parentX = 0; parentX < BLOCK_DIMENSION / 2; parentX++)
					{
						for (int parentY = 0; parentY < BLOCK_DIMENSION / 2; parentY++)
						{
							// Compute coordinates of voxel in parent cell
							int xi = qx*BLOCK_DIMENSION / 2 + parentX;
							int yi = qy*BLOCK_DIMENSION / 2 + parentY;
							int zi = qz*BLOCK_DIMENSION / 2 + parentZ;

							// Get address of parent voxel in memory buffer
							BlockVoxelData::Index parentIndex(xi, yi, zi);

							MaterialId parentMaterial = VOXEL_ONLY_COORDS;
							bool parentVector = false;
							double pvx, pvy, pvz;
							bool solved = false;

							int childX = 2 * parentX;
							int childY = 2 * parentY;
							int childZ = 2 * parentZ;

							//they signalize the directions in which we should move the scope of the parent voxel
							bool bleedX = false;
							bool bleedY = false;
							bool bleedZ = false;

							double px, py, pz;
							//try to find a solution for the voxel
							bool solutionQEF = calculatePlanesQEF(childrenTriangles,
								this, childData, childId, childX, childY, childZ, 2, 2, 2,
								px, py, pz, parentMaterial, parentVector, bleedX, bleedY, bleedZ);

							if (!parentVector && (xi == BLOCK_DIMENSION - 1 || yi == BLOCK_DIMENSION - 1 || zi == BLOCK_DIMENSION - 1))
							{
								int nx = (xi == BLOCK_DIMENSION - 1) ? 1 : 0;
								int ny = (yi == BLOCK_DIMENSION - 1) ? 1 : 0;
								int nz = (zi == BLOCK_DIMENSION - 1) ? 1 : 0;
								solutionQEF = calculatePlanesQEF(childrenTriangles,
									this, childData, childId, childX, childY, childZ, 2 + nx, 2 + ny, 2 + nz,
									px, py, pz, parentMaterial, parentVector, bleedX, bleedY, bleedZ);
							}

							if (parentVector)
							{
								pvx = std::min(highValue, std::max(lowValue, 0.5*px));
								pvy = std::min(highValue, std::max(lowValue, 0.5*py));
								pvz = std::min(highValue, std::max(lowValue, 0.5*pz));

								solved = solutionQEF;

								//check for material changes
								if (parentMaterial == 0 || parentMaterial == VOXEL_ONLY_COORDS)
								{
									//get if there is any pending operation
									int operation = operations[parentIndex] & 0x07;
									if (operation != 0)
									{
										bool parentChange = false;
										//decode the direction in which we should move the scope of the voxel
										int offX = operation & 0x01;
										int offY = (operation & 0x02) >> 1;
										int offZ = (operation & 0x04) >> 2;

										for (int dir = 0; dir<3 && (parentMaterial == 0 || parentMaterial == VOXEL_ONLY_COORDS); dir++)
										{
											int dirX = VoxelMovementDirections[dir][0];
											int dirY = VoxelMovementDirections[dir][1];
											int dirZ = VoxelMovementDirections[dir][2];

											//align the scope with the suggestion direction
											if (!((offX == 1 && dirX == 1) || (offY == 1 && dirY == 1) || (offZ == 1 && dirZ == 1)))
											{
												continue;
											}

											//get the coordinates for the first child voxel of the new parent scope
											int childX = 2 * parentX + dirX;
											int childY = 2 * parentY + dirY;
											int childZ = 2 * parentZ + dirZ;

											MaterialId childMaterial = VOXEL_ONLY_COORDS;
											BlockVoxelData::Index childIndex(childX, childY, childZ);
											if (childData->hasType(childIndex, VOXEL_HAS_MATERIAL))
											{
												childMaterial = childData->getMaterial(childIndex);
											}

											//take the material even when there is not a better solution, but only if it is solid
											if (parentMaterial != VOXEL_ONLY_COORDS)
											{
												if (parentMaterial == 0)
												{
													if ((childMaterial != VOXEL_ONLY_COORDS) && (childMaterial != 0))
													{
														parentChange = true;
														parentMaterial = childMaterial;
													}
												}
												else
												{
													if ((childMaterial == VOXEL_ONLY_COORDS) || (childMaterial == 0))
													{
														parentChange = true;
													}
												}
											}
											else if (childMaterial != VOXEL_ONLY_COORDS)
											{
												parentChange = true;
												parentMaterial = childMaterial;
											}
										}

										if (parentChange)
										{
											//add new pending operations for the neighbor voxels backward
											for (int ndz = -1; ndz < 1; ndz++)
											{
												for (int ndx = -1; ndx < 1; ndx++)
												{
													for (int ndy = -1; ndy < 1; ndy++)
													{
														//do not mark the current parent
														if ((ndx == 0) && (ndy == 0) && (ndz == 0))
														{
															continue;
														}

														int nxi = xi + ndx;
														int nyi = yi + ndy;
														int nzi = zi + ndz;

														//mark only inside the cell
														if (nxi<0 || nyi<0 || nzi<0)
														{
															continue;
														}

														BlockVoxelData::Index neighborIndex(nxi, nyi, nzi);
														operations[neighborIndex] |= operation;
														pendingOperations[pendingBufferIndex].insert(neighborIndex);
													}
												}
											}
										}
									}
								}
							}
							else
							{
								//if there is no vector the voxel is empty or solid and this is the solution
								solved = true;
								//mark the voxel solution as empty/solid
								operations[parentIndex] |= 0x80;
							}

							if (!solved)
							{
								for (int dir = 0; dir<7 && !solved; dir++)
								{
									int offX = VoxelMovementDirections[dir][0];
									int offY = VoxelMovementDirections[dir][1];
									int offZ = VoxelMovementDirections[dir][2];

									//align the scope with the bleading direction
									if (!((offX == 1 && bleedX) || (offY == 1 && bleedY) || (offZ == 1 && bleedZ)))
									{
										continue;
									}

									// Compute coordinates of the first child voxel in parent voxel scope
									childX = 2 * parentX + offX;
									childY = 2 * parentY + offY;
									childZ = 2 * parentZ + offZ;

									double px, py, pz;
									MaterialId offMat;
									bool hasVector, bx, by, bz;
									//try to find a solution for the voxel
									bool solutionQEF = calculatePlanesQEF(childrenTriangles,
										this, childData, childId, childX, childY, childZ, 2, 2, 2,
										px, py, pz, offMat, hasVector, bx, by, bz);
									if (solutionQEF && hasVector)
									{
										solved = true;

										//translate the coordinates from the child scope to the parent voxel
										pvx = std::min(highValue, std::max(lowValue, 0.5*(px + offX)));
										pvy = std::min(highValue, std::max(lowValue, 0.5*(py + offY)));
										pvz = std::min(highValue, std::max(lowValue, 0.5*(pz + offZ)));

										//assign the material
										if (xi>0 && yi>0 && zi>0)
										{
											if (parentMaterial == VOXEL_ONLY_COORDS)
											{
												parentMaterial = offMat;
											}
											else if (parentMaterial == 0)
											{
												if ((offMat != VOXEL_ONLY_COORDS) && (offMat != 0))
												{
													parentMaterial = offMat;
												}
											}
										}

										int operationCode = offX | (offY << 1) | (offZ << 2);
										//mark the solution for the voxel
										operations[parentIndex] |= (operationCode << 4);
										operations[parentIndex] |= operationCode;

										//mark the neighbor voxels for pending operations
										for (int ndz = -1; ndz < 2; ndz++)
										{
											for (int ndx = -1; ndx < 2; ndx++)
											{
												for (int ndy = -1; ndy < 2; ndy++)
												{
													//do not mark the current parent
													if ((ndx == 0) && (ndy == 0) && (ndz == 0))
													{
														continue;
													}

													//do not mark the current parent
													int nxi = xi + ndx;
													int nyi = yi + ndy;
													int nzi = zi + ndz;

													//mark only inside the cell
													if (nxi<0 || nyi<0 || nzi<0 || nxi >= BLOCK_DIMENSION || nyi >= BLOCK_DIMENSION || nzi >= BLOCK_DIMENSION)
													{
														continue;
													}

													BlockVoxelData::Index neighborIndex(nxi, nyi, nzi);
													operations[neighborIndex] |= operationCode;
													pendingOperations[pendingBufferIndex].insert(neighborIndex);
												}
											}
										}
									}
								}
							}

							//mark the voxel as solved
							if (solved)
							{
								operations[parentIndex] |= 0x08;
							}

							// Set parent voxel state
							if (parentMaterial != VOXEL_ONLY_COORDS)
							{
								blockData->setMaterial(parentIndex, parentMaterial);
							}

							if (parentVector)
							{
								blockData->setVector(parentIndex, pvx, pvy, pvz);
							}
						}
					}
				}
			}
		}
	}

	//second round, analize the voxels with pending operations
	//we have two buffers one for adding new pending operations and another one for processing previous operations
	int steps = 0;
	CellId prevCell = 0;
	BlockVoxelData* childData = NULL;

	while (!pendingOperations[pendingBufferIndex].empty() && (steps < 10))
	{
		processingBufferIndex = pendingBufferIndex;
		pendingBufferIndex = (pendingBufferIndex + 1) % 2;
		pendingOperations[pendingBufferIndex].clear();

		//try to find a better solution for voxels with pending operations
		for (TVFSet<int>::iterator i = pendingOperations[processingBufferIndex].begin(); i != pendingOperations[processingBufferIndex].end(); ++i)
		{
			int parentIndex = *i;

			//decode the coordinates for the parent voxel
			int zi = (int)(parentIndex / (BLOCK_DIMENSION*BLOCK_DIMENSION));
			int rest = parentIndex % (BLOCK_DIMENSION*BLOCK_DIMENSION);
			int xi = (int)(rest / BLOCK_DIMENSION);
			int yi = rest % BLOCK_DIMENSION;

			int cellX = 2 * xc + (int)(xi * 2 / BLOCK_DIMENSION);
			int cellY = 2 * yc + (int)(yi * 2 / BLOCK_DIMENSION);
			int cellZ = 2 * zc + (int)(zi * 2 / BLOCK_DIMENSION);

			CellId childId = packCellId(level - 1, cellX, cellY, cellZ);
			BlockVoxelData* childData = fetchData(childId, false);

			if ((prevCell != childId) || (childData == NULL))
			{
				prevCell = childId;
				childData = fetchData(childId, false);
			}

			if (childData == NULL)
			{
				continue;
			}

			//get the previous solution and the pending operations for the voxel
			int solution = (operations[parentIndex] >> 4) & 0x07;
			int operation = operations[parentIndex] & 0x07;

			//check if the operation was already applied to the voxel
			if ((operation & solution) != operation)
			{
				//figure out if the voxel is empty/solid
				bool voxelEmpty = (operations[parentIndex] & 0x80) != 0;
				//figure out if the voxel was solved
				bool voxelSolved = (operations[parentIndex] & 0x08) != 0;

				//decode the direction in which we should move the scope of the voxel
				int offX = operation & 0x01;
				int offY = (operation & 0x02) >> 1;
				int offZ = (operation & 0x04) >> 2;

				if (voxelEmpty)
				{
					//get the coordinates for the first child voxel of the new parent scope
					int childX = (2 * xi % BLOCK_DIMENSION);
					int childY = (2 * yi % BLOCK_DIMENSION);
					int childZ = (2 * zi % BLOCK_DIMENSION);

					MaterialId parentMaterial;
					double px, py, pz;
					bool vector, bx, by, bz;
					bool solved = false;
					bool solutionQEF = calculatePlanesQEF(childrenTriangles,
						this, childData, childId, childX, childY, childZ, 2 + offX, 2 + offY, 2 + offZ,
						px, py, pz, parentMaterial, vector, bx, by, bz);
					if (!vector)
					{
						solutionQEF = calculatePlanesQEF(childrenTriangles,
							this, childData, childId, childX, childY, childZ, 3, 3, 3,
							px, py, pz, parentMaterial, vector, bx, by, bz);
					}

					if (vector)
					{
						//an empty/solid vector takes a vector from the neighbor voxels
						double pvx = std::min(highValue, std::max(lowValue, 0.5*px));
						double pvy = std::min(highValue, std::max(lowValue, 0.5*py));
						double pvz = std::min(highValue, std::max(lowValue, 0.5*pz));
						blockData->setVector(parentIndex, pvx, pvy, pvz);
						operations[parentIndex] |= 0x80;
					}
				}
				else if (voxelSolved)
				{
					int offX = solution & 0x01;
					int offY = (solution & 0x02) >> 1;
					int offZ = (solution & 0x04) >> 2;

					int newOperation = solution ^ operation;
					int dirX = newOperation & 0x01;
					int dirY = (newOperation & 0x02) >> 1;
					int dirZ = (newOperation & 0x04) >> 2;

					//get the coordinates for the first child voxel of the new parent scope
					int childX = (2 * xi % BLOCK_DIMENSION) + offX;
					int childY = (2 * yi % BLOCK_DIMENSION) + offY;
					int childZ = (2 * zi % BLOCK_DIMENSION) + offZ;

					MaterialId offMaterial;
					double px, py, pz;
					bool vector, bx, by, bz;
					bool solutionQEF = calculatePlanesQEF(childrenTriangles,
						this, childData, childId, childX, childY, childZ, 2 + dirX, 2 + dirY, 2 + dirZ,
						px, py, pz, offMaterial, vector, bx, by, bz);
					if (solutionQEF)
					{
						//an empty/solid vector takes a vector from the neighbor voxels
						double pvx = std::min(highValue, std::max(lowValue, 0.5*(px + offX)));
						double pvy = std::min(highValue, std::max(lowValue, 0.5*(py + offY)));
						double pvz = std::min(highValue, std::max(lowValue, 0.5*(pz + offZ)));

						blockData->setVector(parentIndex, pvx, pvy, pvz);
					}
				}
			}
		}
	}

	//free the buffers
	VF_FREE(operations);
	VF_FREE(childrenTriangles);
}

/// Copies voxel data from a block and its neighbors
void BlockDataLayer::blockContourData(CellId cell, ContourVoxelData* data, bool& empty, BlockVoxelData* blockKernel[3][3][3])
{
	// nx, ny and nz specifies which BlockVoxelData
	// dx, dy, dz specifies which voxel inside BlockVoxelData[nz][ny][nx] 

	// The getNeighborCoords function is called to covert from the contouring voxel
	// coordinate into the corresponding voxel coordinate within a cell

	if (data == NULL)
	{
		return;
	}

	int nx, ny, nz, dx, dy, dz;

	ContourVoxelData* originalData = VF_NEW ContourVoxelData();
	originalData->copy(*data);

	HyVoxel::Vector* userDataPoints = VF_ALLOC(HyVoxel::Vector, 50);
	HyVoxel::Vector* contourDataPoints = VF_ALLOC(HyVoxel::Vector, 50);

	for (int z = 0; z < BLOCK_SIZE; z++)
	{
		getNeighborCoords(z, nz, dz);
		for (int x = 0; x < BLOCK_SIZE; x++)
		{
			getNeighborCoords(x, nx, dx);
			for (int y = 0; y < BLOCK_SIZE; y++)
			{
				getNeighborCoords(y, ny, dy);
				BlockVoxelData* kernelData = blockKernel[nz][ny][nx];
				if (kernelData != NULL)
				{
					int bx = x + dx;
					int by = y + dy;
					int bz = z + dz;
					// There is data for this cell
					//Voxel value = kernelData[BLOCK_DIMENSION*BLOCK_DIMENSION*bz + BLOCK_DIMENSION*bx + by];
					BlockVoxelData::Index blockIndex(bx, by, bz);
					ContourVoxelData::Index index(x, y, z);

					double lowValue = VECTOR_MIN_VALUE;
					double highValue = VECTOR_MAX_VALUE;

					bool dataHasVector = originalData->hasType(index, VOXEL_HAS_VECTOR);

					bool hasMaterial = kernelData->hasType(blockIndex, VOXEL_HAS_MATERIAL);
					bool hasVector = kernelData->hasType(blockIndex, VOXEL_HAS_VECTOR);
					if (hasMaterial || hasVector)
					{
						if (hasVector)
						{
							double vx, vy, vz;
							kernelData->getVector(blockIndex, vx, vy, vz);

							if (dataHasVector)
							{
								double dvx, dvy, dvz;
								originalData->getVector(index, dvx, dvy, dvz);

								if (ContourKernelBoundaryConflict(originalData, blockKernel, x, y, z))
								{
									int contourCount = getContourIntersections(contourDataPoints, originalData, x, y, z);

									if (contourCount > 0)
									{
										int userCount = getKernelIntersections(userDataPoints, blockKernel, x, y, z);
										if (userCount > 0)
										{
											const int SOLVER_MAX_POINTS = 64;
											double m[SOLVER_MAX_POINTS][3];
											double v[SOLVER_MAX_POINTS];
											int count = 0;

											// Mass point for the voxel intersections
											double originOffset[3] = { (float)VECTOR_DEFAULT_VALUE, (float)VECTOR_DEFAULT_VALUE, (float)VECTOR_DEFAULT_VALUE };

											HyVoxel::Vector userA = HyVoxel::Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);
											HyVoxel::Vector contourA = HyVoxel::Algebra::Vector_withValues((float)dvx, (float)dvy, (float)dvz);
											HyVoxel::Vector massUserA = HyVoxel::Algebra::Vector_withValues((float)(vx - originOffset[0]), (float)(vy - originOffset[1]),
												(float)(vz - originOffset[2]));
											HyVoxel::Vector massContourA = HyVoxel::Algebra::Vector_withValues((float)(dvx - originOffset[0]), (float)(dvy - originOffset[1]),
												(float)(dvz - originOffset[2]));

											// Load normals and points into arrays for solver
											for (int userIndex = 0; userIndex < userCount; userIndex += 2)
											{
												HyVoxel::Vector B = userDataPoints[userIndex];
												HyVoxel::Vector C = userDataPoints[userIndex + 1];

												HyVoxel::Vector ab = HyVoxel::Algebra::Vector_withValues(B.x - userA.x, B.y - userA.y, B.z - userA.z);
												HyVoxel::Vector ac = HyVoxel::Algebra::Vector_withValues(C.x - userA.x, C.y - userA.y, C.z - userA.z);

												HyVoxel::Vector normal = HyVoxel::Algebra::Vector_cross(ab, ac);
												HyVoxel::Algebra::Vector_normalize(&normal);
												m[count][0] = normal.x;
												m[count][1] = normal.y;
												m[count][2] = normal.z;
												v[count] = HyVoxel::Algebra::Vector_dot(normal, massUserA);
												count++;
											}

											for (int contourIndex = 0; contourIndex < contourCount; contourIndex += 2)
											{
												HyVoxel::Vector B = contourDataPoints[contourIndex];
												HyVoxel::Vector C = contourDataPoints[contourIndex + 1];

												HyVoxel::Vector ab = HyVoxel::Algebra::Vector_withValues(B.x - contourA.x, B.y - contourA.y, B.z - contourA.z);
												HyVoxel::Vector ac = HyVoxel::Algebra::Vector_withValues(C.x - contourA.x, C.y - contourA.y, C.z - contourA.z);

												HyVoxel::Vector normal = HyVoxel::Algebra::Vector_cross(ab, ac);
												HyVoxel::Algebra::Vector_normalize(&normal);
												m[count][0] = normal.x;
												m[count][1] = normal.y;
												m[count][2] = normal.z;
												v[count] = HyVoxel::Algebra::Vector_dot(normal, massContourA);
												count++;
											}

											double px, py, pz;
											HyVoxel::Algebra::QEF::evaluate(m, v, count, px, py, pz);
											px += originOffset[0];
											py += originOffset[1];
											pz += originOffset[2];

											vx = std::min(highValue, std::max(lowValue, px));
											vy = std::min(highValue, std::max(lowValue, py));
											vz = std::min(highValue, std::max(lowValue, pz));
										}
									}
								}
							}

							data->setVector(index, vx, vy, vz);
						}

						if (hasMaterial)
						{
							int material = kernelData->getMaterial(blockIndex);
							data->setMaterial(index, material);
						}

						empty = false;
					}
				}
			}
		}
	}

	VF_FREE(userDataPoints);
	VF_FREE(contourDataPoints);
	VF_DELETE originalData;
}

BlockVoxelData* BlockDataLayer::fetchData(
	/// Cell ID for the cell to be retrieved
	CellId cell,
	/// If set to TRUE will create an empty cell buffer if no buffer is found for the cell
	bool create)
{
	BlockVoxelData* blockData = fetchCacheData(cell, create);
	if (blockData == NULL)
	{
		// Run callbacks
		blockData = loadCell(cell);
	}

	if (blockData == NULL && create)
	{
		blockData = VF_NEW BlockVoxelData();
		blockData->clear();
	}

	if (blockData != NULL)
	{
		lock.Lock();
		blockCache[cell] = blockData;
		lock.Unlock();
	}

	return blockData;
}

BlockVoxelData* BlockDataLayer::loadCell(CellId cell)
{
	if (blockIO != NULL)
	{
		return blockIO->loadCell(cell);
	}
	else
	{
		return NULL;
	}
}

BlockVoxelData* BlockDataLayer::fetchCacheData(
	/// Cell ID for the cell to be retrieved
	CellId cell,
	/// If set to true will create an empty cell buffer if no buffer is found for the cell
	bool create)
{
	BlockVoxelData* blockData = NULL;

	lock.Lock();
	TVFMap<CellId, BlockVoxelData*>::iterator i = blockCache.find(cell);
	if (i != blockCache.end())
	{
		blockData = i->second;
	}
	lock.Unlock();

	if (blockData == NULL && create)
	{
		blockData = VF_NEW BlockVoxelData();
		blockData->clear();
	}

	if (blockData != NULL)
	{
		lock.Lock();
		blockCache[cell] = blockData;
		lock.Unlock();
	}

	return blockData;
}

int getContourIntersections(HyVoxel::Vector* result, ContourVoxelData* data, int x, int y, int z)
{
	int pointCount = 0;

	for (int edge = 0; edge < 12; edge++)
	{
		int nsx = x + VoxelEdgeNeighbors[edge][0][0];
		int nsy = y + VoxelEdgeNeighbors[edge][0][1];
		int nsz = z + VoxelEdgeNeighbors[edge][0][2];

		int ndx = x + VoxelEdgeNeighbors[edge][1][0];
		int ndy = y + VoxelEdgeNeighbors[edge][1][1];
		int ndz = z + VoxelEdgeNeighbors[edge][1][2];

		int nsMaterial = 0;
		if ((nsx >= 0) && (nsx < BLOCK_SIZE) && (nsy >= 0) && (nsy < BLOCK_SIZE) && (nsz >= 0) && (nsz < BLOCK_SIZE))
		{
			ContourVoxelData::Index nsIndex(nsx, nsy, nsz);
			if (data->hasType(nsIndex, VOXEL_HAS_MATERIAL))
			{
				nsMaterial = data->getMaterial(nsIndex);
			}
		}

		int ndMaterial = 0;
		if ((ndx >= 0) && (ndx < BLOCK_SIZE) && (ndy >= 0) && (ndy < BLOCK_SIZE) && (ndz >= 0) && (ndz < BLOCK_SIZE))
		{
			ContourVoxelData::Index ndIndex(ndx, ndy, ndz);
			if (data->hasType(ndIndex, VOXEL_HAS_MATERIAL))
			{
				ndMaterial = data->getMaterial(ndIndex);
			}
		}

		if ((nsMaterial != ndMaterial) && ((nsMaterial == 0) || (ndMaterial == 0)))
		{
			double points[3][3];
			int count = 0;
			for (int neighbor = 0; neighbor < 2; neighbor++)
			{
				int px = x + VoxelTriangle[edge][neighbor][0];
				int py = y + VoxelTriangle[edge][neighbor][1];
				int pz = z + VoxelTriangle[edge][neighbor][2];

				if ((px < 0) || (px >= BLOCK_SIZE) || (py < 0) || (py >= BLOCK_SIZE) || (pz < 0) || (pz >= BLOCK_SIZE))
				{
					break;
				}

				ContourVoxelData::Index pIndex(px, py, pz);
				if (data->hasType(pIndex, VOXEL_HAS_VECTOR))
				{
					data->getVector(pIndex, points[neighbor][0], points[neighbor][1], points[neighbor][2]);

					points[neighbor][0] += VoxelTriangle[edge][neighbor][0];
					points[neighbor][1] += VoxelTriangle[edge][neighbor][1];
					points[neighbor][2] += VoxelTriangle[edge][neighbor][2];

					count++;
				}
			}
			if (count == 2)
			{
				HyVoxel::Vector v0 = HyVoxel::Algebra::Vector_withValues((float)points[0][0], (float)points[0][1], (float)points[0][2]);
				HyVoxel::Vector v1 = HyVoxel::Algebra::Vector_withValues((float)points[1][0], (float)points[1][1], (float)points[1][2]);

				result[pointCount++] = v0;
				result[pointCount++] = v1;
			}
		}
	}

	return pointCount;
}

int getVoxelIntersections(
	Vector* result,
	BlockDataLayer* blocks,
	BlockVoxelData* data,
	CellId cell,
	int x, int y, int z,
	bool air)
{
	if (data == NULL)
	{
		return 0;
	}

	int pointCount = 0;
	int airMaterial = air ? VOXEL_ONLY_COORDS : 0;
	CellId nsLastCell = 0;
	CellId ndLastCell = 0;
	CellId pLastCell = 0;
	BlockVoxelData* nsData = NULL;
	BlockVoxelData* ndData = NULL;
	BlockVoxelData* pData = NULL;

	for (int edge = 0; edge < 12; edge++)
	{
		CellId nsCell = cell;
		int nsx = x + VoxelEdgeNeighbors[edge][0][0];
		int nsy = y + VoxelEdgeNeighbors[edge][0][1];
		int nsz = z + VoxelEdgeNeighbors[edge][0][2];
		adjustCellCoordinates(nsCell, nsx, nsy, nsz);

		if ((nsLastCell != nsCell) || (nsData == NULL))
		{
			nsData = (cell == nsCell) ? data : blocks->fetchData(nsCell, false);
			nsLastCell = nsCell;
		}

		CellId ndCell = cell;
		int ndx = x + VoxelEdgeNeighbors[edge][1][0];
		int ndy = y + VoxelEdgeNeighbors[edge][1][1];
		int ndz = z + VoxelEdgeNeighbors[edge][1][2];
		adjustCellCoordinates(ndCell, ndx, ndy, ndz);

		if ((ndLastCell != ndCell) || (ndData == NULL))
		{
			if (ndCell == cell)
			{
				ndData = data;
			}
			else if (ndCell == nsCell)
			{
				ndData = nsData;
			}
			else
			{
				ndData = blocks->fetchData(ndCell, false);
			}

			ndLastCell = ndCell;
		}

		int nsMaterial = airMaterial;
		BlockVoxelData::Index nsIndex(nsx, nsy, nsz);
		if ((nsData != NULL) && nsData->hasType(nsIndex, VOXEL_HAS_MATERIAL))
		{
			nsMaterial = nsData->getMaterial(nsIndex);
		}

		int ndMaterial = airMaterial;
		BlockVoxelData::Index ndIndex(ndx, ndy, ndz);
		if ((ndData != NULL) && ndData->hasType(ndIndex, VOXEL_HAS_MATERIAL))
		{
			ndMaterial = ndData->getMaterial(ndIndex);
		}

		if (nsMaterial != ndMaterial)
		{
			double points[3][3];
			int count = 0;
			for (int neighbor = 0; neighbor < 3; neighbor++)
			{
				CellId pCell = cell;
				int px = x + VoxelTriangle[edge][neighbor][0];
				int py = y + VoxelTriangle[edge][neighbor][1];
				int pz = z + VoxelTriangle[edge][neighbor][2];
				adjustCellCoordinates(pCell, px, py, pz);

				if ((pLastCell != pCell) || (pData == NULL))
				{
					if (pCell == cell)
					{
						pData = data;
					}
					else if (pCell == nsCell)
					{
						pData = nsData;
					}
					else if (pCell == ndCell)
					{
						pData = ndData;
					}
					else
					{
						pData = blocks->fetchData(pCell, false);
					}

					pLastCell = pCell;
				}

				if (pData != NULL)
				{
					BlockVoxelData::Index pIndex(px, py, pz);
					if (pData->hasType(pIndex, VOXEL_HAS_VECTOR))
					{
						pData->getVector(pIndex, points[neighbor][0], points[neighbor][1], points[neighbor][2]);
					}
					else
					{
						points[neighbor][0] = VECTOR_DEFAULT_VALUE;
						points[neighbor][1] = VECTOR_DEFAULT_VALUE;
						points[neighbor][2] = VECTOR_DEFAULT_VALUE;
					}

					points[neighbor][0] += VoxelTriangle[edge][neighbor][0];
					points[neighbor][1] += VoxelTriangle[edge][neighbor][1];
					points[neighbor][2] += VoxelTriangle[edge][neighbor][2];

					count++;
				}
			}

			if (count == 3)
			{
				HyVoxel::Vector v0 = HyVoxel::Algebra::Vector_withValues((float)points[0][0], (float)points[0][1], (float)points[0][2]);
				HyVoxel::Vector v1 = HyVoxel::Algebra::Vector_withValues((float)points[1][0], (float)points[1][1], (float)points[1][2]);
				HyVoxel::Vector v2 = HyVoxel::Algebra::Vector_withValues((float)points[2][0], (float)points[2][1], (float)points[2][2]);

				if (nsMaterial == airMaterial)
				{
					result[pointCount++] = v2;
					result[pointCount++] = v1;
					result[pointCount++] = v0;
				}
				else
				{
					result[pointCount++] = v0;
					result[pointCount++] = v1;
					result[pointCount++] = v2;
				}
			}
		}
	}

	return pointCount;
}

}