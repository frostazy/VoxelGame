#include "HyVoxelPrivatePCH.h"

#include "Voxelization.h"

#include "Voxel/VoxelData.h"
#include "HyVoxelLib.h"
#include "HyVoxelImpl.h"

#include "Common/Ray.h"
#include "Common/RTree.h"
#include "Common/MatrixAlg.h"
#include "Common/qef.h"
#include "Common/HyContainerImpl.h"
#include "Common/Utility.h"

using namespace HyVoxel::Algebra;

#define CROSS(dest,v1,v2) dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; dest[2]=v1[0]*v2[1]-v1[1]*v2[0];
#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
#define SUB(dest,v1,v2) dest[0]=v1[0]-v2[0]; dest[1]=v1[1]-v2[1]; dest[2]=v1[2]-v2[2];

namespace HyVoxel {

struct StampRayTest
{
	double pos;
	double basePos[2];
	bool facing;
	double normal[3];
	int solid;
	int next;
};

typedef RTree<size_t, double, 2, double, 32> PlaneIndex;

const int MAX_RAYS_SHOT = 30;
const int RAYS_SHOT = 9;
const double RAY_DELTA[RAYS_SHOT][2] =
{
	{ 0.00, 0.00 },
	{ 0.00, 1.00 },
	{ 1.00, 0.00 },
	{ 1.00, 1.00 },
	{ 0.50, 0.50 },
	{ 0.25, 0.25 },
	{ 0.25, 0.75 },
	{ 0.75, 0.25 },
	{ 0.75, 0.75 }
};

const int SMOOTH_RAYS_SHOT = 9;
const double SMOOTH_RAY_DELTA[SMOOTH_RAYS_SHOT][2] =
{
	{ 0.00, 0.00 },
	{ 0.00, 1.00 },
	{ 1.00, 0.00 },
	{ 1.00, 1.00 },
	{ 0.50, 0.50 },
	{ 0.25, 0.25 },
	{ 0.25, 0.75 },
	{ 0.75, 0.25 },
	{ 0.75, 0.75 }
};

struct PlaneRayIntersectContext
{
	/// List of triangles to test
	Vector* triangles;

	/// List of solids that originated the triangles
	int* solids;

	/// Ray origin
	double origin[3];

	const double* worldPos;

	/// Ray direction
	double dir[3];

	/// Main axis for the ray (X = 0, Y = 1, z = 2)
	int axis0;

	/// Secondary axis
	int axis1;

	/// Secondary axis
	int axis2;

	/// Ray origin point along the main axis
	double blockSize;
	double axisMin;

	int rayTestCount[MAX_RAYS_SHOT];

	StampRayTest** rayTests;

	int rayTestHead[MAX_RAYS_SHOT];

	int maxTestIntersections[MAX_RAYS_SHOT];

	/// initializes the maxTestIntersections buffer
	void setMaxTestIntersections(int MTI)
	{
		for (int i = 0; i < MAX_RAYS_SHOT; i++)
			maxTestIntersections[i] = MTI;
	};
	/// initialize default MTI buffer
	PlaneRayIntersectContext()
	{
		setMaxTestIntersections(MAX_TEST_INTERSECTIONS);
	};
};

// This function adds an intersection point to the list
void addIntersectionToIndex(
	int& intersectionCount,
	int& intersectionListCount,
	int* intersectionIndex,
	RayIntersection* intersections,
	RayIntersectionList* intersectionLists,
	double pvx, double pvy, double pvz,
	double nx, double ny, double nz,
	int x, int y, int z, int material)
{
	// Test that the maximum allowed intersection count was not reached
	if (intersectionCount < MAX_INTERSECTIONS - 1)
	{
		// Compute voxel index in buffer
		int voxelIndex = z*BLOCK_DIMENSION*BLOCK_DIMENSION + x*BLOCK_DIMENSION + y;

		// Find index to list of intersection for that voxel
		RayIntersection* listItem = NULL;
		int listIndex = intersectionIndex[voxelIndex];
		if (listIndex == -1)
		{
			// Voxel had no previous interesections, start a new list
			if (intersectionListCount < MAX_INTERSECTIONS - 1)
			{
				// Create list record
				listIndex = intersectionListCount;
				RayIntersectionList& list = intersectionLists[listIndex];

				// Initialize list
				list.count = 1;
				list.x = x;
				list.y = y;
				list.z = z;
				list.head = intersectionCount;
				list.tail = intersectionCount;
				listItem = &intersections[intersectionCount];
				intersectionListCount++;
				intersectionIndex[voxelIndex] = listIndex;
			}
		}
		else
		{
			// Get the existing list
			RayIntersectionList& list = intersectionLists[listIndex];

			// Increment number of elements in the list
			list.count++;
			RayIntersection& tail = intersections[list.tail];

			// Create a new record and insert it to list
			tail.next = intersectionCount;
			list.tail = intersectionCount;
			listItem = &intersections[intersectionCount];
		}
		if (listItem != NULL)
		{
			// Set intersection data
			listItem->next = -1;
			listItem->normal.x = (float)nx;
			listItem->normal.y = (float)ny;
			listItem->normal.z = (float)nz;
			listItem->position.x = (float)pvx;
			listItem->position.y = (float)pvy;
			listItem->position.z = (float)pvz;
			listItem->material = material;
			intersectionCount++;
		}
	}
}

void addIntersectionToIndex(
	int xsize, int ysize, int zsize,
	int& intersectionCount,
	int& intersectionListCount,
	int* intersectionIndex,
	RayIntersection* intersections,
	RayIntersectionList* intersectionLists,
	double pvx, double pvy, double pvz,
	double nx, double ny, double nz,
	int x, int y, int z, int material)
{
	// Test that the maximum allowed intersection count was not reached
	if (intersectionCount < MAX_INTERSECTIONS - 1)
	{
		// Compute voxel index in buffer
		int voxelIndex = z*ysize*xsize + x*ysize + y;

		// Find index to list of intersection for that voxel
		RayIntersection* listItem = NULL;
		int listIndex = intersectionIndex[voxelIndex];
		if (listIndex == -1)
		{
			// Voxel had no previous intersections, start a new list
			if (intersectionListCount < MAX_INTERSECTIONS - 1)
			{
				// Create list record
				listIndex = intersectionListCount;
				RayIntersectionList& list = intersectionLists[listIndex];

				// Initialize list
				list.count = 1;
				list.x = x;
				list.y = y;
				list.z = z;
				list.head = intersectionCount;
				list.tail = intersectionCount;
				listItem = &intersections[intersectionCount];
				intersectionListCount++;
				intersectionIndex[voxelIndex] = listIndex;
			}
		}
		else
		{
			if (listIndex < MAX_INTERSECTIONS - 1)
			{
				// Get the existing list
				RayIntersectionList& list = intersectionLists[listIndex];

				// Increment number of elements in the list
				list.count++;
				RayIntersection& tail = intersections[list.tail];

				// Create a new record and insert it to list
				tail.next = intersectionCount;
				list.tail = intersectionCount;
				listItem = &intersections[intersectionCount];
			}
		}

		if (listItem != NULL)
		{
			// Set intersection data
			listItem->next = -1;
			listItem->normal.x = (float)nx;
			listItem->normal.y = (float)ny;
			listItem->normal.z = (float)nz;
			listItem->position.x = (float)pvx;
			listItem->position.y = (float)pvy;
			listItem->position.z = (float)pvz;
			listItem->material = material;
			intersectionCount++;
		}
	}
}

inline double TriangleArea(double x1, double y1, double x2, double y2, double x3, double y3)
{
	// no need to div2; just mult2 the constant EPS below.
	return abs(((x1 - x3) * (y2 - y3)) - ((x2 - x3) * (y1 - y3)));
}

bool __cdecl planeIntersectCallback(size_t triangle, void* context)
{
	// Get the search context
	PlaneRayIntersectContext* ctx = (PlaneRayIntersectContext*)context;

	// Get which instance originated the of the matching triangle
	int solid = ctx->solids[triangle / 3];

	// Get triangle vertices
	Vector& tv0 = ctx->triangles[triangle];
	Vector& tv1 = ctx->triangles[triangle + 1];
	Vector& tv2 = ctx->triangles[triangle + 2];
	double v0[3] = { tv0.x, tv0.y, tv0.z };
	double v1[3] = { tv1.x, tv1.y, tv1.z };
	double v2[3] = { tv2.x, tv2.y, tv2.z };

	for (int i = 0; i < 3; i++)
	{
		v0[i] += ctx->worldPos[i];
		v1[i] += ctx->worldPos[i];
		v2[i] += ctx->worldPos[i];
	}

	// compute normal
	double v10[3] = { v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2] };
	double v20[3] = { v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2] };
	double normal[3];
	CROSS(normal, v20, v10);
	if (abs(normal[ctx->axis0]) < 0.000000000001)
	{
		return true;
	}

	double p1x, p1y, p2x, p2y, p3x, p3y;
	p1x = v0[ctx->axis1];
	p1y = v0[ctx->axis2];
	p2x = v1[ctx->axis1];
	p2y = v1[ctx->axis2];
	p3x = v2[ctx->axis1];
	p3y = v2[ctx->axis2];

	for (int i = 0; i < RAYS_SHOT; i++)
	{
		// project triangle to main axis plane and test in 2D
		double p2Dx, p2Dy;
		p2Dx = ctx->origin[ctx->axis1] + RAY_DELTA[i][0] * ctx->blockSize;
		p2Dy = ctx->origin[ctx->axis2] + RAY_DELTA[i][1] * ctx->blockSize;

		// See if projected point for ray is inside triangle
		static const double eps = 0.00002;	// old eps * 2
		double TotalArea = TriangleArea(p1x, p1y, p2x, p2y, p3x, p3y);
		double Area1 = TriangleArea(p2Dx, p2Dy, p2x, p2y, p3x, p3y);
		double Area2 = TriangleArea(p2Dx, p2Dy, p1x, p1y, p3x, p3y);
		double Area3 = TriangleArea(p2Dx, p2Dy, p1x, p1y, p2x, p2y);
		bool pointInside = (Area1 + Area2 + Area3) - TotalArea <= eps;

		if (pointInside)
		{
			// Ray actually intersects triangle.
			// Add aah n intersection record to the list
			int rayIndex = ctx->rayTestCount[i];
			StampRayTest* ri = ctx->rayTests[i] + rayIndex;
			ri->next = -1;

			// Compute intersection position
			double planeD = -(normal[0] * v0[0] + normal[1] * v0[1] + normal[2] * v0[2]);
			double d = -normal[ctx->axis1] * p2Dx - normal[ctx->axis2] * p2Dy - planeD;
			d /= normal[ctx->axis0];
			ri->pos = d - ctx->axisMin;
			ri->basePos[0] = p2Dx;
			ri->basePos[1] = p2Dy;

			// Store normal data
			ri->normal[0] = normal[0];
			ri->normal[1] = normal[1];
			ri->normal[2] = normal[2];

			// Compute if triangle faces ray
			double angle = DOT(normal, ctx->dir);
			ri->facing = angle < 0.0;
			ri->solid = solid;

			// Add intersection to list
			if (rayIndex == 0)
			{
				// It is the firt hit
				ctx->rayTestHead[i] = 0;
			}
			else
			{
				// Insert sorted
				StampRayTest* prev = NULL;
				StampRayTest* next = ctx->rayTests[i] + ctx->rayTestHead[i];
				while (next->next != -1 && next->pos < ri->pos)
				{
					prev = next;
					next = ctx->rayTests[i] + next->next;
				}
				if (next->next == -1 && next->pos <= ri->pos)
				{
					// It goes at the end of the list
					next->next = rayIndex;
				}
				else
				{
					// Insert before next
					if (prev == NULL)
					{
						// It is a new head for the list
						ri->next = ctx->rayTestHead[i];
						ctx->rayTestHead[i] = rayIndex;
					}
					else
					{
						// Insert between nodes
						ri->next = prev->next;
						prev->next = rayIndex;
					}
				}
			}

			ctx->rayTestCount[i]++;
			if (ctx->rayTestCount[i] == ctx->maxTestIntersections[i])
			{
				ctx->maxTestIntersections[i] *= 2;
				ctx->rayTests[i] = (StampRayTest*)VF_RAWREALLOC(ctx->rayTests[i], ctx->maxTestIntersections[i] * sizeof(StampRayTest));
			}
		}
	}

	// Return true so RTree search continues
	return true;
}

inline bool hasVoxelOccupancy(unsigned char* occupancyData, int hitCell[3])
{
	static const unsigned char VOXEL_OCCUPANCY_BLOCK = 22;	// ensure this matches!!
															// Read voxel occupancy for that position
	int x = (hitCell[0] + BLOCK_MARGIN) / 2;
	int y = (hitCell[1] + BLOCK_MARGIN) / 2;
	int z = (hitCell[2] + BLOCK_MARGIN) / 2;
	unsigned char v = occupancyData[z*VOXEL_OCCUPANCY_BLOCK*VOXEL_OCCUPANCY_BLOCK + x*VOXEL_OCCUPANCY_BLOCK + y];
	int vdx = hitCell[0] % 2;
	int vdy = hitCell[1] % 2;
	int vdz = hitCell[2] % 2;
	int vsub = 4 * vdz + 2 * vdx + vdy;
	unsigned char bit = 1 << vsub;
	v &= bit;
	return v != 0;
}

/// A compact 32 bit representation for a voxel
typedef uint32_t CompactVoxel;

void stampMeshQEF(
	BlockDataLayer* blockData,
	IMeshStampSource* mesh,
	IMeshStampMaterialSource* materials,
	const double worldPos[3],
	const Matrix& transform,
	std::set<CellId>* changedCells,
	Physics::CNoiseBrush* noise,
	Physics::CVoxelBuffer* nVoxBuffer,
	const TVFMap<CellId, unsigned char*>* voxelCache,
	const int hitMaterial)
{
	bool inserting = false;
	const double margin2D = 0.0001;
	double scale = (1 << CELL_LOD_MIN)*CELL_WIDTH;
	double voxelSize = scale / BLOCK_DIMENSION;

	CompactVoxel* overflow = VF_ALLOC(CompactVoxel, BLOCK_DATA_SIZE);

	// determine whether or not this is a physics request
	bool physicsRequest = !(noise == NULL);
	int noiseDim, noiseCells, ncellSpanned = -1;
	unsigned short* ncellCount = NULL;
	if (physicsRequest)
	{
		noiseDim = noise->dimension[0];
		noiseCells = noise->dimension[1];
		ncellCount = (unsigned short*)VF_RAWCALLOC(noiseDim*noiseDim*noiseDim, sizeof(unsigned short));
	}

	unsigned long* voxelMaterials = VF_ALLOC(unsigned long, BLOCK_DATA_SIZE);
	HyVoxel::Vector* userDataPoints = VF_ALLOC(HyVoxel::Vector, 50);

	BlockVoxelData* buffer = VF_NEW BlockVoxelData();
	RayIntersection* intersections = VF_ALLOC(RayIntersection, MAX_INTERSECTIONS);
	RayIntersectionList* intersectionLists = VF_ALLOC(RayIntersectionList, MAX_INTERSECTIONS);
	int* intersectionIndex = VF_ALLOC(int, BLOCK_DIMENSION*BLOCK_DIMENSION*BLOCK_DIMENSION);
	StampRayTest* rayTests[RAYS_SHOT];
	for (int i = 0; i < RAYS_SHOT; i++)
	{
		rayTests[i] = VF_ALLOC(StampRayTest, MAX_TEST_INTERSECTIONS);
	}

	int solidCount = mesh->getSolidCount();
	for (int solid = 0; solid < solidCount; solid++)
	{
		double meshMax[3] = { 0.0, 0.0, 0.0 };
		double meshMin[3] = { 0.0, 0.0, 0.0 };
		int triangleCount = mesh->getFaceCount(solid);

		int* solids = VF_ALLOC(int, triangleCount);
		Vector* triangles = VF_ALLOC(Vector, 3 * triangleCount);

		int index = 0;
		for (int face = 0; face < triangleCount; face++)
		{
			solids[face] = solid;
			Vector& v0 = triangles[index++];
			Vector& v1 = triangles[index++];
			Vector& v2 = triangles[index++];
			mesh->getFace(solid, face, v0, v1, v2);
			v0 = Matrix_multiplyVector(transform, v0);
			v1 = Matrix_multiplyVector(transform, v1);
			v2 = Matrix_multiplyVector(transform, v2);
			if (face == 0)
			{
				// First iteration, initialize min and max
				meshMin[0] = (double)std::min(v0.x, std::min(v1.x, v2.x));
				meshMin[1] = (double)std::min(v0.y, std::min(v1.y, v2.y));
				meshMin[2] = (double)std::min(v0.z, std::min(v1.z, v2.z));
				meshMax[0] = (double)std::max(v0.x, std::max(v1.x, v2.x));
				meshMax[1] = (double)std::max(v0.y, std::max(v1.y, v2.y));
				meshMax[2] = (double)std::max(v0.z, std::max(v1.z, v2.z));
			}
			else
			{
				meshMin[0] = std::min(meshMin[0], (double)std::min(v0.x, std::min(v1.x, v2.x)));
				meshMin[1] = std::min(meshMin[1], (double)std::min(v0.y, std::min(v1.y, v2.y)));
				meshMin[2] = std::min(meshMin[2], (double)std::min(v0.z, std::min(v1.z, v2.z)));
				meshMax[0] = std::max(meshMax[0], (double)std::max(v0.x, std::max(v1.x, v2.x)));
				meshMax[1] = std::max(meshMax[1], (double)std::max(v0.y, std::max(v1.y, v2.y)));
				meshMax[2] = std::max(meshMax[2], (double)std::max(v0.z, std::max(v1.z, v2.z)));
			}
		}

		for (int i = 0; i < 3; i++)
		{
			meshMin[i] += worldPos[i];
			meshMax[i] += worldPos[i];
		}

		int cellDX = (int)((meshMax[0] - meshMin[0]) / scale) + 1;
		int cellDY = (int)((meshMax[1] - meshMin[1]) / scale) + 1;
		int cellDZ = (int)((meshMax[2] - meshMin[2]) / scale) + 1;

		PlaneRayIntersectContext context;
		context.blockSize = voxelSize;
		context.triangles = triangles;
		context.solids = solids;
		context.rayTests = rayTests;
		context.worldPos = worldPos;

		// Define three 2D RTree indices, one for each main coordinate plane.
		PlaneIndex xy, xz, yz;

		// Collect all its triangles and insert their projections
		// into the corresponding plane RTree
		for (size_t i = 0; i < (size_t)3 * triangleCount; i += 3)
		{
			Vector& v0 = triangles[i];
			Vector& v1 = triangles[i + 1];
			Vector& v2 = triangles[i + 2];

			// Compute triangle AABB and insert into corresponding plane index
			double xymin[2] =
			{
				worldPos[0] + (double)std::min(std::min(v0.x, v1.x), v2.x) - margin2D,
				worldPos[1] + (double)std::min(std::min(v0.y, v1.y), v2.y) - margin2D
			};
			double xymax[2] =
			{
				worldPos[0] + (double)std::max(std::max(v0.x, v1.x), v2.x) + margin2D,
				worldPos[1] + (double)std::max(std::max(v0.y, v1.y), v2.y) + margin2D
			};
			xy.Insert(xymin, xymax, i);

			double zxmin[2] =
			{
				worldPos[0] + (double)std::min(std::min(v0.x, v1.x), v2.x) - margin2D,
				worldPos[2] + (double)std::min(std::min(v0.z, v1.z), v2.z) - margin2D
			};
			double zxmax[2] =
			{
				worldPos[0] + (double)std::max(std::max(v0.x, v1.x), v2.x) + margin2D,
				worldPos[2] + (double)std::max(std::max(v0.z, v1.z), v2.z) + margin2D
			};
			xz.Insert(zxmin, zxmax, i);

			double zymin[2] =
			{
				worldPos[1] + (double)std::min(std::min(v0.y, v1.y), v2.y) - margin2D,
				worldPos[2] + (double)std::min(std::min(v0.z, v1.z), v2.z) - margin2D
			};
			double zymax[2] =
			{
				worldPos[1] + (double)std::max(std::max(v0.y, v1.y), v2.y) + margin2D,
				worldPos[2] + (double)std::max(std::max(v0.z, v1.z), v2.z) + margin2D
			};
			yz.Insert(zymin, zymax, i);
		}

		// Voxelize over each affected cell
		const double cellMargin = 0.001;
		int xc = (int)(meshMin[0] / scale);
		int yc = (int)(meshMin[1] / scale);
		int zc = (int)(meshMin[2] / scale);

		if (meshMin[0] < 0)
			xc--;
		if (meshMin[1] < 0)
			yc--;
		if (meshMin[2] < 0)
			zc--;

		if (xc != ((int)((meshMax[0]) / scale)))
		{
			cellDX++;
		}
		if (yc != ((int)((meshMax[1]) / scale)))
		{
			cellDY++;
		}
		if (zc != ((int)((meshMax[2]) / scale)))
		{
			cellDZ++;
		}

		for (int zi = 0; zi < cellDZ; zi++)
		{
			for (int xi = 0; xi < cellDX; xi++)
			{
				for (int yi = 0; yi < cellDY; yi++)
				{
					ncellSpanned++;

					memset(intersectionIndex, 0xFF, BLOCK_DIMENSION*BLOCK_DIMENSION*BLOCK_DIMENSION*sizeof(int));
					int intersectionCount = 0;
					int intersectionListCount = 0;

					CellId currentCell = packCellId(CELL_LOD_MIN, xc + xi, yc + yi, zc + zi);
					BlockVoxelData* voxelData = blockData->fetchData(currentCell, true);

					// TODO: liang
					// For undo
					// blockData->trackCellChanges(currentCell, voxelData);
					if (changedCells != NULL)
					{
						changedCells->insert(currentCell);
					}
					if (physicsRequest && nVoxBuffer->material == -1)
						nVoxBuffer->material = hitMaterial;

					buffer->copy(*voxelData);
					memset(voxelMaterials, 0, BLOCK_DATA_SIZE*sizeof(unsigned long));

					double minX = meshMin[0];
					double minY = meshMin[1];
					double minZ = meshMin[2];

					// Step back a couple of voxels before voxelization
					// This is to make sure an intersection will happen from air to solid
					minX -= 2.0*voxelSize;
					minY -= 2.0*voxelSize;
					minZ -= 2.0*voxelSize;

					bool changed = false;
					bool changedNeighbor[3][2] = { { false, false },{ false, false },{ false, false } };

					// X rays
					for (int z = 0; z < BLOCK_DIMENSION; z++)
					{
						double rayPosZ = (zc + zi)*scale + z*voxelSize;
						for (int y = 0; y < BLOCK_DIMENSION; y++)
						{
							// Set up ray information
							double rayPosY = (yc + yi)*scale + y*voxelSize;

							double zymin[2] =
							{
								rayPosY - cellMargin,
								rayPosZ - cellMargin
							};
							double zymax[2] =
							{
								rayPosY + cellMargin + voxelSize,
								rayPosZ + cellMargin + voxelSize
							};

							context.dir[0] = 1.0;
							context.dir[1] = 0.0;
							context.dir[2] = 0.0;
							context.origin[0] = minX;
							context.origin[1] = rayPosY;
							context.origin[2] = rayPosZ;
							context.axis0 = 0;
							context.axis1 = 1;
							context.axis2 = 2;
							context.axisMin = minX;
							for (int i = 0; i < RAYS_SHOT; i++)
							{
								context.rayTestCount[i] = 0;
							}

							// Search the plane index
							yz.Search(zymin, zymax, planeIntersectCallback, (void*)&context);

							for (int r = 0; r < RAYS_SHOT; r++)
							{
								// If no intersections were found, skip to the next ray
								if (context.rayTestCount[r] == 0)
								{
									continue;
								}

								// Counts how many times the ray has entered a volume
								// When a volume is exited, it counts backwards
								bool inside = false;

								// Tracks the position for the last intersection detected
								int lastX = 0;
								bool firstIntersectionX = true;
								bool prevFacingX = false;
								bool initialFacingX = true;

								// Iterate over each intersection
								for (StampRayTest* ray = &rayTests[r][context.rayTestHead[r]];
								ray != NULL;
									ray = (ray->next != -1) ? &rayTests[r][ray->next] : NULL)
								{
									StampRayTest& intersection = *ray;

									MaterialId materialId = mesh->getSolidMaterial(solid);
									if (materials != NULL)
									{
										double voxelWorldPos[3] = { minX + intersection.pos, context.origin[1], context.origin[2] };
										materialId = materials->translateMaterial(voxelWorldPos, materialId);
									}
									if (materialId != 0)
									{
										inserting = true;
									}

									// Compute intersection point
									double pvx = (minX + intersection.pos - (xc + xi)*scale) / voxelSize;
									int nextX = (pvx >= 0) ? (int)pvx : (int)pvx - 1;
									bool difRayX = firstIntersectionX || (prevFacingX != intersection.facing);
									if (firstIntersectionX)
									{
										initialFacingX = intersection.facing;
									}
									firstIntersectionX = false;

									if ((nextX >= 0) && (nextX < BLOCK_DIMENSION))
									{
										double pvy = RAY_DELTA[r][0];
										double pvz = RAY_DELTA[r][1];

										pvx = pvx - nextX;
										addIntersectionToIndex(
											intersectionCount,
											intersectionListCount,
											intersectionIndex,
											intersections,
											intersectionLists,
											pvx, pvy, pvz,
											intersection.normal[0],
											intersection.normal[1],
											intersection.normal[2],
											nextX, y, z,
											materialId);
									}

									if (r == 0)
									{
										if (difRayX)
										{
											inside = intersection.facing == initialFacingX;
											prevFacingX = intersection.facing;
											if (inside)
											{
												lastX = nextX + 1;
											}
										}
										if (!inside)
										{
											// If inside, iterate to the next intersection and set all
											// voxels in-between as solid
											for (int x = std::max(0, lastX); x <= std::min(BLOCK_DIMENSION - 1, nextX); x++)
											{
												BlockVoxelData::Index voxelIndex2 = BlockVoxelData::Index(x, y, z);
												unsigned char matCount = (voxelMaterials[voxelIndex2] & 0xFF) | 1;

												voxelMaterials[voxelIndex2] = (materialId << 8) | matCount;
											}

											lastX = std::min(BLOCK_DIMENSION - 1, nextX) + 1;
										}
									}
								}
							}
						}
					}

					// Y rays. Same as X rays, but no need to set voxels to solid
					for (int z = 0; z < BLOCK_DIMENSION; z++)
					{
						double rayPosZ = (zc + zi)*scale + z*voxelSize;
						for (int x = 0; x < BLOCK_DIMENSION; x++)
						{
							double rayPosX = (xc + xi)*scale + x*voxelSize;

							double zxmin[2] =
							{
								rayPosX - cellMargin,
								rayPosZ - cellMargin
							};
							double zxmax[2] =
							{
								rayPosX + cellMargin + voxelSize,
								rayPosZ + cellMargin + voxelSize
							};

							context.dir[0] = 0.0;
							context.dir[1] = 1.0;
							context.dir[2] = 0.0;
							context.origin[0] = rayPosX;
							context.origin[1] = minY;
							context.origin[2] = rayPosZ;
							context.axis0 = 1;
							context.axis1 = 0;
							context.axis2 = 2;
							context.axisMin = minY;
							for (int i = 0; i < RAYS_SHOT; i++)
							{
								context.rayTestCount[i] = 0;
							}

							xz.Search(zxmin, zxmax, planeIntersectCallback, (void*)&context);

							for (int r = 0; r < RAYS_SHOT; r++)
							{
								if (context.rayTestCount[r] == 0)
								{
									continue;
								}

								bool inside = false;
								int lastY = 0;
								bool firstIntersectionY = true;
								bool prevFacingY = false;
								bool initialFacingY = true;

								for (StampRayTest* ray = &rayTests[r][context.rayTestHead[r]];
								ray != NULL;
									ray = (ray->next != -1) ? &rayTests[r][ray->next] : NULL)
								{
									StampRayTest& intersection = *ray;
									MaterialId materialId = mesh->getSolidMaterial(solid);
									if (materials != NULL)
									{
										double voxelWorldPos[3] = { context.origin[0], minY + intersection.pos, context.origin[2] };
										materialId = materials->translateMaterial(voxelWorldPos, materialId);
									}
									if (materialId != 0)
									{
										inserting = true;
									}

									double pvy = (minY + intersection.pos - (yc + yi)*scale) / voxelSize;
									int nextY = (pvy >= 0) ? (int)pvy : (int)pvy - 1;
									bool difRayY = firstIntersectionY || (prevFacingY != intersection.facing);
									if (firstIntersectionY)
									{
										initialFacingY = intersection.facing;
									}
									firstIntersectionY = false;

									if (nextY >= 0 && nextY < BLOCK_DIMENSION)
									{
										double pvx = RAY_DELTA[r][0];
										double pvz = RAY_DELTA[r][1];
										pvy = pvy - nextY;

										addIntersectionToIndex(
											intersectionCount,
											intersectionListCount,
											intersectionIndex,
											intersections,
											intersectionLists,
											pvx, pvy, pvz,
											intersection.normal[0],
											intersection.normal[1],
											intersection.normal[2],
											x, nextY, z,
											materialId);
									}

									if (r == 0)
									{
										if (difRayY)
										{
											inside = (intersection.facing == initialFacingY);
											prevFacingY = intersection.facing;
											if (inside)
											{
												lastY = nextY + 1;
											}
										}
										if (!inside)
										{
											// If inside, iterate to the next intersection and set all
											// voxels in-between as solid
											for (int y = std::max(0, lastY); y <= std::min(BLOCK_DIMENSION - 1, nextY); y++)
											{
												BlockVoxelData::Index voxelIndex2 = BlockVoxelData::Index(x, y, z);
												unsigned char matCount = (voxelMaterials[voxelIndex2] & 0xFF) | 2;

												voxelMaterials[voxelIndex2] = (materialId << 8) | matCount;
											}

											lastY = std::min(BLOCK_DIMENSION - 1, nextY) + 1;
										}
									}
								}
							}
						}
					}

					// Z rays. Same as X rays, but no need to set voxels to solid
					for (int x = 0; x < BLOCK_DIMENSION; x++)
					{
						double rayPosX = (xc + xi)*scale + x*voxelSize;
						for (int y = 0; y < BLOCK_DIMENSION; y++)
						{
							double rayPosY = (yc + yi)*scale + y*voxelSize;

							double xymin[2] =
							{
								rayPosX - cellMargin,
								rayPosY - cellMargin
							};
							double xymax[2] =
							{
								rayPosX + cellMargin + voxelSize,
								rayPosY + cellMargin + voxelSize
							};
							context.dir[0] = 0.0;
							context.dir[1] = 0.0;
							context.dir[2] = 1.0;
							context.origin[0] = rayPosX;
							context.origin[1] = rayPosY;
							context.origin[2] = minZ;
							context.axis0 = 2;
							context.axis1 = 0;
							context.axis2 = 1;
							context.axisMin = minZ;
							for (int i = 0; i < RAYS_SHOT; i++)
							{
								context.rayTestCount[i] = 0;
							}

							xy.Search(xymin, xymax, planeIntersectCallback, (void*)&context);

							for (int r = 0; r < RAYS_SHOT; r++)
							{
								if (context.rayTestCount[r] == 0)
								{
									continue;
								}

								bool inside = false;
								int lastZ = 0;
								bool firstIntersectionZ = true;
								bool prevFacingZ = false;
								bool initialFacingZ = true;

								for (StampRayTest* ray = &rayTests[r][context.rayTestHead[r]];
								ray != NULL;
									ray = (ray->next != -1) ? &rayTests[r][ray->next] : NULL)
								{
									StampRayTest& intersection = *ray;
									MaterialId materialId = mesh->getSolidMaterial(solid);
									if (materials != NULL)
									{
										double voxelWorldPos[3] = { context.origin[0], context.origin[1], minZ + intersection.pos };
										materialId = materials->translateMaterial(voxelWorldPos, materialId);
									}
									if (materialId != 0)
									{
										inserting = true;
									}

									double pvz = (minZ + intersection.pos - (zc + zi)*scale) / voxelSize;
									int nextZ = (pvz >= 0) ? (int)pvz : (int)pvz - 1;
									bool difRayZ = firstIntersectionZ || (prevFacingZ != intersection.facing);
									if (firstIntersectionZ)
									{
										initialFacingZ = intersection.facing;
									}
									firstIntersectionZ = false;

									if ((nextZ >= 0) && (nextZ < BLOCK_DIMENSION))
									{
										double pvx = RAY_DELTA[r][0];
										double pvy = RAY_DELTA[r][1];
										pvz = pvz - nextZ;

										addIntersectionToIndex(
											intersectionCount,
											intersectionListCount,
											intersectionIndex,
											intersections,
											intersectionLists,
											pvx, pvy, pvz,
											intersection.normal[0],
											intersection.normal[1],
											intersection.normal[2],
											x, y, nextZ,
											materialId);
									}

									if (r == 0)
									{
										if (difRayZ)
										{
											inside = intersection.facing == initialFacingZ;
											prevFacingZ = intersection.facing;
											if (inside)
											{
												lastZ = nextZ + 1;
											}
										}
										if (!inside)
										{
											// If inside, iterate to the next intersection and set all
											// voxels in-between as solid
											for (int z = std::max(0, lastZ); z <= std::min(BLOCK_DIMENSION - 1, nextZ); z++)
											{
												BlockVoxelData::Index voxelIndex2 = BlockVoxelData::Index(x, y, z);
												unsigned char matCount = (voxelMaterials[voxelIndex2] & 0xFF) | 4;

												voxelMaterials[voxelIndex2] = (materialId << 8) | matCount;
											}

											lastZ = std::min(BLOCK_DIMENSION - 1, nextZ) + 1;
										}
									}
								}
							}
						}
					}

					//set vectors
					// Solve the QEF for each voxel
					// Iterate over each list of intersections
					const int SOLVER_MAX_POINTS = 64;
					double matrix[SOLVER_MAX_POINTS][3];
					double vector[SOLVER_MAX_POINTS];
					bool validTriangles[SOLVER_MAX_POINTS];
					double lowValue = VECTOR_MIN_VALUE;
					double highValue = VECTOR_MAX_VALUE;

					memset(overflow, 0, BLOCK_DATA_SIZE*sizeof(CompactVoxel));

					for (int idx = 0; idx < intersectionListCount; idx++)
					{
						// Retrieve list
						RayIntersectionList& list = intersectionLists[idx];

						// Get number of points in the list.
						// Make sure they do not exceed the number of points the QEF solver can take.
						int count = std::min(SOLVER_MAX_POINTS, list.count);

						// If no points were found, bail out
						if (count == 0)
						{
							continue;
						}

						// Compute voxel index
						int userDataCount = 0;

						BlockVoxelData::Index voxelIndex = BlockVoxelData::Index(list.x, list.y, list.z);
						if (count < SOLVER_MAX_POINTS)
						{
							userDataCount = getVoxelIntersections(userDataPoints, blockData, buffer, currentCell, list.x, list.y, list.z);
						}

						int userDataTriangles = (int)(userDataCount * 2 / 3);
						if (count + userDataTriangles > SOLVER_MAX_POINTS)
						{
							userDataTriangles = 0;
						}

						if (count + userDataCount == 1)
						{
							// If only one intersection, no need for QEF (also it cannot be solved)
							// Get intersection and use it for voxel data
							RayIntersection& ri = intersections[list.head];

							// Get intersection coordinates within voxel
							float px = (float)std::min(highValue, std::max(lowValue, (double)ri.position.x));
							float py = (float)std::min(highValue, std::max(lowValue, (double)ri.position.y));
							float pz = (float)std::min(highValue, std::max(lowValue, (double)ri.position.z));

							voxelData->setVector(voxelIndex, px, py, pz);
							overflow[voxelIndex] = 2;

							changed = true;
							if (list.x < BLOCK_MARGIN)
							{
								changedNeighbor[0][0] |= true;
							}
							else if (list.x >= BLOCK_DIMENSION - BLOCK_MARGIN)
							{
								changedNeighbor[0][1] |= true;
							}
							if (list.y < BLOCK_MARGIN)
							{
								changedNeighbor[1][0] |= true;
							}
							else if (list.y >= BLOCK_DIMENSION - BLOCK_MARGIN)
							{
								changedNeighbor[1][1] |= true;
							}
							if (list.z < BLOCK_MARGIN)
							{
								changedNeighbor[2][0] |= true;
							}
							else if (list.z >= BLOCK_DIMENSION - BLOCK_MARGIN)
							{
								changedNeighbor[2][1] |= true;
							}

							// Continue to the next list
							continue;
						}

						//normalize normals
						RayIntersection* ri = &intersections[list.head];
						for (int i = 0; i < count; i++)
						{
							Vector& normal = ri->normal;
							Vector_normalize(&normal);

							matrix[i][0] = normal.x;
							matrix[i][1] = normal.y;
							matrix[i][2] = normal.z;

							validTriangles[i] = true;

							if (ri->next != -1)
							{
								ri = &intersections[ri->next];
							}
						}

						// Load normals into arrays for solver
						int total = count;

						if (userDataTriangles > 0)
						{
							double vx, vy, vz;
							buffer->getVector(voxelIndex, vx, vy, vz);
							HyVoxel::Vector A = HyVoxel::Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);

							for (int userIndex = 0; userIndex < userDataCount; userIndex += 3)
							{
								HyVoxel::Vector B = userDataPoints[userIndex];
								HyVoxel::Vector C = userDataPoints[userIndex + 1];
								HyVoxel::Vector D = userDataPoints[userIndex + 2];

								HyVoxel::Vector ab = HyVoxel::Algebra::Vector_subtract(B, A);
								HyVoxel::Vector ac = HyVoxel::Algebra::Vector_subtract(C, A);
								HyVoxel::Vector ad = HyVoxel::Algebra::Vector_subtract(D, A);

								HyVoxel::Vector n = HyVoxel::Algebra::Vector_cross(ab, ac);
								validTriangles[total] = (abs(n.x) > 0.0001) || (abs(n.y) > 0.0001) || (abs(n.z) > 0.0001);

								HyVoxel::Algebra::Vector_normalize(&n);
								matrix[total][0] = n.x;
								matrix[total][1] = n.y;
								matrix[total][2] = n.z;
								total++;

								n = HyVoxel::Algebra::Vector_cross(ac, ad);
								validTriangles[total] = (abs(n.x) > 0.0001) || (abs(n.y) > 0.0001) || (abs(n.z) > 0.0001);

								HyVoxel::Algebra::Vector_normalize(&n);
								matrix[total][0] = n.x;
								matrix[total][1] = n.y;
								matrix[total][2] = n.z;
								total++;
							}
						}

						double QEF_error = 0.05;

						// Mass point for the voxel intersections
						double generalMassPoint[3] = { 0.0, 0.0, 0.0 };
						float originOffset[3] = { 0.5f, 0.5f, 0.5f };

						// Load normals and points into arrays for solver
						ri = &intersections[list.head];
						for (int i = 0; i < count; i++)
						{
							Vector pt = ri->position;
							Vector point = HyVoxel::Algebra::Vector_withValues(pt.x - originOffset[0], pt.y - originOffset[1], pt.z - originOffset[2]);
							Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
							vector[i] = Vector_dot(normal, point);

							if (ri->next != -1)
							{
								ri = &intersections[ri->next];
							}
						}

						int planeCount = count;
						if (userDataTriangles > 0)
						{
							double vx, vy, vz;
							buffer->getVector(voxelIndex, vx, vy, vz);
							vx -= originOffset[0];
							vy -= originOffset[1];
							vz -= originOffset[2];
							Vector point = HyVoxel::Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);

							for (int i = count; i < userDataTriangles + count; i++)
							{
								Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
								vector[i] = Vector_dot(normal, point);
							}

							//remove triangles with normal 0
							int idx = count;
							while (idx < total)
							{
								if (validTriangles[idx])
								{
									if (planeCount != idx)
									{
										matrix[planeCount][0] = matrix[idx][0];
										matrix[planeCount][1] = matrix[idx][1];
										matrix[planeCount][2] = matrix[idx][2];
										vector[planeCount] = vector[idx];
									}

									planeCount++;
								}

								idx++;
							}
						}
						bool resolvedQEF = false;

						double px, py, pz;
						QEF::evaluate(matrix, vector, planeCount, px, py, pz);
						px += originOffset[0];
						py += originOffset[1];
						pz += originOffset[2];

						if ((px < 0.0 && px >= -0.3) || (px > 1.0 && px <= 1.3) ||
							(py < 0.0 && py >= -0.3) || (py > 1.0 && py <= 1.3) ||
							(pz < 0.0 && pz >= -0.3) || (pz > 1.0 && pz <= 1.3))
						{
							/*
							try to correct the problem
							take the closest corner to the bad QEF solution point as a new mass point and calculate the QEF again.
							*/
							float massPt[3] = { (px > 0.5) ? 1.0f : 0.0f, (py > 0.5) ? 1.0f : 0.0f, (pz > 0.5) ? 1.0f : 0.0f };

							//run the QEF again
							// Load normals and points into arrays for solver
							ri = &intersections[list.head];
							for (int i = 0; i < count; i++)
							{
								Vector pt = ri->position;
								Vector point = HyVoxel::Algebra::Vector_withValues(pt.x - massPt[0], pt.y - massPt[1], pt.z - massPt[2]);
								Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
								vector[i] = Vector_dot(normal, point);

								if (ri->next != -1)
								{
									ri = &intersections[ri->next];
								}
							}

							int planeCount = count;
							if (userDataTriangles > 0)
							{
								double vx, vy, vz;
								buffer->getVector(voxelIndex, vx, vy, vz);
								vx -= massPt[0];
								vy -= massPt[1];
								vz -= massPt[2];
								Vector point = HyVoxel::Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);

								for (int i = count; i < userDataTriangles + count; i++)
								{
									Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
									vector[i] = Vector_dot(normal, point);
								}

								//remove triangles with normal 0
								int idx = count;
								while (idx < total)
								{
									if (validTriangles[idx])
									{
										if (planeCount != idx)
										{
											matrix[planeCount][0] = matrix[idx][0];
											matrix[planeCount][1] = matrix[idx][1];
											matrix[planeCount][2] = matrix[idx][2];
											vector[planeCount] = vector[idx];
										}

										planeCount++;
									}

									idx++;
								}
							}

							QEF::evaluate(matrix, vector, total, px, py, pz);
							px += massPt[0];
							py += massPt[1];
							pz += massPt[2];
						}

						if (px >= lowValue - QEF_error && px <= highValue + QEF_error &&
							py >= lowValue - QEF_error && py <= highValue + QEF_error &&
							pz >= lowValue - QEF_error && pz <= highValue + QEF_error)
						{
							resolvedQEF = true;
						}

						// Compute mass point
						ri = &intersections[list.head];
						for (int i = 0; i < count; i++)
						{
							Vector& point = ri->position;

							generalMassPoint[0] += point.x;
							generalMassPoint[1] += point.y;
							generalMassPoint[2] += point.z;

							if (ri->next != -1)
							{
								ri = &intersections[ri->next];
							}
						}

						generalMassPoint[0] /= count;
						generalMassPoint[1] /= count;
						generalMassPoint[2] /= count;

						if (!resolvedQEF)
						{
							//calculate the vector without the user data, only the mesh
							if (userDataTriangles > 0)
							{
								MaterialId airMat = VOXEL_ONLY_COORDS;
								if (inserting)
								{
									airMat = 0;
									userDataCount = getVoxelIntersections(userDataPoints, blockData, buffer, currentCell, list.x, list.y, list.z, false);
								}

								int voxelOverflow = 0;
								double npx = 0.5;
								double npy = 0.5;
								double npz = 0.5;
								bool resolved = false;

								if (count == 1)
								{
									// If only one intersection, no need for QEF (also it cannot be solved)
									// Get intersection and use it for voxel data
									RayIntersection& ri = intersections[list.head];

									// Get intersection coordinates within voxel
									npx = (float)std::min(highValue, std::max(lowValue, (double)ri.position.x));
									npy = (float)std::min(highValue, std::max(lowValue, (double)ri.position.y));
									npz = (float)std::min(highValue, std::max(lowValue, (double)ri.position.z));
									resolved = true;
									voxelOverflow = 2;
								}
								else
								{
									// Load normals and points into arrays for solver
									ri = &intersections[list.head];
									for (int i = 0; i < count; i++)
									{
										Vector pt = ri->position;
										Vector point = HyVoxel::Algebra::Vector_withValues(pt.x - originOffset[0], pt.y - originOffset[1], pt.z - originOffset[2]);
										Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
										vector[i] = Vector_dot(normal, point);

										if (ri->next != -1)
										{
											ri = &intersections[ri->next];
										}
									}

									QEF::evaluate(matrix, vector, count, npx, npy, npz);
									npx += originOffset[0];
									npy += originOffset[1];
									npz += originOffset[2];

									if ((npx < 0.0 && npx >= -0.3) || (npx > 1.0 && npx <= 1.3) ||
										(npy < 0.0 && npy >= -0.3) || (npy > 1.0 && npy <= 1.3) ||
										(npz < 0.0 && npz >= -0.3) || (npz > 1.0 && npz <= 1.3))
									{
										/*
										try to correct the problem
										take the closest corner to the bad QEF solution point as a new mass point and calculate the QEF again.
										*/
										float massPt[3] = { (npx > 0.5) ? 1.0f : 0.0f, (npy > 0.5) ? 1.0f : 0.0f, (npz > 0.5) ? 1.0f : 0.0f };

										//run the QEF again
										// Load normals and points into arrays for solver
										ri = &intersections[list.head];
										for (int i = 0; i < count; i++)
										{
											Vector pt = ri->position;
											Vector point = HyVoxel::Algebra::Vector_withValues(pt.x - massPt[0], pt.y - massPt[1], pt.z - massPt[2]);
											Vector normal = HyVoxel::Algebra::Vector_withValues((float)matrix[i][0], (float)matrix[i][1], (float)matrix[i][2]);
											vector[i] = Vector_dot(normal, point);

											if (ri->next != -1)
											{
												ri = &intersections[ri->next];
											}
										}

										QEF::evaluate(matrix, vector, count, npx, npy, npz);
										npx += massPt[0];
										npy += massPt[1];
										npz += massPt[2];
									}

									if (npx >= lowValue - QEF_error && npx <= highValue + QEF_error &&
										npy >= lowValue - QEF_error && npy <= highValue + QEF_error &&
										npz >= lowValue - QEF_error && npz <= highValue + QEF_error)
									{
										resolved = true;
									}
								}

								if (!resolved)
								{
									voxelOverflow = 2;

									npx = generalMassPoint[0];
									npy = generalMassPoint[1];
									npz = generalMassPoint[2];
								}
								else
								{
									if (npx < -QEF_error || npx > 1.0 + QEF_error ||
										npy < -QEF_error || npy > 1.0 + QEF_error ||
										npz < -QEF_error || npz > 1.0 + QEF_error)
									{
										npx = std::min(highValue, std::max(lowValue, npx));
										npy = std::min(highValue, std::max(lowValue, npy));
										npz = std::min(highValue, std::max(lowValue, npz));

										voxelOverflow = 1;
										int ddx = encodeVectorCoord(generalMassPoint[0]);
										int ddy = encodeVectorCoord(generalMassPoint[1]);
										int ddz = encodeVectorCoord(generalMassPoint[2]);
										overflow[voxelIndex] |= (ddx << 8) | (ddy << 16) | (ddz << 24);
									}
									else
									{
										voxelOverflow = 2;
										npx = std::min(1.0, std::max(0.0, npx));
										npy = std::min(1.0, std::max(0.0, npy));
										npz = std::min(1.0, std::max(0.0, npz));
									}
								}

								if (userDataCount > 0)
								{
									//find out if the dot is inside the volumen
									int position = -1;

									double vx, vy, vz;
									buffer->getVector(voxelIndex, vx, vy, vz);

									HyVoxel::Vector p0 = HyVoxel::Algebra::Vector_withValues((float)npx, (float)npy, (float)npz);
									HyVoxel::Vector v0 = HyVoxel::Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);
									HyVoxel::Vector point;

									bool positionResolved = false;
									for (int corner = 0; (corner < 8) && !positionResolved; corner++)
									{
										HyVoxel::Vector p1 = HyVoxel::Algebra::Vector_withValues((float)VoxelCorners[corner][0], (float)VoxelCorners[corner][1],
											(float)VoxelCorners[corner][2]);

										//check is the point is close to a corner
										Vector p1p0 = HyVoxel::Algebra::Vector_subtract(p1, p0);
										float mag = HyVoxel::Algebra::Vector_magnitudeSquared(p1p0);
										if (mag > 0.0025)
										{
											bool intersect = false;
											int triang = -1;
											for (int userIndex = 0; userIndex < userDataCount; userIndex += 3)
											{
												HyVoxel::Vector v1 = userDataPoints[userIndex];
												HyVoxel::Vector v2 = userDataPoints[userIndex + 1];

												int result = HyVoxel::Algebra::segmentIntersectTriangle3D(p0, p1, v0, v1, v2, &point);
												if (result == 1)
												{
													Vector p0Point = HyVoxel::Algebra::Vector_subtract(point, p0);
													float mag = HyVoxel::Algebra::Vector_magnitudeSquared(p0Point);
													if (mag > 0.00001)
													{
														triang = userIndex;
														intersect = true;
													}
													else
													{
														position = 0;
														positionResolved = true;
													}

													break;
												}
												if ((result == 2) || (result == 3))
												{
													position = 0;
													positionResolved = true;
													intersect = true;
													break;
												}

												HyVoxel::Vector v3 = userDataPoints[userIndex + 2];
												result = HyVoxel::Algebra::segmentIntersectTriangle3D(p0, p1, v0, v2, v3, &point);
												if (result == 1)
												{
													Vector p0Point = HyVoxel::Algebra::Vector_subtract(point, p0);
													float mag = HyVoxel::Algebra::Vector_magnitudeSquared(p0Point);
													if (mag > 0.00001)
													{
														triang = userIndex + 1;
														intersect = true;
													}
													else
													{
														position = 0;
														positionResolved = true;
													}

													break;
												}
												if ((result == 2) || (result == 3))
												{
													position = 0;
													positionResolved = true;
													intersect = true;
													break;
												}
											}

											if (intersect && !positionResolved)
											{
												//calculate normal
												Vector A = Algebra::Vector_withValues((float)vx, (float)vy, (float)vz);
												Vector B = userDataPoints[triang];
												Vector C = userDataPoints[triang + 1];

												Vector ab = Algebra::Vector_subtract(B, A);
												Vector ac = Algebra::Vector_subtract(C, A);

												Vector normal = Algebra::Vector_cross(ac, ab);

												double ang = Algebra::Vector_dot(normal, p1p0);
												if (ang < 0.0)
												{
													position = -1;
												}
												else
												{
													position = 1;
												}

												positionResolved = true;
											}
										}
										else
										{
											//the point is very close to a corner, so take the material in the corner
											//if there is material the point is inside the volumen
											//calculate material
											CellId nCell = currentCell;

											// Compute voxel coordinates
											int nx = list.x + VoxelCorners[corner][0];
											int ny = list.y + VoxelCorners[corner][1];
											int nz = list.z + VoxelCorners[corner][2];

											// Adjust cell ID and voxel coordinates so they are consistent
											adjustCellCoordinates(nCell, nx, ny, nz);
											BlockVoxelData::Index nVoxelIndex = BlockVoxelData::Index(nx, ny, nz);

											MaterialId nMat = airMat;
											BlockVoxelData* nVoxelData = buffer;
											if (nCell != currentCell)
											{
												nVoxelData = blockData->fetchData(nCell, false);
											}
											if ((nVoxelData != NULL) && nVoxelData->hasType(nVoxelIndex, VOXEL_HAS_MATERIAL))
											{
												nMat = nVoxelData->getMaterial(nVoxelIndex);
											}

											position = (nMat == airMat) ? -1 : 1;
											positionResolved = true;
										}
									}

									if (!positionResolved)
									{
										//take material of current voxel
										//if there is material the point is inside the volumen
										MaterialId mat = airMat;
										if (buffer->hasType(voxelIndex, VOXEL_HAS_MATERIAL))
										{
											mat = buffer->getMaterial(voxelIndex);
										}

										position = (mat == airMat) ? -1 : 1;
										positionResolved = true;
									}

									if (position <= 0) //outside
									{
										if (inserting) //inserting
										{
											px = npx;
											py = npy;
											pz = npz;
											overflow[voxelIndex] = voxelOverflow;
										}
										else //deleting
										{
											px = vx;
											py = vy;
											pz = vz;
											overflow[voxelIndex] = 2;
										}
									}
									else //inside
									{
										if (inserting) //inserting
										{
											px = vx;
											py = vy;
											pz = vz;
											overflow[voxelIndex] = 2;
										}
										else //deleting
										{
											px = npx;
											py = npy;
											pz = npz;
											overflow[voxelIndex] = voxelOverflow;
										}
									}
								}
								else
								{
									px = npx;
									py = npy;
									pz = npz;
									overflow[voxelIndex] = voxelOverflow;
								}
							}
							else
							{
								overflow[voxelIndex] = 3;
								px = generalMassPoint[0];
								py = generalMassPoint[1];
								pz = generalMassPoint[2];
							}
						}
						else
						{
							if (px < -QEF_error || px > 1.0 + QEF_error ||
								py < -QEF_error || py > 1.0 + QEF_error ||
								pz < -QEF_error || pz > 1.0 + QEF_error)
							{
								px = std::min(highValue, std::max(lowValue, px));
								py = std::min(highValue, std::max(lowValue, py));
								pz = std::min(highValue, std::max(lowValue, pz));

								overflow[voxelIndex] = 1;
								int ddx = encodeVectorCoord(generalMassPoint[0]);
								int ddy = encodeVectorCoord(generalMassPoint[1]);
								int ddz = encodeVectorCoord(generalMassPoint[2]);
								overflow[voxelIndex] |= (ddx << 8) | (ddy << 16) | (ddz << 24);
							}
							else
							{
								overflow[voxelIndex] = 2;
								px = std::min(1.0, std::max(0.0, px));
								py = std::min(1.0, std::max(0.0, py));
								pz = std::min(1.0, std::max(0.0, pz));
							}
						}

						voxelData->setVector(voxelIndex, px, py, pz);
						changed = true;
						if (list.x < BLOCK_MARGIN)
						{
							changedNeighbor[0][0] = true;
						}
						else if (list.x >= BLOCK_DIMENSION - BLOCK_MARGIN)
						{
							changedNeighbor[0][1] = true;
						}
						if (list.y < BLOCK_MARGIN)
						{
							changedNeighbor[1][0] = true;
						}
						else if (list.y >= BLOCK_DIMENSION - BLOCK_MARGIN)
						{
							changedNeighbor[1][1] = true;
						}
						if (list.z < BLOCK_MARGIN)
						{
							changedNeighbor[2][0] = true;
						}
						else if (list.z >= BLOCK_DIMENSION - BLOCK_MARGIN)
						{
							changedNeighbor[2][1] = true;
						}
					}

#ifdef VECTOR_OVERFLOW
					//detect errors with vector overflow
					for (int idx = 0; idx < intersectionListCount; idx++)
					{
						// Retrieve list
						RayIntersectionList& list = intersectionLists[idx];
						BlockVoxelData::Index vIndex = BlockVoxelData::Index(list.x, list.y, list.z);

						if ((overflow[vIndex] & 0xFF) == 1)
						{
							double vx, vy, vz;
							voxelData->getVector(vIndex, vx, vy, vz);

							int nx = list.x;
							if (vx < 0.0)
							{
								nx--;
							}
							if (vx > 1.0)
							{
								nx++;
							}

							int ny = list.y;
							if (vy < 0.0)
							{
								ny--;
							}
							if (vy > 1.0)
							{
								ny++;
							}

							int nz = list.z;
							if (vz < 0.0)
							{
								nz--;
							}
							if (vz > 1.0)
							{
								nz++;
							}

							if (nx >= 0 && nx < BLOCK_DIMENSION && ny >= 0 && ny < BLOCK_DIMENSION && nz >= 0 && nz < BLOCK_DIMENSION)
							{
								BlockVoxelData::Index destVoxelIndex = BlockVoxelData::Index(nx, ny, nz);
								if ((overflow[destVoxelIndex] & 0xFF) == 0)
								{
									double ptx = decodeVectorCoord((overflow[vIndex] & 0x0000FF00) >> 8);
									double pty = decodeVectorCoord((overflow[vIndex] & 0x00FF0000) >> 16);
									double ptz = decodeVectorCoord((overflow[vIndex] & 0xFF000000) >> 24);

									voxelData->setVector(vIndex, ptx, pty, ptz);
								}
							}
						}
					}
#endif

					//set materials
					for (int x = 0; x < BLOCK_DIMENSION; x++)
					{
						for (int y = 0; y < BLOCK_DIMENSION; y++)
						{
							for (int z = 0; z < BLOCK_DIMENSION; z++)
							{
								BlockVoxelData::Index voxelIndex2 = BlockVoxelData::Index(x, y, z);
								unsigned char matCount = (voxelMaterials[voxelIndex2] & 0xFF);

								if (matCount != 0 && matCount != 1 && matCount != 2 && matCount != 4)
								{
									unsigned short mat = ((voxelMaterials[voxelIndex2] >> 8) & 0xFFFF);
									if (mat != VOXEL_ONLY_COORDS)
										voxelData->setMaterial(voxelIndex2, mat);
									else
									{
										voxelData->removeType(voxelIndex2, VOXEL_HAS_MATERIAL);
									}

									if (intersectionIndex[voxelIndex2] == -1)
									{
										voxelData->removeType(voxelIndex2, VOXEL_HAS_VECTOR);
									}

									// buffer these voxels to compute physics meshes
									if (physicsRequest && voxelCache->count(currentCell) != 0)
									{
										int nidx = (z % noiseDim)*noiseDim*noiseDim + (x % noiseDim)*noiseDim + (y % noiseDim);
										unsigned short nmask = noise->field[nidx];
										int voxPos[3] = { x, y, z };

										if (nmask > 0 && hasVoxelOccupancy(voxelCache->at(currentCell), voxPos))
										{
											Voxel curVox;
											voxelData->getVoxel(voxelIndex2, curVox);
											Physics::CVoxelBuffer::nVoxel curnVox;

											// voxel vector
											if (curVox.vx > 0 || curVox.vy > 0 || curVox.vz > 0)
											{
												#ifdef VECTOR_OVERFLOW
												curnVox.v[0] = (unsigned char)(255.f * decodeVectorCoord(curVox.vx));
												curnVox.v[1] = (unsigned char)(255.f * decodeVectorCoord(curVox.vy));
												curnVox.v[2] = (unsigned char)(255.f * decodeVectorCoord(curVox.vz));
												#else
												curnVox.v[0] = curVox.vx;
												curnVox.v[1] = curVox.vy;
												curnVox.v[2] = curVox.vz;
												#endif
											}
											else
											{
												curnVox.v[0] = 128;
												curnVox.v[1] = 128;
												curnVox.v[2] = 128;
											}

											// voxel type
											curnVox.type = nmask + ncellCount[nidx] * noiseCells;

											// update voxel range buffer
											if (nVoxBuffer->voxelRanges[curnVox.type].cellId < 0)
											{
												nVoxBuffer->voxelRanges[curnVox.type].cellId = ncellSpanned * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE;
												nVoxBuffer->voxelRanges[curnVox.type].pos[0] = (xc + xi)*scale;
												nVoxBuffer->voxelRanges[curnVox.type].pos[1] = (yc + yi)*scale;
												nVoxBuffer->voxelRanges[curnVox.type].pos[2] = (zc + zi)*scale;
											}
											nVoxBuffer->voxelRanges[curnVox.type].minx = std::min(x, nVoxBuffer->voxelRanges[curnVox.type].minx);
											nVoxBuffer->voxelRanges[curnVox.type].maxx = std::max(x, nVoxBuffer->voxelRanges[curnVox.type].maxx);
											nVoxBuffer->voxelRanges[curnVox.type].miny = std::min(y, nVoxBuffer->voxelRanges[curnVox.type].miny);
											nVoxBuffer->voxelRanges[curnVox.type].maxy = std::max(y, nVoxBuffer->voxelRanges[curnVox.type].maxy);
											nVoxBuffer->voxelRanges[curnVox.type].minz = std::min(z, nVoxBuffer->voxelRanges[curnVox.type].minz);
											nVoxBuffer->voxelRanges[curnVox.type].maxz = std::max(z, nVoxBuffer->voxelRanges[curnVox.type].maxz);

											// update voxel type buffer
											nVoxBuffer->voxelIds[(ncellSpanned * BLOCK_SIZE * BLOCK_SIZE * BLOCK_SIZE + z* BLOCK_SIZE * BLOCK_SIZE + x * BLOCK_SIZE + y)] = curnVox;
											ncellCount[nidx]++;
										}
									}

									changed = true;
									if (x < BLOCK_MARGIN)
									{
										changedNeighbor[0][0] = true;
									}
									else if (x >= BLOCK_DIMENSION - BLOCK_MARGIN)
									{
										changedNeighbor[0][1] = true;
									}
									if (y < BLOCK_MARGIN)
									{
										changedNeighbor[1][0] = true;
									}
									else if (y >= BLOCK_DIMENSION - BLOCK_MARGIN)
									{
										changedNeighbor[1][1] = true;
									}
									if (z < BLOCK_MARGIN)
									{
										changedNeighbor[2][0] = true;
									}
									else if (z >= BLOCK_DIMENSION - BLOCK_MARGIN)
									{
										changedNeighbor[2][1] = true;
									}
								}
							}
						}
					}

					if (changed && (changedCells != NULL))
					{
						for (int nx = -1; nx < 2; nx++)
						{
							bool changeX = true;
							if (nx == -1)
							{
								changeX = changedNeighbor[0][0];
							}
							else if (nx == 1)
							{
								changeX = changedNeighbor[0][1];
							}

							int nCellX = xc + xi + nx;

							for (int ny = -1; (ny < 2) && changeX; ny++)
							{
								bool changeY = true;
								if (ny == -1)
								{
									changeY = changedNeighbor[1][0];
								}
								else if (ny == 1)
								{
									changeY = changedNeighbor[1][1];
								}

								int nCellY = yc + yi + ny;

								for (int nz = -1; (nz < 2) && changeY; nz++)
								{
									bool changeZ = true;
									if (nz == -1)
									{
										changeZ = changedNeighbor[2][0];
									}
									else if (nz == 1)
									{
										changeZ = changedNeighbor[2][1];
									}

									int nCellZ = zc + zi + nz;
									if (changeZ)
									{
										changedCells->insert(packCellId(CELL_LOD_MIN, nCellX, nCellY, nCellZ));
										// TODO: liang
										// For undo
										// blockData->trackCellChanges(packCellId(CELL_LOD_MIN, nCellX, nCellY, nCellZ), NULL);
									}
								}
							}
						}
					}
				}
			}
		}

		VF_FREE(solids);
		VF_FREE(triangles);
	}

	VF_DELETE buffer;
	VF_FREE(overflow);
	VF_FREE(voxelMaterials);
	VF_FREE(userDataPoints);
	VF_FREE(intersections);
	VF_FREE(intersectionLists);
	VF_FREE(intersectionIndex);
	for (int i = 0; i < RAYS_SHOT; i++)
	{
		VF_FREE(rayTests[i]);
	}
}

class VoxelizationSource : public IMeshStampSource
{
public:
	VoxelizationSource(const VoxelizationInputMesh& inMesh)
		: mesh(inMesh)
	{}

	virtual int getSolidCount() override { return 1; }
	virtual MaterialId getSolidMaterial(int /*solid*/) override { return VOXEL_ONLY_COORDS; }
	virtual int getFaceCount(int /*solid*/) override { return mesh.indexNum / 3; }

	virtual void getFace(int /*solid*/,
		int triIndex,
		Vector& v0,
		Vector& v1,
		Vector& v2) override
	{
		int vi = triIndex * 3;
		v0 = mesh.verts[mesh.indices[vi]];
		v1 = mesh.verts[mesh.indices[vi+1]];
		v2 = mesh.verts[mesh.indices[vi+2]];
	}

private:
	const VoxelizationInputMesh& mesh;
};

class CTestStampMeshMaterials2 : public IMeshStampMaterialSource
{
public:
	MaterialId id;
	virtual MaterialId translateMaterial(const double worldPos[3], MaterialId meshMaterial)
	{
		if (meshMaterial == VOXEL_ONLY_COORDS)
		{
			return id;
		}
		else
		{
			return meshMaterial;
		}
	}
};

void HyVoxelInterfaceImpl::Voxelize(IVoxelLayer* voxelLayer,
	const VoxelizationInputMesh& inputMesh,
	const Matrix& transformMtx,
	bool isAdd,
	HyVoxel::OutVector<HyVoxel::CellId>*& changedCells)
{
	static const Matrix flipYZMtx = Matrix_withValues(
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f);

	VoxelizationSource voxelizationSource(inputMesh);

	CTestStampMeshMaterials2 stampMat;
	// TODO: liang
	// Just hardcode it to be 1.
	stampMat.id = isAdd ? 1 : 0;

	double worldPos[3] = { 0.0, 0.0, 0.0 };

	// TODO: liang
	#if 0
	FMatrix scaleCMToDMMtx = FScaleMatrix::Make(HyVoxel::Vector(0.1f));
	FMatrix adjustOffset = FTranslationMatrix::Make(gAdjustOffset);
	Matrix mtx(transform.ToMatrixWithScale() * adjustOffset *scaleCMToDMMtx * flipYZMtx);
	#else
	Matrix scaleCMToDMMtx = Matrix_identity();
	Matrix_scale(&scaleCMToDMMtx, 0.1f, 0.1f, 0.1f);
	Matrix mtx = Matrix_multiplied(scaleCMToDMMtx, transformMtx);
	mtx = Matrix_multiplied(flipYZMtx, mtx);
	#endif

	BlockDataLayer* blockLayer = static_cast<BlockDataLayer*>(voxelLayer);

	std::set<CellId> changedCellSet;
	stampMeshQEF(blockLayer,
		&voxelizationSource,
		&stampMat,
		worldPos,
		mtx,
		&changedCellSet);

	if (changedCellSet.size() > 0)
	{
		OutVectorImpl<CellId>* ocells = new OutVectorImpl<CellId>();

		for (auto& cellId : changedCellSet)
			ocells->push_back(cellId);

		ocells->SetOutVectorValues();
		changedCells = ocells;
	}
}

struct VoxelHFImpl : public VoxelHF
{
    int nSegmentsX;
    int nSegmentsY;
    int nSegmentsZ;
    HyVoxel::Vector minPos;
    HyVoxel::Vector maxPos;
    HyVoxel::Vector voxDimReciprocal;
    typedef std::vector<int> Indices;
    std::vector<Indices> indicesArray;      // Per voxel
    typedef std::vector<Vector> Vectors;
    std::vector<Vectors> newVerticesArrayPerTriangle; // Per triangle
    std::vector<Vectors> newVerticesArrayPerVoxel;    // Per Voxel

    virtual void Release() override { delete this; }

    virtual void GetIndices(int idx, int*& outData, int& outLength) override
    {
        Indices& tmp = indicesArray[idx];
        outLength = (int)tmp.size();
        outData = (outLength > 0) ? &tmp[0] : nullptr;
    }

    virtual void GetSegments(int& _nSegmentsX, int& _nSegmentsY, int& _nSegmentsZ) override
    {
        _nSegmentsX = nSegmentsX;
        _nSegmentsY = nSegmentsY;
        _nSegmentsZ = nSegmentsZ;
    }

    virtual int GetVoxelIndex(int IndexX, int IndexY, int IndexZ) override
    {
        return IndexX * nSegmentsY * nSegmentsZ + IndexY * nSegmentsZ + IndexZ;
    }

    virtual int GetVoxelNumber() override
    {
        return nSegmentsX * nSegmentsY * nSegmentsZ;
    }

    virtual void GetMinMaxPos(HyVoxel::Vector& _minPos, HyVoxel::Vector& _maxPos) override
    {
        _minPos = minPos;
        _maxPos = maxPos;
    }

    virtual void GetVoxDimReciprocal(HyVoxel::Vector& _voxDimReciprocal) override
    {
        _voxDimReciprocal = voxDimReciprocal;
    }

    virtual void GetNewVerticesPerTriangle(int idx, HyVoxel::Vector*& outData, int& outLength) override
    {
        Vectors& tmp = newVerticesArrayPerTriangle[idx];
        outLength = (int)tmp.size();
        outData = (outLength > 0) ? &tmp[0] : nullptr;
    }

    virtual void GetNewVerticesPerVoxel(int idx, HyVoxel::Vector*& outData, int& outLength) override
    {
        Vectors& tmp = newVerticesArrayPerVoxel[idx];
        outLength = (int)tmp.size();
        outData = (outLength > 0) ? &tmp[0] : nullptr;
    }
};

static bool AlmostEqual(const float a, const float b, float epsilon = 1.0e-5f)
{
    return fabs(a - b) <= epsilon;
}

static bool AlmostEqual(const Vector & v1, const Vector & v2, float epsilon = 1.0e-5f)
{
    return AlmostEqual(Vector::Dist(v1, v2), 0.0f, epsilon);
}

static bool LineSegmentIntersectsPlane(const Vector & P0, const Vector & P1, const Vector & N, const Vector & V0, float & t, Vector & P)
{
    t = Vector::DotProduct(N, V0 - P0) / Vector::DotProduct(N, P1 - P0);

    if (t >= 0.0f && t <= 1.0f) // Intersect
    {
        P = P0 + t * (P1 - P0);
        return true;
    }

    return false;
}

// Only two intersection points need further processing
static bool CheckResults(const int nNewVertices, std::vector<Vector> & newVertices)
{
    if (nNewVertices == 1)  // Something wrong
        newVertices.pop_back();
    else if (nNewVertices == 2)
    {
        // @TODO: Double check if duplicate (intersects right with one of the triangle points)
        return true;
    }
    else
    {
        assert(nNewVertices == 3);
        newVertices.pop_back();
        newVertices.pop_back();
        newVertices.pop_back();
        // @TODO: Something wrong
    }
    return false;
}

static int GetVoxelIndex(const float pos, const float base, const float delta_reciprocal, const int upper_bound, const int lower_bound)
{
    int index = (int)((pos - base) * delta_reciprocal);
    HyVoxel::Clamp(index, upper_bound, lower_bound);

    return index;
}

static int GetVoxelIndex(const int x, const int y, const int z, const int nSegmentsY, const int nSegmentsZ)
{
    return x * (nSegmentsY * nSegmentsZ) + y * (nSegmentsZ) + z;
}

static Vector CalcPolygonCenter(const std::vector<Vector>& polygon)
{
    int nVertices = polygon.size();
    Vector center(0.0f, 0.0f, 0.0f);
    for (int i = 0; i < nVertices; i++)
        center += polygon[i];
    center /= (float)nVertices;

    return center;
}

static float GetVoxelPos(const float base, const float delta_reciprocal, const int index)
{
    return base + (float)index / delta_reciprocal;
}

static float GetRatio(const float begin, const float end, const float pos)
{
    if (AlmostEqual(begin, end))    // Specially handle begin and end are equal
        return -1.0f;
    else
        return (pos - begin) / (end - begin);
}

static bool IsValidRatio(const float t)
{
    return t >= 0.0f && t <= 1.0f;
}

static bool IsValidVoxelIndex(const int voxelID, const int minID, const int maxID)
{
    return (voxelID >= minID) && (voxelID <= maxID);
}

static int CalcDeltaBasedOnDirection(const int begIdx, const int endIdx)
{
    if (begIdx < endIdx)
        return 1;
    else if (begIdx > endIdx)
        return -1;
    else
        return 0;
}

static int CalcVoxelIndexForIntersection(const int delta, const int voxelIdx)
{
    return (delta > 0) ? (voxelIdx + 1) : voxelIdx;
}

//
// Triangle v[3] intersects with Plane (V0, N).
// Returns the number of intersection points, and the first intersection point position IP1 (when the number is 1), 
//     and the second intersection point position IP2 (when the number is 2).
// Returns 0 when no intersection points is available.
// See http://geomalgorithms.com/a06-_intersect-2.html for more details.
//
static int TrianglePlaneIntersection(const Vector v[3], const Vector& V0, const Vector& N, Vector& IP1, int& edgeID1, float& t1, Vector& IP2, int& edgeID2, float& t2)
{
    // Check how many vertices are on the plane
    int nOnPlaneVertices = 0;
    for (int i = 0; i < 3; ++i)
    {
        const Vector & P0 = v[i];
        float denominator = Vector::DotProduct(N, P0 - V0);
        if (AlmostEqual(denominator, 0.0f)) // Vertex on plane
            nOnPlaneVertices++;
        if (nOnPlaneVertices >= 2)
            return 0;
    }

    // Perform the actual intersection
    int nIPVertices = 0;
    for (int i = 0; i < 3; ++i)
    {
        const Vector & P0 = v[i];
        const Vector & P1 = v[(i + 1) % 3];

        float denominator = Vector::DotProduct(N, P1 - P0);
        if (!AlmostEqual(denominator, 0.0f)) // Edge does not parallel to the plane
        {
            float t = Vector::DotProduct(N, V0 - P0) / denominator;
            if (t >= 0.0f && t <= 1.0f) // Intersect
            {
                Vector P = P0 + t * (P1 - P0); // Intersection point

                // Cache IP1 and IP2
                if (nIPVertices == 0)
                {
                    IP1 = P;
                    edgeID1 = i;
                    t1 = t;
                    nIPVertices++;
                }
                else if (nIPVertices == 1)
                {
                    if (!AlmostEqual(P, IP1)) // Think about the intersection line acrosses one triangle vertex.
                    {
                        IP2 = P;
                        edgeID2 = i;
                        t2 = t;
                        nIPVertices++;
                        break;
                    }
                }
            }
        }
    }

    assert(nIPVertices >= 0 && nIPVertices <= 2);
    return nIPVertices;
}

//
// Both end points of the segment are on the plane: no intersection.
// One point above (positive sideness), one point on the plane: no intersection.
// One point below (negative sideness), one point on the plane: has intersection.
//
// Return true if the polygon intersects the plane, else return false.
//
static bool PolygonPlaneIntersection(const std::vector<Vector>& polygon, const Vector& V0, const Vector& N, 
    std::vector<Vector>& newPoly1, std::vector<Vector>& newPoly2)
{
    // Determine the sideness of each vertex.
    int nEdge = polygon.size();
    std::vector<int> sideness;
    sideness.resize(nEdge);
    for (int i = 0; i < nEdge; i++)
    {
        const Vector & P = polygon[i];
        float dotProduct = Vector::DotProduct(N, P - V0);
        if (AlmostEqual(dotProduct, 0.0f)) // On the plane
            sideness[i] = 0;
        else if (dotProduct > 0.0f) // Above the plane
            sideness[i] = 1;
        else // Below the plane
            sideness[i] = -1;
    }

    // Perform the intersection calculation.
    int nIPVertices = 0;
    Vector IP1, IP2;
    int edgeID1, edgeID2;
    float t1, t2;
    for (int i = 0; i < nEdge; i++)
    {
        if (sideness[i] == sideness[(i + 1) % nEdge]) // Segment on the same side of the plane
            continue;

        const Vector & P0 = polygon[i];
        const Vector & P1 = polygon[(i + 1) % nEdge];

        float denominator = Vector::DotProduct(N, P1 - P0);
        if (!AlmostEqual(denominator, 0.0f)) // Edge does not parallel to the plane
        {
            float t = Vector::DotProduct(N, V0 - P0) / denominator;
            if (t >= 0.0f && t <= 1.0f) // Intersect
            {
                Vector P = P0 + t * (P1 - P0); // Intersection point

                // Cache IP1 and IP2
                if (nIPVertices == 0)
                {
                    IP1 = P;
                    edgeID1 = i;
                    t1 = t;
                    nIPVertices++;
                }
                else if (nIPVertices == 1)
                {
                    if (!AlmostEqual(P, IP1)) // Think about the intersection line acrosses one vertex.
                    {
                        IP2 = P;
                        edgeID2 = i;
                        t2 = t;
                        nIPVertices++;
                        break;
                    }
                }
            }
        }
    }

    // Split the polygon.
    if (nIPVertices == 2)
    {
        // Generate newPoly1
        newPoly1.clear();
        for (int i = 0; i <= edgeID1; i++)
        {
            const Vector & P = polygon[i];
            newPoly1.push_back(P);
        }
        if (!AlmostEqual(t1, 0.0f))
            newPoly1.push_back(IP1);
        if (!AlmostEqual(t2, 1.0f))
            newPoly1.push_back(IP2);
        for (int i = edgeID2 + 1; i < nEdge; i++)
        {
            const Vector & P = polygon[i];
            newPoly1.push_back(P);
        }

        // Generate newPoly2
        newPoly2.clear();
        if (!AlmostEqual(t1, 1.0f))
            newPoly2.push_back(IP1);
        for (int i = edgeID1 + 1; i <= edgeID2; i++)
        {
            const Vector & P = polygon[i];
            newPoly2.push_back(P);
        }
        if (!AlmostEqual(t2, 0.0f))
            newPoly2.push_back(IP2);

        // If intersects right on an edge, one of newPolys will have just 2 vertices.
        return newPoly1.size() > 2 && newPoly2.size() > 2;
    }

    return false;
}

//
// Line segment intersects with a set of axis-aligned parallel planes.
//
static void LineSegmentIntersectsParallelPlanes(const float IP1, const float IP2, 
    const float base, const float delta_reciprocal, const int nSegment,
    std::vector<float>& ParamArray)
{
    ParamArray.clear();

    // Ensure consistency with GetRatio().
    if (AlmostEqual(IP1, IP2)) // Degenerate, no need to go further.
        return;

    int begIdx = GetVoxelIndex(IP1, base, delta_reciprocal, nSegment - 1, 0);
    int endIdx = GetVoxelIndex(IP2, base, delta_reciprocal, nSegment - 1, 0);
    int delta = CalcDeltaBasedOnDirection(begIdx, endIdx);
    for (int i = begIdx; i != endIdx; i += delta)
    {
        assert(IsValidVoxelIndex(i, 0, nSegment - 1));
        float pos = GetVoxelPos(base, delta_reciprocal, CalcVoxelIndexForIntersection(delta, i));
        float t = GetRatio(IP1, IP2, pos);
        if (IsValidRatio(t))
            ParamArray.push_back(t);
    }
}

static void GenerateIntersectionPoints(const Vector& IP1, const Vector& IP2, const std::vector<float>& ParamArray, std::vector<Vector>& newVertices)
{
    newVertices.push_back(IP1);
    int nIP = ParamArray.size();
    for (int idx = 0; idx < nIP; idx++)
    {
        float t = ParamArray[idx];
        Vector IP = IP1 + (IP2 - IP1) * t;

        bool bDuplicate = false;
        if (idx == 0)   // First vertex, check if very close to IP1
        {
            if (AlmostEqual(IP1, IP))
                bDuplicate = true;
        }
        if (idx == nIP - 1) // Last vertex, check if very close to IP2 (Cover only one case)
        {
            if (AlmostEqual(IP2, IP))
                bDuplicate = true;
        }
        if (!bDuplicate)
        {
            newVertices.push_back(IP);
            newVertices.push_back(IP);
        }
    }
    newVertices.push_back(IP2);
}

//
// Merge two param arrays.
// Input: ParamArray1 and ParamArray2, two ascending float arrays within the range of [0.0, 1.0].
// Output: ParamArray, ascending float array that is the merged results.
//
static void SortIntersectionPointParams(const std::vector<float>& ParamArray1, const std::vector<float>& ParamArray2, std::vector<float>& ParamArray)
{
    int idx1 = 0, idx2 = 0;
    int nIP1 = ParamArray1.size(), nIP2 = ParamArray2.size();
    while (idx1 < nIP1 && idx2 < nIP2)
    {
        float t1 = ParamArray1[idx1];
        float t2 = ParamArray2[idx2];
        if (AlmostEqual(t1, t2)) // Move along 1 and 2 simultaneously
        {
            ParamArray.push_back(t1);
            idx1++;
            idx2++;
        }
        else if (t1 < t2) // Move along 1
        {
            ParamArray.push_back(t1);
            idx1++;
        }
        else // Move along 2
        {
            ParamArray.push_back(t2);
            idx2++;
        }
    }
    while (idx1 < nIP1) // Continue moving along 1
    {
        float t1 = ParamArray1[idx1];
        ParamArray.push_back(t1);
        idx1++;
    }
    while (idx2 < nIP2) // Continue moving along 2
    {
        float t2 = ParamArray2[idx2];
        ParamArray.push_back(t2);
        idx2++;
    }
}

void HyVoxelInterfaceImpl::VoxelizeHF(const VoxelizationInputMesh& inputMesh, const Matrix& transformMtx, 
    const int nSegmentsX, const int nSegmentsY, const int nSegmentsZ, 
    VoxelHF*& voxelHF)
{
    // Allocated in HyVoxel (one dll), read-only used in VoxelHF (another dll).
    // Application should call VoxelizeHF_Release() to release the allocated momery.
    VoxelHFImpl* voxelHFImpl = new VoxelHFImpl();
    if (voxelHFImpl == NULL)
        return;
    voxelHF = static_cast<VoxelHF *>(voxelHFImpl);

    int nVertices = inputMesh.vertNum;
    int nIndices = inputMesh.indexNum;

    voxelHFImpl->nSegmentsX = nSegmentsX;
    voxelHFImpl->nSegmentsY = nSegmentsY;
    voxelHFImpl->nSegmentsZ = nSegmentsZ;

    // Calculate the min/max position of the mesh
    HyVoxel::Vector minPos( FLT_MAX,  FLT_MAX,  FLT_MAX);
    HyVoxel::Vector maxPos(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    for (int i = 0; i < nVertices; ++i)
    {
        HyVoxel::Vector & pos = inputMesh.verts[i];
        minPos.x = (pos.x < minPos.x) ? pos.x : minPos.x;
        minPos.y = (pos.y < minPos.y) ? pos.y : minPos.y;
        minPos.z = (pos.z < minPos.z) ? pos.z : minPos.z;
        maxPos.x = (maxPos.x < pos.x) ? pos.x : maxPos.x;
        maxPos.y = (maxPos.y < pos.y) ? pos.y : maxPos.y;
        maxPos.z = (maxPos.z < pos.z) ? pos.z : maxPos.z;
    }

    voxelHFImpl->minPos = minPos;
    voxelHFImpl->maxPos = maxPos;

    // Calculate the reciprocal of the voxel dimension
    HyVoxel::Vector voxDimReciprocal;
    voxDimReciprocal.x = (float)nSegmentsX / (maxPos.x - minPos.x);
    voxDimReciprocal.y = (float)nSegmentsY / (maxPos.y - minPos.y);
    voxDimReciprocal.z = (float)nSegmentsZ / (maxPos.z - minPos.z);

    voxelHFImpl->voxDimReciprocal = voxDimReciprocal;

    // Go though each triangle, put it into overlapping voxels based on its bounding box.
    HyVoxel::Vector v[3];
    HyVoxel::Vector minVtx, maxVtx;
    typedef std::vector<int> Indices;
    std::vector<Indices>& IndicesArray = voxelHFImpl->indicesArray;
    IndicesArray.resize(nSegmentsX*nSegmentsY*nSegmentsZ);
    typedef std::vector<Vector> Vectors;
    Vectors newVertices;
    std::vector<Vectors>& newVerticesArrayPerTriangle = voxelHFImpl->newVerticesArrayPerTriangle;
    newVerticesArrayPerTriangle.resize(nIndices / 3);
    std::vector<Vectors>& newVerticesArrayPerVoxel = voxelHFImpl->newVerticesArrayPerVoxel;
    newVerticesArrayPerVoxel.resize(nSegmentsX*nSegmentsY*nSegmentsZ);
    for (int i = 0; i < nIndices; i += 3)
    {
        v[0] = inputMesh.verts[inputMesh.indices[i]];
        v[1] = inputMesh.verts[inputMesh.indices[i+1]];
        v[2] = inputMesh.verts[inputMesh.indices[i+2]];

        // Calculate the bounding box of the triangle
        minVtx = v[0];
        maxVtx = v[0];
        for (int j = 1; j < 3; ++j)
        {
            minVtx.x = (v[j].x < minVtx.x) ? v[j].x : minVtx.x;
            minVtx.y = (v[j].y < minVtx.y) ? v[j].y : minVtx.y;
            minVtx.z = (v[j].z < minVtx.z) ? v[j].z : minVtx.z;

            maxVtx.x = (maxVtx.x < v[j].x) ? v[j].x : maxVtx.x;
            maxVtx.y = (maxVtx.y < v[j].y) ? v[j].y : maxVtx.y;
            maxVtx.z = (maxVtx.z < v[j].z) ? v[j].z : maxVtx.z;
        }

        // Calculate overlapping voxels of the triangle
        int begIdxX = GetVoxelIndex(minVtx.x, minPos.x, voxDimReciprocal.x, nSegmentsX - 1, 0);
        int endIdxX = GetVoxelIndex(maxVtx.x, minPos.x, voxDimReciprocal.x, nSegmentsX - 1, 0);
        int begIdxY = GetVoxelIndex(minVtx.y, minPos.y, voxDimReciprocal.y, nSegmentsY - 1, 0);
        int endIdxY = GetVoxelIndex(maxVtx.y, minPos.y, voxDimReciprocal.y, nSegmentsY - 1, 0);
        int begIdxZ = GetVoxelIndex(minVtx.z, minPos.z, voxDimReciprocal.z, nSegmentsZ - 1, 0);
        int endIdxZ = GetVoxelIndex(maxVtx.z, minPos.z, voxDimReciprocal.z, nSegmentsZ - 1, 0);

        for (int x = begIdxX; x <= endIdxX; ++x)
        {
            for (int y = begIdxY; y <= endIdxY; ++y)
            {
                for (int z = begIdxZ; z <= endIdxZ; ++z)
                {
                    int idx = GetVoxelIndex(x, y, z, nSegmentsY, nSegmentsZ);
                    IndicesArray[idx].push_back(i);
                }
            }
        }

        // Calculate intersection lines of the triangle (by XOY, YOZ and ZOX).
        Vector V0, N;
        typedef std::vector<Vector> Polygon;
        typedef std::vector<Polygon> PolygonArray;
        PolygonArray Polygons, NewPolygons;
        Polygon poly, newPoly1, newPoly2;
        poly.push_back(v[0]);
        poly.push_back(v[1]);
        poly.push_back(v[2]);
        Polygons.push_back(poly);
        PolygonArray& TraversePolygons = Polygons;
        PolygonArray& GeneratePolygons = NewPolygons;
        for (int x = begIdxX + 1; x <= endIdxX; ++x) // No need to traverse a single voxel (begIdxX == endIdxX).
        {
            // Construct the intersection plane (V0, N).
            N = Vector(-1.0f, 0.0f, 0.0f);
            V0.x = GetVoxelPos(minPos.x, voxDimReciprocal.x, x);
            V0.y = 0.0f;
            V0.z = 0.0f;

            // For each polygon to be traversed
            int nPoly = TraversePolygons.size();
            for (int n = 0; n < nPoly; n++)
            {
                Polygon& polygon = TraversePolygons[n];
                if (PolygonPlaneIntersection(polygon, V0, N, newPoly1, newPoly2))
                {
                    GeneratePolygons.push_back(newPoly1);
                    GeneratePolygons.push_back(newPoly2);
                }
                else
                {
                    GeneratePolygons.push_back(polygon);
                }
            }

            // Swith polygon buffers
            PolygonArray& tmp = TraversePolygons;
            TraversePolygons = GeneratePolygons;
            GeneratePolygons = tmp;
            GeneratePolygons.clear();
        }

        for (int y = begIdxY+1; y <= endIdxY; ++y) // No need to traverse a single voxel (begIdxY == endIdxY).
        {
            N = Vector(0.0f, -1.0f, 0.0f);
            V0.x = 0.0f;
            V0.y = GetVoxelPos(minPos.y, voxDimReciprocal.y, y);
            V0.z = 0.0f;

            // For each polygon to be traversed
            int nPoly = TraversePolygons.size();
            for (int n = 0; n < nPoly; n++)
            {
                Polygon& polygon = TraversePolygons[n];
                if (PolygonPlaneIntersection(polygon, V0, N, newPoly1, newPoly2))
                {
                    GeneratePolygons.push_back(newPoly1);
                    GeneratePolygons.push_back(newPoly2);
                }
                else
                {
                    GeneratePolygons.push_back(polygon);
                }
            }

            // Swith polygon buffers
            PolygonArray& tmp = TraversePolygons;
            TraversePolygons = GeneratePolygons;
            GeneratePolygons = tmp;
            GeneratePolygons.clear();
        }
        
        for (int z = begIdxZ+1; z <= endIdxZ; ++z) // No need to traverse a single voxel (begIdxZ == endIdxZ).
        {
            N = Vector(0.0f, 0.0f, -1.0f);
            V0.x = 0.0f;
            V0.y = 0.0f;
            V0.z = GetVoxelPos(minPos.z, voxDimReciprocal.z, z);

            // For each polygon to be traversed
            int nPoly = TraversePolygons.size();
            for (int n = 0; n < nPoly; n++)
            {
                Polygon& polygon = TraversePolygons[n];
                if (PolygonPlaneIntersection(polygon, V0, N, newPoly1, newPoly2))
                {
                    GeneratePolygons.push_back(newPoly1);
                    GeneratePolygons.push_back(newPoly2);
                }
                else
                {
                    GeneratePolygons.push_back(polygon);
                }
            }

            // Swith polygon buffers
            PolygonArray& tmp = TraversePolygons;
            TraversePolygons = GeneratePolygons;
            GeneratePolygons = tmp;
            GeneratePolygons.clear();
        }

        // Generate triangles based on resulted polygons.
        newVertices.clear();
        int nPoly = TraversePolygons.size();
        int triangleID = i / 3;
        for (int n = 0; n < nPoly; n++)
        {
            Polygon& polygon = TraversePolygons[n];
            int nVertices = polygon.size();
            assert(nVertices >= 3);
            for (int v = 0; v < nVertices; v++)
            {
                newVerticesArrayPerTriangle[triangleID].push_back(polygon[v]);
                newVerticesArrayPerTriangle[triangleID].push_back(polygon[(v+1)%nVertices]);
            }

            // Calculate polygon center, assuming that the entire polygon is within a single voxel.
            Vector center = CalcPolygonCenter(polygon);
            int idxX = GetVoxelIndex(center.x, minPos.x, voxDimReciprocal.x, nSegmentsX - 1, 0);
            int idxY = GetVoxelIndex(center.y, minPos.y, voxDimReciprocal.y, nSegmentsY - 1, 0);
            int idxZ = GetVoxelIndex(center.z, minPos.z, voxDimReciprocal.z, nSegmentsZ - 1, 0);
            int voxelID = GetVoxelIndex(idxX, idxY, idxZ, nSegmentsY, nSegmentsZ);
            /*
            // Generate edges based on the polygon
            for (int v = 0; v < nVertices; v++)
            {
                newVerticesArrayPerVoxel[voxelID].push_back(polygon[v]);
                newVerticesArrayPerVoxel[voxelID].push_back(polygon[(v + 1) % nVertices]);
            }
            */

            // Generate triangles based on the polygon
            for (int v = 1; v < nVertices - 1; v++)
            {
                newVerticesArrayPerVoxel[voxelID].push_back(polygon[0]);
                newVerticesArrayPerVoxel[voxelID].push_back(polygon[v]);
                newVerticesArrayPerVoxel[voxelID].push_back(polygon[(v + 1) % nVertices]);
            }
        }
    }
}

}