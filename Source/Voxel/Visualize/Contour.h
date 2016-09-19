#pragma once

#include "Voxel/VoxelData.h"
#include "HyVoxelLib.h"
#include "Common/MatrixSSE.h"
#include "Common/Matrix.h"
#include "Voxel/Material/MaterialLibrary.h"
#include "Voxel/Visualize/CellMesh.h"

namespace HyVoxel {

/// How many times the contouring octree can subdivide
const int CONTOUR_OCTREE_LEVELS = 7;
/// Number of voxels along one dimension in the contouring octree
const int CONTOUR_OCTREE_SIZE = 64;

/// A node in the contouring octree. This octree is used to compress the polygonal output of the contouring stage.
struct OctreeNode
{
	/// Links the node to another node in the octree. This is set when the octree is compressed
	unsigned int map;
	/// Pointer to the vertex information
	unsigned int index;
	/// Pointers to the eight children
	unsigned int children[8];
	/// Octree can collapse below this specified level
	short targetLevel;
	/// A point inside the cell where the isosurface crosses
	Vector p;
	/// Number of isosurface crossings
	int count;
	/// Indicates whether it is a end node in the octree
	bool leaf;
	/// Indicates whether the node is in the world cell boundary
	bool boundary;
	/// Maximum simplification error allowed for the node
	double maxError;
	/// Quadratic Error collected for the cell
	Algebra::QEFMatrix qef;
	/// Material id for the node
	int material;
	/// Node is made of same material inside
	bool homogeneus;
};

/// Generation of voxel data and then meshes out if it can be time consuming. The engine allows multiple threads to call the contouring functions so the workload can be split among them. Each calling thread must create its own ContourThreadContext object.
class ContourThreadContext
{
public:
	ContourThreadContext();
	virtual ~ContourThreadContext();
	ContourVoxelData* voxelData;
	Vector* points;
	int* pointIndex;
	OctreeNode* nodes;
	int* quadIndex;
	int* quadTypes;
	int* vertexTypes;
	bool* edgeSign;
	MaterialId* materials;

	Algebra::QEFMatrixSSE* comp_quadrics;
	bool* comp_fVert;
	bool* comp_sVert;
	bool* comp_lVert;
	bool* comp_vPlanes;
	double* comp_qPlanes;
	unsigned short* comp_mMap;
};

int GenerateCellMesh(
	/// A ContourThreadContext instance containing the work buffers for the contouring operation
	ContourThreadContext* tc,
	/// Material library to be used in contouring
	CMaterialLibrary* materialLibrary,
	/// The CellMesh object that will receive the polygonal output
	CellMesh* cellData,
	/// Signals whether simplify the meshes
	bool compress,
	/// Signals whether to simplify only the cell boundaries
	bool compressOnlyBoundaries = false
	);

extern CMaterialLibrary gMaterialLibrary;


}

