#pragma once

#include "HyVoxelLib.h"

#include "Common/Vector.h"
#include "Voxel/Layer/BlockDataLayer.h"
#include "Physics/PhysicsData.h"

namespace HyVoxel {

// TODO: liang
// Using virtual function for each face is too luxury.
// Use template instead!
/// It allows to access to the faces and materials of a list of solids
class IMeshStampSource
{
public:
	/// Returns the number of solids
	virtual int getSolidCount() = 0;
	/// Returns the material of a solid
	virtual MaterialId getSolidMaterial(int solid) = 0;
	/// Returns the number of faces in a solid
	virtual int getFaceCount(int solid) = 0;
	/// Get the vertices of a face of a solid
	virtual void getFace(int solid,
		int index,
		HyVoxel::Vector& v0,
		HyVoxel::Vector& v1,
		HyVoxel::Vector& v2) = 0;
};

/// It translates a given material depending of the position in the world
class IMeshStampMaterialSource
{
public:
	virtual MaterialId translateMaterial(const double worldPos[3], MaterialId meshMaterial) = 0;
};

/// It voxelizes a mesh and stamps it into a block layer and minimizes the error with the existing data
void stampMeshQEF(
	/// A reference to the block layer
	BlockDataLayer* blockData,
	/// The mesh to be stamped
	IMeshStampSource* mesh,
	/// It translates a given material depending of the position in the world
	IMeshStampMaterialSource* materials,
	/// Position in world to stamp the mesh
	const double worldPos[3],
	/// A matrix for transforming the mesh before process it
	const Matrix& transform,
	/// Returns a list of all affected cells during the process
	std::set<CellId>* changedCells,

	/// Physics ----
	Physics::CNoiseBrush* noise = NULL,
	Physics::CVoxelBuffer* nVoxBuffer = NULL,
	const TVFMap<CellId, unsigned char*>* voxelCache = NULL,
	const int hitMaterial = 0
	);


}