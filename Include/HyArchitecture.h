#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

namespace HyVoxel {

	#pragma pack(push, 8)
	/// Enumerates possible scope types
	enum ScopeType
	{
		SCOPE_BOX,
		SCOPE_PRISM,
		SCOPE_PRISM_LEFT,
		SCOPE_PRISM_RIGHT
	};

	/// Tracks a single instance of a mesh. An architecture entity will be composed of many mesh instances.
	struct InstancedMesh
	{
		/// A reference to a polygonal mesh for the instance
		void* mesh;

		/// Material for the instance
		int material;
		/// Reference to the entity to which the instance belongs
		void* entity;
		/// Position and rotation of the instance
		Matrix transform;
		/// Scale for the instance
		Vector size;
		/// Type of scope for the instance (Box or Prism)
		ScopeType scopeType;
	};

	class IHyArchitectureInterface : public IHyVoxelInterfaceBase
	{
	public:
		// Load Dat format meshes.
		// Please refer to:
		//		http://docs.voxelfarm.com/mesh-asset-format
		// This feature is not used in our products yet.
		virtual int LoadDatMeshFromFile(const char* meshName, const char* fullPath) = 0;
		virtual int LoadDatMeshFromString(const char* meshName, const char* buf) = 0;
		virtual bool IsDatMeshLoaded(const char* meshName) = 0;

		// Load prefabs.
		// Please refer to:
		//		http://docs.voxelfarm.com/introduction-to-prefabs
		// Please note that, the prefab has to follow ".mod" grammar.
		// Please DO NOT USE the ".code" file.
		// This feature is used in HyBuilder.
		virtual void SetLSystemPath(const char* path) = 0;
		virtual int LoadPrefabFromFile(const char* prefabName, const char* fullPath) = 0;
		virtual int LoadPrefabFromString(const char* prefabName, const char* buf) = 0;
		virtual bool IsPrefabLoaded(const char* prefabName) = 0;
		virtual void UnloadPrefab(const char* prefabName) = 0;

		// Position user meshes according to the input prefab
		// This feature is used in HyBuilder.
		virtual void AddNameToUserMeshMapping(const char* name, void* userMesh) = 0;
		virtual void ClearNameToUserMeshMappings() = 0;
		virtual bool ProcedurallyPositionUserMeshes(const char* prefabName, double boxSize[3], OutVector<InstancedMesh>*& outMeshes) = 0;

		// Generate a mesh according to
		//		- prefab
		//		- loaded Dat format meshes
		// This feature is not used in our products yet.
		virtual int GenerateMeshFromPrefab(const char* prefabName, const Vector& boxSize, OutputMesh& outMesh) = 0;
	};
	#pragma pack(pop)
}
