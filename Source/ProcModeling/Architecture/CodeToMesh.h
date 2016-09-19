#pragma once

#include "HyVoxelLib.h"

#include "Common/HyContainerImpl.h"

namespace HyVoxel {

typedef std::map<std::string, void*> Name2Mesh;

class HyArchitectureInterfaceImpl : public IHyArchitectureInterface
{
public:
	virtual void Release() override { delete this; }

	// Load Dat format meshes.
	// Please refer to:
	//		http://docs.voxelfarm.com/mesh-asset-format
	// This feature is not used in our products yet.
	virtual int LoadDatMeshFromFile(const char* meshName, const char* fullPath) override;
	virtual int LoadDatMeshFromString(const char* meshName, const char* buf) override;
	virtual bool IsDatMeshLoaded(const char* meshName) override;

	// Load prefabs.
	// Please refer to:
	//		http://docs.voxelfarm.com/introduction-to-prefabs
	// This feature is used in HyBuilder.
	virtual void SetLSystemPath(const char* path) override;
	virtual int LoadPrefabFromFile(const char* prefabName, const char* fullPath) override;
	virtual int LoadPrefabFromString(const char* prefabName, const char* buf) override;
	virtual bool IsPrefabLoaded(const char* prefabName) override;
	virtual void UnloadPrefab(const char* prefabName) override;

	// Position user meshes according to the input prefab
	// This feature is used in HyBuilder.
	virtual void AddNameToUserMeshMapping(const char* name, void* userMesh) override;
	virtual void ClearNameToUserMeshMappings() override;
	virtual bool ProcedurallyPositionUserMeshes(const char* prefabName, double boxSize[3], OutVector<InstancedMesh>*& outMeshes) override;

	// Generate a mesh according to
	//		- prefab
	//		- loaded Dat format meshes
	// This feature is not used in our products yet.
	virtual int GenerateMeshFromPrefab(const char* prefabName, const Vector& boxSize, OutputMesh& outMesh) override;

private:
	bool ProcedurallyPositionUserMeshesImpl(const Name2Mesh& name2Mesh, const char* prefabName, double boxSize[3], hyvector<InstancedMesh>& outMeshes);
};

}

// TODO: liang
// to be removed
void TestVoxelization(const char* datPath, const HyVoxel::Vector& boxSize, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, bool isCubeContour);
void TestMeshToVoxel(std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, const HyVoxel::Vector& boxSize, bool isCubeContour);
