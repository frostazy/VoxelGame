#include "HyVoxelPrivatePCH.h"

#include "CodeToMesh.h"

#include "HyVoxelConfig.h"
#include "ProcModeling/Architecture/Grammar.h"
#include "Common/FastQuadrics.h"
#include "Common/fileutils.h"
#include "GoldEngine/ModCompiler.h"
#include "Voxel/VoxelData.h"
#include "Voxel/Voxelization.h"
#include "Voxel/ClipmapView.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <map>

using namespace HyVoxel;
using namespace HyVoxel::Architecture;
using namespace HyVoxel::Algebra;

namespace HyVoxel{


// POD
struct ArchPrefab
{
	CGrammar		*grammar;
};
std::map<std::string, ArchPrefab>		gArchPrefabs;

#if 1
// TODO:
// The code in these functions are copied from ArchitectureManager.cpp.
// Currently, we do not draw the whole ArchitectureManager.cpp in, because we don't need the voxel related logic for now.
template <typename TStream>
static CFastQuadrics* LoadMeshFromStream(TStream& file)
{

	CFastQuadrics* result = VF_NEW CFastQuadrics();

	float minX = (float)INT_MAX;
	float minY = (float)INT_MAX;
	float minZ = (float)INT_MAX;
	float maxX = 0.0f;
	float maxY = 0.0f;
	float maxZ = 0.0f;

	float mesh_minX = (float)INT_MAX;
	float mesh_minY = (float)INT_MAX;
	float mesh_minZ = (float)INT_MAX;
	float mesh_maxX = 0.0f;
	float mesh_maxY = 0.0f;
	float mesh_maxZ = 0.0f;

	// Read number of submeshes
	int count;
	file >> count;

	int totalVertexCount = 0;
	int totalFaceCount = 0;

	bool hasReference = false;

	// Store submeshes in a list
	TVFVector<CFastQuadrics*> meshes;

	// Repeat for each submesh in the file
	for (int iMesh = 0; iMesh < count; iMesh++)
	{
		// Read header values
		int isReference;
		int vertexCount;
		int faceCount;
		file >> isReference;
		file >> vertexCount;
		file >> faceCount;

		hasReference = hasReference || isReference;

		// Is this the reference mesh?
		if (!isReference)
		{
			// The mesh has polygons that should be used
			totalVertexCount += vertexCount;
			totalFaceCount += faceCount;

			// Create mesh object
			CFastQuadrics* tmpMesh = VF_NEW CFastQuadrics();
			tmpMesh->allocate(vertexCount, faceCount);
			meshes.push_back(tmpMesh);

			// Read vertices
			for (int i = 0; i < vertexCount; i++)
			{
				float x, y, z;
				file >> x;
				file >> z;
				file >> y;

				mesh_minX = std::min(mesh_minX, x);
				mesh_minY = std::min(mesh_minY, y);
				mesh_minZ = std::min(mesh_minZ, z);
				mesh_maxX = std::max(mesh_maxX, x);
				mesh_maxY = std::max(mesh_maxY, y);
				mesh_maxZ = std::max(mesh_maxZ, z);

				tmpMesh->addVertex(x, y, z, 0);
			}

			// Read faces
			for (int i = 0; i < faceCount; i++)
			{
				int id1, id2, id3, group;
				file >> id1;
				file >> id2;
				file >> id3;
				file >> group;
				tmpMesh->addFace(id1, id2, id3, 0);
			}
		}
		else
		{
			// It is the reference mesh
			// Go over its vertices to compute boundaries (min, max)
			for (int i = 0; i < vertexCount; i++)
			{
				float x, y, z;
				file >> x;
				file >> z;
				file >> y;
				minX = std::min(minX, x);
				minY = std::min(minY, y);
				minZ = std::min(minZ, z);
				maxX = std::max(maxX, x);
				maxY = std::max(maxY, y);
				maxZ = std::max(maxZ, z);
			}
			for (int i = 0; i < faceCount; i++)
			{
				int id1, id2, id3, group;
				file >> id1;
				file >> id2;
				file >> id3;
				file >> group;
			}
		}
	}

	if (!hasReference)
	{
		minX = mesh_minX;
		minY = mesh_minY;
		minZ = mesh_minZ;
		maxX = mesh_maxX;
		maxY = mesh_maxY;
		maxZ = mesh_maxZ;
	}


	// Allocated number of read vertices and faces in the result mesh
	result->allocate(totalVertexCount, totalFaceCount);

	// Add vertices and faces from all submeshes
	int vertIdBase = 0;
	for (TVFVector<CFastQuadrics*>::iterator imesh = meshes.begin();
		imesh != meshes.end(); ++imesh)
	{
		CFastQuadrics* mesh = *imesh;
		if (mesh == NULL)
		{
			continue;
		}
		for (int i = 0; i < mesh->vertexCount; i++)
		{
			FQ_Vertex& v = mesh->vertices[i];
			result->addVertex((v.x - minX) / (maxX - minX), (v.y - minY) / (maxY - minY), (v.z - minZ) / (maxZ - minZ), 0);
		}
		for (int i = 0; i < mesh->faceCount; i++)
		{
			FQ_Face& f = mesh->faces[i];
			result->addFace(vertIdBase + f[0], vertIdBase + f[1], vertIdBase + f[2], f[3]);
		}
		vertIdBase += mesh->vertexCount;
		VF_DELETE mesh;
	}
	return result;
}


/// An RTree index to accelerate spatial instance queries
typedef RTree<InstancedMesh*, double, 3, double, 32> InstanceIndex;

/// An Entity is an unique architecture element. The entity struct tracks the location and type of the entity, along with some other data.
struct Entity
{
	/// Class identifier for the entity. Grammar rules to generate the entity are defined at the class level, as many different entities can be generated from the same rules.
	int classId;
	/// Position for the entity in world coordinates
	double position[3];
	/// Orientation for the entity in degrees over each main coordinate axis.
	double orientation[3];
	/// Scale for the entity. Determines the size of the initial scope used by the generation system.
	double size[3];
	/// Index of all the instances that were produced for the entity, grouped by architecture LOD
	HyVoxel::TVFMap<int, InstanceIndex*> lod;
	/// A map that resolves alphanumeric material IDs as defined in the grammar source code to numeric material IDs as defined by the material library.
	HyVoxel::TVFMap<VFString, int>* materials;
	/// Height of the ground level for the instance.
	double groundLevel;
	/// An entropy value to simulate decay and breackage for the entity.
	double decay;
};

enum PaletteType
{
	/// do not use a palette...take the material from the entity material list
	PT_NONE,
	/// take as material the user material in use
	PT_DEFAULT,
	/// set material to 0 for deleting with the instance
	PT_ERASE,
	/// choose the material from a palette of materials
	PT_PALETTE
};

/// Implements IInstanceCreator to collect all instances created during the program run
class CInstanceCollector : public IInstanceCreator
{
public:
	CInstanceCollector(const Name2Mesh& name2Mesh, hyvector<InstancedMesh>& outInstances)
		: meshes(name2Mesh)
		, instances(outInstances)
	{
		type = PT_NONE;
		entity = NULL;
	}

	Entity* entity;
	PaletteType type;
	virtual void createInstance(VFString id, Matrix& matrix, Vector& size, ScopeType scopeType, VFString material) override;

public:
	const Name2Mesh& meshes;
	hyvector<InstancedMesh>& instances;
};

void CInstanceCollector::createInstance(VFString id, Matrix& matrix, Vector& size, ScopeType scopeType, VFString material)
{
	// Test if there is a mesh registered for the specified ID
	Name2Mesh::const_iterator imesh = meshes.find(id);
	if (imesh == meshes.end())
	{
		return;
	}

	// Create a new record for the instance and store it
	InstancedMesh instance;
	instance.mesh = imesh->second;
	instance.entity = entity;
	instance.transform = matrix;
	instance.scopeType = scopeType;
	instance.size = size;

	if (type == PT_PALETTE)
	{
		if (entity != NULL && material != "")
		{
			TVFMap<VFString, int>::iterator i = entity->materials->find(material);
			if (i != entity->materials->end())
			{
				instance.material = i->second;
			}
			else
			{
				instance.material = VOXEL_ONLY_COORDS;
			}
		}
		else
		{
			instance.material = VOXEL_ONLY_COORDS;
		}
	}
	else if (type == PT_DEFAULT)
	{
		instance.material = VOXEL_ONLY_COORDS;
	}
	else if (type == PT_ERASE)
	{
		instance.material = 0;
	}
	else
	{
		if (entity != NULL)
		{
			instance.material = (*(entity->materials))[material];
		}
		else
		{
			instance.material = 1;
		}
	}

	instances.push_back(instance);
}

#endif

bool ReadFileIntoString(const char* fullPath, std::string& outStr)
{
	if (!HyVoxel::fileExists(fullPath))
		// failed to load
		return false;

	std::ifstream t(fullPath);
	std::stringstream buffer;
	buffer << t.rdbuf();

	outStr = buffer.str();

	return true;
}

bool HyArchitectureInterfaceImpl::IsPrefabLoaded(const char* prefabName)
{
	return gArchPrefabs.find(prefabName) != gArchPrefabs.end();
}

void HyArchitectureInterfaceImpl::UnloadPrefab(const char* prefabName)
{
	std::map<std::string, ArchPrefab>::iterator it = gArchPrefabs.find(prefabName);

	if (it == gArchPrefabs.end())
		return;

	ArchPrefab oldPrefab = it->second;
	gArchPrefabs.erase(it);

	delete oldPrefab.grammar;
}

static std::string sLSystemPath;
void HyArchitectureInterfaceImpl::SetLSystemPath(const char* path)
{
	sLSystemPath = path;
}

int HyArchitectureInterfaceImpl::LoadPrefabFromString(const char* prefabName, const char* buf)
{
	if (IsPrefabLoaded(prefabName))
		return 0;

	std::vector<GoldCPP::SourceFile> files;
	files.push_back(GoldCPP::SourceFile());
	files[0].modSrc = buf;

	std::string outCode, outDbg;
	std::vector<GoldCPP::GrammarError> errors;

	GoldCPP::ModCompiler compiler(sLSystemPath.c_str());
	compiler.Compile(files, outCode, outDbg, errors);

	CGrammar* grammar = new CGrammar();

	CModule* mainModule = &grammar->define("main");
	mainModule->baseAxiom = "main";
	mainModule->loadFromString(outCode);

	ArchPrefab newPrefab;
	newPrefab.grammar = grammar;

	gArchPrefabs[prefabName] = newPrefab;
	return 0;
}

int HyArchitectureInterfaceImpl::LoadPrefabFromFile(const char* prefabName, const char* fullPath)
{
	// early out
	if (IsPrefabLoaded(prefabName))
		return 0;

	std::string outStr;
	if (!ReadFileIntoString(fullPath, outStr))
		return -1;

	return LoadPrefabFromString(prefabName, outStr.c_str());
}

#if 1

// TODO:
// Copied from ArchitectureManager.h
/// This class allows to access meshes in a list of instances.
class CArchitectureMesh : public IMeshStampSource
{
public:
	CArchitectureMesh(hyvector<InstancedMesh>* instances, Matrix);
	virtual int getSolidCount();
	virtual MaterialId getSolidMaterial(int solid);
	virtual int getFaceCount(int solid);
	virtual void getFace(int solid, int index, Vector& v0, Vector& v1, Vector& v2);
private:
	/// List of instances
	hyvector<InstancedMesh>* meshInstances;
	/// A rotation matrix
	Matrix transform;
};

CArchitectureMesh::CArchitectureMesh(hyvector<InstancedMesh>* instances, Matrix transform)
	: transform(transform), meshInstances(instances)
{
}

int CArchitectureMesh::getSolidCount()
{
	return (int)meshInstances->size();
}

MaterialId CArchitectureMesh::getSolidMaterial(int solid)
{
	hyvector<InstancedMesh>::iterator imesh = meshInstances->begin() + solid;
	InstancedMesh& instance = *imesh;
	return instance.material;
}

int CArchitectureMesh::getFaceCount(int solid)
{
	hyvector<InstancedMesh>::iterator imesh = meshInstances->begin() + solid;
	InstancedMesh& instance = *imesh;
	return reinterpret_cast<CFastQuadrics*>(instance.mesh)->faceCount;
}

void CArchitectureMesh::getFace(int solid, int index, Vector& v0, Vector& v1, Vector& v2)
{
	hyvector<InstancedMesh>::iterator imesh = meshInstances->begin() + solid;
	InstancedMesh& instance = *imesh;

	CFastQuadrics* mesh = reinterpret_cast<CFastQuadrics*>(instance.mesh);
	FQ_Face& f = mesh->faces[index];
	FQ_Vertex& mv0 = mesh->vertices[f[0]];
	FQ_Vertex& mv1 = mesh->vertices[f[2]];
	FQ_Vertex& mv2 = mesh->vertices[f[1]];

	float isx = instance.size.x;
	float isy = instance.size.y;
	float isz = instance.size.z;
	if (instance.scopeType == SCOPE_BOX)
	{
		v0 = Vector_withValues(isx*(float)mv0.x, isy*(float)mv0.y, isz*(float)mv0.z);
		v1 = Vector_withValues(isx*(float)mv1.x, isy*(float)mv1.y, isz*(float)mv1.z);
		v2 = Vector_withValues(isx*(float)mv2.x, isy*(float)mv2.y, isz*(float)mv2.z);
	}
	else if (instance.scopeType == SCOPE_PRISM)
	{
		v0 = Vector_withValues(
			instance.size.x*(float)mv0.x,
			instance.size.y*(float)mv0.y,
			(1.0f - (float)mv0.y)*instance.size.z*(float)mv0.z + 0.5f*(float)mv0.y*instance.size.z);

		v1 = Vector_withValues(
			instance.size.x*(float)mv1.x,
			instance.size.y*(float)mv1.y,
			(1.0f - (float)mv1.y)*instance.size.z*(float)mv1.z + 0.5f*(float)mv1.y*instance.size.z);

		v2 = Vector_withValues(
			instance.size.x*(float)mv2.x,
			instance.size.y*(float)mv2.y,
			(1.0f - (float)mv2.y)*instance.size.z*(float)mv2.z + 0.5f*(float)mv2.y*instance.size.z);
	}
	else if (instance.scopeType == SCOPE_PRISM_LEFT)
	{
		v0 = Vector_withValues(
			instance.size.x*(float)mv0.x,
			instance.size.y*(float)mv0.y,
			(1.0f - (float)mv0.y)*instance.size.z*(float)mv0.z);

		v1 = Vector_withValues(
			instance.size.x*(float)mv1.x,
			instance.size.y*(float)mv1.y,
			(1.0f - (float)mv1.y)*instance.size.z*(float)mv1.z);

		v2 = Vector_withValues(
			instance.size.x*(float)mv2.x,
			instance.size.y*(float)mv2.y,
			(1.0f - (float)mv2.y)*instance.size.z*(float)mv2.z);
	}
	else if (instance.scopeType == SCOPE_PRISM_RIGHT)
	{
		v0 = Vector_withValues(
			instance.size.x*(float)mv0.x,
			instance.size.y*(float)mv0.y,
			(1.0f - (float)mv0.y)*instance.size.z*(float)mv0.z + (float)mv0.y*instance.size.z);

		v1 = Vector_withValues(
			instance.size.x*(float)mv1.x,
			instance.size.y*(float)mv1.y,
			(1.0f - (float)mv1.y)*instance.size.z*(float)mv1.z + (float)mv1.y*instance.size.z);

		v2 = Vector_withValues(
			instance.size.x*(float)mv2.x,
			instance.size.y*(float)mv2.y,
			(1.0f - (float)mv2.y)*instance.size.z*(float)mv2.z + (float)mv2.y*instance.size.z);
	}

	// Translate and rotate vertex by instance matrix
	v0 = Matrix_multiplyVector(instance.transform, v0);
	v1 = Matrix_multiplyVector(instance.transform, v1);
	v2 = Matrix_multiplyVector(instance.transform, v2);

	// Architecture is in meters while HyVoxel are in decimeters.
	// Scale by 10

	v0.x *= 10.0f;
	v0.y *= 10.0f;
	v0.z *= 10.0f;
	v1.x *= 10.0f;
	v1.y *= 10.0f;
	v1.z *= 10.0f;
	v2.x *= 10.0f;
	v2.y *= 10.0f;
	v2.z *= 10.0f;

	v0 = Matrix_multiplyVector(transform, v0);
	v1 = Matrix_multiplyVector(transform, v1);
	v2 = Matrix_multiplyVector(transform, v2);
}

#endif


static Name2Mesh sName2UserMesh;

void HyArchitectureInterfaceImpl::AddNameToUserMeshMapping(const char* name, void* userMesh)
{
	sName2UserMesh[name] = userMesh;
}

void HyArchitectureInterfaceImpl::ClearNameToUserMeshMappings()
{
	sName2UserMesh.clear();
}

bool HyArchitectureInterfaceImpl::ProcedurallyPositionUserMeshes(const char* prefabName, double boxSize[3], OutVector<InstancedMesh>*& outMeshes)
{
	OutVectorImpl<InstancedMesh>* outMeshesImpl = new OutVectorImpl<InstancedMesh>();

	bool ret = ProcedurallyPositionUserMeshesImpl(sName2UserMesh, prefabName, boxSize, *outMeshesImpl);

	outMeshesImpl->SetOutVectorValues();
	outMeshes = outMeshesImpl;

	return ret;
}

bool HyArchitectureInterfaceImpl::ProcedurallyPositionUserMeshesImpl(const Name2Mesh& name2Mesh, const char* prefabName, double boxSize[3], hyvector<InstancedMesh>& outMeshes)
{
	std::map<std::string, ArchPrefab>::iterator it = gArchPrefabs.find(prefabName);
	if (it == gArchPrefabs.end())
		return false;

	ArchPrefab& archPrefab = it->second;
	CGrammar* grammar = archPrefab.grammar;

	CInstanceCollector collector(name2Mesh, outMeshes);
	collector.type = PT_DEFAULT;
	collector.entity = NULL;
	CEvaluator evaluator;
	evaluator.instanceCreator = &collector;

	evaluator.root.scope.sx = (float)boxSize[0];
	evaluator.root.scope.sy = (float)boxSize[1];
	evaluator.root.scope.sz = (float)boxSize[2];

	evaluator.root.scope.ax = 0;
	evaluator.root.scope.ay = 0;
	evaluator.root.scope.az = 0;

	evaluator.root.scope.x = 0.0f;
	evaluator.root.scope.y = 0.0f;
	evaluator.root.scope.z = 0.0f;

	// Run evaluation
	// This will produce multiple calls into the CInstanceCollector callback interface
	evaluator.run(*grammar, "main", NULL, -1);

	return true;
}


static Name2Mesh sName2FastQuadrics;

bool HyArchitectureInterfaceImpl::IsDatMeshLoaded(const char* meshName)
{
	return sName2FastQuadrics.find(meshName) != sName2FastQuadrics.end();
}

int HyArchitectureInterfaceImpl::LoadDatMeshFromString(const char* meshName, const char* buf)
{
	if (IsDatMeshLoaded(meshName))
		return 0;

	std::istringstream iss(buf);
	CFastQuadrics* mesh = LoadMeshFromStream<>(iss);
	if (mesh)
	{
		sName2FastQuadrics[meshName] = mesh;
		return 0;
	}

	return -1;
}

int HyArchitectureInterfaceImpl::LoadDatMeshFromFile(const char* meshName, const char* fullPath)
{
	// early out
	if (IsDatMeshLoaded(meshName))
		return 0;

	std::string outStr;
	if (!ReadFileIntoString(fullPath, outStr))
		return -1;

	return LoadDatMeshFromString(meshName, outStr.c_str());
}

int HyArchitectureInterfaceImpl::GenerateMeshFromPrefab(const char* prefabName, const Vector& boxSize, OutputMesh& outMesh)
{
	double dBoxSize[] = { boxSize.x, boxSize.y, boxSize.z };

	hyvector<InstancedMesh> instancedMeshes;
	if (!ProcedurallyPositionUserMeshesImpl(sName2FastQuadrics, prefabName, dBoxSize, instancedMeshes))
		return -1;

	CArchitectureMesh archMesh(&instancedMeshes, Matrix_identity());

	int faceCount = 0;
	int instanceCount = archMesh.getSolidCount();
	for (int ins = 0; ins < instanceCount; ins++)
		faceCount += archMesh.getFaceCount(ins);

	int vertCount = faceCount * 3;

	OutVectorImpl<Vector>* overts = new OutVectorImpl<Vector>();
	OutVectorImpl<Vector>* onormals = new OutVectorImpl<Vector>();
	OutVectorImpl<int>* otris = new OutVectorImpl<int>();

	hyvector<Vector>& verts = *overts; verts.resize(vertCount);
	hyvector<Vector>& normals = *onormals; normals.resize(vertCount);
	hyvector<int>& tris = *otris; tris.resize(faceCount * 3);

	Vector tmpv0, tmpv1, tmpv2;

	int index = 0;
	for (int ins = 0; ins < instanceCount; ins++)
	{
		int faceCount = archMesh.getFaceCount(ins);
		for (int face = 0; face < faceCount; face++)
		{
			tris[index] = index;
			HyVoxel::Vector& v0 = verts[index];
			HyVoxel::Vector& n0 = normals[index];
			index++;

			tris[index] = index;
			HyVoxel::Vector& v1 = verts[index];
			HyVoxel::Vector& n1 = normals[index];
			index++;

			tris[index] = index;
			HyVoxel::Vector& v2 = verts[index];
			HyVoxel::Vector& n2 = normals[index];
			index++;

			// correct: CW for front face
			archMesh.getFace(ins, face, tmpv0, tmpv2, tmpv1);

			v0.x = tmpv0.x; v0.y = tmpv0.y; v0.z = tmpv0.z;
			v1.x = tmpv1.x; v1.y = tmpv1.y; v1.z = tmpv1.z;
			v2.x = tmpv2.x; v2.y = tmpv2.y; v2.z = tmpv2.z;

			n0 = HyVoxel::Vector::CrossProduct(v2 - v0, v1 - v0);
			n2 = n1 = n0;
		}
	}

	overts->SetOutVectorValues(); outMesh.verts = overts;
	onormals->SetOutVectorValues(); outMesh.normals = onormals;
	otris->SetOutVectorValues(); outMesh.tris = otris;

	return 0;
}

}

// TODO: liang
// to be removed
#if 1

#include "Voxel/ClipmapView.h"

class CTestStampMeshMaterials : public IMeshStampMaterialSource
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

void AddVoxelTriangleAt(float x, float y, float z, float cubeWidth, HyVoxel::Vector v0, HyVoxel::Vector v1, HyVoxel::Vector v2, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals)
{
	// float scale = cubeWidth * 0.5f;
	float scale = cubeWidth;

	// TODO: liang
	// flip y and z
	HyVoxel::Vector posOffset(x, z, y);
	// HyVoxel::Vector posOffset(x, y, z);

	int triOffset = verts.size();

	v0 = v0 * scale + posOffset;
	v1 = v1 * scale + posOffset;
	v2 = v2 * scale + posOffset;

	HyVoxel::Vector n = HyVoxel::Vector::CrossProduct(v2 - v0, v1 - v0);

	verts.push_back(v0); normals.push_back(n);
	tris.push_back(verts.size() - 1);
	verts.push_back(v1); normals.push_back(n);
	tris.push_back(verts.size() - 1);
	verts.push_back(v2); normals.push_back(n);
	tris.push_back(verts.size() - 1);
}

static void AddVoxelAt(float x, float y, float z, float cubeWidth, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals)
{
	// float scale = cubeWidth * 0.5f;
	float scale = cubeWidth;

	// TODO: liang
	// flip y and z
	HyVoxel::Vector posOffset(x, z, y);
	// HyVoxel::Vector posOffset(x, y, z);

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

	// CCW
	const int TRIS[] = {
		3,1,0,
		3,2,1,
		2,5,1,
		2,6,5,
		6,4,5,
		6,7,4,
		7,0,4,
		7,3,0,
		7,2,3,
		7,6,2,
		0,5,4,
		0,1,5,
	};

	int triOffset = verts.size();

	int triCount = sizeof(TRIS) / sizeof(int) / 3;
	for (int iTri = 0; iTri < triCount; ++iTri)
	{
		int i0 = TRIS[iTri * 3 + 0],
			i1 = TRIS[iTri * 3 + 1],
			i2 = TRIS[iTri * 3 + 2];

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

int GetMaterial(int x, int y, int z, int margin, BlockVoxelData* buffer[3][3][3])
{
	int i, j, k;
	if (x < margin) {
		x += 40 - margin;
		i = 0;
	}
	else if (x >= 40 + margin) {
		x += -40 - margin;
		i = 2;
	}
	else {
		x += -margin;
		i = 1;
	}
	if (y < margin) {
		y += 40 - margin;
		j = 0;
	}
	else if (y >= 40 + margin) {
		y += -40 - margin;
		j = 2;
	}
	else {
		y += -margin;
		j = 1;
	}
	if (z < margin) {
		z += 40 - margin;
		k = 0;
	}
	else if (z >= 40 + margin) {
		z += -40 - margin;
		k = 2;
	}
	else {
		z += -margin;
		k = 1;
	}
	Voxel vx;
	BlockVoxelData::Index vxIdx(x, y, z);
	if (buffer[i][j][k] == nullptr)
		return 0;
	buffer[i][j][k]->getVoxel(vxIdx, vx);
	int material = vx.material;
	return material;
}

void CalculateVisible(BlockDataLayer& blockLayer, int lod, int cellX, int cellY, int cellZ, int(&visible)[40][40][40], int(&material)[40][40][40])
{
	BlockVoxelData* voxel_data[3][3][3];
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			for (int k = 0; k < 3; ++k)
				voxel_data[i][j][k] = blockLayer.fetchData(packCellId(lod, cellX + i - 1, cellY + j - 1, cellZ + k - 1), false);

	bool is_block[42][42][42];
	int m[42][42][42];
	for (int x = 0; x < 42; ++x)
		for (int y = 0; y < 42; ++y)
			for (int z = 0; z < 42; ++z) {
				m[x][y][z] = GetMaterial(x, y, z, 1, voxel_data);
				is_block[x][y][z] = (m[x][y][z] != 0);
			}

	for (int x = 1; x < 41; x++) {
		for (int y = 1; y < 41; y++) {
			for (int z = 1; z < 41; z++) {
				// x+ x- y+ y- z+ z-
				int v = 0x00;
				if (is_block[x][y][z]) {
					v = 0x3f;
					if (is_block[x + 1][y][z])
						v &= 0x1f;
					if (is_block[x - 1][y][z])
						v &= 0x2f;
					if (is_block[x][y + 1][z])
						v &= 0x37;
					if (is_block[x][y - 1][z])
						v &= 0x3b;
					if (is_block[x][y][z + 1])
						v &= 0x3d;
					if (is_block[x][y][z - 1])
						v &= 0x3e;
				}
				visible[x - 1][y - 1][z - 1] = v;
				material[x - 1][y - 1][z - 1] = m[x][y][z];
			}
		}
	}
}

void VisualizeBlock(BlockDataLayer& blockLayer, int lod, int cellX, int cellY, int cellZ, bool showBlockBoundary, BlockVoxelData* blockData, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, bool hideInner)
{
	const float oneOnDiv = 1.0f / 40.0f;

	double scale = (1 << lod) * CELL_WIDTH * 10.0f; // unit: centimeter

	float lodXOffset = (lod - CELL_LOD_MIN) * 3000.0f;

	int visible[40][40][40];
	int material[40][40][40];

	CalculateVisible(blockLayer, lod, cellX, cellY, cellZ, visible, material);

	for (int x = 0; x < blockData->cDimX; ++x)
	{
		for (int y = 0; y < blockData->cDimY; ++y)
		{
			for (int z = 0; z < blockData->cDimZ; ++z)
			{
				BlockVoxelData::Index vxIdx(x, y, z);
				Voxel vx;
				blockData->getVoxel(vxIdx, vx);

				float xoff = (cellX + x * oneOnDiv) * scale + lodXOffset;
				float yoff = (cellY + y * oneOnDiv) * scale;
				float zoff = (cellZ + z * oneOnDiv) * scale;
				float cubeWidth = scale * oneOnDiv;

				if (!hideInner) {
					if (material[x][y][z])
						AddVoxelAt(xoff, yoff, zoff, cubeWidth / 2.0, verts, tris, normals);
				}
				else if (visible[x][y][z]) {
					if (visible[x][y][z] & 0x20) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(1, 0, 0), HyVoxel::Vector(1, 0, 1), HyVoxel::Vector(1, 1, 1), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(1, 0, 0), HyVoxel::Vector(1, 1, 1), HyVoxel::Vector(1, 1, 0), verts, tris, normals);
					}
					if (visible[x][y][z] & 0x10) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(0, 1, 1), HyVoxel::Vector(0, 0, 1), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(0, 1, 0), HyVoxel::Vector(0, 1, 1), verts, tris, normals);
					}
					if (visible[x][y][z] & 0x02) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 1, 0), HyVoxel::Vector(1, 1, 1), HyVoxel::Vector(0, 1, 1), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 1, 0), HyVoxel::Vector(1, 1, 0), HyVoxel::Vector(1, 1, 1), verts, tris, normals);
					}
					if (visible[x][y][z] & 0x01) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(0, 0, 1), HyVoxel::Vector(1, 0, 1), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(1, 0, 1), HyVoxel::Vector(1, 0, 0), verts, tris, normals);
					}
					if (visible[x][y][z] & 0x08) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 1), HyVoxel::Vector(0, 1, 1), HyVoxel::Vector(1, 1, 1), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 1), HyVoxel::Vector(1, 1, 1), HyVoxel::Vector(1, 0, 1), verts, tris, normals);
					}
					if (visible[x][y][z] & 0x04) {
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(1, 1, 0), HyVoxel::Vector(0, 1, 0), verts, tris, normals);
						AddVoxelTriangleAt(xoff, yoff, zoff, cubeWidth, HyVoxel::Vector(0, 0, 0), HyVoxel::Vector(1, 0, 0), HyVoxel::Vector(1, 1, 0), verts, tris, normals);
					}
				}
			}
		}
	}
}

void VisualizeBlockLayer(BlockDataLayer& blockLayer, int layerMax, bool showBlockBoundary, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals)
{
	for (auto& kv : blockLayer.blockCache)
	{
		CellId cellId = kv.first;
		BlockVoxelData* blockData = kv.second;

		int level, cellX, cellY, cellZ;
		unpackCellId(cellId, level, cellX, cellY, cellZ);

		if (level <= layerMax)
			VisualizeBlock(blockLayer, level, cellX, cellY, cellZ, showBlockBoundary, blockData, verts, tris, normals, false);
	}
}

void TestVoxelization(const char* datPath, const HyVoxel::Vector& boxSize, std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, bool isCubeContour)
{
	ClipmapView clipmapView;

	HyArchitectureInterfaceImpl impl;
	impl.LoadDatMeshFromFile(datPath, datPath);

	InstancedMesh insMesh;
	insMesh.mesh = sName2FastQuadrics[datPath];
	insMesh.material = VOXEL_ONLY_COORDS;
	insMesh.entity = NULL;
	insMesh.transform = Matrix_identity();
	Matrix_translate(&(insMesh.transform), 0.0f, 1.0f, 5.0f);
	Vector sizeVec = { boxSize.x, boxSize.y, boxSize.z };
	insMesh.size = sizeVec;
	insMesh.scopeType = SCOPE_BOX;

	hyvector<InstancedMesh> insMeshes;
	insMeshes.push_back(insMesh);

	CArchitectureMesh archMesh(&insMeshes, Matrix_identity());

	CTestStampMeshMaterials stampMat;
	// TODO: liang
	// Just hardcode it to be 1.
	stampMat.id = 1;

	double worldPos[3] = { 0.0, 0.0, 0.0 };

	std::set<CellId> changedCellSet;
	stampMeshQEF(&clipmapView.blockLayer,
		&archMesh,
		&stampMat,
		worldPos,
		Matrix_identity(),
		&changedCellSet);

	std::vector<CellId> changedCells(changedCellSet.begin(), changedCellSet.end());
	clipmapView.blockLayer.UpdateLod(changedCells.data(), changedCells.size());

	if (isCubeContour)
	{
	
		VisualizeBlockLayer(clipmapView.blockLayer, 2, false, verts, tris, normals);
	}
	else
	{
		clipmapView.CalculateNextSceneIfNeeded();
		clipmapView.ContourLoop(verts, tris, normals);
	}
}

void AddTestCube(CFastQuadrics* mesh, float x, float y, float z, float s)
{
	CFastQuadrics newMesh = CFastQuadrics();
	newMesh.allocate(8, 12);
	newMesh.addVertex(x, y, z, 0);
	newMesh.addVertex(x + s, y, z, 0);
	newMesh.addVertex(x, y + s, z, 0);
	newMesh.addVertex(x + s, y + s, z, 0);
	newMesh.addVertex(x, y, z + s, 0);
	newMesh.addVertex(x + s, y, z + s, 0);
	newMesh.addVertex(x, y + s, z + s, 0);
	newMesh.addVertex(x + s, y + s, z + s, 0);
	newMesh.addFace(0, 2, 3, 0);
	newMesh.addFace(3, 1, 0, 0);
	newMesh.addFace(4, 5, 7, 0);
	newMesh.addFace(7, 6, 4, 0);
	newMesh.addFace(0, 1, 5, 0);
	newMesh.addFace(5, 4, 0, 0);
	newMesh.addFace(1, 3, 7, 0);
	newMesh.addFace(7, 5, 1, 0);
	newMesh.addFace(3, 2, 6, 0);
	newMesh.addFace(6, 7, 3, 0);
	newMesh.addFace(2, 0, 4, 0);
	newMesh.addFace(4, 6, 2, 0);
	mesh->append(newMesh);
}

void TestMeshToVoxel(std::vector<HyVoxel::Vector>& verts, std::vector<int32_t>& tris, std::vector<HyVoxel::Vector>& normals, const HyVoxel::Vector& boxSize, bool isCubeContour)
{
	ClipmapView clipmapView;

	CFastQuadrics* mesh = VF_NEW CFastQuadrics();
	mesh->allocate(verts.size(), tris.size() / 3);
	for (int i = 0; i < verts.size(); i++) {
		mesh->addVertex(verts[i].x, verts[i].z, verts[i].y, 0);
	}
	for (int i = 0; i < tris.size(); i += 3) {
		mesh->addFace(tris[i], tris[i + 1], tris[i + 2], 0);
	}

	InstancedMesh insMesh;
	insMesh.mesh = mesh;
	insMesh.material = VOXEL_ONLY_COORDS;
	insMesh.entity = NULL;
	insMesh.transform = Matrix_identity();
	Matrix_translate(&(insMesh.transform), 0.0f, 1.0f, 5.0f);
	Vector sizeVec = { boxSize.x, boxSize.y, boxSize.z };
	insMesh.size = sizeVec;
	insMesh.scopeType = SCOPE_BOX;

	hyvector<InstancedMesh> insMeshes;
	insMeshes.push_back(insMesh);

	CArchitectureMesh archMesh(&insMeshes, Matrix_identity());

	CTestStampMeshMaterials stampMat;
	// TODO: liang
	// Just hardcode it to be 1.
	stampMat.id = 1;

	double worldPos[3] = { 0.0, 0.0, 0.0 };

	TVFSet<CellId> changedCells;

	//stampMeshQEF(&clipmapView.blockLayer,
	//	&archMesh,
	//	&stampMat,
	//	worldPos,
	//	Matrix_identity(),
	//	&changedCells);



	verts.clear();
	tris.clear();
	normals.clear();

	if (isCubeContour)
	{
		VisualizeBlockLayer(clipmapView.blockLayer, 2, false, verts, tris, normals);
	}
	else
	{
		clipmapView.CalculateNextSceneIfNeeded();
		clipmapView.ContourLoop(verts, tris, normals);
	}

	VF_DELETE mesh;
}

#endif

