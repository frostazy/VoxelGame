/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#include "HyVoxelPrivatePCH.h"

#include "HyVoxelConfig.h"
#include "MaterialLibrary.h"
#include "Voxel/Visualize/CellMesh.h"
#include "Voxel/Physics/PhysicsMaterials.h"

#include <windows.h>
#include <algorithm>

#undef min

#define _USE_MATH_DEFINES
#include <math.h>

using namespace HyVoxel;
using namespace std;

namespace HyVoxel {

CMaterialLibrary::CMaterialLibrary()
	: mapSize(0)
	, materialCount(0)
	, materialIndex(NULL)
	, billboardPack()
	, maxHomogeneousLOD(CELL_LOD_MAX)
{
}

CMaterialLibrary::~CMaterialLibrary(void)
{
	for (int i = 0; i < materialCount; ++i)
	{
		//VF_DELETE materialIndex[i].subMaterials;
	}
	VF_FREE(materialIndex);
}

void CMaterialLibrary::init(int mapSize, int materialCount)
{
	this->mapSize = mapSize;
	this->materialCount = materialCount;
	materialIndex = VF_ALLOC(CMaterial, materialCount);
	::memset(materialIndex, 0, sizeof(materialIndex[0]) * materialCount);
}

//.ini tool
double GetPrivateProfileFloat(LPCSTR lpAppName, LPCSTR lpKeyName, LPCSTR lpDefault, LPCSTR lpFileName)
{
	char str[MAX_PATH];
	GetPrivateProfileStringA(lpAppName, lpKeyName, lpDefault, str, MAX_PATH, lpFileName);
	return atof(str);
};

bool readMaterialDefinitions(const char* matDefPath, CMaterialLibrary& materialLibrary)
{
	char projectfilename[MAX_PATH];
	sprintf_s(projectfilename, MAX_PATH, "%s", matDefPath);
	materialLibrary.materialCount = GetPrivateProfileIntA("Project", "materials", -1, projectfilename);
	materialLibrary.materialIndex = VF_ALLOC(CMaterial, materialLibrary.materialCount);
	char defultSimplificationError[MAX_PATH];
	GetPrivateProfileStringA("Project", "simplificationError", "0.01", defultSimplificationError, MAX_PATH, projectfilename);
	int defaultResolution = GetPrivateProfileIntA("Project", "resolution", 2, projectfilename);
	materialLibrary.mapSize = GetPrivateProfileIntA("Project", "mapsize", 512, projectfilename);
	char defaultSnowLine[MAX_PATH];
	GetPrivateProfileStringA("Project", "snowLine", "60000.0", defaultSnowLine, MAX_PATH, projectfilename);

	int lastMeshId = 0;
	int lastInstanceSlice = 0;

	for (int i = 0; i < materialLibrary.materialCount; i++)
	{
		CMaterial& m = materialLibrary.materialIndex[i];
		char sectionId[MAX_PATH];
		sprintf_s(sectionId, MAX_PATH, "Material%d", i);
		m.id = GetPrivateProfileIntA(sectionId, "id", 0, projectfilename);
		m.medium = GetPrivateProfileIntA(sectionId, "medium", 0, projectfilename);
		m.billboard = GetPrivateProfileIntA(sectionId, "billboard", -1, projectfilename);
		m.billboardAngleMin = GetPrivateProfileFloat(sectionId, "billboardAngleMin", "-1.0", projectfilename);
		m.billboardAngleMax = GetPrivateProfileFloat(sectionId, "billboardAngleMax", "1.0", projectfilename);
		m.billboardType = GetPrivateProfileIntA(sectionId, "billboardType", 0, projectfilename);
		m.applyColor = GetPrivateProfileIntA(sectionId, "applyColor", 0, projectfilename);
		m.billboardSize = GetPrivateProfileFloat(sectionId, "billboardSize", "1.0", projectfilename);
		m.billboardRigidity = GetPrivateProfileFloat(sectionId, "billboardRigidity", "1.0", projectfilename);
		m.billboardDensity = GetPrivateProfileFloat(sectionId, "billboardDensity", "0.05", projectfilename);
		m.angleMin = GetPrivateProfileFloat(sectionId, "angleMin", "-1.0", projectfilename);
		m.angleMax = GetPrivateProfileFloat(sectionId, "angleMax", "1.0", projectfilename);
		m.resolution = GetPrivateProfileIntA(sectionId, "resolution", defaultResolution, projectfilename);
		m.maxlevel = GetPrivateProfileIntA(sectionId, "maxlevel", 20, projectfilename);
		m.selfIllumination = GetPrivateProfileFloat(sectionId, "selfIllumination", "0.0", projectfilename);
		m.nearFreq = GetPrivateProfileFloat(sectionId, "nearFreq", "1.0", projectfilename);
		m.farFreq = GetPrivateProfileFloat(sectionId, "farFreq", "1.0", projectfilename);
		m.macroFreq = GetPrivateProfileFloat(sectionId, "macroFreq", "1.0", projectfilename);
		m.faceSmoothAngle = cos(GetPrivateProfileFloat(sectionId, "faceSmoothAngle", "90", projectfilename)*M_PI / 180.0);
		m.simplificationError = GetPrivateProfileFloat(sectionId, "simplificationError", defultSimplificationError, projectfilename);
		m.carved = GetPrivateProfileIntA(sectionId, "carved", i, projectfilename);
		m.displacementSize = GetPrivateProfileFloat(sectionId, "displacementSize", "1.0", projectfilename);
		m.displacementShift = GetPrivateProfileFloat(sectionId, "displacementShift", "0.0", projectfilename);
		m.displacementFreq = GetPrivateProfileFloat(sectionId, "displacementFreq", "1.0", projectfilename);
		m.placementType = GetPrivateProfileIntA(sectionId, "placementType", 0, projectfilename);
		m.placementOctaves = GetPrivateProfileIntA(sectionId, "placementOctaves", 1, projectfilename);
		m.placementFreq = GetPrivateProfileFloat(sectionId, "placementFreq", "0.0", projectfilename);
		m.placementStep = GetPrivateProfileFloat(sectionId, "placementStep", "0.5", projectfilename);
		m.placementLacunarity = GetPrivateProfileFloat(sectionId, "placementLacunarity", "2.0", projectfilename);
		m.placementPhaseX = GetPrivateProfileFloat(sectionId, "placementPhaseX", "0.0", projectfilename);
		m.placementPhaseY = GetPrivateProfileFloat(sectionId, "placementPhaseY", "0.0", projectfilename);
		m.placementPhaseZ = GetPrivateProfileFloat(sectionId, "placementPhaseZ", "0.0", projectfilename);
		m.placementScaleX = GetPrivateProfileFloat(sectionId, "placementScaleX", "1.0", projectfilename);
		m.placementScaleY = GetPrivateProfileFloat(sectionId, "placementScaleY", "1.0", projectfilename);
		m.placementScaleZ = GetPrivateProfileFloat(sectionId, "placementScaleZ", "1.0", projectfilename);
		m.placementClampMin = GetPrivateProfileFloat(sectionId, "placementClampMin", "0.0", projectfilename);
		m.placementClampMax = GetPrivateProfileFloat(sectionId, "placementClampMax", "1.0", projectfilename);
		m.instanceMaskMode = GetPrivateProfileIntA(sectionId, "instanceMaskMode", 0, projectfilename);
		m.heightMin = GetPrivateProfileFloat(sectionId, "heightMin", "0.0", projectfilename);
		m.heightMax = GetPrivateProfileFloat(sectionId, "heightMax", "-1.0", projectfilename);
		m.snowLine = GetPrivateProfileFloat(sectionId, "snowLine", defaultSnowLine, projectfilename);
		m.wearMaterial = GetPrivateProfileIntA(sectionId, "wear", i + 1, projectfilename);
		m.transparent = GetPrivateProfileIntA(sectionId, "transparent", 0, projectfilename) == 1;
		m.layeredMaterial = GetPrivateProfileIntA(sectionId, "layeredMaterial", m.id, projectfilename);
		m.layeredMaterialAngle = GetPrivateProfileFloat(sectionId, "layeredMaterialAngle", "1.0f", projectfilename);
		GetPrivateProfileStringA(sectionId, "stepSound", "", m.stepSoundId, MAX_PATH, projectfilename);
		GetPrivateProfileStringA(sectionId, "digSound", "", m.digSoundId, MAX_PATH, projectfilename);
		m.stepSoundCount = GetPrivateProfileIntA(sectionId, "stepSoundCount", 0, projectfilename);
		m.digSoundCount = GetPrivateProfileIntA(sectionId, "digSoundCount", 0, projectfilename);
		m.colorH = GetPrivateProfileFloat(sectionId, "colorH", "0.0", projectfilename);
		m.colorS = GetPrivateProfileFloat(sectionId, "colorS", "0.0", projectfilename);
		m.colorV = GetPrivateProfileFloat(sectionId, "colorV", "1.0", projectfilename);
		m.colorFreqX = GetPrivateProfileFloat(sectionId, "colorFreqX", "0.0001", projectfilename);
		m.colorFreqY = GetPrivateProfileFloat(sectionId, "colorFreqY", "0.001", projectfilename);
		m.colorFreqZ = GetPrivateProfileFloat(sectionId, "colorFreqZ", "0.0001", projectfilename);
		m.colorDeltaH = GetPrivateProfileFloat(sectionId, "colorDeltaH", "0.0", projectfilename);
		m.colorDeltaS = GetPrivateProfileFloat(sectionId, "colorDeltaS", "0.0", projectfilename);
		m.colorDeltaV = GetPrivateProfileFloat(sectionId, "colorDeltaV", "0.0", projectfilename);
		m.useMacroColor = GetPrivateProfileIntA(sectionId, "useMacroColor", 0, projectfilename) == 1;
		m.displacementMap = NULL;
		m.macroColor = NULL;
		m.homogeneous = GetPrivateProfileIntA(sectionId, "homogeneous", 0, projectfilename) == 1;
		m.blend = GetPrivateProfileIntA(sectionId, "blend", 0, projectfilename) == 1;
		m.maxEdgeLength = GetPrivateProfileFloat(sectionId, "maxEdgeLength", "1.0", projectfilename);
		m.materialType = GetPrivateProfileIntA(sectionId, "materialType", Physics::MAT_DEFAULT, projectfilename);

		char colorList[256];
		GetPrivateProfileStringA(sectionId, "submaterial", "", colorList, 256, projectfilename);
		char* nextMat = colorList;
		if (*nextMat != 0)
		{
			m.subMaterials = new TVFVector<int>();
			while (*nextMat != 0)
			{
				char* separator = nextMat;
				while (*separator != ',' && *separator != 0)
				{
					separator++;
				}
				bool ended = *separator == 0;
				*separator = 0;
				int submat = atoi(nextMat);
				if (!ended)
				{
					nextMat = separator + 1;
				}
				else
				{
					nextMat = separator;
				}
				if (submat != 0)
				{
					m.subMaterials->push_back(submat);
				}
			}
		}
		else
		{
			m.subMaterials = NULL;
		}

		// read material instances
		m.instanceCount = std::min(CMaterial::MAX_INSTANCES, (int)GetPrivateProfileIntA(sectionId, "instanceCount", 0, projectfilename));
		for (int j = 0; j < m.instanceCount; j++)
		{
			char propertyName[MAX_PATH];
			sprintf_s(propertyName, MAX_PATH, "instance%d", j);
			char instanceId[MAX_PATH];
			GetPrivateProfileStringA(sectionId, propertyName, "", instanceId, MAX_PATH, projectfilename);

			sprintf_s(propertyName, MAX_PATH, "instance%dDensity", j);
			m.instances[j].density = GetPrivateProfileFloat(sectionId, propertyName, "0.05", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dAngleMin", j);
			m.instances[j].minAngle = GetPrivateProfileFloat(sectionId, propertyName, "0.0", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dAngleMax", j);
			m.instances[j].maxAngle = GetPrivateProfileFloat(sectionId, propertyName, "1.0", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dSizeMin", j);
			m.instances[j].minSize = GetPrivateProfileFloat(sectionId, propertyName, "5", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dSizeMax", j);
			m.instances[j].maxSize = GetPrivateProfileFloat(sectionId, propertyName, "5", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dNoiseFreq", j);
			m.instances[j].noiseFreq = GetPrivateProfileFloat(sectionId, propertyName, "0.005", projectfilename);
			sprintf_s(propertyName, MAX_PATH, "instance%dNoiseShift", j);
			m.instances[j].noiseShift = GetPrivateProfileFloat(sectionId, propertyName, "0.0", projectfilename);



			map<string, int>::iterator inst = instanceNameIdMap.find(instanceId);
			if (inst != instanceNameIdMap.end())
			{
				m.instances[j].id = inst->second;
			}
			else
			{
				instanceNameIdMap[instanceId] = lastMeshId;
				instanceIdNameMap[lastMeshId] = instanceId;
				m.instances[j].id = lastMeshId;
				lastMeshId++;
			}
		}

		GetPrivateProfileStringA(sectionId, "diffuse", "", m.diffuse, MAX_PATH, projectfilename);
		GetPrivateProfileStringA(sectionId, "normal", "", m.normal, MAX_PATH, projectfilename);
		GetPrivateProfileStringA(sectionId, "displacement", "", m.displacement, MAX_PATH, projectfilename);
		GetPrivateProfileStringA(sectionId, "displacementNormal", "", m.displacementNormal, MAX_PATH, projectfilename);
	}
	return true;
}


}