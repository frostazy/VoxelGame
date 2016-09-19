/************************************************************
* (C) Voxel Farm Inc. 2015
*/

#pragma once

// 260 was taken from windef.h
#ifndef MAX_PATH
	#define MAX_PATH  260
#endif

#include "HyVoxelConfig.h"

namespace HyVoxel
{
	/// A file needed for a material was not found
	const int ERROR_MATERIAL_FILE_NOT_FOUND = 301;

	/// Describes a material instance
	struct InstanceDescriptor
	{
		/// ID for the mesh to be used
		int id;
		/// Controls how dense the instance distribution over the material is
		double density;
		/// Minimum angle at which the instance appears
		double minAngle;
		/// Maximum angle at which the instance appears
		double maxAngle;
		/// Minimum size for the instance
		double minSize;
		/// Maximum size for the instance
		double maxSize;
		/// Frequency for the instance distribution noise
		double noiseFreq;
		/// Amplitude shift for the instance distribution noise
		double noiseShift;
	};

	class IMacroColorSource
	{
	public:
		virtual void getColor(double x, double y, double z, unsigned char& r, unsigned char& g, unsigned char& b) = 0;
	};

	/// A material definition for the HyVoxel.com engine
	struct CMaterial
	{
		/// Identifies the material family
		int id;
		/// Indentifies to which medium the material belongs to
		int medium;
		/// Path to a BMP image for the diffuse map
		char diffuse[MAX_PATH];
		/// Index in the texture array for diffuse map
		int diffuseMapIndex;
		/// Path to a BMP image for the normal map
		char normal[MAX_PATH];
		/// Index in the texture array for diffuse map
		int normalMapIndex;
		/// Index in the texture array for displacement map intensity
		int displacementMapIndex;
		/// Index in the texture array for displacement map normal
		int displacementNormalMapIndex;
		/// Controls how many polygons the material takes on screen. Zero has the highest polygon density possible.
		int resolution;
		/// Maximum level at which the material will be displayed
		int maxlevel;
		/// Maximum simplification error allowed for the material
		double simplificationError;
		/// Which material will be set when this material is carved
		int carved;
		/// Minimum slope angle for the material to be applied
		double angleMin;
		/// Maximum slope angle for the material to be applied
		double angleMax;
		/// Neighboring faces will have their normals smoothen if they are below this separation angle.
		double faceSmoothAngle;
		/// A value from zero to one specifying how much lightign conditions affect the appearance of the material. A value of one makes the material not affected by external light.
		double selfIllumination;
		/// Texture tiling frequency for the near range of the material
		double nearFreq;
		/// Texture tiling frequency for the mid range of the material
		double farFreq;
		/// Texture tiling frequency for the far range of the material
		double macroFreq;
		/// Tints the material using the blue channel as a mask and the red channel as an intensity
		unsigned int applyColor;
		/// Altitude at which the material begins to show snow
		double snowLine;
		/// Maximum edge length
		double maxEdgeLength;

		bool homogeneous;
		bool blend;

		int wearMaterial;
		int layeredMaterial;
		double layeredMaterialAngle;
		bool transparent;

		/// Identifier for the billboard texture to be placed on top of the material
		int billboard;
		/// Billboard type. Zero is standing up, like grass. One is around the object, like foliage.
		int billboardType;
		/// Minimum slope angle for the billboards to be applied
		double billboardAngleMin;
		/// Maximum slope angle for the billboards to be applied
		double billboardAngleMax;
		/// Scale of the billboard. Default is one.
		double billboardSize;
		/// Controls how much the billboard is affected by wind
		double billboardRigidity;
		/// Controls billboard density. Smaller values will produce higher density;
		double billboardDensity;

		bool useMacroColor;
		IMacroColorSource* macroColor;

		/// Path to a displacement map
		char displacement[MAX_PATH];
		/// Path to a normal map for the material displacement
		char displacementNormal[MAX_PATH];
		/// Size of the displacement effect
		double displacementSize;
		/// Additional value for the displacement map
		double displacementShift;
		/// Frequency of the displacement tiling
		double displacementFreq;
		/// Contains the values read from the displacement map
		unsigned char* displacementMap;

		/// Type of placement mask. Only one type is currently supported, zero, which is a Perlin noise
		int placementType;
		/// Frequency of the placement noise
		double placementFreq;
		/// Amplitude multiplier of the placement noise
		double placementStep;
		/// Frequency multiplier of the placement noise
		double placementLacunarity;
		/// Phase along X axis of placement noise
		double placementPhaseX;
		/// Phase along Y axis of placement noise
		double placementPhaseY;
		/// Phase along Z axis of placement noise
		double placementPhaseZ;
		/// Number of octaves in placement noise
		int placementOctaves;
		/// Scale along X axis of placement noise
		double placementScaleX;
		/// Scale along Y axis of placement noise
		double placementScaleY;
		/// Scale along Z axis of placement noise
		double placementScaleZ;
		/// Placement noise clamp minimum
		double placementClampMin;
		/// Placement noise clamp maximum
		double placementClampMax;

		/// Used to mask the material by using the instance map
		int instanceMaskMode;

		/// Minimum world height at which the material could appear
		double heightMin;
		/// Maximum world height at which the material could appear
		double heightMax;

		/// Used for material HSV colorization, determines the H (hue) component for the base color
		double colorH;
		/// Used for material HSV colorization, determines the S (saturation) component for the base color
		double colorS;
		/// Used for material HSV colorization, determines the V (value) component for the base color
		double colorV;
		/// Used for material HSV colorization, determines the frequency at which components change over the X axis
		double colorFreqX;
		/// Used for material HSV colorization, determines the frequency at which components change over the Y axis
		double colorFreqY;
		/// Used for material HSV colorization, determines the frequency at which components change over the Z axis
		double colorFreqZ;
		/// Used for material HSV colorization, determines the range of variation for the H (hue) component
		double colorDeltaH;
		/// Used for material HSV colorization, determines the range of variation for the S (saturation) component
		double colorDeltaS;
		/// Used for material HSV colorization, determines the range of variation for the V (value) component
		double colorDeltaV;

		TVFVector<int>* subMaterials;

		/// Maximum amount of different material instance types
		static const int MAX_INSTANCES = 20;
		/// Number of material instance types
		int instanceCount;
		/// Material instance definition
		InstanceDescriptor instances[MAX_INSTANCES];

		// sound
		char stepSoundId[MAX_PATH];
		int stepSoundCount;
		char digSoundId[MAX_PATH];
		int digSoundCount;

		/************************************
		*	Permitted Material Types:		*
		*	(see PhysicsMaterials.h)		*
		*									*
		*		MATERIAL_CORK		=	0	*
		*		MATERIAL_WOOD		=	1	*
		*		MATERIAL_ICE		=	2	*
		*		MATERIAL_WATER		=	3	*
		*		MATERIAL_PLASTIC	=	4	*
		*		MATERIAL_CONCRETE	=	5	*
		*		MATERIAL_ROCK		=	6	*
		*		MATERIAL_METAL		=	7	*
		*		MATERIAL_HELIUM		=	8	*
		*		MATERIAL_DEFAULT	=	9	*
		************************************/
		int materialType;
	};

	/// Texture coordinate packing for billboard atlases
	struct CBillboardLODPacking
	{
		float LOD0 [6][4][2];
		float LOD1 [3][4][2];
		float LOD1B[8][4][2]; // for type2 billboards
		float LOD2 [2][4][2];
		float rowHeight;
	};

	/// Contains all materials available in the system
	class CMaterialLibrary
	{
	public:
		CMaterialLibrary();
		virtual ~CMaterialLibrary();

		void init(int mapSize, int materialCount);
        CMaterial& getMaterial( int index ) { return materialIndex[index]; }
        virtual int getMaterialMapping(int material) const { return material; }

		inline CMaterialLibrary& operator=(const CMaterialLibrary& m)
		{
			mapSize = m.mapSize;
			materialCount = m.materialCount;
			materialIndex = VF_ALLOC(CMaterial, materialCount);
			memcpy(materialIndex, m.materialIndex, materialCount*sizeof(CMaterial));
			return *this;
		};

	public:
		/// Defines the map size shared by all the bitmaps in the library. Must be a power of 2
		int mapSize;
		/// An index to all materials
		CMaterial* materialIndex;
		/// Number of materials
		int materialCount;
		/// Billboard atlas packing data
		CBillboardLODPacking billboardPack;
		int maxHomogeneousLOD;
	};

	bool readMaterialDefinitions(const char* matDefPath, CMaterialLibrary& materialLibrary);

}

