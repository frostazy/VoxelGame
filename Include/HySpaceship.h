#pragma once

// Please don't include this file directly.
// Please include "HyVoxelLib.h"

namespace HyVoxel {

	#pragma pack(push, 8)

	struct ShipComponent
	{
		OutVector<Vector>*		verts = NULL;
		OutVector<Vector>*		normals = NULL;
		OutVector<int>*			tris = NULL;
		OutVector<Vector2f>*	uvs = NULL;

		~ShipComponent()
		{
			if (verts)
			{
				verts->Release();
				verts = NULL;
			}
			if (normals)
			{
				normals->Release();
				normals = NULL;
			}
			if (tris)
			{
				tris->Release();
				tris = NULL;
			}
			if (uvs)
			{
				uvs->Release();
				uvs = NULL;
			}
		}
	};

	enum Direction {
		FRONT = 0,
		BACK = 3,
		LEFT = 1,
		RIGHT = 4,
		TOP = 2,
		BOTTOM = 5
	};

	struct Spaceship
	{
		enum ComponentType {
			HULL = 0,
			CABIN,
			ACCESSORY,
			ENGINE_EMIT,
			WINDOW_EMIT,
			HULL_SPHERE,
			LIGHT,
			COMPONENT_COUNT,
		};

		enum CabinType {
			MAIN,			// Possible antenna, sphere, disc, cylinder added
			FRONT_END,		// Possible antenna added
			BACK_END,		// Main engine added
			RIBBED,			// Nothing added
			SIDE_BASE,		// Possible engine added; no sphere
			SIDE,			// Possible engine added; no sphere
			SIDE_END,		// Possible engine added; no sphere, antenna
							//...
		};

		enum FaceType {
			INVALID,
			UNSPECIFIED,
			SIDED,
			ENGINE,
			GRID,
			ANTENNA,
			WEAPON,
			SPHERE,
			DISC,
			CYLINDER
		};

		static const int CABIN_COUNT = 100;
		ShipComponent components[COMPONENT_COUNT];

	};

	class IHySpaceshipInterface : public IHyVoxelInterfaceBase
	{
	public:
		virtual void GenerateSpaceship(Vector center, Vector rotation, float spaceship_size, unsigned int seed, Spaceship& outship) = 0;

		// TODO: liang
		//  to be moved into somewhere like help function interface.
		virtual void GenerateSphere(const Vector& center, float radius, int num, OutputMesh& omesh) = 0;
        virtual void GenerateTriangle(const Vector& v0, const Vector& v1, const Vector& v2, OutputMesh& omesh) = 0;
	};

	#pragma pack(pop)
}